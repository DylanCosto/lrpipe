//! `lrpipe run`'s full-screen display, on a terminal: every stage with its progress, the run's
//! with its memory, and the end of a stage's log.
//!
//! A stage's progress is what it writes to its progress file ([`lrpipe_reads::progress`]:
//! lrpipe-snv and the read pass do), or else its time so far against what it took when it last
//! succeeded in the output directory, or else only that it runs. The display runs on a thread of
//! its own and reads the stages' files a few times a second; the run's events come over a channel.
//!
//! The run's lines (those `lrpipe run` prints without the display) are kept and printed when the
//! display ends, so the terminal is left with the record a run without it leaves. `q` ends the
//! display and the run goes on printing them; Ctrl-C stops the run, as SIGINT does (the terminal
//! sends no signal while the display has it).

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Once, mpsc};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

use ratatui::Terminal;
use ratatui::backend::CrosstermBackend;
use ratatui::crossterm::event::{self as term_event, Event as TermEvent, KeyCode, KeyEventKind, KeyModifiers};
use ratatui::crossterm::terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode};
use ratatui::crossterm::{cursor, execute};
use ratatui::layout::{Constraint, Layout as Split, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph};

use crate::manifest::{StageRecord, StageStatus};
use crate::plan::{Layout, Stage, StageName};
use crate::run::Event;

/// What the display's first lines say of the run.
#[derive(Clone, Debug)]
pub struct Header {
    /// The sample.
    pub sample: String,
    /// `hifi` or `ont`.
    pub platform: String,
    /// The output directory.
    pub output: PathBuf,
    /// The threads, as the run's first line gives them.
    pub threads: String,
    /// The memory budget, GiB.
    pub memory_gib: f64,
    /// `[snv] device`.
    pub device: String,
    /// lrpipe-snv's model, if the run calls SNVs.
    pub model: Option<String>,
    /// The analyses the run makes, as [`analyses_line`] names them.
    pub analyses: String,
    /// What follows the stages shown, when they are not all of the run: the stages planned once
    /// the reads are aligned.
    pub then: Option<String>,
}

/// The analyses of a run, named for the display's header.
#[must_use]
pub fn analyses_line(analyses: impl IntoIterator<Item = crate::plan::Analysis>) -> String {
    use crate::plan::Analysis;
    analyses
        .into_iter()
        .map(|a| match a {
            Analysis::Snv => "SNVs and indels",
            Analysis::Phase => "phasing",
            Analysis::Sv => "SVs",
            Analysis::Cnv => "CNVs and LoH",
            Analysis::Tr => "tandem repeats",
            Analysis::Depth => "depth",
            Analysis::Methyl => "methylation",
        })
        .collect::<Vec<_>>()
        .join(" · ")
}

/// The display, while it runs.
#[derive(Debug)]
pub struct Display {
    tx: mpsc::Sender<Message>,
    thread: Option<JoinHandle<()>>,
}

#[derive(Debug)]
enum Message {
    /// An event of the run, and the line a run without the display prints for it.
    Event(Event, String),
    /// The run's stages are done: the display ends.
    Finish,
}

/// Whether the display can run: stdin and stderr are a terminal, which is not `dumb`.
#[must_use]
pub fn available() -> bool {
    use std::io::IsTerminal as _;
    io::stderr().is_terminal()
        && io::stdin().is_terminal()
        && std::env::var("TERM").map_or(true, |t| !t.is_empty() && t != "dumb")
}

impl Display {
    /// Takes the terminal and shows `stages`, of `layout`'s run, each with what it took when it
    /// last succeeded among `earlier`; those `ran` records as succeeded (this run aligned the
    /// reads before planning the rest) as done. Only those are listed: an analysis the run does
    /// not make has no row. Ctrl-C sets `stop`.
    ///
    /// The terminal is taken only once a stage starts or waits for memory: a run whose stages
    /// were all done prints its lines as without the display. If it cannot be taken then, the run
    /// prints its lines, saying why.
    ///
    /// # Errors
    ///
    /// If the display's thread cannot be started.
    pub fn start(
        header: Header,
        stages: &[Stage],
        earlier: &[StageRecord],
        ran: &[StageRecord],
        layout: &Layout,
        stop: Arc<AtomicBool>,
    ) -> io::Result<Self> {
        let rows = rows(stages, earlier, ran, layout);
        let (tx, rx) = mpsc::channel();
        let screen =
            Screen { header, rows, selected: 0, started: Instant::now(), stop, lines: Vec::new(), log: Vec::new() };
        let thread = thread::Builder::new().name("display".into()).spawn(move || screen.run(&rx))?;
        Ok(Self { tx, thread: Some(thread) })
    }
}

/// The display's rows: the stages of `stages`, and those `ran` records as succeeded, in the
/// run's order. An analysis the run does not make has none.
fn rows(stages: &[Stage], earlier: &[StageRecord], ran: &[StageRecord], layout: &Layout) -> Vec<Row> {
    StageName::ALL
        .iter()
        .filter_map(|&name| {
            let stage = stages.iter().find(|s| s.name == name);
            let done = ran.iter().find(|r| r.name == name && r.status == StageStatus::Succeeded);
            if stage.is_none() && done.is_none() {
                return None;
            }
            let expected = earlier
                .iter()
                .find(|r| r.name == name && matches!(r.status, StageStatus::Succeeded | StageStatus::Skipped))
                .and_then(|r| r.wall_seconds)
                .filter(|&s| s > 0.0);
            Some(Row {
                name,
                module: stage.map_or_else(|| module(name).to_owned(), module_of),
                after: stage.map(|s| s.after.clone()).unwrap_or_default(),
                state: match done {
                    Some(record) => State::Done(record.wall_seconds.unwrap_or_default()),
                    None => State::Pending,
                },
                expected,
                log: layout.log(name),
                progress_file: layout.progress(name),
                reported: None,
                rss_gib: 0.0,
                peak_gib: 0.0,
            })
        })
        .collect()
}

impl Display {
    /// Shows `event`; `line` is what a run without the display prints for it.
    pub fn report(&self, event: Event, line: String) {
        // The display's thread ends only when told to, so it takes it.
        let _ = self.tx.send(Message::Event(event, line));
    }

    /// Ends the display: the terminal is given back, and the run's lines printed.
    pub fn finish(mut self) {
        self.end();
    }

    fn end(&mut self) {
        if let Some(thread) = self.thread.take() {
            let _ = self.tx.send(Message::Finish);
            let _ = thread.join();
        }
    }
}

impl Drop for Display {
    fn drop(&mut self) {
        self.end();
    }
}

/// Whether the display has the terminal.
enum Mode {
    /// Not yet: no stage has started.
    Before,
    /// It has.
    Shown(Terminal<CrosstermBackend<io::Stderr>>),
    /// Not any more: `q`, or a terminal that could not be drawn on.
    Lines,
}

/// The terminal, for the display; as it was if it cannot be.
fn take_terminal() -> io::Result<Terminal<CrosstermBackend<io::Stderr>>> {
    enter()?;
    Terminal::new(CrosstermBackend::new(io::stderr())).inspect_err(|_| leave())
}

/// Raw mode and the alternate screen, with a panic hook that gives the terminal back: with
/// `panic = "abort"` nothing else would.
fn enter() -> io::Result<()> {
    static HOOK: Once = Once::new();
    HOOK.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            leave();
            previous(info);
        }));
    });
    enable_raw_mode()?;
    if let Err(e) = execute!(io::stderr(), EnterAlternateScreen, cursor::Hide) {
        leave();
        return Err(e);
    }
    Ok(())
}

/// Gives the terminal back as it was; harmless when it already is.
fn leave() {
    let _ = execute!(io::stderr(), LeaveAlternateScreen, cursor::Show);
    let _ = disable_raw_mode();
}

/// How long a run goes before the display says how long is left: before, it has too little to
/// go on.
const ESTIMATE_AFTER: f64 = 10.0;

/// How often the display is drawn, and the stages' files read.
const TICK: Duration = Duration::from_millis(250);

/// What a stage is doing, as far as the display knows.
#[derive(Clone, Debug, PartialEq)]
enum State {
    /// Waits for the stages it follows.
    Pending,
    /// Ready, but waits for memory: it expects this many GiB.
    Memory(f64),
    /// Runs as this process.
    Running(Instant, u32),
    /// Succeeded after this many seconds.
    Done(f64),
    /// Was already done.
    Skipped,
    /// Failed, for this reason.
    Failed(String),
    /// Was killed: the run was stopped, or another stage failed.
    Cancelled,
}

#[derive(Debug)]
struct Row {
    name: StageName,
    /// The tool it runs, and for the read pass what it computes.
    module: String,
    after: Vec<StageName>,
    state: State,
    /// Seconds it took when it last succeeded here.
    expected: Option<f64>,
    log: PathBuf,
    progress_file: PathBuf,
    /// What the stage last wrote to its progress file.
    reported: Option<(f64, String)>,
    rss_gib: f64,
    peak_gib: f64,
}

impl Row {
    fn finished(&self) -> bool {
        matches!(self.state, State::Done(_) | State::Skipped)
    }

    /// How far it has got, if that is known, and whether by an estimate from its last run.
    fn fraction(&self) -> Option<(f64, bool)> {
        match &self.state {
            State::Done(_) | State::Skipped => Some((1.0, false)),
            State::Running(started, _) => match (&self.reported, self.expected) {
                (Some((f, _)), _) => Some((*f, false)),
                (None, Some(expected)) => Some(((started.elapsed().as_secs_f64() / expected).min(0.99), true)),
                (None, None) => None,
            },
            // Where it got to.
            State::Failed(_) | State::Cancelled => Some((self.reported.as_ref().map_or(0.0, |r| r.0), false)),
            _ => Some((0.0, false)),
        }
    }
}

struct Screen {
    header: Header,
    rows: Vec<Row>,
    selected: usize,
    started: Instant,
    stop: Arc<AtomicBool>,
    /// The run's lines, printed when the display ends.
    lines: Vec<String>,
    /// The end of the selected stage's log.
    log: Vec<String>,
}

impl Screen {
    fn run(mut self, rx: &mpsc::Receiver<Message>) {
        let mut mode = Mode::Before;
        loop {
            let Mode::Shown(terminal) = &mut mode else {
                // Lines, as they come, until a stage starts.
                let message = match rx.recv_timeout(TICK) {
                    Ok(message) => message,
                    Err(mpsc::RecvTimeoutError::Timeout) => continue,
                    Err(mpsc::RecvTimeoutError::Disconnected) => Message::Finish,
                };
                if !self.receive(message, &mut mode) {
                    break;
                }
                continue;
            };
            self.refresh();
            if terminal.draw(|frame| self.draw(frame)).is_err() {
                // A terminal that cannot be drawn on: the run's lines instead.
                self.give_back();
                mode = Mode::Lines;
                continue;
            }
            if term_event::poll(TICK).unwrap_or(false)
                && let Ok(TermEvent::Key(key)) = term_event::read()
                && key.kind != KeyEventKind::Release
            {
                match key.code {
                    KeyCode::Char('c') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                        self.stop.store(true, Ordering::Relaxed);
                    }
                    KeyCode::Char('q') => {
                        self.give_back();
                        mode = Mode::Lines;
                        continue;
                    }
                    KeyCode::Up | KeyCode::Char('k') => self.select(-1),
                    KeyCode::Down | KeyCode::Char('j') => self.select(1),
                    _ => {}
                }
            }
            let mut running = true;
            while running {
                match rx.try_recv() {
                    Ok(message) => running = self.receive(message, &mut mode),
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => running = self.receive(Message::Finish, &mut mode),
                }
            }
            if !running {
                break;
            }
        }
        if matches!(mode, Mode::Shown(_)) {
            self.give_back();
        }
    }

    /// Takes in `message`: an event is shown, its line kept while the display is (or printed),
    /// and the first stage to start or wait takes the terminal. Whether the display goes on.
    fn receive(&mut self, message: Message, mode: &mut Mode) -> bool {
        let Message::Event(event, line) = message else { return false };
        let takes = matches!(event, Event::Started(..) | Event::Waiting(..));
        self.apply(event);
        match mode {
            Mode::Shown(_) => self.lines.push(line),
            Mode::Before | Mode::Lines => say(&line),
        }
        if takes && matches!(mode, Mode::Before) {
            *mode = match take_terminal() {
                Ok(terminal) => Mode::Shown(terminal),
                Err(e) => {
                    say(&format!("lrpipe: warning: the display could not start ({e}): printing lines instead"));
                    Mode::Lines
                }
            };
        }
        true
    }

    /// Leaves the display and prints the run's lines so far.
    fn give_back(&mut self) {
        leave();
        for line in self.lines.drain(..) {
            say(&line);
        }
    }

    fn select(&mut self, step: isize) {
        if self.rows.is_empty() {
            return;
        }
        self.selected = self.selected.saturating_add_signed(step).min(self.rows.len() - 1);
    }

    fn row(&mut self, name: StageName) -> &mut Row {
        let i = self.rows.iter().position(|r| r.name == name).expect("every stage has a row");
        &mut self.rows[i]
    }

    fn apply(&mut self, event: Event) {
        match event {
            Event::Skipped(s) => self.row(s).state = State::Skipped,
            Event::Waiting(s, needs, _) => self.row(s).state = State::Memory(needs),
            Event::Started(s, pid) => {
                let row = self.row(s);
                row.state = State::Running(Instant::now(), pid);
                row.reported = None;
                // The stage that just started is the one to watch, unless one is failing.
                if !self.rows.iter().any(|r| matches!(r.state, State::Failed(_))) {
                    self.selected = self.rows.iter().position(|r| r.name == s).unwrap_or(self.selected);
                }
            }
            Event::Succeeded(s, secs) => {
                let row = self.row(s);
                row.state = State::Done(secs);
                row.rss_gib = 0.0;
            }
            Event::Failed(s, why) => {
                let row = self.row(s);
                row.state = State::Failed(why);
                row.rss_gib = 0.0;
                self.selected = self.rows.iter().position(|r| r.name == s).unwrap_or(self.selected);
            }
            Event::Cancelled(s) | Event::Abandoned(s, _) => {
                let row = self.row(s);
                row.state = State::Cancelled;
                row.rss_gib = 0.0;
            }
            Event::Interrupted => self.stop.store(true, Ordering::Relaxed),
        }
    }

    /// Reads the running stages' progress and memory, and the selected stage's log.
    fn refresh(&mut self) {
        for row in &mut self.rows {
            if let State::Running(_, pid) = row.state {
                if let Some(reported) = lrpipe_reads::progress::read(&row.progress_file) {
                    row.reported = Some(reported);
                }
                if let Some((rss, peak)) = memory(pid) {
                    row.rss_gib = rss;
                    row.peak_gib = row.peak_gib.max(peak);
                }
            }
        }
        let row = &self.rows[self.selected];
        self.log = match row.state {
            State::Pending | State::Memory(_) => Vec::new(),
            _ => tail(&row.log, LOG_BYTES),
        };
    }

    fn draw(&self, frame: &mut ratatui::Frame<'_>) {
        let area = frame.area();
        let [top, table, log, keys] = Split::vertical([
            Constraint::Length(6),
            Constraint::Length(self.rows.len() as u16 + 1 + u16::from(self.header.then.is_some())),
            Constraint::Min(3),
            Constraint::Length(1),
        ])
        .areas(area);
        frame.render_widget(Paragraph::new(self.top(top.width)), top);
        frame.render_widget(Paragraph::new(self.table(table.width)), table);
        self.draw_log(frame, log);
        let help = "↑↓ select a stage   q leave the display (the run goes on)   Ctrl-C stop the run";
        frame.render_widget(Paragraph::new(Line::styled(help, dim())), keys);
    }

    fn top(&self, width: u16) -> Vec<Line<'static>> {
        let h = &self.header;
        let planned: Vec<&Row> = self.rows.iter().collect();
        // Each stage weighs what it took last time, when every stage's is known; else its share
        // of a whole genome's run.
        let measured = planned.iter().all(|r| r.expected.is_some());
        let (mut done, mut all) = (0.0, 0.0);
        for row in &planned {
            let weight = if measured { row.expected.unwrap_or_default() } else { share(row.name) };
            all += weight;
            done += weight * row.fraction().map_or(0.0, |(f, _)| f);
        }
        let fraction = if all > 0.0 { done / all } else { 0.0 };
        let failed = planned.iter().any(|r| matches!(r.state, State::Failed(_)));
        let stopping = self.stop.load(Ordering::Relaxed);
        let elapsed = self.started.elapsed().as_secs_f64();
        let (state, color) = if failed {
            ("failed: stopping".to_owned(), Color::Red)
        } else if stopping {
            ("stopping".to_owned(), Color::Yellow)
        } else if planned.iter().all(|r| r.finished()) {
            ("done".to_owned(), Color::Green)
        } else if fraction > 0.02 && elapsed >= ESTIMATE_AFTER {
            (format!("about {} left", duration(elapsed * (1.0 - fraction) / fraction)), Color::Cyan)
        } else {
            (String::new(), Color::Cyan)
        };
        // Room for the percentage, the time and what the run is doing.
        let bar_width = usize::from(width).saturating_sub(58).clamp(10, 50);
        let rss: f64 = planned.iter().map(|r| r.rss_gib).sum();
        let memory_bar = bar(rss / h.memory_gib, 20);
        vec![
            Line::from(vec![
                Span::styled(format!("lrpipe {}", env!("CARGO_PKG_VERSION")), bold()),
                Span::raw(format!("  run  {} · {} · ", h.sample, h.platform)),
                Span::styled(h.output.display().to_string(), dim()),
            ]),
            Line::from(vec![
                Span::styled("model ", dim()),
                Span::raw(h.model.clone().unwrap_or_else(|| "none (no SNV calling)".into())),
                Span::styled("  device ", dim()),
                Span::raw(h.device.clone()),
            ]),
            Line::from(vec![Span::styled("analyses ", dim()), Span::raw(h.analyses.clone())]),
            Line::raw(""),
            Line::from(vec![
                Span::styled("Overall  ", bold()),
                Span::styled(bar(fraction, bar_width), Style::new().fg(color)),
                Span::raw(format!(" {:>3.0}%  {} elapsed  ", fraction * 100.0, duration(elapsed))),
                Span::styled(state, Style::new().fg(color)),
            ]),
            Line::from(vec![
                Span::styled("memory   ", dim()),
                Span::styled(memory_bar, Style::new().fg(if rss > h.memory_gib { Color::Red } else { Color::Reset })),
                Span::raw(format!(" {rss:.1} of {:.0} GiB", h.memory_gib)),
                Span::styled("   threads ", dim()),
                Span::raw(h.threads.clone()),
            ]),
        ]
    }

    fn table(&self, width: u16) -> Vec<Line<'static>> {
        // Stage, module, bar, %, time, peak memory, and what it is doing in what is left; on a
        // narrow terminal without the module, then the peak, so that what it is doing shows.
        let width = usize::from(width);
        let module_width = if width >= WIDE { MODULE } else { 0 };
        let peak_width = if width >= MEDIUM { PEAK } else { 0 };
        let fixed = 2 + STAGE + module_width + PERCENT + TIME + peak_width;
        let bar_width = width.saturating_sub(fixed + DOING).clamp(8, 32);
        let mut lines = vec![Line::styled(
            format!(
                "  {:<STAGE$}{:<module_width$}{:<bar_width$}{:>PERCENT$}{:>TIME$}{:>peak_width$}",
                "STAGE",
                if module_width > 0 { "MODULE" } else { "" },
                "PROGRESS",
                "%",
                "TIME  ",
                if peak_width > 0 { "PEAK  " } else { "" },
            ),
            dim(),
        )];
        for (i, row) in self.rows.iter().enumerate() {
            let (icon, color) = match row.state {
                State::Pending => ("·", Color::DarkGray),
                State::Memory(_) => ("…", Color::Yellow),
                State::Running(..) => ("▶", Color::Cyan),
                State::Done(_) => ("✓", Color::Green),
                State::Skipped => ("✓", Color::DarkGray),
                State::Failed(_) => ("✗", Color::Red),
                State::Cancelled => ("■", Color::Yellow),
            };
            let mut spans = vec![
                Span::styled(format!("{icon} "), Style::new().fg(color)),
                Span::styled(
                    format!("{:<STAGE$}", row.name.as_str()),
                    if i == self.selected { bold() } else { plain() },
                ),
                Span::styled(format!("{:<module_width$}", fit(&row.module, module_width)), dim()),
            ];
            let note_width = bar_width + PERCENT;
            let note = |text: String| Span::styled(format!("{:<note_width$}", fit(&text, note_width)), dim());
            match &row.state {
                State::Pending => {
                    let waits: Vec<&str> = row
                        .after
                        .iter()
                        .filter(|a| self.rows.iter().any(|r| r.name == **a && !r.finished()))
                        .map(|a| a.as_str())
                        .collect();
                    spans.push(note(format!("waits on {}", waits.join(", "))));
                }
                State::Memory(needs) => spans.push(note(format!("waits for memory ({needs:.1} GiB)"))),
                State::Skipped => spans.push(note("already done".into())),
                _ => match row.fraction() {
                    Some((f, estimated)) => {
                        let style = Style::new().fg(color);
                        spans.push(Span::styled(bar(f, bar_width), style));
                        let pct = format!("{:>3.0}%", (f * 100.0).floor());
                        let pct = if estimated { pct.replacen(' ', "~", 1) } else { pct };
                        spans.push(Span::raw(format!("{pct:>PERCENT$}")));
                    }
                    None => {
                        let elapsed = match row.state {
                            State::Running(started, _) => started.elapsed().as_secs_f64(),
                            _ => 0.0,
                        };
                        spans.push(Span::styled(bounce(elapsed, bar_width), Style::new().fg(color)));
                        spans.push(Span::raw(" ".repeat(PERCENT)));
                    }
                },
            }
            let time = match row.state {
                State::Running(started, _) => duration(started.elapsed().as_secs_f64()),
                State::Done(secs) => duration(secs),
                _ => String::new(),
            };
            let peak =
                if row.peak_gib > 0.0 && peak_width > 0 { format!("{:.1} GiB", row.peak_gib) } else { String::new() };
            spans.push(Span::raw(format!("{:>TIME$}{:>peak_width$}", format!("{time}  "), format!("{peak}  "))));
            let doing = match &row.state {
                State::Running(..) => match (&row.reported, row.expected) {
                    (Some((_, what)), _) => what.clone(),
                    (None, Some(e)) => format!("took {} last time", duration(e)),
                    (None, None) => "no estimate: first run here".into(),
                },
                State::Failed(why) => why.clone(),
                State::Cancelled => "stopped".into(),
                _ => String::new(),
            };
            spans.push(Span::styled(doing, dim()));
            lines.push(Line::from(spans));
        }
        if let Some(then) = &self.header.then {
            lines.push(Line::styled(format!("  {then}"), dim()));
        }
        lines
    }

    fn draw_log(&self, frame: &mut ratatui::Frame<'_>, area: Rect) {
        let row = &self.rows[self.selected];
        let title = format!(" log: {} ({}) ", row.name, row.log.display());
        let block = Block::new().borders(Borders::TOP).title(title).border_style(dim());
        let height = usize::from(block.inner(area).height);
        let mut lines: Vec<Line<'static>> = Vec::new();
        if let State::Failed(why) = &row.state {
            lines.push(Line::styled(format!("failed: {why}"), Style::new().fg(Color::Red)));
        }
        if self.log.is_empty() {
            let why = match row.state {
                State::Pending | State::Memory(_) => "not started",
                State::Skipped => "already done: this run left its log as it was",
                _ => "nothing written yet",
            };
            lines.push(Line::styled(why, dim()));
        }
        let room = height.saturating_sub(lines.len());
        let skip = self.log.len().saturating_sub(room);
        lines.extend(self.log[skip..].iter().map(|l| Line::raw(l.clone())));
        frame.render_widget(Paragraph::new(lines).block(block), area);
    }
}

/// The stage column's width.
const STAGE: usize = 13;
/// The module column's width, shown from [`WIDE`] columns.
const MODULE: usize = 30;
/// The percentage's width.
const PERCENT: usize = 5;
/// The time's width.
const TIME: usize = 9;
/// The peak memory's width, shown from [`MEDIUM`] columns.
const PEAK: usize = 10;
/// What a stage is doing gets at least this much, where the terminal has it.
const DOING: usize = 32;
/// A terminal this wide shows every column.
const WIDE: usize = 120;
/// A terminal this wide shows the peak memory.
const MEDIUM: usize = 96;

/// How much of a log's end is read for its pane.
const LOG_BYTES: u64 = 32 << 10;

/// The tool `stage` runs, and for the read pass what it computes, from its command: only the
/// analyses the run makes.
fn module_of(stage: &Stage) -> String {
    if stage.name != StageName::Reads {
        return module(stage.name).to_owned();
    }
    let has = |option: &str| stage.args.iter().any(|a| a == option);
    let made: Vec<&str> = [("--vcf", "SVs"), ("--depth", "depth"), ("--repeats", "TR"), ("--bedmethyl", "5mC")]
        .into_iter()
        .filter(|&(option, _)| has(option))
        .map(|(_, what)| what)
        .collect();
    if made.is_empty() { "lrpipe-sv".to_owned() } else { format!("lrpipe-sv: {}", made.join(", ")) }
}

/// The tool a stage runs.
fn module(stage: StageName) -> &'static str {
    match stage {
        StageName::Align => "bonefish",
        StageName::MergeReads => "lrpipe merge-reads",
        StageName::Snv => "lrpipe-snv",
        StageName::Reads => "lrpipe-sv + depth",
        StageName::Haplotag => "lrpipe-phase",
        StageName::Repeats => "lrpipe-tr",
        StageName::SvPhase => "lrpipe-sv phase",
        StageName::Cnv => "lrpipe-cnv",
        StageName::Methylation => "lrpipe-methyl",
        StageName::Merge => "lrpipe merge",
        StageName::Qc => "lrpipe qc",
    }
}

/// A stage's share of a whole genome's run, roughly (HG002 ONT, 32 cores and a GPU: lrpipe-snv
/// 22 min, the read pass beside it 14 min, the stages after it a few minutes in all), for the
/// overall progress of a run without earlier times.
fn share(stage: StageName) -> f64 {
    match stage {
        StageName::Align => 60.0,
        StageName::MergeReads => 5.0,
        StageName::Snv => 22.0,
        StageName::Reads => 14.0,
        StageName::Haplotag => 1.5,
        StageName::Repeats => 1.0,
        StageName::SvPhase => 0.3,
        StageName::Cnv => 1.0,
        StageName::Methylation => 3.0,
        StageName::Merge => 0.5,
        StageName::Qc => 0.3,
    }
}

/// A bar `width` cells wide, `fraction` of it filled.
fn bar(fraction: f64, width: usize) -> String {
    let filled = ((fraction.clamp(0.0, 1.0) * width as f64).round() as usize).min(width);
    format!("{}{}", "█".repeat(filled), "░".repeat(width - filled))
}

/// A bar `width` cells wide with a block moving back and forth: a stage whose progress is not
/// known.
fn bounce(seconds: f64, width: usize) -> String {
    const BLOCK: usize = 4;
    let span = width.saturating_sub(BLOCK).max(1);
    let step = (seconds * 8.0) as usize % (2 * span);
    let at = if step < span { step } else { 2 * span - step };
    let at = at.min(width.saturating_sub(BLOCK));
    format!("{}{}{}", "░".repeat(at), "█".repeat(BLOCK.min(width)), "░".repeat(width.saturating_sub(at + BLOCK)))
}

/// `seconds` as `45s`, `12m03s` or `1h02m`.
fn duration(seconds: f64) -> String {
    let s = seconds.max(0.0).round() as u64;
    match (s / 3600, s / 60 % 60, s % 60) {
        (0, 0, s) => format!("{s}s"),
        (0, m, s) => format!("{m}m{s:02}s"),
        (h, m, _) => format!("{h}h{m:02}m"),
    }
}

/// `text` cut to `width` characters, with `…` where it was cut.
fn fit(text: &str, width: usize) -> String {
    if width == 0 {
        return String::new();
    }
    if text.chars().count() <= width {
        return text.to_owned();
    }
    let mut cut: String = text.chars().take(width.saturating_sub(1)).collect();
    cut.push('…');
    cut
}

/// A process's resident memory now and at its peak, GiB, from `/proc`.
fn memory(pid: u32) -> Option<(f64, f64)> {
    let status = std::fs::read_to_string(format!("/proc/{pid}/status")).ok()?;
    let kib = |key: &str| -> Option<f64> {
        let line = status.lines().find(|l| l.starts_with(key))?;
        line[key.len()..].trim().trim_end_matches("kB").trim().parse::<f64>().ok()
    };
    let rss = kib("VmRSS:")? / 1_048_576.0;
    Some((rss, kib("VmHWM:").map_or(rss, |k| k / 1_048_576.0)))
}

/// The last lines of `path`, of its last `bytes`, but blank ones, with a logger's
/// `[time LEVEL target]` taken off.
fn tail(path: &Path, bytes: u64) -> Vec<String> {
    let Ok(mut file) = File::open(path) else { return Vec::new() };
    let len = file.metadata().map_or(0, |m| m.len());
    // From the byte before, to tell whether the first line read starts there or in its middle.
    let from = len.saturating_sub(bytes).saturating_sub(1);
    if file.seek(SeekFrom::Start(from)).is_err() {
        return Vec::new();
    }
    let mut buf = Vec::new();
    if file.read_to_end(&mut buf).is_err() {
        return Vec::new();
    }
    // Read from the middle: whatever comes before the first line break is part of a line.
    let start = if from > 0 { buf.iter().position(|&b| b == b'\n').map_or(buf.len(), |i| i + 1) } else { 0 };
    let text = String::from_utf8_lossy(&buf[start..]);
    let lines: Vec<&str> = text.lines().collect();
    lines
        .into_iter()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            // Tabs, carriage returns and escapes would move the cursor over the display.
            let l: String = l.chars().map(|c| if c.is_control() { ' ' } else { c }).collect();
            let l = l.trim_end();
            match l.strip_prefix('[').and_then(|rest| rest.split_once("] ")) {
                Some((prefix, message)) if prefix.len() < 80 && prefix.contains(' ') => message.to_owned(),
                _ => l.to_owned(),
            }
        })
        .collect()
}

/// `line` on stderr, as `lrpipe run` prints its lines: a failed write is dropped.
fn say(line: &str) {
    let _ = writeln!(io::stderr(), "{line}");
}

fn bold() -> Style {
    Style::new().add_modifier(Modifier::BOLD)
}

fn dim() -> Style {
    Style::new().fg(Color::DarkGray)
}

fn plain() -> Style {
    Style::new()
}

#[cfg(test)]
mod tests {
    use std::ffi::OsString;

    use super::*;

    /// A display of a run without the align and merge-reads stages, every other stage after snv.
    fn screen() -> Screen {
        let rows = StageName::ALL
            .iter()
            .filter(|&&name| !matches!(name, StageName::Align | StageName::MergeReads))
            .map(|&name| Row {
                name,
                module: module(name).to_owned(),
                after: if name == StageName::Snv { Vec::new() } else { vec![StageName::Snv] },
                state: State::Pending,
                expected: None,
                log: PathBuf::new(),
                progress_file: PathBuf::new(),
                reported: None,
                rss_gib: 0.0,
                peak_gib: 0.0,
            })
            .collect();
        Screen {
            header: Header {
                sample: "S".into(),
                platform: "ont".into(),
                output: "/out".into(),
                threads: "16".into(),
                memory_gib: 32.0,
                device: "cpu".into(),
                model: None,
                analyses: "SNVs and indels".into(),
                then: None,
            },
            rows,
            selected: 0,
            started: Instant::now(),
            stop: Arc::new(AtomicBool::new(false)),
            lines: Vec::new(),
            log: Vec::new(),
        }
    }

    fn state(screen: &mut Screen, name: StageName) -> State {
        screen.row(name).state.clone()
    }

    #[test]
    fn events_move_the_rows_and_the_selection_follows_them() {
        let mut s = screen();
        s.apply(Event::Waiting(StageName::Reads, 7.1, 14.5));
        assert_eq!(state(&mut s, StageName::Reads), State::Memory(7.1));
        s.apply(Event::Started(StageName::Snv, 1));
        assert!(matches!(state(&mut s, StageName::Snv), State::Running(_, 1)));
        assert_eq!(s.rows[s.selected].name, StageName::Snv, "the stage that started is shown");
        s.apply(Event::Started(StageName::Reads, 2));
        assert_eq!(s.rows[s.selected].name, StageName::Reads);
        s.apply(Event::Succeeded(StageName::Snv, 12.0));
        assert_eq!(state(&mut s, StageName::Snv), State::Done(12.0));
        // The QC summary's failure does not fail the run: the others go on, and it stays shown.
        s.apply(Event::Started(StageName::Qc, 3));
        s.apply(Event::Failed(StageName::Qc, "exit code 1".into()));
        assert_eq!(s.rows[s.selected].name, StageName::Qc);
        s.apply(Event::Started(StageName::Merge, 4));
        assert_eq!(s.rows[s.selected].name, StageName::Qc, "a failure keeps the selection");
        s.apply(Event::Interrupted);
        assert!(s.stop.load(Ordering::Relaxed));
        s.apply(Event::Cancelled(StageName::Reads));
        s.apply(Event::Abandoned(StageName::Merge, 4));
        assert_eq!(state(&mut s, StageName::Reads), State::Cancelled);
        assert_eq!(state(&mut s, StageName::Merge), State::Cancelled);
        assert_eq!(s.rows.iter().map(|r| r.rss_gib).sum::<f64>(), 0.0, "no memory is counted for ended stages");
    }

    #[test]
    fn the_selection_stops_at_the_ends() {
        let mut s = screen();
        for _ in 0..3 {
            s.select(-1);
        }
        assert_eq!(s.rows[s.selected].name, StageName::Snv, "the first stage in the run is the top");
        for _ in 0..30 {
            s.select(1);
        }
        assert_eq!(s.rows[s.selected].name, StageName::Qc, "the last is the bottom");
        s.select(-1);
        assert_eq!(s.rows[s.selected].name, StageName::Merge);
    }

    #[test]
    fn a_failed_stage_keeps_how_far_it_got() {
        let mut s = screen();
        s.apply(Event::Started(StageName::Snv, 1));
        s.row(StageName::Snv).reported = Some((0.4, "pileup".into()));
        s.apply(Event::Failed(StageName::Snv, "killed by signal 9".into()));
        assert_eq!(s.row(StageName::Snv).fraction(), Some((0.4, false)));
    }

    #[test]
    fn bars_and_times_read_as_meant() {
        assert_eq!(bar(0.5, 4), "██░░");
        assert_eq!(bar(1.7, 3), "███");
        assert_eq!(bar(-1.0, 2), "░░");
        for t in 0..200 {
            assert_eq!(bounce(f64::from(t) / 8.0, 12).chars().count(), 12);
        }
        assert_eq!(duration(45.4), "45s");
        assert_eq!(duration(723.0), "12m03s");
        assert_eq!(duration(3720.0), "1h02m");
        assert_eq!(fit("abcdef", 4), "abc…");
        assert_eq!(fit("abc", 4), "abc");
    }

    #[test]
    fn a_logs_tail_drops_blank_lines_and_logger_prefixes() {
        let path = std::env::temp_dir().join(format!("lrpipe-tui-tail-{}", std::process::id()));
        std::fs::write(
            &path,
            "first\n\n[2026-10-30T00:41:28Z INFO  pellucid::pipeline] [1/7] pileup calling\n[1/7] kept\n",
        )
        .expect("written");
        assert_eq!(tail(&path, 1 << 20), ["first", "[1/7] pileup calling", "[1/7] kept"]);
        assert_eq!(tail(&path, 11), ["[1/7] kept"], "a line read from its start is kept");
        assert_eq!(tail(&path, 9), Vec::<String>::new(), "a partial first line is dropped");
        std::fs::write(&path, b"bell\x07 tab\tescape\x1b[2J\r\n\xff\xfe bytes\n").expect("written");
        assert_eq!(tail(&path, 1 << 20), ["bell  tab escape [2J", "\u{fffd}\u{fffd} bytes"]);
        assert_eq!(tail(Path::new("/nonexistent/log"), 100), Vec::<String>::new());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn the_display_draws_on_small_and_large_terminals() {
        use ratatui::backend::TestBackend;
        let rows = StageName::ALL
            .iter()
            .filter(|&&name| name != StageName::Align)
            .map(|&name| Row {
                name,
                module: module(name).to_owned(),
                after: vec![StageName::Snv],
                state: match name {
                    StageName::Snv => State::Running(Instant::now(), std::process::id()),
                    StageName::Reads => State::Memory(7.1),
                    StageName::Haplotag => State::Failed("exit code 1".into()),
                    StageName::Qc => State::Done(3.0),
                    _ => State::Pending,
                },
                expected: None,
                log: PathBuf::from("/nonexistent/snv.log"),
                progress_file: PathBuf::new(),
                reported: (name == StageName::Snv).then(|| (0.4, "pileup calling: 5 of 13 chunks".to_owned())),
                rss_gib: 1.0,
                peak_gib: 2.0,
            })
            .collect();
        let screen = Screen {
            header: Header {
                sample: "HG002".into(),
                platform: "ont".into(),
                output: "/out".into(),
                threads: "16".into(),
                memory_gib: 32.0,
                device: "cpu".into(),
                model: None,
                analyses: "SNVs and indels · SVs".into(),
                then: Some("then, once the reads are aligned: more".into()),
            },
            rows,
            selected: 1,
            started: Instant::now(),
            stop: Arc::new(AtomicBool::new(false)),
            lines: Vec::new(),
            log: vec!["[1/7] pileup calling".into()],
        };
        for (w, h) in [(20, 6), (80, 24), (140, 30), (200, 60)] {
            let mut terminal = Terminal::new(TestBackend::new(w, h)).expect("a test terminal");
            terminal.draw(|f| screen.draw(f)).expect("drawn");
            let buffer = terminal.backend().buffer();
            let text: String = buffer.content().iter().map(|c| c.symbol()).collect();
            if std::env::var_os("LRPIPE_SHOW_TUI").is_some() {
                for y in 0..h {
                    let line: String = (0..w).map(|x| buffer[(x, y)].symbol()).collect();
                    eprintln!("{}", line.trim_end());
                }
            }
            if w >= 80 && h >= 24 {
                assert!(text.contains("analyses SNVs and indels · SVs"), "{w}x{h}");
                assert!(text.contains("then, once the reads are aligned"), "{w}x{h}");
            }
            if w >= 80 {
                assert!(text.contains("pileup calling: 5 of 13"), "{w}x{h}");
                assert!(text.contains("waits for memory"), "{w}x{h}");
                assert_eq!(text.contains("lrpipe-snv"), usize::from(w) >= WIDE, "the module column, {w}x{h}");
            }
        }
    }

    #[test]
    fn only_the_runs_stages_get_rows_named_for_what_they_make() {
        // The analyses asked for were SNVs, SVs and depth: no CNV, repeats or methylation rows,
        // and the read pass named for its SVs and depth alone. The align stage, which this run
        // ran before planning the rest, is listed as done.
        let stage = |name, args: &[&str]| Stage {
            name,
            program: PathBuf::from("lrpipe"),
            args: args.iter().map(OsString::from).collect(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            after: Vec::new(),
            memory_gib: 1.0,
            priority: crate::plan::Priority::Critical,
            scratch: None,
        };
        let stages = [
            stage(StageName::Snv, &[]),
            stage(StageName::Reads, &["sv", "--vcf", "sv.vcf", "--depth", "depth/S"]),
            stage(StageName::Merge, &[]),
            stage(StageName::Qc, &[]),
        ];
        let aligned: StageRecord = serde_json::from_value(serde_json::json!({
            "name": "align", "status": "succeeded", "command": [], "log": "logs/align.log", "wall_seconds": 840.0
        }))
        .expect("a stage record");
        let rows = rows(&stages, &[], &[aligned], &Layout::new(Path::new("/out"), "S"));
        let names: Vec<StageName> = rows.iter().map(|r| r.name).collect();
        assert_eq!(names, [StageName::Align, StageName::Snv, StageName::Reads, StageName::Merge, StageName::Qc]);
        assert_eq!(rows[0].state, State::Done(840.0));
        assert_eq!(rows[2].module, "lrpipe-sv: SVs, depth");
        let all = stage(StageName::Reads, &["--vcf", "v", "--depth", "d", "--repeats", "c", "--bedmethyl", "b"]);
        assert_eq!(module_of(&all), "lrpipe-sv: SVs, depth, TR, 5mC");
        use crate::plan::Analysis;
        assert_eq!(analyses_line([Analysis::Snv, Analysis::Sv, Analysis::Depth]), "SNVs and indels · SVs · depth");
    }

    #[test]
    fn a_stage_without_a_report_is_estimated_from_its_last_run() {
        let row = |reported, expected| Row {
            name: StageName::Snv,
            module: String::new(),
            after: Vec::new(),
            state: State::Running(Instant::now() - Duration::from_secs(30), 1),
            expected,
            log: PathBuf::new(),
            progress_file: PathBuf::new(),
            reported,
            rss_gib: 0.0,
            peak_gib: 0.0,
        };
        assert_eq!(row(Some((0.4, "pileup".into())), Some(60.0)).fraction(), Some((0.4, false)));
        let (f, estimated) = row(None, Some(60.0)).fraction().expect("estimated");
        assert!(estimated && (0.45..0.6).contains(&f), "{f}");
        assert_eq!(row(None, Some(10.0)).fraction(), Some((0.99, true)), "never 100% before it ends");
        assert_eq!(row(None, None).fraction(), None);
    }
}
