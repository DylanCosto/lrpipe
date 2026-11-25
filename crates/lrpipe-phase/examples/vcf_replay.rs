//! The lrpipe-phase side of the VCF I/O oracle: replays DIR/ops.txt written by
//! tests/vcf/gen.py and prints the same lines as tests/vcf/replay.py, writing
//! the output VCF to DIR/rs.vcf.
use std::collections::HashMap;

use anyhow::Result;
use lrpipe_phase::hts::{
    BCF_HL_FMT, BCF_HL_INFO, FormatScalar, Header, IndexedVcf, Record, VcfFetchError, VcfFile, VcfWriter,
};
use lrpipe_phase::vcf::{PhaseTag, open_for_phased_output};

fn hexf(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn show_number(v: Option<f64>, typ: &str) -> String {
    match v {
        None => ".".into(),
        Some(v) if typ == "Integer" => format!("{}", v as i64),
        Some(v) => hexf(v),
    }
}

fn show_gt(gt: &Option<Vec<Option<u32>>>) -> String {
    match gt {
        None => "-".into(),
        Some(a) => {
            let items: Vec<String> = a.iter().map(|x| x.map_or(".".into(), |v| v.to_string())).collect();
            format!("[{}]", items.join(","))
        }
    }
}

fn dump_header(header: &Header) {
    println!("H samples={}", header.samples().join(","));
    println!("H contigs={}", header.contigs().join(","));
    for (name, hl) in [("INFO", BCF_HL_INFO), ("FORMAT", BCF_HL_FMT)] {
        let items: Vec<String> =
            header.defined(hl).into_iter().map(|(k, v)| format!("{k}:{}:{}", v.number, v.typ)).collect();
        println!("H {name}={}", items.join(","));
    }
    println!("H keys={}", header.record_keys().join(","));
}

fn dump_sample(rec: &Record, header: &Header, j: usize) -> Result<String> {
    let mut out =
        vec![format!("gt={}", show_gt(&rec.gt(header, j)?)), format!("ph={}", rec.sample_phased(header, j)? as u8)];
    let formats: HashMap<String, _> = header.defined(BCF_HL_FMT).into_iter().collect();
    for key in ["GL", "PL", "HP", "PS", "PQ"] {
        let Some(meta) = formats.get(key) else { continue };
        let (number, typ) = (meta.number.as_str(), meta.typ.as_str());
        let numeric = typ == "Integer" || typ == "Float";
        if numeric && number == "G" {
            match rec.format_numbers(header, j, key)? {
                None => out.push(format!("{key}=-")),
                Some(v) => {
                    let items: Vec<String> = v.into_iter().map(|x| show_number(x, typ)).collect();
                    out.push(format!("{key}={}", items.join(",")));
                }
            }
        } else if numeric && number == "1" {
            out.push(match rec.format_scalar(header, j, key)? {
                None => format!("{key}=-"),
                Some(None) => format!("{key}=None"),
                Some(Some(FormatScalar::Int(v))) => format!("{key}=i:{v}"),
                Some(Some(FormatScalar::Float(v))) => format!("{key}=f:{}", hexf(v as f64)),
            });
            if key == "PS" {
                let has = rec.has_format(header, "PS")?;
                out.push(format!("hasPS={}", has as u8));
                out.push(format!("presentPS={}", (has && rec.format_value_present(header, j, "PS")?) as u8));
            }
        } else if typ == "String" && number != "1" {
            match rec.format_strings(header, j, key)? {
                None => out.push(format!("{key}=-")),
                Some(v) => {
                    let items: Vec<String> = v.into_iter().map(|x| x.unwrap_or_else(|| "<None>".into())).collect();
                    out.push(format!("{key}={}", items.join(",")));
                }
            }
        }
    }
    Ok(out.join(" "))
}

fn dump_record(prefix: &str, i: usize, rec: &Record, header: &Header) -> Result<()> {
    let alts = rec.alts()?.map_or("-".into(), |a| a.join(","));
    println!(
        "{prefix} {i} {} {} {} {alts} info={} format={}",
        rec.chrom(header),
        rec.pos(),
        rec.reference()?.unwrap_or_else(|| "None".into()),
        rec.info_keys(header)?.join(","),
        rec.format_keys(header)?.join(",")
    );
    for j in 0..rec.n_sample() {
        println!("  S {j} {}", dump_sample(rec, header, j)?);
    }
    Ok(())
}

fn read_pass(dir: &str) -> Result<()> {
    let mut f = VcfFile::open(&format!("{dir}/in.vcf"))?;
    dump_header(&f.header);
    let mut rec = Record::new();
    let mut i = 0;
    while f.read(&mut rec)? {
        dump_record("R", i, &rec, &f.header)?;
        i += 1;
    }
    Ok(())
}

fn alleles(s: &str) -> Vec<u32> {
    s.split(',').map(|a| a.parse().unwrap()).collect()
}

fn write_pass(
    dir: &str,
    mode: PhaseTag,
    cmdline: Option<&str>,
    edits: &HashMap<usize, Vec<Vec<String>>>,
) -> Result<()> {
    let mut reader = open_for_phased_output(&format!("{dir}/in.vcf"), cmdline, mode)?;
    let mut writer = VcfWriter::create(&format!("{dir}/rs.vcf"), &reader.header)?;
    dump_header(&reader.header);
    let mut i = 0;
    loop {
        let mut rec = Record::new();
        if !reader.read(&mut rec)? {
            break;
        }
        let header = &reader.header;
        for e in edits.get(&i).into_iter().flatten() {
            let s: usize = e[1].parse()?;
            match e[0].as_str() {
                "unphase" => {
                    if !rec.has_format(header, "GT")? {
                        continue;
                    }
                    rec.set_phased(header, s, false)?;
                    if let Some(gt) = rec.gt(header, s)?
                        && gt.iter().all(|a| a.is_some())
                    {
                        let mut a: Vec<u32> = gt.iter().map(|a| a.unwrap()).collect();
                        a.sort_unstable();
                        rec.set_gt(header, s, &a)?;
                    }
                }
                "gt" => rec.set_gt(header, s, &alleles(&e[2]))?,
                "setps" => {
                    let (comp, phase) = e[2].split_once(' ').unwrap();
                    rec.format_set_int(header, s, "PS", &[Some(comp.parse()?)])?;
                    rec.set_gt(header, s, &alleles(phase))?;
                    rec.set_phased(header, s, true)?;
                }
                "ps" => rec.format_set_int(header, s, "PS", &[None])?,
                "hp" => rec.format_set_str(header, s, "HP", (e[2] != ".").then_some(e[2].as_str()))?,
                op => panic!("unknown op {op}"),
            }
        }
        if edits.contains_key(&i) {
            dump_record("W", i, &rec, header)?;
        }
        writer.write(&mut rec, header)?;
        i += 1;
    }
    writer.close()
}

fn main() -> Result<()> {
    let dir = std::env::args().nth(1).expect("usage: vcf_replay DIR");
    let text = std::fs::read_to_string(format!("{dir}/ops.txt"))?;
    let ops: Vec<Vec<String>> = text.lines().map(|l| l.splitn(5, ' ').map(String::from).collect()).collect();
    let mode = match ops.iter().find(|o| o[0] == "mode").map(|o| o[1].as_str()) {
        Some("HP") => PhaseTag::HP,
        _ => PhaseTag::PS,
    };
    let cmdline = ops.iter().find(|o| o[0] == "cmdline").map(|o| o[1..].join(" ")).filter(|c| c != "-");
    let mut edits: HashMap<usize, Vec<Vec<String>>> = HashMap::new();
    for o in &ops {
        if o[0] == "rec" {
            edits.entry(o[1].parse()?).or_default().push(o[2..].to_vec());
        }
    }

    if read_pass(&dir).is_err() {
        println!("ERROR read");
    }
    if write_pass(&dir, mode, cmdline.as_deref(), &edits).is_err() {
        println!("ERROR write");
    }

    let mut f = IndexedVcf::open(&format!("{dir}/in.vcf.gz"))?;
    for o in ops.iter().filter(|o| o[0] == "fetch") {
        let (contig, start) = (&o[1], o[2].parse::<i64>()?);
        let stop = if o[3] == "-" { None } else { Some(o[3].parse::<i64>()?) };
        println!("fetch {contig} {start} {}", o[3]);
        match f.fetch(contig, start, stop) {
            Ok(recs) => {
                for rec in recs {
                    println!("  F {} {} {}", rec.chrom(f.header()), rec.pos(), rec.reference()?.unwrap_or_default());
                }
            }
            Err(VcfFetchError::InvalidContig) => println!("  ERR invalid contig"),
            Err(_) => println!("  ERR other"),
        }
    }
    Ok(())
}
