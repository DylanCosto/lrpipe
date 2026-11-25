//! Clair3's models, by their standard directory names, and the reads each was trained on: what
//! [`super::choose`] picks a model from and checks a given one against. The names are those of
//! the directories Clair3's Bioconda package and Docker image install (HKU's models and ONT's
//! bundled ones) and of ONT's Rerio models converted for Clair3 2
//! (<https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/>, as listed 2026-11-03; its
//! R10.4 E8.1 models are left out).
//!
//! A model for a newer basecaller is a new row; nothing else changes.

use std::fmt;

/// A Clair3 model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Model {
    /// Its directory name.
    pub name: &'static str,
    /// The reads it was trained on.
    pub reads: Reads,
    /// Whether lrpipe picks it from a header.
    pub pick: Pick,
    /// Where it comes from.
    pub source: Source,
}

/// Where a Clair3 model comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// Clair3's own models, which its Bioconda package and Docker image install in `models/`.
    Bundled,
    /// ONT's Rerio models converted for Clair3 2, downloaded on their own.
    Rerio,
}

impl fmt::Display for Source {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Bundled => "one of Clair3's own models, which its Bioconda package and Docker image install",
            Self::Rerio => {
                "one of ONT's Rerio models converted for Clair3 2, not among Clair3's own: \
                 https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/"
            }
        })
    }
}

/// Whether lrpipe picks a model from the reads' header.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pick {
    /// The model for its reads.
    Yes,
    /// Only when given, for the reason held.
    No(&'static str),
}

/// The reads a model was trained on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reads {
    /// Oxford Nanopore.
    Ont {
        /// The flow cell's pore, as dorado names it (`r10.4.1`).
        pore: &'static str,
        /// The motor protein (`e8.2`).
        motor: &'static str,
        /// Bases per second through the pore.
        bps: u16,
        /// The signal's sampling rate in kHz: 4 for Guppy's and dorado's models before v4.2.0,
        /// 5 from then.
        khz: u8,
        /// The basecall model's accuracy.
        tier: Tier,
        /// What basecalled them.
        basecaller: Basecaller,
    },
    /// PacBio HiFi, from the instruments whose `@RG PM` is one of these.
    Hifi {
        /// The instruments.
        instruments: &'static [&'static str],
        /// Instruments it was not trained on but is Clair3's nearest model for.
        nearest_for: &'static [&'static str],
    },
    /// Illumina.
    Illumina,
}

/// The accuracy of an Oxford Nanopore basecall model, from the fastest to the most accurate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tier {
    /// `fast`.
    Fast,
    /// `hac`, high accuracy.
    Hac,
    /// `sup`, super accuracy.
    Sup,
}

impl Tier {
    /// The name dorado's and Clair3's model names spell it with.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fast => "fast",
            Self::Hac => "hac",
            Self::Sup => "sup",
        }
    }

    /// The tier named `name` (`fast`, `hac`, `sup`).
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        [Self::Fast, Self::Hac, Self::Sup].into_iter().find(|t| t.as_str() == name)
    }
}

/// A basecaller's version: major, minor and patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Version(pub [u16; 3]);

impl Version {
    /// `5.0.0` or `3.3` (the patch then 0), with or without a leading `v`.
    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        let mut parts = text.strip_prefix('v').unwrap_or(text).split('.');
        let mut v = [0; 3];
        for (i, slot) in v.iter_mut().enumerate() {
            match parts.next() {
                Some(p) => *slot = p.parse().ok()?,
                None if i >= 2 => {}
                None => return None,
            }
        }
        parts.next().is_none().then_some(Self(v))
    }

    /// The version as one number, for how far apart two are.
    #[must_use]
    pub fn ordinal(self) -> u64 {
        let [major, minor, patch] = self.0.map(u64::from);
        major * 1_000_000 + minor * 1000 + patch
    }
}

impl fmt::Display for Version {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [major, minor, patch] = self.0;
        write!(f, "v{major}.{minor}.{patch}")
    }
}

/// What basecalled Oxford Nanopore reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Basecaller {
    /// dorado, with the version of its basecall model (`@v5.0.0`).
    Dorado(Version),
    /// Guppy, at this version.
    Guppy(Version),
    /// Not known.
    Unknown,
}

impl fmt::Display for Basecaller {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dorado(v) => write!(f, "dorado {v}"),
            Self::Guppy(v) => write!(f, "Guppy {v}"),
            Self::Unknown => f.write_str("an unknown basecaller"),
        }
    }
}

const fn dorado(major: u16, minor: u16, patch: u16) -> Basecaller {
    Basecaller::Dorado(Version([major, minor, patch]))
}

const fn guppy(major: u16, minor: u16, patch: u16) -> Basecaller {
    Basecaller::Guppy(Version([major, minor, patch]))
}

/// R10.4.1 E8.2 at `bps`, sampled at `khz`.
const fn r1041(bps: u16, khz: u8, tier: Tier, basecaller: Basecaller) -> Reads {
    Reads::Ont { pore: "r10.4.1", motor: "e8.2", bps, khz, tier, basecaller }
}

const fn r941(tier: Tier, basecaller: Basecaller) -> Reads {
    Reads::Ont { pore: "r9.4.1", motor: "e8", bps: 450, khz: 4, tier, basecaller }
}

const fn pick(name: &'static str, reads: Reads, source: Source) -> Model {
    Model { name, reads, pick: Pick::Yes, source }
}

const fn given(name: &'static str, reads: Reads, why: &'static str) -> Model {
    Model { name, reads, pick: Pick::No(why), source: Source::Bundled }
}

/// The move-table models need `--enable_dwell_time` and the reads' `mv` tags.
const DWELL_TIME: &str = "a move-table model, which needs --enable_dwell_time, which lrpipe run does not pass";

const SEQUEL2: Reads = Reads::Hifi { instruments: &["SEQUELII", "SEQUELIIE"], nearest_for: &["SEQUEL"] };

/// Every standard Clair3 model. The first HiFi model lrpipe picks is the one for an instrument
/// the table does not know.
pub const MODELS: &[Model] = {
    use Source::{Bundled, Rerio};
    use Tier::{Fast, Hac, Sup};
    &[
        // ONT's models for dorado: bundled with Clair3, or converted from Rerio.
        pick("r1041_e82_400bps_hac_v600", r1041(400, 5, Hac, dorado(6, 0, 0)), Bundled),
        pick("r1041_e82_400bps_sup_v520", r1041(400, 5, Sup, dorado(5, 2, 0)), Bundled),
        pick("r1041_e82_400bps_hac_v520", r1041(400, 5, Hac, dorado(5, 2, 0)), Bundled),
        pick("r1041_e82_400bps_sup_v500", r1041(400, 5, Sup, dorado(5, 0, 0)), Bundled),
        pick("r1041_e82_400bps_hac_v500", r1041(400, 5, Hac, dorado(5, 0, 0)), Bundled),
        pick("r1041_e82_400bps_sup_v430", r1041(400, 5, Sup, dorado(4, 3, 0)), Rerio),
        pick("r1041_e82_400bps_hac_v430", r1041(400, 5, Hac, dorado(4, 3, 0)), Rerio),
        pick("r1041_e82_400bps_sup_v420", r1041(400, 5, Sup, dorado(4, 2, 0)), Rerio),
        pick("r1041_e82_400bps_hac_v420", r1041(400, 5, Hac, dorado(4, 2, 0)), Rerio),
        pick("r1041_e82_400bps_sup_v410", r1041(400, 4, Sup, dorado(4, 1, 0)), Bundled),
        pick("r1041_e82_400bps_hac_v410", r1041(400, 4, Hac, dorado(4, 1, 0)), Bundled),
        pick("r1041_e82_400bps_sup_v400", r1041(400, 4, Sup, dorado(4, 0, 0)), Rerio),
        pick("r1041_e82_400bps_hac_v400", r1041(400, 4, Hac, dorado(4, 0, 0)), Rerio),
        pick("r1041_e82_260bps_sup_v410", r1041(260, 4, Sup, dorado(4, 1, 0)), Rerio),
        pick("r1041_e82_260bps_hac_v410", r1041(260, 4, Hac, dorado(4, 1, 0)), Rerio),
        pick("r1041_e82_260bps_sup_v400", r1041(260, 4, Sup, dorado(4, 0, 0)), Rerio),
        pick("r1041_e82_260bps_hac_v400", r1041(260, 4, Hac, dorado(4, 0, 0)), Rerio),
        // ONT's for Guppy's R10.4.1 models, converted from Rerio.
        pick("r1041_e82_400bps_sup_g615", r1041(400, 4, Sup, guppy(6, 1, 5)), Rerio),
        pick("r1041_e82_400bps_hac_g632", r1041(400, 4, Hac, guppy(6, 3, 2)), Rerio),
        pick("r1041_e82_400bps_hac_g615", r1041(400, 4, Hac, guppy(6, 1, 5)), Rerio),
        pick("r1041_e82_400bps_fast_g632", r1041(400, 4, Fast, guppy(6, 3, 2)), Rerio),
        pick("r1041_e82_400bps_fast_g615", r1041(400, 4, Fast, guppy(6, 1, 5)), Rerio),
        pick("r1041_e82_260bps_sup_g632", r1041(260, 4, Sup, guppy(6, 3, 2)), Rerio),
        pick("r1041_e82_260bps_hac_g632", r1041(260, 4, Hac, guppy(6, 3, 2)), Rerio),
        pick("r1041_e82_260bps_fast_g632", r1041(260, 4, Fast, guppy(6, 3, 2)), Rerio),
        // HKU's.
        given("r1041_e82_400bps_hac_v600_with_mv", r1041(400, 5, Hac, dorado(6, 0, 0)), DWELL_TIME),
        given("r1041_e82_400bps_hac_v520_with_mv", r1041(400, 5, Hac, dorado(5, 2, 0)), DWELL_TIME),
        given("r1041_e82_400bps_sup_v520_with_mv", r1041(400, 5, Sup, dorado(5, 2, 0)), DWELL_TIME),
        given("r1041_e82_400bps_hac_with_mv", r1041(400, 5, Hac, Basecaller::Unknown), DWELL_TIME),
        given("r1041_e82_400bps_sup_with_mv", r1041(400, 5, Sup, Basecaller::Unknown), DWELL_TIME),
        given(
            "r1041_e82_400bps_sup_v430_bacteria_finetuned",
            r1041(400, 5, Sup, dorado(4, 3, 0)),
            "fine-tuned on bacterial genomes",
        ),
        // R9.4.1, basecalled with Guppy: Guppy 5's SUP model is also Clair3's for HAC reads.
        pick("r941_prom_sup_g5014", r941(Sup, guppy(5, 0, 14)), Bundled),
        pick("r941_prom_hac_g360+g422", r941(Hac, guppy(4, 2, 2)), Bundled),
        given("ont_guppy5", r941(Sup, guppy(5, 0, 14)), "the old name of r941_prom_sup_g5014"),
        given("ont", r941(Hac, guppy(4, 2, 2)), "the old name of r941_prom_hac_g360+g422"),
        // PacBio HiFi; Sequel II's model is the nearest for the first Sequel's HiFi reads.
        pick("hifi_revio", Reads::Hifi { instruments: &["REVIO"], nearest_for: &[] }, Bundled),
        pick("hifi_sequel2", SEQUEL2, Bundled),
        given("hifi", SEQUEL2, "the old name of hifi_sequel2"),
        given("ilmn", Reads::Illumina, "Illumina's"),
    ]
};

/// The model named `name`, if it is a standard Clair3 model.
#[must_use]
pub fn find(name: &str) -> Option<&'static Model> {
    MODELS.iter().find(|m| m.name == name)
}

impl fmt::Display for Reads {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ont { pore, motor, bps, khz, tier, basecaller } => write!(
                f,
                "Oxford Nanopore {} {} {bps} bps {khz} kHz reads basecalled with {basecaller} {}",
                pore.to_uppercase(),
                motor.to_uppercase(),
                tier.as_str().to_uppercase()
            ),
            Self::Hifi { instruments, .. } => write!(f, "PacBio HiFi reads from {}", instruments.join(" or ")),
            Self::Illumina => f.write_str("Illumina reads"),
        }
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;

    /// The directory name Clair3 gives a model for `reads`, before any suffix (`_with_mv`).
    fn clair3_name(reads: &Reads) -> Option<String> {
        let Reads::Ont { pore, motor, bps, tier, basecaller, .. } = reads else { return None };
        let digits = |s: &str| s.chars().filter(char::is_ascii_digit).collect::<String>();
        let tier = tier.as_str();
        // R9.4.1's are named after the PromethION, the others after the chemistry.
        let chemistry = if *pore == "r9.4.1" {
            format!("r{}_prom", digits(pore))
        } else {
            format!("r{}_{}_{bps}bps", digits(pore), motor.replace('.', ""))
        };
        Some(match basecaller {
            Basecaller::Dorado(Version([a, b, c])) => format!("{chemistry}_{tier}_v{a}{b}{c}"),
            Basecaller::Unknown => format!("{chemistry}_{tier}"),
            Basecaller::Guppy(Version([a, b, c])) => format!("{chemistry}_{tier}_g{a}{b}{c}"),
        })
    }

    #[test]
    fn every_name_is_a_clair3_directory_name_for_its_reads() {
        // A row whose name and reads disagree (a typo in either) would choose a directory that
        // is not there, or a model for other reads: fail here instead.
        let mut names = std::collections::HashSet::new();
        for model in MODELS {
            assert!(names.insert(model.name), "{} is listed twice", model.name);
            let well_formed =
                model.name.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || "_+".contains(c));
            assert!(well_formed, "{:?} is not a Clair3 directory name", model.name);
            match &model.reads {
                Reads::Ont { pore, motor, basecaller, khz, .. } => {
                    // dorado's models sample at 5 kHz from v4.2.0, Guppy's (here) at 4.
                    let five = matches!(basecaller, Basecaller::Dorado(v) if *v >= Version([4, 2, 0]))
                        || *basecaller == Basecaller::Unknown;
                    assert_eq!(*khz, if five { 5 } else { 4 }, "{}: sampled at {khz} kHz", model.name);
                    // `r10.4.1`, `e8.2`, `e8`: a letter, then numbers joined by dots.
                    let dotted = |s: &str, letter: char| {
                        s.strip_prefix(letter).is_some_and(|n| {
                            n.split('.').all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
                        })
                    };
                    assert!(dotted(pore, 'r'), "{}: pore {pore}", model.name);
                    assert!(dotted(motor, 'e'), "{}: motor {motor}", model.name);
                    let stem = clair3_name(&model.reads).expect("ONT");
                    let alias = matches!(model.pick, Pick::No(why) if why.starts_with("the old name of "));
                    let suffixed = model.name.strip_prefix(stem.as_str()).is_some_and(|rest| {
                        rest.is_empty() || rest == "_with_mv" || rest.starts_with('_') && model.pick != Pick::Yes
                    });
                    // Guppy 3.6.0 and 4.2.2's HAC model is named after both.
                    let both = model.name == "r941_prom_hac_g360+g422" && *basecaller == guppy(4, 2, 2);
                    assert!(
                        alias || suffixed || both,
                        "{} is not named as Clair3 names a model for {}",
                        model.name,
                        model.reads
                    );
                    if model.pick == Pick::Yes {
                        assert!(model.name == stem || both, "{}: a model lrpipe picks is {stem}", model.name);
                    }
                }
                Reads::Hifi { instruments, nearest_for } => {
                    assert!(model.name.starts_with("hifi"), "{}", model.name);
                    assert!(!instruments.is_empty());
                    for i in instruments.iter().chain(*nearest_for) {
                        assert!(i.chars().all(|c| c.is_ascii_uppercase()), "{}: {i}", model.name);
                    }
                }
                Reads::Illumina => assert_eq!(model.name, "ilmn"),
            }
            if let Pick::No(why) = model.pick
                && let Some(old) = why.strip_prefix("the old name of ")
            {
                let new =
                    find(old).unwrap_or_else(|| panic!("{} is the old name of {old}, which is not listed", model.name));
                assert_eq!(new.reads, model.reads, "{} and {old} are the same model", model.name);
                assert_eq!(new.pick, Pick::Yes);
            }
        }
    }

    #[test]
    fn no_two_models_lrpipe_picks_are_for_the_same_reads() {
        for (i, a) in MODELS.iter().enumerate().filter(|(_, m)| m.pick == Pick::Yes) {
            for b in MODELS[i + 1..].iter().filter(|m| m.pick == Pick::Yes) {
                let same = match (&a.reads, &b.reads) {
                    (
                        Reads::Hifi { instruments: x, nearest_for: n },
                        Reads::Hifi { instruments: y, nearest_for: m },
                    ) => x.iter().chain(*n).any(|i| y.contains(i) || m.contains(i)),
                    (x, y) => x == y,
                };
                assert!(!same, "{} and {} are both for {}", a.name, b.name, a.reads);
            }
        }
    }

    /// What Clair3 2.0.3's Bioconda package installs in `bin/models/`.
    pub(crate) const BUNDLED: &[&str] = &[
        "hifi",
        "hifi_revio",
        "hifi_sequel2",
        "ilmn",
        "ont",
        "ont_guppy5",
        "r1041_e82_400bps_hac_v410",
        "r1041_e82_400bps_hac_v500",
        "r1041_e82_400bps_hac_v520",
        "r1041_e82_400bps_hac_v520_with_mv",
        "r1041_e82_400bps_hac_v600",
        "r1041_e82_400bps_hac_v600_with_mv",
        "r1041_e82_400bps_hac_with_mv",
        "r1041_e82_400bps_sup_v410",
        "r1041_e82_400bps_sup_v430_bacteria_finetuned",
        "r1041_e82_400bps_sup_v500",
        "r1041_e82_400bps_sup_v520",
        "r1041_e82_400bps_sup_v520_with_mv",
        "r1041_e82_400bps_sup_with_mv",
        "r941_prom_hac_g360+g422",
        "r941_prom_sup_g5014",
    ];

    #[test]
    fn the_bundled_models_are_clair3s() {
        let mut bundled: Vec<&str> = MODELS.iter().filter(|m| m.source == Source::Bundled).map(|m| m.name).collect();
        bundled.sort_unstable();
        assert_eq!(bundled, BUNDLED, "the models said to be Clair3's own are those it installs");
    }

    #[test]
    fn versions_parse() {
        assert_eq!(Version::parse("v5.0.0"), Some(Version([5, 0, 0])));
        assert_eq!(Version::parse("3.3"), Some(Version([3, 3, 0])));
        assert_eq!(Version::parse("v5"), None);
        assert_eq!(Version::parse("v5.0.0.1"), None);
        assert_eq!(Version::parse("vx.0"), None);
        assert!(Version([4, 3, 0]).ordinal() < Version([5, 0, 0]).ordinal());
        assert_eq!(Version([5, 2, 0]).to_string(), "v5.2.0");
    }
}
