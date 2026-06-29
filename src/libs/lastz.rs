//! Lastz aligner presets and scoring matrices ported from UCSC.

/// Default scoring matrix for lastz (Human vs Mouse / Macaque / Cow).
pub const MATRIX_DEFAULT: &str = "   A    C    G    T
A  91 -114  -31 -123
C -114  100 -125  -31
G  -31 -125  100 -114
T -123  -31 -114   91
";

/// Distant-species scoring matrix (Human vs Zebrafish / Opossum).
pub const MATRIX_DISTANT: &str = "   A    C    G    T
A  91  -90  -25 -100
C -90  100 -100  -25
G -25 -100  100  -90
T -100  -25  -90   91
";

/// Close-species scoring matrix (Human vs Chimp).
pub const MATRIX_SIMILAR: &str = "   A    C    G    T
A  100 -300 -150 -300
C -300  100 -300 -150
G -150 -300  100 -300
T -300 -150 -300  100
";

/// Close-species scoring matrix variant (Human vs Primate, more sensitive).
#[allow(dead_code)]
pub const MATRIX_SIMILAR2: &str = "   A    C    G    T
A  90 -330 -236 -356
C -330  100 -318 -236
G -236 -318  100 -330
T -356 -236 -330   90
";

/// A predefined lastz parameter set with optional scoring matrix.
#[derive(Debug)]
pub struct Preset {
    pub name: &'static str,
    pub desc: &'static str,
    pub params: &'static str,
    pub matrix: Option<&'static str>,
}

/// UCSC-derived lastz presets for common pairwise vertebrate alignments.
pub const PRESETS: &[Preset] = &[
    Preset {
        name: "set01",
        desc: "Hg17vsPanTro1 (Human vs Chimp)",
        params: "C=0 E=30 K=3000 L=2200 O=400 Y=3400 Q=similar",
        matrix: Some(MATRIX_SIMILAR),
    },
    Preset {
        name: "set02",
        desc: "Hg19vsPanTro2 (Human vs Primate, more sensitive)",
        params: "C=0 E=150 H=2000 K=4500 L=2200 M=254 O=600 T=2 Y=15000 Q=similar2",
        matrix: Some(MATRIX_SIMILAR2),
    },
    Preset {
        name: "set03",
        desc: "Hg17vsMm5 (Human vs Mouse)",
        params: "C=0 E=30 K=3000 L=2200 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set04",
        desc: "Hg17vsRheMac2 (Human vs Macaque)",
        params: "C=0 E=30 H=2000 K=3000 L=2200 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set05",
        desc: "Hg17vsBosTau2 (Human vs Cow)",
        params: "C=0 E=30 H=2000 K=3000 L=2200 M=50 O=400 Q=default",
        matrix: Some(MATRIX_DEFAULT),
    },
    Preset {
        name: "set06",
        desc: "Hg17vsDanRer3 (Human vs Zebrafish)",
        params: "C=0 E=30 H=2000 K=2200 L=6000 O=400 Y=3400 Q=distant",
        matrix: Some(MATRIX_DISTANT),
    },
    Preset {
        name: "set07",
        desc: "Hg17vsMonDom1 (Human vs Opossum)",
        params: "C=0 E=30 H=2000 K=2200 L=10000 O=400 Y=3400 Q=distant",
        matrix: Some(MATRIX_DISTANT),
    },
];

/// Look up a preset by name.
pub fn find_preset(name: &str) -> Option<&'static Preset> {
    PRESETS.iter().find(|p| p.name == name)
}

/// Collect all preset names (for clap PossibleValuesParser).
pub fn preset_names() -> Vec<&'static str> {
    PRESETS.iter().map(|p| p.name).collect()
}

/// Build the preset help string used in `--help` output.
pub fn preset_help() -> String {
    let mut help = String::from("Presets from UCSC:\n");
    for p in PRESETS {
        help.push_str(&format!(
            "    {}: {}\n           {}\n",
            p.name, p.desc, p.params
        ));
    }
    help
}
