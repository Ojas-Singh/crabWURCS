use crate::{CoreError, CoreResult, Modification, Monosaccharide, parse_wurcs};
use std::sync::OnceLock;

macro_rules! residue_kinds {
    ($(($variant:ident, $name:literal, $unique:expr, $generic:expr, [$($alias:literal),*])),+ $(,)?) => {
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub enum ResidueKind {
            $($variant),+
        }

        impl ResidueKind {
            pub const ALL: &'static [Self] = &[$(Self::$variant),+];

            pub const fn canonical_name(self) -> &'static str {
                match self {
                    $(Self::$variant => $name),+
                }
            }

            /// Canonical WURCS 2.0 UniqueRES body, without square brackets.
            ///
            /// `None` is used for SNFG display categories whose identity cannot
            /// be preserved by a standard WURCS UniqueRES (`Sia`, `Unknown`,
            /// and `Assigned`).
            pub const fn unique_residue(self) -> Option<&'static str> {
                match self {
                    $(Self::$variant => $unique),+
                }
            }

            pub const fn is_generic(self) -> bool {
                match self {
                    $(Self::$variant => $generic),+
                }
            }

            /// Return whether this registry kind, when used as a motif query,
            /// accepts `candidate`.
            ///
            /// Concrete kinds match only themselves. Official SNFG generic
            /// classes match themselves and every named residue in their row.
            pub fn matches_family(self, candidate: Self) -> bool {
                use ResidueKind::*;
                match self {
                    Hex => matches!(candidate, Hex | Glc | Man | Gal | Gul | Alt | All | Tal | Ido),
                    HexNAc => matches!(
                        candidate,
                        HexNAc | GlcNAc | ManNAc | GalNAc | GulNAc | AltNAc | AllNAc | TalNAc | IdoNAc
                    ),
                    HexN => matches!(
                        candidate,
                        HexN | GlcN | ManN | GalN | GulN | AltN | AllN | TalN | IdoN
                    ),
                    HexA => matches!(
                        candidate,
                        HexA | GlcA | ManA | GalA | GulA | AltA | AllA | TalA | IdoA
                    ),
                    DHex => matches!(candidate, DHex | Qui | Rha | SixDGul | SixDAlt | SixDTal | Fuc),
                    DHexNAc => matches!(
                        candidate,
                        DHexNAc | QuiNAc | RhaNAc | SixDAltNAc | SixDTalNAc | FucNAc
                    ),
                    DDHex => matches!(candidate, DDHex | Oli | Tyv | Abe | Par | Dig | Col),
                    Pen => matches!(candidate, Pen | Ara | Lyx | Xyl | Rib),
                    NulO => matches!(candidate, NulO | Kdn | Neu5Ac | Neu5Gc | Neu | Sia),
                    Sia => matches!(candidate, Sia | Kdn | Neu5Ac | Neu5Gc | Neu),
                    DDNulO => matches!(candidate, DDNulO | Pse | Leg | Aci | FourELeg),
                    _ => self == candidate,
                }
            }

            pub fn aliases(self) -> &'static [&'static str] {
                match self {
                    $(Self::$variant => &[$($alias),*]),+
                }
            }

            pub fn from_name(name: &str) -> Option<Self> {
                Self::ALL.iter().copied().find(|kind| {
                    kind.canonical_name().eq_ignore_ascii_case(name)
                        || kind.aliases()
                            .iter()
                            .any(|alias| alias.eq_ignore_ascii_case(name))
                })
            }
        }
    };
}

// UniqueRES definitions are the canonical GlyCosmos Composition-to-WURCS
// templates for the SNFG 2.0.4 table.
residue_kinds! {
    (Hex, "Hex", Some("axxxxh-1x_1-5"), true, ["Hexose"]),
    (Glc, "Glc", Some("a2122h-1x_1-5"), false, ["Glucose"]),
    (Man, "Man", Some("a1122h-1x_1-5"), false, ["Mannose"]),
    (Gal, "Gal", Some("a2112h-1x_1-5"), false, ["Galactose"]),
    (Gul, "Gul", Some("a2212h-1x_1-5"), false, ["Gulose"]),
    (Alt, "Alt", Some("a1222h-1x_1-5"), false, ["Altrose"]),
    (All, "All", Some("a1111h-1x_1-5"), false, ["Allose"]),
    (Tal, "Tal", Some("a1112h-1x_1-5"), false, ["Talose"]),
    (Ido, "Ido", Some("a2121h-1x_1-5"), false, ["Idose"]),

    (HexNAc, "HexNAc", Some("axxxxh-1x_1-5_2*NCC/3=O"), true, ["N-Acetylhexosamine"]),
    (GlcNAc, "GlcNAc", Some("a2122h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylglucosamine"]),
    (ManNAc, "ManNAc", Some("a1122h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylmannosamine"]),
    (GalNAc, "GalNAc", Some("a2112h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylgalactosamine"]),
    (GulNAc, "GulNAc", Some("a2212h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylgulosamine"]),
    (AltNAc, "AltNAc", Some("a1222h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylaltrosamine"]),
    (AllNAc, "AllNAc", Some("a1111h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylallosamine"]),
    (TalNAc, "TalNAc", Some("a1112h-1x_1-5_2*NCC/3=O"), false, ["N-Acetyltalosamine"]),
    (IdoNAc, "IdoNAc", Some("a2121h-1x_1-5_2*NCC/3=O"), false, ["N-Acetylidosamine"]),

    (HexN, "HexN", Some("axxxxh-1x_1-5_2*N"), true, ["Hexosamine"]),
    (GlcN, "GlcN", Some("a2122h-1x_1-5_2*N"), false, ["Glucosamine"]),
    (ManN, "ManN", Some("a1122h-1x_1-5_2*N"), false, ["Mannosamine"]),
    (GalN, "GalN", Some("a2112h-1x_1-5_2*N"), false, ["Galactosamine"]),
    (GulN, "GulN", Some("a2212h-1x_1-5_2*N"), false, ["Gulosamine"]),
    (AltN, "AltN", Some("a1222h-1x_1-5_2*N"), false, ["Altrosamine"]),
    (AllN, "AllN", Some("a1111h-1x_1-5_2*N"), false, ["Allosamine"]),
    (TalN, "TalN", Some("a1112h-1x_1-5_2*N"), false, ["Talosamine"]),
    (IdoN, "IdoN", Some("a2121h-1x_1-5_2*N"), false, ["Idosamine"]),

    (HexA, "HexA", Some("axxxxA-1x_1-5"), true, ["Hexuronate", "HexuronicAcid"]),
    (GlcA, "GlcA", Some("a2122A-1x_1-5"), false, ["GlucuronicAcid"]),
    (ManA, "ManA", Some("a1122A-1x_1-5"), false, ["MannuronicAcid"]),
    (GalA, "GalA", Some("a2112A-1x_1-5"), false, ["GalacturonicAcid"]),
    (GulA, "GulA", Some("a2212A-1x_1-5"), false, ["GuluronicAcid"]),
    (AltA, "AltA", Some("a1222A-1x_1-5"), false, ["AltruronicAcid"]),
    (AllA, "AllA", Some("a1111A-1x_1-5"), false, ["AlluronicAcid"]),
    (TalA, "TalA", Some("a1112A-1x_1-5"), false, ["TaluronicAcid"]),
    (IdoA, "IdoA", Some("a2121A-1x_1-5"), false, ["IduronicAcid"]),

    (DHex, "dHex", Some("axxxxm-1x_1-5"), true, ["Deoxyhexose"]),
    (Qui, "Qui", Some("a2122m-1x_1-5"), false, ["Quinovose"]),
    (Rha, "Rha", Some("a2211m-1x_1-5"), false, ["Rhamnose"]),
    (SixDGul, "6dGul", Some("a2212m-1x_1-5"), false, ["6-DeoxyGulose"]),
    (SixDAlt, "6dAlt", Some("a2111m-1x_1-5"), false, ["6-DeoxyAltrose"]),
    (SixDTal, "6dTal", Some("a1112m-1x_1-5"), false, ["6-DeoxyTalose"]),
    (Fuc, "Fuc", Some("a1221m-1x_1-5"), false, ["Fucose"]),

    (DHexNAc, "dHexNAc", Some("axxxxm-1x_1-5_2*NCC/3=O"), true, ["DeoxyhexNAc"]),
    (QuiNAc, "QuiNAc", Some("a2122m-1x_1-5_2*NCC/3=O"), false, ["N-Acetylquinovosamine"]),
    (RhaNAc, "RhaNAc", Some("a2211m-1x_1-5_2*NCC/3=O"), false, ["N-Acetylrhamnosamine"]),
    (SixDAltNAc, "6dAltNAc", Some("a2111m-1x_1-5_2*NCC/3=O"), false, ["N-Acetyl6-DeoxyAltrosamine"]),
    (SixDTalNAc, "6dTalNAc", Some("a1112m-1x_1-5_2*NCC/3=O"), false, ["N-Acetyl6-DeoxyTalosamine"]),
    (FucNAc, "FucNAc", Some("a1221m-1x_1-5_2*NCC/3=O"), false, ["N-Acetylfucosamine"]),

    (DDHex, "ddHex", Some("adxxxm-1x_1-5"), true, ["Di-deoxyhexose", "Dideoxyhexose"]),
    (Oli, "Oli", Some("ad122m-1x_1-5"), false, ["Olivose"]),
    (Tyv, "Tyv", Some("a1d22m-1x_1-5"), false, ["Tyvelose"]),
    (Abe, "Abe", Some("a2d12m-1x_1-5"), false, ["Abequose"]),
    (Par, "Par", Some("a2d22m-1x_1-5"), false, ["Paratose"]),
    (Dig, "Dig", Some("ad222m-1x_1-5"), false, ["Digitoxose"]),
    (Col, "Col", Some("a1d21m-1x_1-5"), false, ["Colitose"]),

    (Pen, "Pen", Some("axxxh-1x_1-5"), true, ["Pentose"]),
    (Ara, "Ara", Some("a211h-1x_1-5"), false, ["Arabinose"]),
    (Lyx, "Lyx", Some("a112h-1x_1-5"), false, ["Lyxose"]),
    (Xyl, "Xyl", Some("a212h-1x_1-5"), false, ["Xylose"]),
    (Rib, "Rib", Some("a222h-1x_1-5"), false, ["Ribose"]),

    (NulO, "NulO", Some("Aadxxxxxh-2x_2-6"), true, ["3-deoxy-nonulosonicAcid", "Nonulosonate"]),
    (Kdn, "Kdn", Some("Aad21122h-2x_2-6"), false, ["KDN"]),
    (Neu5Ac, "Neu5Ac", Some("Aad21122h-2x_2-6_5*NCC/3=O"), false, ["NeuAc", "Neup5Ac"]),
    (Neu5Gc, "Neu5Gc", Some("Aad21122h-2x_2-6_5*NCCO/3=O"), false, ["NeuGc", "Neup5Gc"]),
    (Neu, "Neu", Some("Aad21122h-2x_2-6_5*N"), false, ["NeuraminicAcid"]),
    (Sia, "Sia", None, true, ["SialicAcid"]),

    (DDNulO, "ddNulO", Some("Aadxxxxxm-2x_2-6_5*N_7*N"), true, ["3,9-dideoxy-nonulosonicAcid"]),
    (Pse, "Pse", Some("Aad22111m-2x_2-6_5*N_7*N"), false, ["PseudaminicAcid"]),
    (Leg, "Leg", Some("Aad21122m-2x_2-6_5*N_7*N"), false, ["LegionaminicAcid"]),
    (Aci, "Aci", Some("Aad21111m-2x_2-6_5*N_7*N"), false, ["AcinetaminicAcid"]),
    (FourELeg, "4eLeg", Some("Aad11122m-2x_2-6_5*N_7*N"), false, ["4-epi-Leg"]),

    (Unknown, "Unknown", None, true, ["UnknownSaccharide"]),
    (Bac, "Bac", Some("a2122m-1x_1-5_2*N_4*N"), false, ["Bacillosamine"]),
    (LDManHep, "LDmanHep", Some("a11221h-1x_1-5"), false, ["L-glycero-D-manno-Heptose"]),
    (Kdo, "Kdo", Some("Aad1122h-2x_2-6"), false, ["KDO"]),
    (Dha, "Dha", Some("Aad112A-2x_2-6"), false, []),
    (DDManHep, "DDmanHep", Some("a11222h-1x_1-5"), false, ["D-glycero-D-manno-Heptose"]),
    (MurNAc, "MurNAc", Some("a2122h-1x_1-5_2*NCC/3=O_3*OC^RCO/4=O/3C"), false, []),
    (MurNGc, "MurNGc", Some("a2122h-1x_1-5_2*NCCO/3=O_3*OC^RCO/4=O/3C"), false, []),
    (Mur, "Mur", Some("a2122h-1x_1-5_2*N_3*OC^RCO/4=O/3C"), false, ["MuramicAcid"]),

    (Assigned, "Assigned", None, true, []),
    (Api, "Api", Some("a15h-1x_1-4_3*CO"), false, ["Apiose"]),
    (Fru, "Fru", Some("ha122h-2x_2-6"), false, ["Fructose"]),
    (Tag, "Tag", Some("ha112h-2x_2-6"), false, ["Tagatose"]),
    (Sor, "Sor", Some("ha121h-2x_2-6"), false, ["Sorbose"]),
    (Psi, "Psi", Some("ha222h-2x_2-6"), false, ["Psicose"])
}

pub fn residue_from_kind(kind: ResidueKind) -> CoreResult<Monosaccharide> {
    let unique = kind
        .unique_residue()
        .ok_or_else(|| CoreError::UnrepresentableResidue(kind.canonical_name().into()))?;
    let graph = parse_wurcs(&format!("WURCS=2.0/1,1,0/[{unique}]/1/"))?;
    let mut residue = graph
        .root()
        .and_then(|root| graph.residue(root))
        .cloned()
        .ok_or_else(|| CoreError::MalformedGraph("registry residue has no root".into()))?;
    residue.residue_kind = Some(kind);
    Ok(residue)
}

fn prefix_class(value: &str) -> char {
    value
        .chars()
        .next()
        .filter(|character| matches!(character, 'A' | 'h'))
        .unwrap_or('a')
}

fn modification_matches(
    required: &Modification,
    actual: &Modification,
    allow_n_derivatives: bool,
) -> bool {
    if required.position != actual.position {
        return false;
    }
    if required.descriptor == "N" {
        return actual.descriptor == "N"
            || (allow_n_derivatives && actual.descriptor.starts_with('N'));
    }
    if required.descriptor.starts_with("NCCO") {
        return actual.descriptor.starts_with("NCCO");
    }
    if required.descriptor.starts_with("NCC") {
        return actual.descriptor.starts_with("NCC") && !actual.descriptor.starts_with("NCCO");
    }
    required.descriptor == actual.descriptor
}

/// Classify a WURCS residue as an official SNFG 2.0.4 residue.
///
/// Extra substituents such as sulfation or methylation do not change the
/// underlying monosaccharide identity. The most chemically specific matching
/// registry entry wins.
pub fn classify_residue(residue: &Monosaccharide) -> Option<ResidueKind> {
    if let Some(kind) = residue.residue_kind {
        return Some(kind);
    }
    if residue.display_name.is_some() {
        return Some(ResidueKind::Assigned);
    }

    static TEMPLATES: OnceLock<Vec<(ResidueKind, Monosaccharide)>> = OnceLock::new();
    let templates = TEMPLATES.get_or_init(|| {
        ResidueKind::ALL
            .iter()
            .copied()
            .filter_map(|kind| {
                let mut residue = residue_from_kind(kind).ok()?;
                residue.residue_kind = None;
                Some((kind, residue))
            })
            .collect()
    });
    let mut best: Option<(usize, ResidueKind)> = None;
    for (kind, template) in templates {
        if template.skeleton_code != residue.skeleton_code
            || prefix_class(&template.anomeric_prefix) != prefix_class(&residue.anomeric_prefix)
        {
            continue;
        }
        let allow_n_derivatives = matches!(
            kind,
            ResidueKind::Bac
                | ResidueKind::DDNulO
                | ResidueKind::Pse
                | ResidueKind::Leg
                | ResidueKind::Aci
                | ResidueKind::FourELeg
        );
        if !template.modifications.iter().all(|required| {
            residue
                .modifications
                .iter()
                .any(|actual| modification_matches(required, actual, allow_n_derivatives))
        }) {
            continue;
        }
        let score = template.modifications.len();
        if best.is_none_or(|(best_score, _)| score > best_score) {
            best = Some((score, *kind));
        }
    }
    best.map(|(_, kind)| kind)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_serializable_registry_entry_constructs_and_classifies() {
        for kind in ResidueKind::ALL {
            let Some(unique_residue) = kind.unique_residue() else {
                continue;
            };
            let mut residue = residue_from_kind(*kind).unwrap();
            let mut graph = crate::ResidueGraph::new();
            graph.add_residue(residue.clone());
            assert_eq!(
                crate::write_wurcs(&graph).unwrap(),
                format!("WURCS=2.0/1,1,0/[{unique_residue}]/1/"),
                "{kind:?}"
            );
            residue.residue_kind = None;
            assert_eq!(classify_residue(&residue), Some(*kind), "{kind:?}");
        }
    }
}
