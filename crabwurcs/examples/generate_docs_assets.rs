//! Regenerates the documentation image assets under `crabwurcs/docs/img/`:
//!   - one tightly-cropped SNFG symbol SVG per registered `ResidueKind`
//!   - example figures used by `rendering.md` and `motif-highlighting.md`
//!
//! Run from anywhere in the workspace:
//!
//! ```text
//! cargo run -p crabwurcs --example generate_docs_assets
//! ```
//!
//! The output directory is resolved from `CARGO_MANIFEST_DIR`, so the command
//! works regardless of the current working directory. Re-run it whenever the
//! symbol table, palette, or layout in `crabwurcs-snfg` changes, then commit
//! the generated files so they ship with the published crate and render on
//! docs.rs.

use std::fs;
use std::path::PathBuf;

use crabwurcs::{
    Adduct, Format, IonKind, MassOptions, Orientation, RenderOptions, ResidueKind, cleavage_mark,
    glycan_fragments, parse_notation, render_svg_with_motifs, render_svg_with_options,
    render_symbol_svg,
};

const TARGET: &str = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc";
const MOTIF: &str = "Gal(b1-4)[Fuc(a1-3)]GlcNAc";
const COMPOSITION: &str = "{Hex}3,{HexNAc}2,{dHex}1";
const N_GLYCAN: &str =
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc";
const REPEAT: &str = "WURCS=2.0/4,5,5/[a2112h-1b_1-5_2*NCC/3=O][a2122A-1b_1-5][a2112h-1a_1-5][a2122h-1b_1-5]/1-2-1-3-4/a3-b1_b4-c1_c4-d1_d3-e1_a1-e4~2-5";
const FRAGMENT: &str = "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}";

fn docs_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("docs")
        .join("img")
}

fn residue_for_legend(kind: ResidueKind) -> crabwurcs::core::Monosaccharide {
    match crabwurcs::residue_from_kind(kind) {
        Ok(residue) => residue,
        Err(_) => {
            // Generic classes without a canonical UniqueRES (Sia, Unknown,
            // Assigned) have no skeleton of their own. Carry them on a Hex
            // backbone; `classify_residue` consults `residue_kind` first, so
            // the skeleton never affects the rendered symbol.
            let mut residue = crabwurcs::residue_from_kind(ResidueKind::Hex).unwrap();
            residue.residue_kind = Some(kind);
            if matches!(kind, ResidueKind::Assigned) {
                residue.display_name = Some("Custom".into());
            }
            residue
        }
    }
}

fn main() -> std::io::Result<()> {
    let opts = RenderOptions::default();
    let root = docs_root();
    let symbols_dir = root.join("symbols");
    let examples_dir = root.join("examples");
    fs::create_dir_all(&symbols_dir)?;
    fs::create_dir_all(&examples_dir)?;

    let mut symbol_count = 0usize;
    for &kind in ResidueKind::ALL {
        let residue = residue_for_legend(kind);
        let svg = render_symbol_svg(&residue, &opts)
            .unwrap_or_else(|err| panic!("symbol for {kind:?}: {err}"));
        let path = symbols_dir.join(format!("{}.svg", kind.canonical_name()));
        fs::write(&path, svg)?;
        symbol_count += 1;
    }

    // Motif-highlight before/after pair.
    let target = crabwurcs::iupac::parse_iupac_condensed(TARGET)
        .unwrap_or_else(|err| panic!("target {TARGET:?}: {err}"));
    let motif = crabwurcs::iupac::parse_iupac_condensed(MOTIF)
        .unwrap_or_else(|err| panic!("motif {MOTIF:?}: {err}"));
    fs::write(
        examples_dir.join("motif-before.svg"),
        render_svg_with_options(&target, &opts).unwrap(),
    )?;
    fs::write(
        examples_dir.join("motif-after.svg"),
        render_svg_with_motifs(&target, &[motif], &opts).unwrap(),
    )?;

    // Composition layout figure.
    let composition = crabwurcs::iupac::parse_iupac_condensed(COMPOSITION)
        .unwrap_or_else(|err| panic!("composition {COMPOSITION:?}: {err}"));
    fs::write(
        examples_dir.join("composition.svg"),
        render_svg_with_options(&composition, &opts).unwrap(),
    )?;

    // The opt-in layouts of rendering.md.
    let glycan = |notation| parse_notation(notation, Format::Auto).unwrap();
    let write = |name: &str, notation, options: &RenderOptions| {
        fs::write(
            examples_dir.join(name),
            render_svg_with_options(&glycan(notation), options).unwrap(),
        )
    };
    let vertical = RenderOptions {
        orientation: Orientation::RootBottom,
        reducing_end: Some("Asn".into()),
        ..opts.clone()
    };
    write("vertical.svg", N_GLYCAN, &vertical)?;
    let compact = RenderOptions {
        compact: true,
        ..opts.clone()
    };
    write("compact.svg", N_GLYCAN, &compact)?;
    let brackets = RenderOptions {
        repeat_brackets: true,
        ..opts.clone()
    };
    write("repeat.svg", REPEAT, &brackets)?;
    let floating = RenderOptions {
        floating_fragments: true,
        ..opts.clone()
    };
    write("floating.svg", FRAGMENT, &floating)?;
    let lacnac = "Gal(b1-4)GlcNAc(b1-2)Man";
    let cleavages = RenderOptions {
        cleavages: glycan_fragments(&glycan(lacnac), &MassOptions::default(), Adduct::Proton)
            .unwrap()
            .iter()
            .filter(|ion| matches!(ion.kind, IonKind::B | IonKind::Y))
            .filter(|ion| ion.edge == 1)
            .map(cleavage_mark)
            .collect(),
        ..opts.clone()
    };
    write("cleavage.svg", lacnac, &cleavages)?;

    println!(
        "Generated {symbol_count} symbol SVGs in {}",
        symbols_dir.display()
    );
    println!("Generated example figures in {}", examples_dir.display());
    Ok(())
}
