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
    RenderOptions, ResidueKind, render_svg_with_motifs, render_svg_with_options, render_symbol_svg,
};

const TARGET: &str = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc";
const MOTIF: &str = "Gal(b1-4)[Fuc(a1-3)]GlcNAc";
const COMPOSITION: &str = "{Hex}3,{HexNAc}2,{dHex}1";

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

    println!(
        "Generated {symbol_count} symbol SVGs in {}",
        symbols_dir.display()
    );
    println!("Generated example figures in {}", examples_dir.display());
    Ok(())
}
