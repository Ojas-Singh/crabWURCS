#[test]
fn glycoshape_smiles_parse_and_canonicalize_stably() {
    let mut failures = Vec::new();
    let mut tested = 0usize;
    let records = include_str!("../../crabwurcs/data/glycoshape_notations.tsv")
        .lines()
        .chain(include_str!("../../crabwurcs/data/glycoshape_derived_notations.tsv").lines());
    for (line_number, line) in records.enumerate() {
        let Some(smiles) = line.split('\t').nth(4) else {
            failures.push(format!("line {} has no SMILES field", line_number + 1));
            continue;
        };
        tested += 1;
        let result = chematic::smiles::parse(smiles).and_then(|molecule| {
            let canonical = chematic::smiles::canonical_smiles(&molecule);
            chematic::smiles::parse(&canonical).map(|reparsed| {
                let second = chematic::smiles::canonical_smiles(&reparsed);
                (canonical, second)
            })
        });
        match result {
            Ok((first, second)) if first == second => {}
            Ok((first, second)) => failures.push(format!(
                "line {} canonical output is unstable: {first} != {second}",
                line_number + 1
            )),
            Err(error) => failures.push(format!("line {}: {error}", line_number + 1)),
        }
    }
    assert_eq!(tested, 938);
    assert!(
        failures.is_empty(),
        "{} of {tested} GlycoShape SMILES failed: {:?}",
        failures.len(),
        &failures[..failures.len().min(20)]
    );
}
