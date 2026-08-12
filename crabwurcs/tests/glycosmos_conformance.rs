//! Live, opt-in conformance checks against the GlyCosmos
//! GlycanFormatConverter WURCS normalizer.
//!
//! This test is deliberately ignored for ordinary/offline CI.  Run it with
//! `--ignored` and `CRABWURCS_RUN_GLYCOSMOS=1`; set
//! `CRABWURCS_GLYCOSMOS_LIMIT=0` to exercise every distinct bundled WURCS.

use crabwurcs::core;
use serde_json::Value;
use std::collections::BTreeMap;
use std::env;
use std::process::Command;
use std::thread;
use std::time::Duration;

fn bundled_wurcs() -> BTreeMap<String, String> {
    const DATA: &str = include_str!("../data/glycoshape_notations.tsv");
    const DERIVED: &str = include_str!("../data/glycoshape_derived_notations.tsv");
    const NOTATION_ONLY: &str = include_str!("../data/glycoshape_notation_only.tsv");

    DATA.lines()
        .chain(DERIVED.lines())
        .chain(NOTATION_ONLY.lines())
        .filter_map(|line| {
            let mut fields = line.splitn(5, '\t');
            let source = fields.next()?.trim();
            if source.is_empty() {
                return None;
            }
            let expected = core::normalize_wurcs(source).ok()?;
            Some((source.to_owned(), expected))
        })
        .collect()
}

fn api_wurcs(endpoint: &str, input: &str) -> Result<String, String> {
    let body = serde_json::json!({ "input": input }).to_string();
    let mut last_error = String::new();
    for attempt in 1..=3 {
        let output = Command::new("curl")
            .args([
                "--fail",
                "--silent",
                "--show-error",
                "--max-time",
                "60",
                "-H",
                "Content-Type: application/json",
                "--data",
                &body,
                endpoint,
            ])
            .output()
            .map_err(|error| format!("could not execute curl: {error}"))?;
        if output.status.success() {
            let response: Value = serde_json::from_slice(&output.stdout)
                .map_err(|error| format!("invalid JSON response: {error}"))?;
            return response
                .get("wurcs")
                .or_else(|| response.get("WURCS"))
                .and_then(Value::as_str)
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(ToOwned::to_owned)
                .ok_or_else(|| "response did not contain a non-empty `wurcs` field".into());
        }
        last_error = format!(
            "curl exited with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        );
        if attempt < 3 {
            // GlyCosmos occasionally accepts a request and then stalls before
            // sending its JSON body. A short retry avoids treating that
            // transient service failure as a normalization mismatch.
            thread::sleep(Duration::from_secs(1));
        }
    }
    Err(format!("{last_error} (after 3 attempts)"))
}

fn requested_limit(total: usize) -> usize {
    match env::var("CRABWURCS_GLYCOSMOS_LIMIT") {
        Ok(value) if value.eq_ignore_ascii_case("all") || value == "0" => total,
        Ok(value) => value.parse::<usize>().unwrap_or_else(|_| {
            panic!("CRABWURCS_GLYCOSMOS_LIMIT must be a positive integer or 0/all")
        }),
        Err(_) => total.min(25),
    }
}

#[test]
#[ignore = "requires network/curl; set CRABWURCS_RUN_GLYCOSMOS=1"]
fn glycosmos_matches_local_normalized_wurcs() {
    if env::var("CRABWURCS_RUN_GLYCOSMOS").as_deref() != Ok("1") {
        println!("skipping live GlyCosmos check; set CRABWURCS_RUN_GLYCOSMOS=1 to enable it");
        return;
    }

    let endpoint = env::var("CRABWURCS_GLYCOSMOS_API").unwrap_or_else(|_| {
        "https://api.glycosmos.org/glycanformatconverter/2.8.2/wurcs2wurcs".into()
    });
    let records = bundled_wurcs();
    let limit = requested_limit(records.len()).min(records.len());
    assert!(
        limit > 0,
        "GlyCosmos conformance limit selected zero records"
    );

    let mut checked = 0usize;
    for (source, expected) in records.iter().take(limit) {
        let remote = api_wurcs(&endpoint, source)
            .unwrap_or_else(|error| panic!("GlyCosmos request failed for {source}: {error}"));
        let remote_normalized = core::normalize_wurcs(&remote).unwrap_or_else(|error| {
            panic!("GlyCosmos returned unsupported WURCS for {source}: {error}")
        });
        assert_eq!(
            remote, remote_normalized,
            "GlyCosmos response was not normalized for source {source}"
        );
        assert_eq!(
            remote_normalized, *expected,
            "local and GlyCosmos normalized WURCS differ for source {source}"
        );
        checked += 1;
    }

    println!(
        "GlyCosmos WURCS normalization conformance: {checked}/{total} records",
        total = records.len()
    );
}
