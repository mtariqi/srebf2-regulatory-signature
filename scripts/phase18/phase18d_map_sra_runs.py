#!/usr/bin/env python3
import csv
import io
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

samples = [
    ("GSM8366913", "SRX25139996", "shControl", "1"),
    ("GSM8366914", "SRX25139997", "shControl", "2"),
    ("GSM8366917", "SRX25140000", "shSREBP2_664", "1"),
    ("GSM8366918", "SRX25140001", "shSREBP2_664", "2"),
    ("GSM8366919", "SRX25140002", "shSREBP2_665", "1"),
    ("GSM8366920", "SRX25140003", "shSREBP2_665", "2"),
    ("GSM8366921", "SRX25140004", "shSREBP2_667", "1"),
    ("GSM8366922", "SRX25140005", "shSREBP2_667", "2"),
]

fastq_dir = Path("data/raw/validation/GSE271000/fastq")
records = []
errors = []

for gsm, srx, condition, replicate in samples:
    params = urlencode({
        "accession": srx,
        "result": "read_run",
        "fields": "run_accession,experiment_accession",
        "format": "tsv"
    })

    url = (
        "https://www.ebi.ac.uk/ena/portal/api/"
        "filereport?" + params
    )

    try:
        with urlopen(url, timeout=30) as response:
            data = response.read().decode("utf-8")

        rows = list(csv.DictReader(
            io.StringIO(data), delimiter="\t"
        ))

        if not rows:
            raise ValueError("No run returned")

        for row in rows:
            srr = row["run_accession"]
            if row["experiment_accession"] != srx:
                raise ValueError("Experiment accession mismatch")

            r1 = fastq_dir / f"{srr}_1.fastq.gz"
            r2 = fastq_dir / f"{srr}_2.fastq.gz"

            records.append({
                "gsm": gsm,
                "srx": srx,
                "srr": srr,
                "condition": condition,
                "replicate": replicate,
                "fastq_r1_exists": r1.exists(),
                "fastq_r2_exists": r2.exists()
            })

    except Exception as exc:
        errors.append(f"{gsm}: {exc}")

    time.sleep(0.2)

out = Path("metadata/phase18/GSE271000_sample_run_map.tsv")
out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "gsm", "srx", "srr", "condition",
            "replicate", "fastq_r1_exists",
            "fastq_r2_exists"
        ],
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(records)

print("=== PHASE 18D: SAMPLE-TO-RUN MAPPING ===")

for item in records:
    status = (
        "READY" if
        item["fastq_r1_exists"] and
        item["fastq_r2_exists"]
        else "MISSING_FASTQ"
    )
    print(
        item["gsm"], item["srr"],
        item["condition"], status
    )

if errors:
    print("\nLookup errors:")
    for err in errors:
        print(err)

ready = sum(
    r["fastq_r1_exists"] and r["fastq_r2_exists"]
    for r in records
)

print(f"\nReady runs: {ready}/{len(records)}")
print("Saved:", out)

if len(records) == 8 and ready == 8 and not errors:
    print("PASS: All expected samples mapped")
else:
    print("INCOMPLETE: Review missing runs or lookup errors")
