cd ~/srebf2-regulatory-signature

python3 -c '
import gzip
from pathlib import Path

root = Path("data/raw/fastq")
count = 0

def read_record(handle):
    header = handle.readline()
    if not header:
        return None
    seq, plus, qual = (handle.readline() for _ in range(3))
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError("Invalid FASTQ record")
    if not seq or not qual or len(seq.rstrip()) != len(qual.rstrip()):
        raise ValueError("Incomplete record or sequence/quality length mismatch")
    identifier = header.split()[0]
    if identifier.endswith((b"/1", b"/2")):
        identifier = identifier[:-2]
    return identifier

with gzip.open(root / "SRR28966297_1.fastq.gz", "rb") as a, gzip.open(root / "SRR28966297_2.fastq.gz", "rb") as b:
    while True:
        r1, r2 = read_record(a), read_record(b)
        if r1 is None and r2 is None:
            break
        if r1 is None or r2 is None or r1 != r2:
            raise ValueError(f"Mate mismatch at pair {count + 1}: {r1!r}, {r2!r}")
        count += 1

assert count == 28254500, f"Unexpected pair count: {count}"
message = f"PASS: {count:,} matched read pairs; FASTQ records valid.\n"
Path("results/qc/pilot/PAIR_VALIDATION.txt").write_text(message)
print(message)
'
