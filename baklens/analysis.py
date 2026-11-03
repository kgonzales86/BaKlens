"""Sequence calculations. Coordinates exposed to users are 1-based, inclusive."""
from collections import Counter
from io import StringIO
import math
from Bio import SeqIO
import pandas as pd

IUPAC = set("ACGTRYSWKMBDHVN")
MAX_BASES = 25_000_000
MAX_RECORDS = 10_000


def parse_fasta(text):
    if not text.strip().startswith(">"):
        raise ValueError("FASTA must begin with a > header.")
    records = {}
    total = 0
    for record in SeqIO.parse(StringIO(text.strip()), "fasta"):
        name, seq = record.id, str(record.seq).upper()
        if name in records:
            raise ValueError(f"Duplicate sequence identifier: {name}")
        if not name or not seq:
            raise ValueError("Each FASTA record needs an identifier and a non-empty sequence.")
        invalid = set(seq) - IUPAC
        if invalid:
            raise ValueError(f"{name}: invalid DNA characters: {', '.join(sorted(invalid))}")
        total += len(seq)
        if total > MAX_BASES or len(records) >= MAX_RECORDS:
            raise ValueError("This release supports up to 25 million bases and 10,000 records.")
        records[name] = seq
    if not records:
        raise ValueError("No sequences found.")
    return records


def composition(seq):
    c = Counter(seq)
    canonical = sum(c[x] for x in "ACGT")
    return {
        "length_bp": len(seq),
        "gc_pct": 100 * (c['G'] + c['C']) / canonical if canonical else math.nan,
        "ambiguous_bp": len(seq) - canonical,
        "ambiguous_pct": 100 * (len(seq) - canonical) / len(seq),
        "n_bp": c['N'],
        "gc_skew": (c['G'] - c['C']) / (c['G'] + c['C']) if c['G'] + c['C'] else math.nan,
    }


def sequence_table(records):
    return pd.DataFrame([{"sequence": name, **composition(seq)} for name, seq in records.items()])


def assembly_stats(records):
    lengths = sorted((len(s) for s in records.values()), reverse=True)
    total, running, n50, l50 = sum(lengths), 0, 0, 0
    for i, length in enumerate(lengths, 1):
        running += length
        if running >= total / 2:
            n50, l50 = length, i
            break
    counts = Counter()
    for seq in records.values():
        counts.update(seq)
    canonical = sum(counts[x] for x in "ACGT")
    return {"sequences": len(records), "total_bp": total, "n50_bp": n50, "l50": l50,
            "gc_pct": 100 * (counts['G'] + counts['C']) / canonical if canonical else math.nan,
            "ambiguous_bp": total - canonical}


def windows(seq, size=1000, step=1000):
    if size < 1 or step < 1:
        raise ValueError("Window and step must be positive.")
    if math.ceil(len(seq) / step) > 50_000:
        raise ValueError("Too many windows; increase the step to keep at most 50,000 windows.")
    rows = []
    for start in range(0, len(seq), step):
        end = min(start + size, len(seq))
        part = seq[start:end]
        rows.append({"start": start + 1, "end": end, "position": (start + 1 + end) / 2,
                     **composition(part)})
        if end == len(seq):
            break
    return pd.DataFrame(rows)


def summary_markdown(records, window_size, step, features=None, provenance=None):
    stats = assembly_stats(records)
    gc = f"{stats['gc_pct']:.2f}%" if math.isfinite(stats['gc_pct']) else "undefined (no A/C/G/T)"
    lines = ["# BaKlens analysis report", "", f"- Sequences: {stats['sequences']:,}",
             f"- Total length: {stats['total_bp']:,} bp", f"- GC content: {gc}",
             f"- Ambiguous bases: {stats['ambiguous_bp']:,}",
             f"- N50: {stats['n50_bp']:,} bp; L50: {stats['l50']:,}",
             f"- Window: {window_size:,} bp; step: {step:,} bp", "",
             "GC% = 100 × (G+C)/(A+C+G+T). Ambiguous IUPAC bases are excluded from this denominator.",
             "GC skew = (G−C)/(G+C); undefined denominators are missing, not zero.",
             "Windows are computed independently per sequence; the final window may be shorter.", ""]
    if features is not None:
        lines += [f"Annotation features: {len(features):,}."]
        for kind, count in features['type'].value_counts().items():
            lines.append(f"- {kind}: {count:,}")
        if provenance:
            lines += ["", "Annotation provenance:", "```json", provenance, "```"]
    lines += ["", "## Interpretation limits",
              "Sequence composition does not establish contamination, replication origins, gene function, or pathway activity.",
              "N50 describes assembly contiguity, not correctness. A circular plot does not establish circular topology.",
              "Annotations are predictions supplied by Bakta. This report uses calculated statistics, not generative AI."]
    return "\n".join(lines)
