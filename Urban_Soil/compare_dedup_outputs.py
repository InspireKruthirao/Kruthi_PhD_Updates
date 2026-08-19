#!/usr/bin/env python3
"""
Compares the representative sequence sets produced by:
  1. The original Python dedup+containment script
  2. cd-hit-est (-c 1.0 -aS 1.0)

Reports how many representative sequences match, and flags any differences
so you can decide whether cd-hit-est is a safe drop-in replacement.

Usage:
    python3 compare_dedup_outputs.py <python_output.fna> <cdhit_output.fna> <report.txt>
"""
import sys


def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, ''.join(seq_chunks)


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <python_output.fna> <cdhit_output.fna> <report.txt>")
        sys.exit(1)

    py_path, cdhit_path, report_path = sys.argv[1:4]

    py_seqs = set(seq for _, seq in fasta_iter(py_path))
    cdhit_seqs = set(seq for _, seq in fasta_iter(cdhit_path))

    only_py = py_seqs - cdhit_seqs
    only_cdhit = cdhit_seqs - py_seqs
    shared = py_seqs & cdhit_seqs

    lines = []
    lines.append("=== Dedup Comparison Report: Python vs cd-hit-est ===\n")
    lines.append(f"Python script representatives:    {len(py_seqs):,}")
    lines.append(f"cd-hit-est representatives:        {len(cdhit_seqs):,}")
    lines.append(f"Shared (identical sequence, in both): {len(shared):,}")
    lines.append(f"Only in Python output (missed by cd-hit-est): {len(only_py):,}")
    lines.append(f"Only in cd-hit-est output (missed by Python): {len(only_cdhit):,}")
    lines.append("")

    if len(py_seqs) > 0:
        pct_shared = 100 * len(shared) / len(py_seqs)
        lines.append(f"Agreement rate (shared / Python total): {pct_shared:.2f}%")

    if only_py:
        lines.append("\n--- Sample sequences ONLY in Python output (first 5) ---")
        for s in list(only_py)[:5]:
            lines.append(f"  len={len(s)}  {s[:80]}{'...' if len(s) > 80 else ''}")

    if only_cdhit:
        lines.append("\n--- Sample sequences ONLY in cd-hit-est output (first 5) ---")
        for s in list(only_cdhit)[:5]:
            lines.append(f"  len={len(s)}  {s[:80]}{'...' if len(s) > 80 else ''}")

    lines.append("\n=== Interpretation guide ===")
    lines.append("- High agreement (>99%) with only a handful of edge-case diffs: cd-hit-est")
    lines.append("  is very likely safe to use for the full 265M-sequence run.")
    lines.append("- Large disagreement: containment-matching logic differs meaningfully between")
    lines.append("  the two tools and needs further investigation before scaling up.")

    report = '\n'.join(lines)
    print(report)
    with open(report_path, 'w') as f:
        f.write(report + '\n')


if __name__ == "__main__":
    main()
