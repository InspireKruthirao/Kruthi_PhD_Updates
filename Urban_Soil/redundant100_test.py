#!/usr/bin/env python3
"""
CLI-friendly version of the original redundant100.py dedup+containment script,
for use in the Python-vs-cd-hit-est comparison test.

Usage:
    python3 redundant100_test.py <input.fna> <output_100NT.fna> <output_matches.tsv>
"""
import sys
from collections import defaultdict


def fasta_iter(path):
    """Minimal fasta reader: yields (header, sequence) tuples."""
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
        print(f"Usage: {sys.argv[0]} <input.fna> <output_100NT.fna> <output_matches.tsv>")
        sys.exit(1)

    in_path, out_100_path, out_matches_path = sys.argv[1:4]

    seqs = list(fasta_iter(in_path))

    def k(h_s):
        h, s = h_s
        return (len(s), s, h)

    print(f"Loaded {len(seqs)}", flush=True)
    seqs.sort(key=k)

    deduped = []
    matches = []
    prev = ''
    h_prev = None
    for h, seq in seqs:
        if seq == prev:
            matches.append((h, '=', h_prev))
        else:
            deduped.append((h, seq))
            prev = seq
            h_prev = h
    del seqs

    min_len = len(deduped[0][1])

    def rolling_hashes(seq, window_size):
        for i in range(len(seq) - window_size + 1):
            yield hash(seq[i:i + window_size])

    hash_matches = defaultdict(list)
    eliminated = set()
    print(f"Deduped: {len(deduped)}", flush=True)

    for ix, (h, seq) in enumerate(deduped):
        if ix % 100_000 == 0:
            print(f"Dedupe iteration {ix}", flush=True)
        candidates = set()
        for ha in set(rolling_hashes(seq, min_len)):
            candidates.update(hash_matches[ha])
            hash_matches[ha].append(ix)
        for ix2 in candidates:
            if ix2 in eliminated:
                continue
            h2, seq2 = deduped[ix2]
            assert ix2 < ix
            if seq2 in seq:
                matches.append((h2, 'C', h))
                assert h2 != h
                eliminated.add(ix2)

    with open(out_100_path, 'w') as f:
        n_ix = 0
        for ix, (h, seq) in enumerate(deduped):
            if ix in eliminated:
                continue
            n_h = f"US.100NT.TEST.{n_ix:09d}"
            f.write(f">{n_h} {h}\n{seq}\n")
            n_ix += 1

    with open(out_matches_path, 'w') as f:
        for h1, rel, h2 in matches:
            f.write(f"{h1}\t{rel}\t{h2}\n")

    print(f"100NT representatives: {n_ix}", flush=True)
    print(f"Created: {out_100_path}", flush=True)
    print(f"Created: {out_matches_path}", flush=True)


if __name__ == "__main__":
    main()
