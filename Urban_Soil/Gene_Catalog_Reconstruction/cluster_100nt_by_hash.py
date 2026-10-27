#!/usr/bin/env python3
"""
100% nucleotide clustering without MMseqs2.

100% identity and 100% coverage means identical sequences, so ORFs are
grouped by a hash of the sequence (MMseqs2 also matches the reverse
complement, so the hash is of whichever of the two comes first). Replaces
run_mmseqs_linclust_splitmem.sh + Run_US_100NT_catalog.py.

    python3 cluster_100nt_by_hash.py ORF_FASTA OUT_DIR
        -> OUT_DIR/US.100NT.fna.zst      representatives (first ORF of each group)
           OUT_DIR/US.100NT.matches.zst  ORF  =/C  US.100NT ID (as before)

    python3 cluster_100nt_by_hash.py --compare OLD_MATCHES NEW_MATCHES
        checks that both group the ORFs in exactly the same way

ORF_FASTA may be plain, .gz, .xz or .zst. Hash: 128-bit BLAKE2b (chance of
two different sequences sharing one is negligible). ORF and group IDs end in
a number (US.ORF.000_000_123). Memory: about 64 GB for 265M ORFs.
"""

import hashlib
from array import array
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

PREFIX = "US.100NT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reader(path):
    """Text lines of a plain or compressed file (decompressed by a fast tool)."""
    path = str(path)
    for ext, cmd in ((".xz", ["xz", "-dc", "-T0"]), (".zst", ["zstd", "-dc"]),
                     (".gz", ["pigz", "-dc"] if shutil.which("pigz") else ["gzip", "-dc"])):
        if path.endswith(ext):
            proc = subprocess.Popen(cmd + [path], stdout=subprocess.PIPE, bufsize=1 << 24)
            yield from (line.decode() for line in proc.stdout)
            if proc.wait():
                sys.exit(f"{cmd[0]} failed on {path}")
            return
    with open(path) as fh:
        yield from fh


def writer(path):
    """zstd-compressed output (multi-threaded)."""
    out = open(path, "wb")
    proc = subprocess.Popen(["zstd", "-q", "-T0", "-c"], stdin=subprocess.PIPE, stdout=out)
    return proc, out


def fasta(path):
    """(ORF ID, sequence)"""
    name, seq = None, []
    for line in reader(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def pad9(i):
    s = f"{i:09d}"
    return f"{PREFIX}.{s[:3]}_{s[3:6]}_{s[6:]}"


def cluster(fasta_path, out_dir):
    t0 = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not shutil.which("zstd"):
        sys.exit("zstd not found")

    # pass 1: hash every ORF (two 64-bit numbers per ORF)
    chunks, block = [], []
    for n, (_, seq) in enumerate(fasta(fasta_path), 1):
        seq = seq.upper()
        rc = seq.translate(COMPLEMENT)[::-1]
        block.append(hashlib.blake2b(min(seq, rc).encode(), digest_size=16).digest())
        if len(block) == 10_000_000:
            chunks.append(np.frombuffer(b"".join(block), dtype=np.uint64))
            block = []
            print(f"  hashed {n:,} ORFs ({time.time() - t0:.0f}s)", flush=True)
    if block:
        chunks.append(np.frombuffer(b"".join(block), dtype=np.uint64))
    h = (np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.uint64)).reshape(-1, 2)
    del chunks, block
    print(f"pass 1: {len(h):,} ORFs hashed ({time.time() - t0:.0f}s)", flush=True)

    # group equal hashes; groups numbered in order of their first ORF
    order = np.lexsort((h[:, 1], h[:, 0]))  # stable: equal hashes stay in input order
    hs = h[order]
    new = np.r_[True, (hs[1:] != hs[:-1]).any(axis=1)]
    del h, hs
    run = np.cumsum(new) - 1  # group of each sorted ORF
    first = order[new]  # first ORF of each group
    rank = np.empty(len(first), dtype=np.int64)
    rank[np.argsort(first, kind="stable")] = np.arange(len(first))
    group = np.empty(len(order), dtype=np.int64)
    group[order] = rank[run]
    rep_of = np.sort(first)  # representative ORF of group 0, 1, ...
    del order, run, new, first, rank
    print(f"groups: {len(rep_of):,} ({time.time() - t0:.0f}s)", flush=True)

    # pass 2: write the representatives and the matches
    fna, fna_file = writer(out_dir / f"{PREFIX}.fna.zst.part")
    mat, mat_file = writer(out_dir / f"{PREFIX}.matches.zst.part")
    n_rep = 0
    for i, (name, seq) in enumerate(fasta(fasta_path)):
        g = int(group[i])
        gid = pad9(g)
        if rep_of[g] == i:
            fna.stdin.write(f">{gid} {name}\n{seq}\n".encode())
            mat.stdin.write(f"{name}\t=\t{gid}\n".encode())
            n_rep += 1
        else:
            mat.stdin.write(f"{name}\tC\t{gid}\n".encode())
        if (i + 1) % 20_000_000 == 0:
            print(f"  wrote {i + 1:,} ORFs ({time.time() - t0:.0f}s)", flush=True)
    for proc, fh in ((fna, fna_file), (mat, mat_file)):
        proc.stdin.close()
        if proc.wait():
            sys.exit("zstd failed")
        fh.close()
    if i + 1 != len(group) or n_rep != len(rep_of):
        sys.exit("the FASTA changed between the two passes")
    for name in (f"{PREFIX}.fna.zst", f"{PREFIX}.matches.zst"):
        (out_dir / f"{name}.part").replace(out_dir / name)

    print(f"\nORFs: {len(group):,}\n100NT representatives: {len(rep_of):,}\n"
          f"collapsed: {len(group) - len(rep_of):,}\nDone in {(time.time() - t0) / 60:.1f} min")


def number(name):
    """US.ORF.000_000_123 -> 123"""
    return int(name.rsplit(".", 1)[1].replace("_", ""))


def groups(path):
    """ORF numbers (sorted) and their group numbers, from a matches file"""
    orf, grp = array("q"), array("q")
    for n, line in enumerate(reader(path), 1):
        o, _, g = line.rstrip("\n").split("\t")
        orf.append(number(o))
        grp.append(number(g))
        if n % 50_000_000 == 0:
            print(f"  read {n:,} lines of {path}", flush=True)
    orf, grp = np.frombuffer(orf, dtype=np.int64), np.frombuffer(grp, dtype=np.int64)
    order = np.argsort(orf, kind="stable")
    return orf[order], grp[order]


def compare(old_path, new_path):
    """Same ORFs, and the same ORFs grouped together?"""
    old_orf, old = groups(old_path)
    new_orf, new = groups(new_path)
    if not np.array_equal(old_orf, new_orf):
        sys.exit(f"different ORFs: {len(old_orf):,} in old, {len(new_orf):,} in new")
    n_old, n_new = len(np.unique(old)), len(np.unique(new))
    n_pairs = len(np.unique(old * (new.max() + 1) + new))
    print(f"ORFs: {len(old_orf):,}\ngroups: old {n_old:,}, new {n_new:,}, (old, new) pairs {n_pairs:,}")
    if n_pairs == n_old == n_new:
        print("SAME: every old group is exactly one new group")
    else:
        print(f"DIFFERENT: {n_pairs - n_new:,} extra splits of new groups, {n_pairs - n_old:,} of old groups")


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 3:
        cluster(sys.argv[1], sys.argv[2])
    else:
        sys.exit(__doc__)
