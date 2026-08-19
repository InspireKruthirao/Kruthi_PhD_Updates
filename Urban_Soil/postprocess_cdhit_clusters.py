#!/usr/bin/env python3
"""
Converts cd-hit-est output (representative fasta + .clstr file) into the
pipeline-standard US.100NT.fna.xz / US.100NT.matches.xz format, matching
the header/ID conventions used in the Shanghai Dogs pipeline.

Usage:
    python3 postprocess_cdhit_clusters.py <cdhit_fasta> <cdhit_clstr> <out_fna_xz> <out_matches_xz>
"""
import sys
import lzma


def pad9(prefix, ix):
    return f"{prefix}.{ix:09d}"


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


def iter_clusters(clstr_path):
    """Yields lists of (header, is_representative) tuples per cluster."""
    cur = None
    with open(clstr_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>Cluster'):
                if cur:
                    yield cur
                cur = []
            else:
                parts = line.split()
                gid_token = parts[2]
                is_rep = parts[-1] == '*'
                gid = gid_token[1:]
                if gid.endswith('...'):
                    gid = gid[:-3]
                cur.append((gid, is_rep))
    if cur:
        yield cur


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <cdhit_fasta> <cdhit_clstr> <out_fna_xz> <out_matches_xz>")
        sys.exit(1)

    cdhit_fasta, cdhit_clstr, out_fna, out_matches = sys.argv[1:5]

    print("Reading cd-hit representative sequences...", flush=True)
    reps_by_header = dict(fasta_iter(cdhit_fasta))
    print(f"Loaded {len(reps_by_header)} representative sequences", flush=True)

    print("Parsing cluster file and writing outputs...", flush=True)
    n_written = 0
    with lzma.open(out_fna, 'wt') as fna_out, lzma.open(out_matches, 'wt') as matches_out:
        for cluster in iter_clusters(cdhit_clstr):
            rep_headers = [h for h, is_rep in cluster if is_rep]
            if len(rep_headers) != 1:
                print(f"WARNING: cluster does not have exactly one representative: {cluster[:3]}...", flush=True)
                continue
            rep_header = rep_headers[0]
            new_id = pad9("US.100NT", n_written)

            seq = reps_by_header.get(rep_header)
            if seq is None:
                print(f"WARNING: representative header not found in fasta: {rep_header}", flush=True)
                continue

            fna_out.write(f">{new_id} {rep_header}\n{seq}\n")

            for h, is_rep in cluster:
                if is_rep:
                    matches_out.write(f"{h}\t=\t{new_id}\n")
                else:
                    matches_out.write(f"{h}\tR\t{new_id}\n")

            n_written += 1
            if n_written % 1_000_000 == 0:
                print(f"  ...written {n_written:,} clusters", flush=True)

    print(f"DONE. Total representative sequences written: {n_written:,}", flush=True)
    print(f"Created: {out_fna}", flush=True)
    print(f"Created: {out_matches}", flush=True)


if __name__ == "__main__":
    main()
EOFcat > /work/microbiome/users/kruthi/intermediate_results/postprocess_cdhit_clusters.py << 'EOF'
#!/usr/bin/env python3
"""
Converts cd-hit-est output (representative fasta + .clstr file) into the
pipeline-standard US.100NT.fna.xz / US.100NT.matches.xz format, matching
the header/ID conventions used in the Shanghai Dogs pipeline.

Usage:
    python3 postprocess_cdhit_clusters.py <cdhit_fasta> <cdhit_clstr> <out_fna_xz> <out_matches_xz>
"""
import sys
import lzma


def pad9(prefix, ix):
    return f"{prefix}.{ix:09d}"


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


def iter_clusters(clstr_path):
    """Yields lists of (header, is_representative) tuples per cluster."""
    cur = None
    with open(clstr_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>Cluster'):
                if cur:
                    yield cur
                cur = []
            else:
                parts = line.split()
                gid_token = parts[2]
                is_rep = parts[-1] == '*'
                gid = gid_token[1:]
                if gid.endswith('...'):
                    gid = gid[:-3]
                cur.append((gid, is_rep))
    if cur:
        yield cur


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <cdhit_fasta> <cdhit_clstr> <out_fna_xz> <out_matches_xz>")
        sys.exit(1)

    cdhit_fasta, cdhit_clstr, out_fna, out_matches = sys.argv[1:5]

    print("Reading cd-hit representative sequences...", flush=True)
    reps_by_header = dict(fasta_iter(cdhit_fasta))
    print(f"Loaded {len(reps_by_header)} representative sequences", flush=True)

    print("Parsing cluster file and writing outputs...", flush=True)
    n_written = 0
    with lzma.open(out_fna, 'wt') as fna_out, lzma.open(out_matches, 'wt') as matches_out:
        for cluster in iter_clusters(cdhit_clstr):
            rep_headers = [h for h, is_rep in cluster if is_rep]
            if len(rep_headers) != 1:
                print(f"WARNING: cluster does not have exactly one representative: {cluster[:3]}...", flush=True)
                continue
            rep_header = rep_headers[0]
            new_id = pad9("US.100NT", n_written)

            seq = reps_by_header.get(rep_header)
            if seq is None:
                print(f"WARNING: representative header not found in fasta: {rep_header}", flush=True)
                continue

            fna_out.write(f">{new_id} {rep_header}\n{seq}\n")

            for h, is_rep in cluster:
                if is_rep:
                    matches_out.write(f"{h}\t=\t{new_id}\n")
                else:
                    matches_out.write(f"{h}\tR\t{new_id}\n")

            n_written += 1
            if n_written % 1_000_000 == 0:
                print(f"  ...written {n_written:,} clusters", flush=True)

    print(f"DONE. Total representative sequences written: {n_written:,}", flush=True)
    print(f"Created: {out_fna}", flush=True)
    print(f"Created: {out_matches}", flush=True)


if __name__ == "__main__":
    main()
