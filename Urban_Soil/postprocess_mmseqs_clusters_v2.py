cat > /work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal/postprocess_mmseqs_clusters_v2.py << 'PYEOF'
import lzma

BASE = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"
MMSEQS_DIR = f"{BASE}/mmseqs_test"

REP_FASTA   = f"{MMSEQS_DIR}/US.mmseqs.100NT.splitmem_rep_seq.fasta"
CLUSTER_TSV = f"{MMSEQS_DIR}/US.mmseqs.100NT.splitmem_cluster.tsv"

OUT_FNA     = f"{BASE}/US.100NT.fna.xz"
OUT_MATCHES = f"{BASE}/US.100NT.matches.xz"


def pad9(prefix, ix):
    return f"{prefix}.{ix:09d}"


def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def xz_out(path):
    return lzma.open(path, "wb")


def main():
    print(f"Pass 1: assigning new IDs from {REP_FASTA} (streaming, no seqs held in memory) ...", flush=True)
    rep_to_newid = {}
    n = 0
    with open(REP_FASTA) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].split()[0].strip()
                rep_to_newid[header] = pad9("US.100NT", n)
                n += 1
                if n % 20_000_000 == 0:
                    print(f"  ...assigned {n:,} ids", flush=True)
    print(f"Total representatives: {n:,}", flush=True)

    print(f"\nPass 2: writing {OUT_FNA} (streaming) ...", flush=True)
    written = 0
    with xz_out(OUT_FNA) as out:
        for header, seq in fasta_iter(REP_FASTA):
            new_id = rep_to_newid[header]
            out.write(f">{new_id} {header}\n{seq}\n".encode("utf-8"))
            written += 1
            if written % 20_000_000 == 0:
                print(f"  ...wrote {written:,} sequences", flush=True)
    print(f"Wrote {written:,} sequences", flush=True)

    print(f"\nPass 3: writing {OUT_MATCHES} (streaming) ...", flush=True)
    n_rows = 0
    n_self = 0
    n_members = 0
    with open(CLUSTER_TSV) as f, xz_out(OUT_MATCHES) as out:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            new_id = rep_to_newid[rep]
            if member == rep:
                out.write(f"{member}\t=\t{new_id}\n".encode("utf-8"))
                n_self += 1
            else:
                out.write(f"{member}\tC\t{new_id}\n".encode("utf-8"))
                n_members += 1
            n_rows += 1
            if n_rows % 20_000_000 == 0:
                print(f"  ...processed {n_rows:,} rows", flush=True)

    print(f"\nTotal rows written: {n_rows:,}", flush=True)
    print(f"  Representative rows ('='): {n_self:,}", flush=True)
    print(f"  Collapsed member rows ('C'): {n_members:,}", flush=True)

    print("\nDONE.")
    print(f"Created: {OUT_FNA}")
    print(f"Created: {OUT_MATCHES}")


if __name__ == "__main__":
    main()
PYEOF
