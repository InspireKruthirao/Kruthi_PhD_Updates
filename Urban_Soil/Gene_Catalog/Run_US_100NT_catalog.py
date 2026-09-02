#!/usr/bin/env python3

import os
import subprocess


# ============================================================
# PATHS
# ============================================================

BASE = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"

MMSEQS_DIR = f"{BASE}/mmseqs_test"

OUTDIR = (
    "/work/microbiome/users/kruthi/intermediate_results/"
    "urban_soil/Gene_Catalog_Urban_Soil"
)


# Successful 100NT MMseqs2 outputs
REP_FASTA = (
    f"{MMSEQS_DIR}/"
    "US.mmseqs.100NT.splitmem_rep_seq.fasta"
)

CLUSTER_TSV = (
    f"{MMSEQS_DIR}/"
    "US.mmseqs.100NT.splitmem_cluster.tsv"
)


# Final catalogue outputs
OUT_FNA = f"{OUTDIR}/US.100NT.fna.xz"
OUT_MATCHES = f"{OUTDIR}/US.100NT.matches.xz"


# Temporary outputs
# Final files are replaced only after successful completion.
OUT_FNA_TMP = OUT_FNA + ".part"
OUT_MATCHES_TMP = OUT_MATCHES + ".part"


# Use the 100 CPUs requested by the PBS job for xz compression
XZ_THREADS = "100"


# Expected catalogue sizes from the successful MMseqs2 run
EXPECTED_REPRESENTATIVES = 256_829_048
EXPECTED_ORFS = 265_392_705


# ============================================================
# ID FORMATTING
# ============================================================

def pad9(prefix, ix):
    """
    Generate catalogue IDs using grouped 9-digit formatting.

    Examples:
        0         -> US.100NT.000_000_000
        1         -> US.100NT.000_000_001
        256829047 -> US.100NT.256_829_047
    """

    number = f"{ix:09d}"

    return (
        f"{prefix}."
        f"{number[:3]}_"
        f"{number[3:6]}_"
        f"{number[6:]}"
    )


# ============================================================
# FASTA READER
# ============================================================

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


# ============================================================
# MULTI-THREADED XZ WRITER
# ============================================================

def xz_writer(path):

    outfile = open(path, "wb")

    proc = subprocess.Popen(
        [
            "xz",
            f"-T{XZ_THREADS}",
            "-c",
        ],
        stdin=subprocess.PIPE,
        stdout=outfile,
    )

    return proc, outfile


# ============================================================
# MAIN
# ============================================================

def main():

    os.makedirs(
        OUTDIR,
        exist_ok=True,
    )


    # Remove incomplete temporary files from a previous failed run
    for path in [
        OUT_FNA_TMP,
        OUT_MATCHES_TMP,
    ]:

        if os.path.exists(path):

            print(
                f"Removing incomplete temporary file: {path}",
                flush=True,
            )

            os.remove(path)


    # ========================================================
    # PASS 1
    # Assign US.100NT IDs
    # ========================================================

    print("=" * 70)
    print("URBAN SOIL 100NT CATALOGUE")
    print("=" * 70)

    print(
        f"\nPass 1: assigning US.100NT IDs from:\n"
        f"{REP_FASTA}\n",
        flush=True,
    )


    rep_to_newid = {}

    n = 0


    with open(REP_FASTA) as f:

        for line in f:

            if line.startswith(">"):

                header = (
                    line[1:]
                    .split()[0]
                    .strip()
                )

                rep_to_newid[header] = pad9(
                    "US.100NT",
                    n,
                )

                n += 1


                if n % 20_000_000 == 0:

                    print(
                        f"  Assigned {n:,} IDs",
                        flush=True,
                    )


    print(
        f"\nTotal representatives: {n:,}",
        flush=True,
    )


    if n != EXPECTED_REPRESENTATIVES:

        raise RuntimeError(
            f"Expected {EXPECTED_REPRESENTATIVES:,} "
            f"representatives but found {n:,}"
        )


    # ========================================================
    # PASS 2
    # Create US.100NT.fna.xz
    # ========================================================

    print()
    print("=" * 70)
    print("PASS 2: CREATE 100NT FASTA")
    print("=" * 70)

    print(
        f"\nCreating:\n{OUT_FNA}\n",
        flush=True,
    )


    written = 0


    xz_proc, xz_file = xz_writer(
        OUT_FNA_TMP
    )


    for header, seq in fasta_iter(
        REP_FASTA
    ):

        new_id = rep_to_newid[header]


        xz_proc.stdin.write(
            (
                f">{new_id} {header}\n"
                f"{seq}\n"
            ).encode("utf-8")
        )


        written += 1


        if written % 20_000_000 == 0:

            print(
                f"  Wrote {written:,} sequences",
                flush=True,
            )


    xz_proc.stdin.close()

    xz_status = xz_proc.wait()

    xz_file.close()


    if xz_status != 0:

        raise RuntimeError(
            "xz failed while creating "
            "US.100NT.fna.xz"
        )


    if written != EXPECTED_REPRESENTATIVES:

        raise RuntimeError(
            f"Expected {EXPECTED_REPRESENTATIVES:,} "
            f"sequences but wrote {written:,}"
        )


    # Replace existing file/symlink only after successful creation
    os.replace(
        OUT_FNA_TMP,
        OUT_FNA,
    )


    print(
        f"\nUS.100NT FASTA sequences written: "
        f"{written:,}",
        flush=True,
    )


    # ========================================================
    # PASS 3
    # Create US.100NT.matches.xz
    # ========================================================

    print()
    print("=" * 70)
    print("PASS 3: CREATE 100NT MAPPING TABLE")
    print("=" * 70)

    print(
        f"\nCreating:\n{OUT_MATCHES}\n",
        flush=True,
    )


    n_rows = 0
    n_self = 0
    n_members = 0


    xz_proc2, xz_file2 = xz_writer(
        OUT_MATCHES_TMP
    )


    with open(CLUSTER_TSV) as f:

        for line in f:

            rep, member = (
                line.rstrip("\n")
                .split("\t")
            )


            new_id = rep_to_newid[rep]


            if member == rep:

                xz_proc2.stdin.write(
                    (
                        f"{member}\t=\t"
                        f"{new_id}\n"
                    ).encode("utf-8")
                )

                n_self += 1


            else:

                xz_proc2.stdin.write(
                    (
                        f"{member}\tC\t"
                        f"{new_id}\n"
                    ).encode("utf-8")
                )

                n_members += 1


            n_rows += 1


            if n_rows % 20_000_000 == 0:

                print(
                    f"  Processed "
                    f"{n_rows:,} mappings",
                    flush=True,
                )


    xz_proc2.stdin.close()

    xz_status2 = xz_proc2.wait()

    xz_file2.close()


    if xz_status2 != 0:

        raise RuntimeError(
            "xz failed while creating "
            "US.100NT.matches.xz"
        )


    if n_rows != EXPECTED_ORFS:

        raise RuntimeError(
            f"Expected {EXPECTED_ORFS:,} mappings "
            f"but found {n_rows:,}"
        )


    if n_self != EXPECTED_REPRESENTATIVES:

        raise RuntimeError(
            f"Expected {EXPECTED_REPRESENTATIVES:,} "
            f"representative mappings but found "
            f"{n_self:,}"
        )


    if n_rows != n_self + n_members:

        raise RuntimeError(
            "Mapping counts are inconsistent."
        )


    # Replace final mapping file only after successful completion
    os.replace(
        OUT_MATCHES_TMP,
        OUT_MATCHES,
    )


    # ========================================================
    # DONE
    # ========================================================

    print()
    print("=" * 70)
    print("100NT CATALOGUE COMPLETE")
    print("=" * 70)

    print(
        f"100NT representatives : "
        f"{n_self:,}"
    )

    print(
        f"Collapsed ORFs        : "
        f"{n_members:,}"
    )

    print(
        f"Total ORF mappings    : "
        f"{n_rows:,}"
    )

    print()

    print(
        f"Created: {OUT_FNA}"
    )

    print(
        f"Created: {OUT_MATCHES}"
    )


if __name__ == "__main__":
    main()
PYEOF
