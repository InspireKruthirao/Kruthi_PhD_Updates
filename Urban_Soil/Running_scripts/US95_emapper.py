cat > US95_emapper.py <<'PY'
#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import subprocess
from pathlib import Path
from contextlib import contextmanager

import numpy as np


# ============================================================
# Shanghai Dogs utility functions
# ============================================================

sys.path.insert(
    0,
    "/work/microbiome/shanghai_dogs/resource_generation"
)

import lib


# ============================================================
# PATHS
# ============================================================

BASE = Path(
    "/work/microbiome/users/kruthi/"
    "intermediate_results/urban_soil/Prodigal"
)

EGGNOG_BASE = Path(
    "/work/microbiome/users/kruthi/"
    "intermediate_results/egg_nog/results"
)

ORIGIN = BASE / "UrbanSoil.ORF.orig.tsv.xz"
CLUSTERS = BASE / "US.clusters.tsv.xz"

OUT = BASE / "US.95NT.emapper.annotations.gz"
OUT_TMP = BASE / "US.95NT.emapper.annotations.gz.part"

TMP_BEST = BASE / "US.tmp.best_orf_for_95.uint32"
TMP_ORF95 = BASE / "US.tmp.orf_to_95.uint32"

SORT_TMP = BASE / "US_emapper_sort_tmp"


# ============================================================
# EXPECTED CATALOGUE SIZES
# ============================================================

N_ORFS = 265_392_705
N95 = 199_172_868

SENTINEL = np.uint32(0xFFFFFFFF)

THREADS = int(
    os.environ.get("PBS_NCPUS", "16")
)


# ============================================================
# XZ reader
# ============================================================

@contextmanager
def xz_reader(fname):

    p = subprocess.Popen(
        [
            "xz",
            "-dc",
            f"-T{THREADS}",
            str(fname)
        ],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024
    )

    try:
        yield p.stdout

    finally:

        if p.stdout is not None:
            p.stdout.close()

        ret = p.wait()

        if ret != 0:
            raise RuntimeError(
                f"xz failed reading {fname}"
            )


# ============================================================
# Numeric part of catalogue ID
#
# Handles:
# US.ORF.081_385_823
# US.95NT.044_734_275
# ============================================================

def numeric_id(identifier):

    return int(
        identifier.rsplit(".", 1)[1]
        .replace("_", "")
    )


# ============================================================
# Input checks
# ============================================================

for f in [ORIGIN, CLUSTERS]:

    if not f.exists():
        raise FileNotFoundError(
            f"Missing required file: {f}"
        )


emapper_files = sorted(
    glob.glob(
        str(
            EGGNOG_BASE /
            "*" /
            "*.emapper.annotations"
        )
    )
)


print("=" * 75, flush=True)
print("URBAN SOIL 95NT eggNOG ANNOTATION TRANSFER", flush=True)
print("=" * 75, flush=True)

print(
    f"eggNOG files found : {len(emapper_files)}",
    flush=True
)

if len(emapper_files) != 58:

    raise RuntimeError(
        f"Expected 58 eggNOG files, "
        f"found {len(emapper_files)}"
    )


# ============================================================
# Build sample -> eggNOG file mapping
# ============================================================

sample_to_emapper = {}

common_header = None


for fname in emapper_files:

    p = Path(fname)

    sample = p.parent.name

    expected_name = (
        sample + ".emapper.annotations"
    )

    if p.name != expected_name:

        raise RuntimeError(
            f"Unexpected eggNOG file name: {p}"
        )

    if sample in sample_to_emapper:

        raise RuntimeError(
            f"Duplicate eggNOG sample: {sample}"
        )

    sample_to_emapper[sample] = p


    # Find eggNOG table header
    this_header = None

    with open(p, "rt") as f:

        for line in f:

            if line.startswith("#query\t"):

                this_header = line.rstrip("\n")
                break


    if this_header is None:

        raise RuntimeError(
            f"Could not find #query header in {p}"
        )


    if common_header is None:

        common_header = this_header

    elif this_header != common_header:

        raise RuntimeError(
            f"eggNOG headers are not identical: {p}"
        )


print(
    "All 58 eggNOG files have matching headers.",
    flush=True
)


# ============================================================
# Safety
# ============================================================

if OUT.exists():

    raise FileExistsError(
        f"{OUT} already exists."
    )


for f in [
    OUT_TMP,
    TMP_BEST,
    TMP_ORF95
]:

    if f.exists():
        f.unlink()


if SORT_TMP.exists():

    shutil.rmtree(SORT_TMP)


SORT_TMP.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# STEP 1
#
# Equivalent to SHD:
#
# eqs = clusters.query(
#     'rel1 == "=" & rel2 == "="'
# )
#
# Then keep the first ORF for every 95NT.
#
# We cannot use pandas for 265M rows, so:
#
# best_orf[p95] = lowest qualifying ORF
#
# orf_to_95[orf] = corresponding 95NT
# ============================================================

print(flush=True)
print("=" * 75, flush=True)
print(
    "STEP 1: Selecting one exact ORF for every 95NT",
    flush=True
)
print("=" * 75, flush=True)


best_orf = np.memmap(
    TMP_BEST,
    dtype=np.uint32,
    mode="w+",
    shape=(N95,)
)

orf_to_95 = np.memmap(
    TMP_ORF95,
    dtype=np.uint32,
    mode="w+",
    shape=(N_ORFS,)
)


best_orf[:] = SENTINEL
orf_to_95[:] = SENTINEL


n_cluster_rows = 0
n_exact_candidates = 0


with xz_reader(CLUSTERS) as f:

    for line in f:

        if not line.strip():
            continue


        fields = line.split()

        if len(fields) != 5:

            raise ValueError(
                "Unexpected cluster row:\n"
                + line
            )


        (
            orf,
            rel1,
            rep100,
            rel2,
            rep95
        ) = fields


        n_cluster_rows += 1


        if rel1 == "=" and rel2 == "=":

            iorf = numeric_id(orf)
            i95 = numeric_id(rep95)


            if iorf >= N_ORFS:

                raise RuntimeError(
                    f"ORF index outside range: {orf}"
                )


            if i95 >= N95:

                raise RuntimeError(
                    f"95NT index outside range: {rep95}"
                )


            orf_to_95[iorf] = i95


            if (
                best_orf[i95] == SENTINEL
                or
                iorf < int(best_orf[i95])
            ):

                best_orf[i95] = iorf


            n_exact_candidates += 1


        if n_cluster_rows % 10_000_000 == 0:

            print(
                f"Cluster rows read: "
                f"{n_cluster_rows:,}",
                flush=True
            )


print(flush=True)

print(
    f"Total cluster rows : "
    f"{n_cluster_rows:,}",
    flush=True
)

print(
    f"Exact candidates   : "
    f"{n_exact_candidates:,}",
    flush=True
)


if n_cluster_rows != N_ORFS:

    raise RuntimeError(
        f"Expected {N_ORFS:,} cluster rows, "
        f"found {n_cluster_rows:,}"
    )


best_orf.flush()
orf_to_95.flush()


# ============================================================
# Same check as Shanghai Dogs:
#
# Every 95NT must have an ORF satisfying
# rel1 = "=" and rel2 = "="
# ============================================================

represented_95 = int(
    np.count_nonzero(
        best_orf != SENTINEL
    )
)


print(
    f"95NT representatives with exact ORF: "
    f"{represented_95:,}",
    flush=True
)


if represented_95 != N95:

    raise RuntimeError(
        f"Expected {N95:,} 95NT representatives, "
        f"but only {represented_95:,} have "
        f"an exact ORF."
    )


print(
    "All 95NT representatives verified.",
    flush=True
)


# ============================================================
# STEP 2
#
# Prepare sorted gzip output.
#
# Shanghai Dogs sorts by SHD.95NT before writing.
#
# We feed rows through GNU sort and then gzip.
# ============================================================

print(flush=True)
print("=" * 75, flush=True)
print(
    "STEP 2: Mapping original ORFs to eggNOG annotations",
    flush=True
)
print("=" * 75, flush=True)


raw_out = open(
    OUT_TMP,
    "wb"
)


gzip_proc = subprocess.Popen(
    [
        "gzip",
        "-c"
    ],
    stdin=subprocess.PIPE,
    stdout=raw_out
)


sort_env = os.environ.copy()
sort_env["LC_ALL"] = "C"


sort_proc = subprocess.Popen(
    [
        "sort",
        f"--parallel={THREADS}",
        "-S",
        "48G",
        "-T",
        str(SORT_TMP),
        "-t",
        "\t",
        "-k1,1"
    ],
    stdin=subprocess.PIPE,
    stdout=gzip_proc.stdin,
    env=sort_env
)


# Parent no longer needs its own copy of gzip stdin
gzip_proc.stdin.close()


# Write table header.
#
# Similar final structure to SHD:
#
# US.95NT  #query  seed_ortholog ...
#
sort_proc.stdin.write(
    (
        "US.95NT\t"
        + common_header
        + "\n"
    ).encode("utf-8")
)


# ============================================================
# Process one sample at a time
#
# For each sample:
#
#   selected Original_ID -> US.95NT
#
# Then stream that sample's eggNOG file.
#
# This avoids loading all 265M ORFs or all eggNOG
# annotations into RAM.
# ============================================================

def process_sample(
    sample,
    selected
):

    if sample not in sample_to_emapper:

        raise RuntimeError(
            f"No eggNOG annotation file for sample "
            f"{sample}"
        )


    fname = sample_to_emapper[sample]


    print(
        f"Processing {sample}: "
        f"{len(selected):,} selected 95NT ORFs",
        flush=True
    )


    found = 0


    with open(fname, "rt") as f:

        for line in f:

            if not line.strip():
                continue


            # Skip eggNOG comments/header/footer
            if line.startswith("#"):
                continue


            query = line.split("\t", 1)[0]


            p95 = selected.get(query)


            if p95 is None:
                continue


            rep95_name = lib.pad9(
                "US.95NT",
                p95
            )


            sort_proc.stdin.write(
                (
                    rep95_name
                    + "\t"
                    + line
                ).encode("utf-8")
            )


            found += 1


    print(
        f"  annotated selected 95NTs: "
        f"{found:,}",
        flush=True
    )


    return found


# ============================================================
# Scan UrbanSoil.ORF.orig.tsv.xz
# ============================================================

current_sample = None
selected = {}

seen_samples = set()

n_origin = 0
n_selected = 0
n_annotated = 0


with xz_reader(ORIGIN) as f:

    for line in f:

        if not line.strip():
            continue


        fields = line.rstrip("\n").split("\t")


        # Header
        if fields[0] == "ORF":
            continue


        if len(fields) < 3:

            raise RuntimeError(
                "Unexpected ORF origin row:\n"
                + line
            )


        orf = fields[0]
        sample = fields[1]
        original_id = fields[2]


        # --------------------------------------------
        # Samples are expected to occur in contiguous
        # blocks because ORFs were collated sample-wise.
        # --------------------------------------------

        if sample != current_sample:

            if current_sample is not None:

                n_annotated += process_sample(
                    current_sample,
                    selected
                )

                seen_samples.add(
                    current_sample
                )


            if sample in seen_samples:

                raise RuntimeError(
                    f"Sample {sample} occurs in "
                    f"multiple non-contiguous blocks."
                )


            current_sample = sample
            selected = {}


        iorf = numeric_id(orf)


        if iorf >= N_ORFS:

            raise RuntimeError(
                f"ORF outside expected range: "
                f"{orf}"
            )


        p95 = int(
            orf_to_95[iorf]
        )


        # Not an ORF with rel1 = rel2 = "="
        if p95 != int(SENTINEL):

            # Keep only the same "first ORF"
            # selected by the SHD strategy.

            if int(best_orf[p95]) == iorf:

                if original_id in selected:

                    raise RuntimeError(
                        f"Duplicate Original_ID "
                        f"within sample {sample}: "
                        f"{original_id}"
                    )


                selected[original_id] = p95

                n_selected += 1


        n_origin += 1


        if n_origin % 10_000_000 == 0:

            print(
                f"Origin ORFs read: "
                f"{n_origin:,}",
                flush=True
            )


# Process final sample

if current_sample is not None:

    n_annotated += process_sample(
        current_sample,
        selected
    )

    seen_samples.add(
        current_sample
    )


print(flush=True)

print(
    f"Origin ORFs read             : "
    f"{n_origin:,}",
    flush=True
)

print(
    f"Selected canonical 95NT ORFs : "
    f"{n_selected:,}",
    flush=True
)

print(
    f"95NTs with eggNOG annotation : "
    f"{n_annotated:,}",
    flush=True
)

print(
    f"Samples processed             : "
    f"{len(seen_samples)}",
    flush=True
)


# ============================================================
# Verification
# ============================================================

if n_origin != N_ORFS:

    raise RuntimeError(
        f"Expected {N_ORFS:,} ORFs in origin table, "
        f"found {n_origin:,}"
    )


if n_selected != N95:

    raise RuntimeError(
        f"Expected {N95:,} selected canonical ORFs, "
        f"found {n_selected:,}"
    )


if len(seen_samples) != 58:

    raise RuntimeError(
        f"Expected 58 samples, "
        f"processed {len(seen_samples)}"
    )


if seen_samples != set(sample_to_emapper):

    missing_annotations = (
        seen_samples - set(sample_to_emapper)
    )

    missing_origin = (
        set(sample_to_emapper) - seen_samples
    )

    raise RuntimeError(
        "Sample mismatch.\n"
        f"Missing eggNOG: {missing_annotations}\n"
        f"Missing origin samples: {missing_origin}"
    )


# ============================================================
# Finish sort -> gzip pipeline
# ============================================================

sort_proc.stdin.close()


sort_status = sort_proc.wait()


if sort_status != 0:

    raise RuntimeError(
        f"sort failed with exit code "
        f"{sort_status}"
    )


gzip_status = gzip_proc.wait()

raw_out.close()


if gzip_status != 0:

    raise RuntimeError(
        f"gzip failed with exit code "
        f"{gzip_status}"
    )


# ============================================================
# Accept output only after everything succeeded
# ============================================================

os.replace(
    OUT_TMP,
    OUT
)


# ============================================================
# Cleanup
# ============================================================

del best_orf
del orf_to_95


TMP_BEST.unlink(
    missing_ok=True
)

TMP_ORF95.unlink(
    missing_ok=True
)


if SORT_TMP.exists():

    shutil.rmtree(
        SORT_TMP
    )


# ============================================================
# DONE
# ============================================================

print(flush=True)

print("=" * 75, flush=True)
print("DONE", flush=True)
print("=" * 75, flush=True)

print(
    f"Total 95NT representatives : "
    f"{N95:,}",
    flush=True
)

print(
    f"Annotated 95NTs            : "
    f"{n_annotated:,}",
    flush=True
)

print(
    f"Unannotated 95NTs          : "
    f"{N95 - n_annotated:,}",
    flush=True
)

print(
    f"Created                    : "
    f"{OUT}",
    flush=True
)
PY
