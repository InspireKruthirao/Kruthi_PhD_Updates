#!/usr/bin/env python3

import os
import sys
import subprocess
from contextlib import contextmanager
from pathlib import Path

import numpy as np


# ============================================================
# USE SHANGHAI DOGS lib.py
# Contains:
#   lib.pad9()
#   lib.xz_out()
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

M100 = BASE / "US.100NT_v3.matches.xz"
M95  = BASE / "US.95NT.matches.tsv.xz"

OUT = BASE / "US.clusters.tsv.xz"
OUT_TMP = BASE / "US.clusters.tsv.xz.part"

TMP_REP95 = BASE / "US.tmp.rep95.uint32"
TMP_REL95 = BASE / "US.tmp.rel95.uint8"


# ============================================================
# EXPECTED COUNTS
# ============================================================

N100 = 256_829_048
EXPECTED_ORFS = 265_392_705

THREADS = int(
    os.environ.get("PBS_NCPUS", "16")
)


# ============================================================
# XZ READER
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
                f"xz failed while reading {fname}"
            )


# ============================================================
# Convert catalogue IDs to integer
#
# Handles BOTH:
#
# US.100NT.199144448
# US.100NT.199_144_448
#
# and:
#
# US.95NT.160_974_800
# ============================================================

def numeric_id(identifier):

    number = identifier.rsplit(".", 1)[1]

    return int(
        number.replace("_", "")
    )


# ============================================================
# CHECK INPUTS
# ============================================================

for fname in [M100, M95]:

    if not fname.exists():

        raise FileNotFoundError(
            f"Missing required file: {fname}"
        )


if OUT.exists():

    raise FileExistsError(
        f"{OUT} already exists. "
        "Move/remove it before running."
    )


# Remove remnants from an interrupted previous run

for fname in [
    OUT_TMP,
    TMP_REP95,
    TMP_REL95
]:

    if fname.exists():
        fname.unlink()


print("=" * 70, flush=True)
print("URBAN SOIL CLUSTER TABLE", flush=True)
print("=" * 70, flush=True)

print(f"100NT mappings : {M100}", flush=True)
print(f"95NT mappings  : {M95}", flush=True)
print(f"Output         : {OUT}", flush=True)

print(flush=True)


# ============================================================
# STEP 1
#
# Build:
#
# US.100NT -> US.95NT
#
# There are ~257 million 100NT sequences.
# We therefore use disk-backed numpy arrays instead of a
# gigantic pandas/Python dictionary.
# ============================================================

print("=" * 70, flush=True)
print(
    "STEP 1: Loading 100NT -> 95NT mappings",
    flush=True
)
print("=" * 70, flush=True)


rep95 = np.memmap(
    TMP_REP95,
    dtype=np.uint32,
    mode="w+",
    shape=(N100,)
)

rel95 = np.memmap(
    TMP_REL95,
    dtype=np.uint8,
    mode="w+",
    shape=(N100,)
)

# 0 means that no mapping has been loaded yet
rel95[:] = 0


n95 = 0


with xz_reader(M95) as f:

    for line in f:

        if not line.strip():
            continue

        fields = line.split()

        if len(fields) != 3:

            raise ValueError(
                "Unexpected 95NT row:\n"
                + line
            )


        rep100_raw, relationship, rep95_raw = fields


        if not rep100_raw.startswith("US.100NT."):

            raise ValueError(
                f"Unexpected 100NT ID: {rep100_raw}"
            )


        if not rep95_raw.startswith("US.95NT."):

            raise ValueError(
                f"Unexpected 95NT ID: {rep95_raw}"
            )


        i100 = numeric_id(rep100_raw)
        i95 = numeric_id(rep95_raw)


        if i100 >= N100:

            raise ValueError(
                f"100NT index outside expected range: "
                f"{rep100_raw}"
            )


        if rel95[i100] != 0:

            raise ValueError(
                f"Duplicate 95NT mapping for "
                f"{rep100_raw}"
            )


        rep95[i100] = i95


        if relationship == "=":

            rel95[i100] = 1

        elif relationship == "R":

            rel95[i100] = 2

        else:

            raise ValueError(
                f"Unexpected 95NT relationship: "
                f"{relationship}"
            )


        n95 += 1


        if n95 % 10_000_000 == 0:

            print(
                f"95NT mappings loaded: "
                f"{n95:,}",
                flush=True
            )


rep95.flush()
rel95.flush()


print(flush=True)

print(
    f"Total 95NT mappings loaded: {n95:,}",
    flush=True
)


if n95 != N100:

    raise RuntimeError(
        f"Expected {N100:,} 100NT -> 95NT "
        f"mappings but found {n95:,}"
    )


print(
    "100NT -> 95NT mapping verified.",
    flush=True
)


# ============================================================
# STEP 2
#
# Build:
#
# ORF -> 100NT -> 95NT
#
# US.100NT_v3.matches.xz ALREADY contains every original ORF,
# including the ORFs that are themselves 100NT representatives.
#
# Therefore we read ONLY this file.
# ============================================================

print(flush=True)

print("=" * 70, flush=True)
print(
    "STEP 2: Writing US.clusters.tsv.xz",
    flush=True
)
print("=" * 70, flush=True)


n_orfs = 0


with lib.xz_out(str(OUT_TMP)) as out:

    with xz_reader(M100) as f:

        for line in f:

            if not line.strip():
                continue


            fields = line.split()


            if len(fields) != 3:

                raise ValueError(
                    "Unexpected 100NT row:\n"
                    + line
                )


            orf, rel100, rep100_raw = fields


            if not orf.startswith("US.ORF."):

                raise ValueError(
                    f"Unexpected ORF ID: {orf}"
                )


            if not rep100_raw.startswith(
                "US.100NT."
            ):

                raise ValueError(
                    f"Unexpected 100NT ID: "
                    f"{rep100_raw}"
                )


            if rel100 not in ("=", "C"):

                raise ValueError(
                    f"Unexpected 100NT relationship: "
                    f"{rel100}"
                )


            # -----------------------------------------------
            # Convert current 100NT ID to numeric ID
            # -----------------------------------------------

            i100 = numeric_id(rep100_raw)


            if i100 >= N100:

                raise ValueError(
                    f"100NT index outside range: "
                    f"{rep100_raw}"
                )


            # -----------------------------------------------
            # Get its 95NT relationship
            # -----------------------------------------------

            code = int(
                rel95[i100]
            )


            if code == 1:

                rel95_char = "="

            elif code == 2:

                rel95_char = "R"

            else:

                raise RuntimeError(
                    f"No 95NT mapping found for "
                    f"{rep100_raw}"
                )


            # -----------------------------------------------
            # IMPORTANT:
            #
            # Use Shanghai Dogs pad9() formatting
            # for BOTH catalogue levels.
            #
            # 199144448
            # becomes:
            #
            # US.100NT.199_144_448
            #
            # and:
            #
            # US.95NT.160_974_800
            # -----------------------------------------------

            rep100_name = lib.pad9(
                "US.100NT",
                i100
            )


            rep95_name = lib.pad9(
                "US.95NT",
                int(rep95[i100])
            )


            # -----------------------------------------------
            # Same 5-column structure as Shanghai Dogs:
            #
            # ORF
            # relationship to 100NT
            # 100NT representative
            # relationship to 95NT
            # 95NT representative
            # -----------------------------------------------

            out.write(
                (
                    f"{orf}\t"
                    f"{rel100}\t"
                    f"{rep100_name}\t"
                    f"{rel95_char}\t"
                    f"{rep95_name}\n"
                ).encode("ascii")
            )


            n_orfs += 1


            if n_orfs % 10_000_000 == 0:

                print(
                    f"ORFs written: "
                    f"{n_orfs:,}",
                    flush=True
                )


# ============================================================
# VERIFY FINAL ROW COUNT
# ============================================================

print(flush=True)

print(
    f"Total cluster rows written: "
    f"{n_orfs:,}",
    flush=True
)


if n_orfs != EXPECTED_ORFS:

    raise RuntimeError(
        f"Expected {EXPECTED_ORFS:,} ORFs, "
        f"but wrote {n_orfs:,}. "
        "Output will NOT be accepted."
    )


# ============================================================
# Rename .part only after successful verification
# ============================================================

os.replace(
    OUT_TMP,
    OUT
)


# ============================================================
# CLEANUP TEMPORARY LOOKUP FILES
# ============================================================

del rep95
del rel95


TMP_REP95.unlink(
    missing_ok=True
)

TMP_REL95.unlink(
    missing_ok=True
)


# ============================================================
# DONE
# ============================================================

print(flush=True)

print("=" * 70, flush=True)
print("DONE", flush=True)
print("=" * 70, flush=True)

print(
    f"Correct cluster rows : "
    f"{n_orfs:,}",
    flush=True
)

print(
    f"Created              : "
    f"{OUT}",
    flush=True
)
