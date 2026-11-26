#!/usr/bin/env python3
"""
STAGE 5 - Turn the raw mmseqs 95NT clustering output into the
final 95NT catalogue.

Input : mmseqs_95nt_64bit/US.mmseqs.95NT_rep_seq.fasta / _cluster.tsv
        (produced by 04_run_mmseqs_95nt_64bit.sh, clustering the
        100NT representatives at 95% identity / 90% coverage)
Output: US.95NT.fna.xz          95NT representative sequences,
                                 IDs renumbered to US.95NT.#########
        US.95NT.matches.tsv.xz  100NT -> 95NT mapping, 3 columns:
                                 100NT_id, relationship ('=' rep /
                                 'R' reduced), 95NT_representative_id

Depends on the shared `pad9`/`xz_out` helpers from the Shanghai Dogs
resource-generation library (kept as an external dependency so both
projects stay in sync -- see lib.py in this directory for a vendored
copy/reference if that path is not available).
"""

import sys
from array import array
from pathlib import Path

# Use exactly the same helper functions as Shanghai Dogs
sys.path.insert(
    0,
    "/work/microbiome/shanghai_dogs/resource_generation"
)

from lib import pad9, xz_out


BASE = Path(
    "/work/microbiome/users/kruthi/"
    "intermediate_results/urban_soil/Prodigal"
)

MMSEQS = BASE / "mmseqs_95nt_64bit"

REP_FASTA = MMSEQS / "US.mmseqs.95NT_rep_seq.fasta"
CLUSTER_TSV = MMSEQS / "US.mmseqs.95NT_cluster.tsv"

OUT_FASTA = BASE / "US.95NT.fna.xz"
OUT_MATCHES = BASE / "US.95NT.matches.tsv.xz"

N_100NT = 256_829_048
EXPECTED_95NT = 199_172_868


def get_100nt_index(raw_id):
    raw_id = raw_id.split(None, 1)[0]

    prefix = b"US.100NT."

    if not raw_id.startswith(prefix):
        raise ValueError(
            f"Unexpected 100NT ID: {raw_id!r}"
        )

    number = raw_id[len(prefix):].replace(b"_", b"")

    return int(number)


print("==========================================")
print("Urban Soil 95NT postprocessing")
print("==========================================")
print(f"100NT sequences expected : {N_100NT:,}")
print(f"95NT clusters expected   : {EXPECTED_95NT:,}")


print()
print("Creating US.95NT.fna.xz ...")

rep_to_95 = array("i", [-1]) * N_100NT

n95 = 0

with open(REP_FASTA, "rb") as infile, \
        xz_out(str(OUT_FASTA)) as outfile:

    for line in infile:

        if line.startswith(b">"):

            header = line[1:].strip()

            rep100 = header.split(None, 1)[0]

            index100 = get_100nt_index(rep100)

            if index100 >= N_100NT:
                raise ValueError(
                    f"100NT index out of range: {rep100!r}"
                )

            if rep_to_95[index100] != -1:
                raise ValueError(
                    f"Duplicate representative: {rep100!r}"
                )

            rep_to_95[index100] = n95

            new_id = pad9("US.95NT", n95)

            outfile.write(
                f">{new_id}\n".encode("ascii")
            )

            n95 += 1

            if n95 % 5_000_000 == 0:
                print(
                    f"{n95:,} representatives processed",
                    flush=True
                )

        else:
            outfile.write(line)


print()
print(f"95NT representatives found: {n95:,}")

if n95 != EXPECTED_95NT:
    raise ValueError(
        f"Expected {EXPECTED_95NT:,} representatives "
        f"but found {n95:,}"
    )


print()
print("Creating US.95NT.matches.tsv.xz ...")

rows = 0
representatives = 0
reduced = 0

with open(CLUSTER_TSV, "rb") as infile, \
        xz_out(str(OUT_MATCHES)) as outfile:

    for line in infile:

        line = line.rstrip(b"\n")

        rep_field, member_field = line.split(b"\t", 1)

        rep100 = rep_field.split(None, 1)[0]
        member100 = member_field.split(None, 1)[0]

        rep_index = get_100nt_index(rep100)

        rep95_index = rep_to_95[rep_index]

        if rep95_index == -1:
            raise ValueError(
                f"Representative not found in rep FASTA: "
                f"{rep100!r}"
            )

        rep95 = pad9(
            "US.95NT",
            rep95_index
        )

        if member100 == rep100:
            relation = "="
            representatives += 1
        else:
            relation = "R"
            reduced += 1

        outfile.write(
            member100
            + b"\t"
            + relation.encode("ascii")
            + b"\t"
            + rep95.encode("ascii")
            + b"\n"
        )

        rows += 1

        if rows % 5_000_000 == 0:
            print(
                f"{rows:,} / {N_100NT:,} mappings processed",
                flush=True
            )


print()
print("==========================================")
print("POSTPROCESSING COMPLETE")
print("==========================================")
print(f"Total mappings         : {rows:,}")
print(f"95NT representatives   : {representatives:,}")
print(f"Reduced members        : {reduced:,}")

print()
print(f"Created: {OUT_FASTA}")
print(f"Created: {OUT_MATCHES}")

if rows != N_100NT:
    raise ValueError(
        f"Expected {N_100NT:,} mappings "
        f"but processed {rows:,}"
    )

if representatives != EXPECTED_95NT:
    raise ValueError(
        f"Expected {EXPECTED_95NT:,} representative mappings "
        f"but found {representatives:,}"
    )

print()
print("All checks passed.")
