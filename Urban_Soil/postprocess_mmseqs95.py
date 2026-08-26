#!/usr/bin/env python3

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

REP_FASTA = (
    MMSEQS
    / "US.mmseqs.95NT_rep_seq.fasta"
)

CLUSTER_TSV = (
    MMSEQS
    / "US.mmseqs.95NT_cluster.tsv"
)

OUT_FASTA = (
    BASE
    / "US.95NT.fna.xz"
)

OUT_MATCHES = (
    BASE
    / "US.95NT.matches.tsv.xz"
)


N_100NT = 256_829_048
EXPECTED_95NT = 199_172_868


def get_100nt_index(raw_id):
    raw_id = raw_id.split(None, 1)[0]

    prefix = b"US.100NT."

    if not raw_id.startswith(prefix):
        raise ValueError(
            f"Unexpected 100NT ID: {raw_id!r}"
        )

    number = (
        raw_id[len(prefix):]
        .replace(b"_", b"")
    )

    return int(number)


print("==========================================")
print("Urban Soil 95NT postprocessing")
print("==========================================")
print(
    f"100NT sequences expected : "
    f"{N_100NT:,}"
)
print(
    f"95NT clusters expected   : "
    f"{EXPECTED_95NT:,}"
)


# ==================================================
# Create clean 95NT FASTA
# ==================================================

print()
print("Creating US.95NT.fna.xz ...")


rep_to_95 = (
    array("i", [-1])
    * N_100NT
)

n95 = 0


with open(REP_FASTA, "rb") as infile, \
        xz_out(str(OUT_FASTA)) as outfile:

    for line in infile:

        if line.startswith(b">"):

            header = line[1:].strip()

            rep100 = (
                header.split(None, 1)[0]
            )

            index100 = (
                get_100nt_index(rep100)
            )

            if index100 >= N_100NT:

                raise ValueError(
                    f"100NT index out of range: "
                    f"{rep100!r}"
                )

            if rep_to_95[index100] != -1:

                raise ValueError(
                    f"Duplicate representative: "
                    f"{rep100!r}"
                )

            rep_to_95[index100] = n95

            new_id = (
                pad9(
                    "US.95NT",
                    n95
                )
            )

            outfile.write(
                f">{new_id}\n"
                .encode("ascii")
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
print(
    f"95NT representatives found: "
    f"{n95:,}"
)


if n95 != EXPECTED_95NT:

    raise ValueError(
        f"Expected {EXPECTED_95NT:,} "
        f"representatives but found "
        f"{n95:,}"
    )


# ==================================================
# Create 100NT -> 95NT mapping
# ==================================================

print()
print(
    "Creating US.95NT.matches.tsv.xz ..."
)


rows = 0
representatives = 0
reduced = 0


with open(CLUSTER_TSV, "rb") as infile, \
        xz_out(str(OUT_MATCHES)) as outfile:

    for line in infile:

        line = line.rstrip(b"\n")

        rep_field, member_field = (
            line.split(b"\t", 1)
        )

        rep100 = (
            rep_field.split(None, 1)[0]
        )

        member100 = (
            member_field.split(None, 1)[0]
        )

        rep_index = (
            get_100nt_index(rep100)
        )

        rep95_index = (
            rep_to_95[rep_index]
        )

        if rep95_index == -1:

            raise ValueError(
                "Representative not found "
                f"in rep FASTA: {rep100!r}"
            )

        rep95 = (
            pad9(
                "US.95NT",
                rep95_index
            )
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
                f"{rows:,} / "
                f"{N_100NT:,} "
                "mappings processed",
                flush=True
            )


# ==================================================
# Final checks
# ==================================================

print()
print("==========================================")
print("POSTPROCESSING COMPLETE")
print("==========================================")

print(
    f"Total mappings         : "
    f"{rows:,}"
)

print(
    f"95NT representatives   : "
    f"{representatives:,}"
)

print(
    f"Reduced members        : "
    f"{reduced:,}"
)

print()

print(
    f"Created: {OUT_FASTA}"
)

print(
    f"Created: {OUT_MATCHES}"
)


if rows != N_100NT:

    raise ValueError(
        f"Expected {N_100NT:,} "
        f"mappings but processed "
        f"{rows:,}"
    )


if representatives != EXPECTED_95NT:

    raise ValueError(
        f"Expected "
        f"{EXPECTED_95NT:,} "
        "representative mappings "
        f"but found "
        f"{representatives:,}"
    )


print()
print("All checks passed.")
