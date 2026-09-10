#!/usr/bin/env python3

"""
smORF neighbourhood prediction failure analysis.
"""

from pathlib import Path
from collections import defaultdict, Counter


# Files
IN_FILE = Path(
    "/work/microbiome/users/kruthi/"
    "smorf_onehot_predictions_ALL_no_min_overlap.tsv"
)

OUT_FILE = IN_FILE.parent / "failure_mode_analysis_11june2026.txt"


# COG descriptions
COG_DESC = {
    "C": "Energy production and conversion",
    "D": "Cell cycle control and division",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation and ribosome",
    "K": "Transcription",
    "L": "DNA replication and repair",
    "M": "Cell wall and membrane",
    "N": "Cell motility",
    "O": "Protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis",
    "R": "General function prediction",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
}


# Load TSV
def load_tsv(path):

    rows = []

    with open(path) as f:

        header = f.readline().strip().split("\t")

        for line in f:

            parts = line.strip().split("\t")

            if len(parts) != len(header):
                continue

            rows.append(dict(zip(header, parts)))

    return rows


# Convert data types
def parse_row(row):

    row["n_occurrences"] = int(row["n_occurrences"])
    row["conservation_pct"] = float(row["conservation_pct"])
    row["vote_score_pct"] = float(row["vote_score_pct"])
    row["correct"] = row["correct"] == "True"

    if row["true_cogs"] == "None":
        row["true_cogs"] = None

    return row


# Section heading
def section(title, width=60):

    return f"\n{'=' * width}\n{title}\n{'=' * width}"


# Subsection heading
def subsection(title):

    return f"\n{title}\n{'-' * len(title)}"


# Percentage
def pct(n, total):

    if total == 0:
        return "n/a"

    return f"{n / total * 100:.1f}%"


# Median
def median(values):

    if not values:
        return 0

    values = sorted(values)

    return values[len(values) // 2]


def main():

    print(f"Loading {IN_FILE} ...")

    rows = [
        parse_row(row)
        for row in load_tsv(IN_FILE)
    ]

    validated = [
        row
        for row in rows
        if row["true_cogs"] is not None
    ]

    correct_rows = [
        row
        for row in validated
        if row["correct"]
    ]

    wrong_rows = [
        row
        for row in validated
        if not row["correct"]
    ]

    lines = []
    write = lines.append


    # Summary
    write(section("smORF PREDICTION FAILURE MODE ANALYSIS"))

    write(f"Input file            : {IN_FILE}")
    write(f"Total rows            : {len(rows):>8,}")
    write(f"Validated             : {len(validated):>8,}")

    write(
        f"Correct               : "
        f"{len(correct_rows):>8,}  "
        f"({pct(len(correct_rows), len(validated))})"
    )

    write(
        f"Wrong                 : "
        f"{len(wrong_rows):>8,}  "
        f"({pct(len(wrong_rows), len(validated))})"
    )


    # 1. Confusion
    write(section("1. TOP COG MISPREDICTIONS"))

    confusion = Counter()

    for row in wrong_rows:

        true_cogs = set(
            row["true_cogs"].split(",")
        )

        predicted_cog = row["predicted_cog_cat"]

        for true_cog in true_cogs:

            confusion[
                (true_cog, predicted_cog)
            ] += 1

    write(
        f"\n{'True':>4}  "
        f"{'Pred':>4}  "
        f"{'N':>6}  "
        f"{'% wrong':>8}  "
        f"True function → Predicted function"
    )

    write("-" * 80)

    for (
        true_cog,
        predicted_cog
    ), n in confusion.most_common(20):

        write(
            f"{true_cog:>4}  "
            f"{predicted_cog:>4}  "
            f"{n:>6,}  "
            f"{pct(n, len(wrong_rows)):>8}  "
            f"{COG_DESC.get(true_cog, '?')[:30]:30} → "
            f"{COG_DESC.get(predicted_cog, '?')[:30]}"
        )


    # 2. Confidence tier
    write(section("2. FAILURE BY CONFIDENCE TIER"))

    write(
        f"{'Tier':<12}  "
        f"{'Total':>7}  "
        f"{'Wrong':>7}  "
        f"{'Error rate':>10}  "
        f"{'Median vote%':>13}  "
        f"{'Median N':>10}  "
        f"{'Median cons%':>13}"
    )

    write("-" * 85)

    for tier in [
        "VERY HIGH",
        "HIGH",
        "MEDIUM",
        "LOW",
    ]:

        tier_rows = [
            row
            for row in validated
            if row["confidence_tier"] == tier
        ]

        tier_wrong = [
            row
            for row in tier_rows
            if not row["correct"]
        ]

        if not tier_rows:
            continue

        med_vote = median([
            row["vote_score_pct"]
            for row in tier_wrong
        ])

        med_occ = median([
            row["n_occurrences"]
            for row in tier_wrong
        ])

        med_cons = median([
            row["conservation_pct"]
            for row in tier_wrong
        ])

        write(
            f"{tier:<12}  "
            f"{len(tier_rows):>7,}  "
            f"{len(tier_wrong):>7,}  "
            f"{pct(len(tier_wrong), len(tier_rows)):>10}  "
            f"{med_vote:>13.1f}  "
            f"{med_occ:>10}  "
            f"{med_cons:>13.1f}"
        )


    # 3. Vote score
    write(section("3. VOTE SCORE VS ACCURACY"))

    vote_bins = [
        (0, 30),
        (30, 40),
        (40, 50),
        (50, 60),
        (60, 70),
        (70, 80),
        (80, 90),
        (90, 101),
    ]

    write(
        f"\n{'Vote score':>12}  "
        f"{'Wrong':>8}  "
        f"{'Correct':>8}  "
        f"{'Error rate':>11}"
    )

    write("-" * 50)

    for low, high in vote_bins:

        wrong_subset = [
            row
            for row in wrong_rows
            if low <= row["vote_score_pct"] < high
        ]

        correct_subset = [
            row
            for row in correct_rows
            if low <= row["vote_score_pct"] < high
        ]

        total = (
            len(wrong_subset)
            + len(correct_subset)
        )

        label = f"{low}-{high - 1}%"

        write(
            f"{label:>12}  "
            f"{len(wrong_subset):>8,}  "
            f"{len(correct_subset):>8,}  "
            f"{pct(len(wrong_subset), total):>11}"
        )


    # 4. Occurrence count
    write(section("4. OCCURRENCE COUNT VS ACCURACY"))

    occurrence_bins = [
        (1, 2),
        (2, 5),
        (5, 10),
        (10, 20),
        (20, 50),
        (50, 200),
        (200, 10000),
    ]

    write(
        f"\n{'Occurrences':>14}  "
        f"{'Total':>7}  "
        f"{'Correct':>8}  "
        f"{'Accuracy':>9}"
    )

    write("-" * 50)

    for low, high in occurrence_bins:

        subset = [
            row
            for row in validated
            if low <= row["n_occurrences"] < high
        ]

        if not subset:
            continue

        n_correct = sum(
            1
            for row in subset
            if row["correct"]
        )

        label = f"{low}-{high - 1}"

        write(
            f"{label:>14}  "
            f"{len(subset):>7,}  "
            f"{n_correct:>8,}  "
            f"{pct(n_correct, len(subset)):>9}"
        )


    # 5. Conservation
    write(section("5. NEIGHBOURHOOD CONSERVATION VS ACCURACY"))

    conservation_bins = [
        (0, 25),
        (25, 50),
        (50, 75),
        (75, 90),
        (90, 101),
    ]

    write(
        f"\n{'Conservation':>14}  "
        f"{'Total':>7}  "
        f"{'Correct':>8}  "
        f"{'Accuracy':>9}"
    )

    write("-" * 50)

    for low, high in conservation_bins:

        subset = [
            row
            for row in validated
            if low <= row["conservation_pct"] < high
        ]

        if not subset:
            continue

        n_correct = sum(
            1
            for row in subset
            if row["correct"]
        )

        label = f"{low}-{high - 1}%"

        write(
            f"{label:>14}  "
            f"{len(subset):>7,}  "
            f"{n_correct:>8,}  "
            f"{pct(n_correct, len(subset)):>9}"
        )


    # 6. Multi-COG
    write(section("6. MULTI-COG ANNOTATIONS"))

    single_cog = [
        row
        for row in validated
        if "," not in row["true_cogs"]
    ]

    multi_cog = [
        row
        for row in validated
        if "," in row["true_cogs"]
    ]

    single_correct = sum(
        1
        for row in single_cog
        if row["correct"]
    )

    multi_correct = sum(
        1
        for row in multi_cog
        if row["correct"]
    )

    write(
        f"\nSingle-COG : "
        f"{len(single_cog):>7,}  "
        f"accuracy = "
        f"{pct(single_correct, len(single_cog))}"
    )

    write(
        f"Multi-COG  : "
        f"{len(multi_cog):>7,}  "
        f"accuracy = "
        f"{pct(multi_correct, len(multi_cog))}"
    )

    if multi_cog:

        write(
            subsection(
                "Common multi-COG combinations"
            )
        )

        multi_wrong = [
            row
            for row in multi_cog
            if not row["correct"]
        ]

        combinations = Counter(
            row["true_cogs"]
            for row in multi_wrong
        )

        for combination, n in combinations.most_common(10):

            write(
                f"{combination:10}  "
                f"N={n:>5,}"
            )


    # 7. Accuracy by COG
    write(section("7. ACCURACY BY TRUE COG CATEGORY"))

    cog_stats = defaultdict(
        lambda: {
            "total": 0,
            "correct": 0,
        }
    )

    for row in validated:

        for cog in row["true_cogs"].split(","):

            cog_stats[cog]["total"] += 1

            if row["correct"]:
                cog_stats[cog]["correct"] += 1

    ranked_cogs = sorted(
        cog_stats.items(),
        key=lambda x:
            x[1]["correct"]
            / x[1]["total"],
    )

    write(
        f"\n{'COG':>4}  "
        f"{'Total':>7}  "
        f"{'Correct':>8}  "
        f"{'Accuracy':>9}  "
        f"Function"
    )

    write("-" * 75)

    for cog, stats in ranked_cogs:

        accuracy = (
            stats["correct"]
            / stats["total"]
        )

        write(
            f"{cog:>4}  "
            f"{stats['total']:>7,}  "
            f"{stats['correct']:>8,}  "
            f"{accuracy:>8.1%}  "
            f"{COG_DESC.get(cog, '?')}"
        )


    # 8. Near misses
    write(section("8. NEAR-MISS ANALYSIS"))

    near_misses = []

    for row in wrong_rows:

        top3_string = row.get(
            "top3_candidates",
            "",
        )

        true_cogs = set(
            row["true_cogs"].split(",")
        )

        top3_cogs = set()

        for candidate in top3_string.split("|"):

            if not candidate.strip():
                continue

            cog = (
                candidate
                .split(":")[0]
                .strip()
            )

            top3_cogs.add(cog)

        if true_cogs & top3_cogs:
            near_misses.append(row)

    hard_misses = (
        len(wrong_rows)
        - len(near_misses)
    )

    write(
        f"\nWrong predictions : "
        f"{len(wrong_rows):>7,}"
    )

    write(
        f"Near misses       : "
        f"{len(near_misses):>7,}  "
        f"({pct(len(near_misses), len(wrong_rows))})"
    )

    write(
        f"Hard misses       : "
        f"{hard_misses:>7,}  "
        f"({pct(hard_misses, len(wrong_rows))})"
    )


    # 9. Examples
    write(section("9. EXAMPLE smORFs"))


    def example_rows(subset, n=5):

        examples = []

        for row in subset[:n]:

            examples.append(
                f"{row['smorf_id']:<30}  "
                f"true={row['true_cogs']:<6}  "
                f"pred={row['predicted_cog_cat']}  "
                f"vote={row['vote_score_pct']}%  "
                f"N={row['n_occurrences']}  "
                f"cons={row['conservation_pct']}%  "
                f"tier={row['confidence_tier']}"
            )

        return "\n".join(examples)


    high_confidence_wrong = sorted(
        [
            row
            for row in wrong_rows
            if row["confidence_tier"]
            in ("VERY HIGH", "HIGH")
        ],
        key=lambda row:
            -row["vote_score_pct"],
    )

    write(
        subsection(
            "High-confidence wrong"
        )
    )

    write(
        example_rows(
            high_confidence_wrong
        )
    )


    rare_wrong = sorted(
        [
            row
            for row in wrong_rows
            if row["n_occurrences"] <= 3
        ],
        key=lambda row:
            -row["vote_score_pct"],
    )

    write(
        subsection(
            "Low-occurrence wrong"
        )
    )

    write(
        example_rows(
            rare_wrong
        )
    )


    write(
        subsection(
            "Near misses"
        )
    )

    write(
        example_rows(
            near_misses
        )
    )


    # Failure modes
    write(section("SUMMARY OF FAILURE MODES"))

    write(
        """
1. NEIGHBOURHOOD AMBIGUITY

2. OPERON / GENE CLUSTER BIAS

3. LOW OCCURRENCE COUNT

4. CONSERVATION WITHOUT SPECIFICITY

5. MULTI-COG ANNOTATIONS

6. CATEGORY IMBALANCE
"""
    )


    # Save
    output = "\n".join(lines)

    print(output)

    with OUT_FILE.open("w") as f:
        f.write(output)

    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
