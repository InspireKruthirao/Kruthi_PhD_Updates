#!/usr/bin/env python3
"""
Held-out check and cascade of rules, from the output of sweep_db.py.

The annotated smORFs are split at random into two halves: TRAIN, used to
choose, and TEST, used only to report. Everything chosen on TRAIN is scored
on TEST, so the TEST accuracies are not inflated by picking the best of
many rules.

1. Held-out check: the 15 rules most accurate on TRAIN, with their TRAIN
   and TEST accuracy and coverage.

2. Cascade: every smORF takes the prediction of the first rule in a list
   that predicts something for it. The list is built on TRAIN: at each
   step, the rule most accurate on the smORFs not yet predicted (with at
   least MIN_EVALUATED of them evaluated); then, last, the rule with the
   most correct predictions on what is left, to predict (nearly) all.
   Scored on TEST step by step: how many smORFs each step predicts and how
   accurately. In the cascade, a prediction is correct if it is any of the
   smORF's known categories (none excluded).

Writes cascade_results.tsv and cascade_predictions.tsv (every smORF,
annotated or not: the cascade's prediction, the step and rule that made it,
and that step's TEST accuracy as the expected accuracy).

    python3 cascade.py [folder with predictions_all_smorfs.tsv and combinations_results.csv]
"""

import csv
import random
import sys
from pathlib import Path

SEED = 0
MIN_EVALUATED = 500  # a step needs at least this many TRAIN smORFs evaluated
MAX_STEPS = 10
NO_CALL = {"NO_EVIDENCE", "UNRESOLVED", "LOW_CONFIDENCE", "TOO_FEW_OCCURRENCES", ""}


def excluded(rule_text):
    """The categories a rule excludes, from its description"""
    if "excluding S and R" in rule_text:
        return {"S", "R"}
    if "excluding " in rule_text and rule_text.rstrip().endswith(" only"):
        return {rule_text.rstrip()[:-5].rsplit("excluding ", 1)[1].strip()}
    return set()


def score(rows, rule, ex, smorfs):
    """(evaluated, correct) of rule on smorfs, ignoring the categories in ex"""
    n = ok = 0
    for i in smorfs:
        p = rows[i][rule]
        known = rows[i]["known"] - ex
        if p in NO_CALL or not known:
            continue
        n += 1
        ok += p in known
    return n, ok


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else ""


def main():
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    with open(folder / "combinations_results.csv", newline="") as fh:
        table = list(csv.DictReader(fh))
    rules = {r["ID"]: r["Rule"] for r in table}
    controls = {r["ID"] for r in table if r["Group"].startswith("J.")}  # never used as predictors
    with open(folder / "predictions_all_smorfs.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    ids = [c for c in rows[0] if c in rules and c not in controls]
    for r in rows:
        r["known"] = set(filter(None, r["known_category"].split(",")))
    ex = {c: excluded(rules[c]) for c in ids}

    annotated = [i for i, r in enumerate(rows) if r["known"]]
    random.Random(SEED).shuffle(annotated)
    half = len(annotated) // 2
    train, test = annotated[:half], annotated[half:]
    out = []
    def say(text=""):
        print(text, flush=True)
        out.append(text)
    say(f"annotated smORFs: {len(annotated):,} (TRAIN {len(train):,}, TEST {len(test):,})")

    # 1. held-out check of the best single rules
    say("\n1. Best rules chosen on TRAIN, scored on TEST")
    say("rule\tTRAIN accuracy %\tTEST accuracy %\tTRAIN coverage %\tTEST coverage %\tdescription")
    stats = {}
    for c in ids:
        (n1, k1), (n2, k2) = score(rows, c, ex[c], train), score(rows, c, ex[c], test)
        stats[c] = (n1, k1, n2, k2)
    best = sorted((c for c in ids if stats[c][0] >= MIN_EVALUATED),
                  key=lambda c: stats[c][1] / stats[c][0], reverse=True)[:15]
    for c in best:
        n1, k1, n2, k2 = stats[c]
        say(f"{c}\t{pct(k1, n1)}\t{pct(k2, n2)}\t{pct(n1, len(train))}\t{pct(n2, len(test))}\t{rules[c]}")

    # 2. cascade, built on TRAIN
    steps, left = [], list(train)
    while left and len(steps) < MAX_STEPS - 1:
        sc = {c: score(rows, c, set(), left) for c in ids if c not in steps}
        cand = [(k / n, n, c) for c, (n, k) in sc.items() if n >= MIN_EVALUATED]
        if not cand:
            break
        acc, n, c = max(cand)
        fall = max(sc, key=lambda x: sc[x][1])  # the catch-all: most correct predictions
        if acc <= sc[fall][1] / max(sc[fall][0], 1):
            break  # nothing beats the catch-all any more
        steps.append(c)
        left = [i for i in left if rows[i][c] in NO_CALL]
    if left:
        sc = {c: score(rows, c, set(), left) for c in ids if c not in steps}
        steps.append(max(sc, key=lambda x: sc[x][1]))

    say("\n2. Cascade (built on TRAIN, scored on TEST): each smORF takes the first rule that predicts")
    say("step\trule\tTEST smORFs predicted at this step\tTEST accuracy %\tcumulative TEST coverage %\t"
        "cumulative TEST accuracy %\tdescription")
    left, tot_n, tot_k, step_acc = list(test), 0, 0, {}
    for s, c in enumerate(steps, 1):
        now = [i for i in left if rows[i][c] not in NO_CALL]
        k = sum(rows[i][c] in rows[i]["known"] for i in now)
        tot_n, tot_k = tot_n + len(now), tot_k + k
        step_acc[c] = pct(k, len(now))
        say(f"{s}\t{c}\t{len(now):,}\t{pct(k, len(now))}\t{pct(tot_n, len(test))}\t{pct(tot_k, tot_n)}\t{rules[c]}")
        left = [i for i in left if rows[i][c] in NO_CALL]

    (folder / "cascade_results.tsv").write_text("\n".join(out) + "\n")
    with open(folder / "cascade_predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bin", "smorf_id", "annotated", "known_category", "predicted", "step", "rule",
                    "expected accuracy (%)", "split"])
        split = {i: "TRAIN" for i in train} | {i: "TEST" for i in test}
        for i, r in enumerate(rows):
            for s, c in enumerate(steps, 1):
                if r[c] not in NO_CALL:
                    w.writerow([r["bin"], r["smorf_id"], r["annotated"], r["known_category"], r[c], s, c,
                                step_acc[c], split.get(i, "")])
                    break
            else:
                w.writerow([r["bin"], r["smorf_id"], r["annotated"], r["known_category"], "NO_PREDICTION",
                            "", "", "", split.get(i, "")])
    print(f"\nSaved {folder / 'cascade_results.tsv'} and {folder / 'cascade_predictions.tsv'}")


if __name__ == "__main__":
    main()
