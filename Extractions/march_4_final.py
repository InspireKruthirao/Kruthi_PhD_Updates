#!/usr/bin/env python3
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────
BASE = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_21_25_no_min_overlap"]

MIN_SUPPORT = 3
ACCURACY_TARGET = 80.0

RULE_KEYS = ["R1_all3", "R2_both_min3", "R3_down", "R4_up", "R5_either"]
RULE_LABELS = {
    "R1_all3": "All three match (upstream + downstream + smORF)",
    "R2_both_min3": f"Both upstream & downstream same (min {MIN_SUPPORT})",
    "R3_down": "Downstream neighbour matches",
    "R4_up": "Upstream neighbour matches",
    "R5_either": "Either neighbour matches"
}

FEATURE_TYPE = "CDS"
UNKNOWN = "Unknown"

# ─────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────
def parse_attrs(s: str) -> dict:
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def cog_set(name: str):
    if not name or name == UNKNOWN or not str(name).startswith("COG"):
        return frozenset()
    i = str(name).rfind("-")
    if i == -1:
        return frozenset()
    return frozenset(ch for ch in str(name)[i+1:] if ch.isalpha())

def is_annotated(cog) -> bool:
    return bool(cog)

def majority(counter: Counter):
    if not counter:
        return None
    return counter.most_common(1)[0][0]

def read_gff_ordered(gff_file: Path):
    features = []
    with gff_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = parse_attrs(parts[8])
            strand = parts[6]
            start = int(parts[3])
            end = int(parts[4])
            features.append({
                "start": start,
                "end": end,
                "strand": strand,
                "is_target": a.get("target", "0") == "1",
                "cog": cog_set(a.get("Name", UNKNOWN))
            })
    features.sort(key=lambda x: x["start"])
    target_strands = [f["strand"] for f in features if f["is_target"]]
    flipped = bool(target_strands and target_strands[0] == "-")
    if flipped:
        features = list(reversed(features))
        flip_map = {"+": "-", "-": "+"}
        for f in features:
            f["norm_strand"] = flip_map.get(f["strand"], f["strand"])
    else:
        for f in features:
            f["norm_strand"] = f["strand"]
    return features, flipped

def aggregate_smorf_evidence(smorf_dir: Path):
    target_votes, up_votes, down_votes = Counter(), Counter(), Counter()
    n_annotated = n_unannotated = n_total = n_flipped = 0

    for contig in smorf_dir.iterdir():
        if not contig.is_dir():
            continue
        for gff in contig.glob("*.gff"):
            seq, flipped = read_gff_ordered(gff)
            if flipped:
                n_flipped += 1
            for i, f in enumerate(seq):
                if not f["is_target"]:
                    continue
                n_total += 1
                tc = f["cog"]
                left_feat = seq[i-1] if i-1 >= 0 else None
                right_feat = seq[i+1] if i+1 < len(seq) else None
                left_cog = left_feat["cog"] if left_feat else frozenset()
                right_cog = right_feat["cog"] if right_feat else frozenset()

                if is_annotated(tc):
                    n_annotated += 1
                    target_votes[tc] += 1
                else:
                    n_unannotated += 1

                if is_annotated(left_cog):
                    up_votes[left_cog] += 1
                if is_annotated(right_cog):
                    down_votes[right_cog] += 1

    return {
        "target": majority(target_votes) or frozenset(),
        "upstream": majority(up_votes) or frozenset(),
        "downstream": majority(down_votes) or frozenset(),
        "n_annotated": n_annotated,
        "n_unannotated": n_unannotated,
        "n_total": n_total,
        "n_flipped": n_flipped
    }

def match_category(ev, is_ann):
    if not is_ann:
        return None
    tc, up, down = ev["target"], ev["upstream"], ev["downstream"]
    up_match = is_annotated(up) and tc == up
    down_match = is_annotated(down) and tc == down
    if up_match and down_match:
        return "matches_both"
    elif up_match:
        return "matches_up_only"
    elif down_match:
        return "matches_down_only"
    else:
        return "matches_neither"

def evaluate_rules(ev):
    tc, uc, dc = ev["target"], ev["upstream"], ev["downstream"]
    is_ann = ev["n_annotated"] > ev["n_unannotated"]
    r = {}
    r["R1_all3"] = {"eligible": is_annotated(uc) and is_annotated(dc) and uc == dc and is_annotated(tc) and tc == uc,
                     "correct": is_annotated(uc) and is_annotated(dc) and uc == dc and is_annotated(tc) and tc == uc}
    r["R2_both_min3"] = {"eligible": is_annotated(uc) and is_annotated(dc) and uc == dc,
                         "correct": is_annotated(uc) and is_annotated(dc) and uc == dc and is_ann and tc == uc}
    r["R3_down"] = {"eligible": is_annotated(dc),
                    "correct": is_annotated(dc) and is_ann and tc == dc}
    r["R4_up"] = {"eligible": is_annotated(uc),
                  "correct": is_annotated(uc) and is_ann and tc == uc}
    r["R5_either"] = {"eligible": is_annotated(uc) or is_annotated(dc),
                      "correct": (is_annotated(uc) or is_annotated(dc)) and is_ann and (tc == uc or tc == dc)}
    return r, is_ann

def plot_line_cov_acc(counters, total_ann, title):
    covs, accs, lbls = [], [], []
    for k in RULE_KEYS:
        el = counters[k]["eligible"]
        co = counters[k]["correct"]
        cov = el / total_ann * 100 if total_ann else 0
        acc = co / el * 100 if el else 0
        covs.append(cov)
        accs.append(acc)
        lbls.append(RULE_LABELS[k])

    fig, ax1 = plt.subplots(figsize=(10, 5.5))
    ax1.plot(range(len(RULE_KEYS)), covs, marker='o', color="#1f77b4", label="Coverage %", linewidth=2.2)
    ax1.set_ylabel("Coverage (%)", color="#1f77b4")
    ax1.tick_params(axis='y', labelcolor="#1f77b4")
    ax1.set_ylim(0, 110)

    ax2 = ax1.twinx()
    ax2.plot(range(len(RULE_KEYS)), accs, marker='s', color="#ff7f0e", label="Accuracy %", linewidth=2.2)
    ax2.set_ylabel("Accuracy (%)", color="#ff7f0e")
    ax2.tick_params(axis='y', labelcolor="#ff7f0e")
    ax2.set_ylim(0, 110)

    ax1.set_xticks(range(len(RULE_KEYS)))
    ax1.set_xticklabels(lbls, rotation=25, ha="right", fontsize=9)
    ax1.set_title(f"Coverage vs Accuracy — {title}")
    ax1.grid(True, ls="--", alpha=0.35)

    lines1, lbls1 = ax1.get_legend_handles_labels()
    lines2, lbls2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, lbls1 + lbls2, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)

    plt.tight_layout()
    plt.show()

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    for folder in RANGE_DIRS:
        root = BASE / folder
        print(f"{'═'*69}")
        print(f"Rule performance — {folder}")
        print(f"{'═'*69}")

        counters = {k: {"eligible": 0, "correct": 0} for k in RULE_KEYS}
        total_sm = ann_sm = 0
        unann_any = unann_up = unann_down = unann_both_same = unann_both_diff = 0

        for d in root.iterdir():
            if not d.is_dir() or not d.name.startswith("SHD1_SM.100AA"):
                continue
            total_sm += 1
            ev = aggregate_smorf_evidence(d)
            rules, is_ann = evaluate_rules(ev)

            if is_ann:
                ann_sm += 1
                for k in RULE_KEYS:
                    if rules[k]["eligible"]:
                        counters[k]["eligible"] += 1
                        counters[k]["correct"] += int(rules[k]["correct"])
            else:
                has_up = bool(ev["upstream"])
                has_down = bool(ev["downstream"])
                if has_up or has_down:
                    unann_any += 1
                    if has_up: unann_up += 1
                    if has_down: unann_down += 1
                    if has_up and has_down:
                        if ev["upstream"] == ev["downstream"]:
                            unann_both_same += 1
                        else:
                            unann_both_diff += 1

        # Print rule performance table
        print(f"{'Rule':<45} {'Eligible':>9} {'Correct':>9} {'Cov%':>7} {'Acc%':>7}")
        print("─"*69)
        for k in RULE_KEYS:
            el = counters[k]["eligible"]
            co = counters[k]["correct"]
            cov = el / ann_sm * 100 if ann_sm else 0
            acc = co / el * 100 if el else 0
            flag = "✓" if acc >= ACCURACY_TARGET else ""
            print(f"{RULE_LABELS[k]:<45} {el:>9} {co:>9} {cov:>7.1f} {acc:>7.1f}{flag}")

        print(f"{'═'*69}")
        print(f"Total unique smORFs:      {total_sm:,}")
        print(f"Annotated:                {ann_sm:,} ({ann_sm/total_sm*100:.1f}%)")
        print(f"Unannotated:              {total_sm - ann_sm:,}\n")

        print("Unannotated smORFs with annotated neighbours:")
        print(f"  • Has any annotated neighbour:   {unann_any:,}")
        print(f"  • Upstream neighbour annotated:  {unann_up:,}")
        print(f"  • Downstream neighbour annotated: {unann_down:,}")
        print(f"  • Both neighbours same function: {unann_both_same:,}")
        print(f"  • Both annotated but different:  {unann_both_diff:,}\n")

        # Show the line plot
        plot_line_cov_acc(counters, ann_sm, folder)

if __name__ == "__main__":
    main()
