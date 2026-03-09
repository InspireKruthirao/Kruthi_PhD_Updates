#!/usr/bin/env python3
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── paths ──────────────────────────────────────────────────────────────────────
BASE       = Path("/work/microbiome/users/kruthi")
RANGE_DIRS = ["SmORF_neighbourhoods_26_30_no_min_overlap"]

# ── tuneable parameters ────────────────────────────────────────────────────────
MIN_SUPPORT_SWEEP = [1, 3, 5, 10, 20, 50]   # values to test for R1 support threshold
MIN_SUPPORT       = 3                         # default used in per-rule table
N_NEIGHBOURS      = 5                         # genes to look at each side
K_OF_N            = 3                         # relaxed conservative: require K/N to agree
ACCURACY_TARGET   = 80.0
# ──────────────────────────────────────────────────────────────────────────────

FEATURE_TYPE = "CDS"
UNKNOWN      = "Unknown"

# Tier labels
TIER_LABELS = {
    "tier1": f"Tier 1 — Conservative ≥{K_OF_N}/{N_NEIGHBOURS} agree (both sides)",
    "tier2": f"Tier 2 — Both neighbours same (support ≥ 10)",
    "tier3": f"Tier 3 — Both neighbours same (support ≥ {MIN_SUPPORT})",
    "tier4":  "Tier 4 — Single neighbour match",
}


# ══════════════════════════════════════════════════════════════════════════════
#  GFF / COG helpers
# ══════════════════════════════════════════════════════════════════════════════

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
    return frozenset(ch for ch in str(name)[i + 1:] if ch.isalpha())


def is_annotated(cog) -> bool:
    return bool(cog)


def majority(counter: Counter):
    if not counter:
        return None
    return counter.most_common(1)[0][0]


def read_gff_ordered(gff_file: Path):
    """
    Parse a GFF, sort by genomic start position, then optionally reverse
    so that 'upstream' always means left-of-target regardless of strand.
    Returns (features, flipped).
    """
    features = []
    with gff_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != FEATURE_TYPE:
                continue
            a = parse_attrs(parts[8])
            features.append({
                "start":     int(parts[3]),
                "end":       int(parts[4]),
                "strand":    parts[6],
                "is_target": a.get("target", "0") == "1",
                "cog":       cog_set(a.get("Name", UNKNOWN)),
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


# ══════════════════════════════════════════════════════════════════════════════
#  Evidence aggregation
# ══════════════════════════════════════════════════════════════════════════════

def aggregate_smorf_evidence(smorf_dir: Path) -> dict:
    """
    Walk every contig/GFF for one smORF cluster and return:
      - immediate neighbour votes + support counts
      - multi-distance neighbour votes (for K-of-N conservative rule)
    """
    target_votes     = Counter()
    up_votes         = Counter()
    down_votes       = Counter()
    multi_up_votes   = defaultdict(Counter)   # dist -> COG -> count
    multi_down_votes = defaultdict(Counter)

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

                left_cog  = seq[i - 1]["cog"] if i - 1 >= 0       else frozenset()
                right_cog = seq[i + 1]["cog"] if i + 1 < len(seq) else frozenset()

                if is_annotated(tc):
                    n_annotated += 1
                    target_votes[tc] += 1
                else:
                    n_unannotated += 1

                if is_annotated(left_cog):
                    up_votes[left_cog] += 1
                if is_annotated(right_cog):
                    down_votes[right_cog] += 1

                # multi-neighbour accumulation
                for dist in range(1, N_NEIGHBOURS + 1):
                    ui, di = i - dist, i + dist
                    if ui >= 0 and is_annotated(seq[ui]["cog"]):
                        multi_up_votes[dist][seq[ui]["cog"]] += 1
                    if di < len(seq) and is_annotated(seq[di]["cog"]):
                        multi_down_votes[dist][seq[di]["cog"]] += 1

    dom_up   = majority(up_votes)   or frozenset()
    dom_down = majority(down_votes) or frozenset()

    multi_up_dominant   = {d: majority(c) for d, c in multi_up_votes.items()}
    multi_down_dominant = {d: majority(c) for d, c in multi_down_votes.items()}

    return {
        "target":              majority(target_votes) or frozenset(),
        "upstream":            dom_up,
        "downstream":          dom_down,
        "up_support":          up_votes.get(dom_up,   0) if up_votes   else 0,
        "down_support":        down_votes.get(dom_down, 0) if down_votes else 0,
        "multi_up_dominant":   multi_up_dominant,
        "multi_down_dominant": multi_down_dominant,
        "multi_up_votes":      multi_up_votes,
        "multi_down_votes":    multi_down_votes,
        "n_annotated":         n_annotated,
        "n_unannotated":       n_unannotated,
        "n_total":             n_total,
        "n_flipped":           n_flipped,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  Rule evaluation helpers
# ══════════════════════════════════════════════════════════════════════════════

def _both_agree_with_support(ev, support_threshold: int) -> bool:
    """True if both immediate neighbours share the same annotated COG
    AND each appears at least support_threshold times."""
    uc, dc = ev["upstream"], ev["downstream"]
    return (
        is_annotated(uc) and is_annotated(dc)
        and uc == dc
        and ev["up_support"]   >= support_threshold
        and ev["down_support"] >= support_threshold
    )


def _k_of_n_consensus(ev, side: str, dominant_cog, k: int, n: int) -> bool:
    """
    Returns True if at least k out of the n positions on `side`
    ('up' or 'down') agree with dominant_cog.
    """
    dom_map = ev[f"multi_{side}_dominant"]
    matches = sum(
        1 for dist in range(1, n + 1)
        if dom_map.get(dist) == dominant_cog
    )
    return matches >= k


def evaluate_rules(ev, support_threshold: int = MIN_SUPPORT) -> tuple:
    """
    Evaluate all rules for one smORF cluster.
    Returns (rule_dict, is_ann, tier).
    `support_threshold` lets the sweep override MIN_SUPPORT for R1.
    """
    tc    = ev["target"]
    uc    = ev["upstream"]
    dc    = ev["downstream"]
    is_ann = ev["n_annotated"] > ev["n_unannotated"]

    r = {}

    # ── R1: both neighbours agree with tunable support ────────────────────
    both_trusted = _both_agree_with_support(ev, support_threshold)
    r["R1_both_minsup"] = {
        "eligible": both_trusted,
        "correct":  both_trusted and is_ann and tc == uc,
    }

    # ── R2: downstream only ───────────────────────────────────────────────
    r["R2_down"] = {
        "eligible": is_annotated(dc),
        "correct":  is_annotated(dc) and is_ann and tc == dc,
    }

    # ── R3: upstream only ─────────────────────────────────────────────────
    r["R3_up"] = {
        "eligible": is_annotated(uc),
        "correct":  is_annotated(uc) and is_ann and tc == uc,
    }

    # ── R4: either neighbour ──────────────────────────────────────────────
    r["R4_either"] = {
        "eligible": is_annotated(uc) or is_annotated(dc),
        "correct":  (
            (is_annotated(uc) or is_annotated(dc))
            and is_ann and (tc == uc or tc == dc)
        ),
    }

    # ── R5: relaxed conservative — K-of-N on BOTH sides, and sides agree ─
    up_ok   = is_annotated(uc) and _k_of_n_consensus(ev, "up",   uc, K_OF_N, N_NEIGHBOURS)
    down_ok = is_annotated(dc) and _k_of_n_consensus(ev, "down", dc, K_OF_N, N_NEIGHBOURS)
    conservative_elig = up_ok and down_ok and uc == dc
    r["R5_conservative"] = {
        "eligible": conservative_elig,
        "correct":  conservative_elig and is_ann and tc == uc,
    }

    # ── Tier assignment (for the tiered-confidence analysis) ──────────────
    tier = _assign_tier(ev, is_ann)

    return r, is_ann, tier


def _assign_tier(ev, is_ann: bool) -> str | None:
    """Assign the highest-confidence tier that applies."""
    uc, dc, tc = ev["upstream"], ev["downstream"], ev["target"]

    # Tier 1 — K-of-N conservative both sides agree
    up_ok   = is_annotated(uc) and _k_of_n_consensus(ev, "up",   uc, K_OF_N, N_NEIGHBOURS)
    down_ok = is_annotated(dc) and _k_of_n_consensus(ev, "down", dc, K_OF_N, N_NEIGHBOURS)
    if up_ok and down_ok and uc == dc:
        return "tier1"

    # Tier 2 — both immediate neighbours same, high support (≥10)
    if _both_agree_with_support(ev, support_threshold=10):
        return "tier2"

    # Tier 3 — both immediate neighbours same, lower support (≥ MIN_SUPPORT)
    if _both_agree_with_support(ev, support_threshold=MIN_SUPPORT):
        return "tier3"

    # Tier 4 — at least one annotated neighbour
    if is_annotated(uc) or is_annotated(dc):
        return "tier4"

    return None


# ══════════════════════════════════════════════════════════════════════════════
#  Frequency-weighted prediction for unannotated smORFs
# ══════════════════════════════════════════════════════════════════════════════

def predict_unannotated(ev) -> dict:
    """
    Score each candidate COG across all N_NEIGHBOURS positions using
    a 1/distance weighting so closer neighbours count more.
    """
    score_map: Counter = Counter()
    for dist in range(1, N_NEIGHBOURS + 1):
        w = 1.0 / dist
        for cog, cnt in ev["multi_up_votes"].get(dist, {}).items():
            score_map[cog] += cnt * w
        for cog, cnt in ev["multi_down_votes"].get(dist, {}).items():
            score_map[cog] += cnt * w

    if not score_map:
        return {"prediction": None, "score": 0.0, "method": "none"}

    best_cog, best_score = score_map.most_common(1)[0]

    uc, dc = ev["upstream"], ev["downstream"]
    if   is_annotated(uc) and is_annotated(dc) and uc == dc: method = "both_agree"
    elif is_annotated(uc) and not is_annotated(dc):           method = "upstream_only"
    elif is_annotated(dc) and not is_annotated(uc):           method = "downstream_only"
    else:                                                      method = "weighted_vote"

    return {"prediction": best_cog, "score": best_score, "method": method}


# ══════════════════════════════════════════════════════════════════════════════
#  Plotting helpers
# ══════════════════════════════════════════════════════════════════════════════

RULE_KEYS_MAIN = ["R1_both_minsup", "R2_down", "R3_up", "R4_either", "R5_conservative"]

def _rule_label(k):
    labels = {
        "R1_both_minsup":  f"Both neighbours same\n(support ≥ {MIN_SUPPORT})",
        "R2_down":          "Downstream\nmatches",
        "R3_up":            "Upstream\nmatches",
        "R4_either":        "Either\nmatches",
        "R5_conservative":  f"Conservative\n≥{K_OF_N}/{N_NEIGHBOURS} agree",
    }
    return labels.get(k, k)


def plot_main_rules(counters, total_ann, title):
    """Coverage vs Accuracy line chart for the five main rules."""
    covs, accs, lbls = [], [], []
    for k in RULE_KEYS_MAIN:
        el  = counters[k]["eligible"]
        co  = counters[k]["correct"]
        covs.append(el / total_ann * 100 if total_ann else 0)
        accs.append(co / el * 100        if el        else 0)
        lbls.append(_rule_label(k))

    fig, ax1 = plt.subplots(figsize=(11, 5.5))
    ax1.plot(range(len(RULE_KEYS_MAIN)), covs, marker="o", color="#1f77b4",
             label="Coverage %", linewidth=2.2)
    ax1.set_ylabel("Coverage (%)", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.set_ylim(0, 110)

    ax2 = ax1.twinx()
    ax2.plot(range(len(RULE_KEYS_MAIN)), accs, marker="s", color="#ff7f0e",
             label="Accuracy %", linewidth=2.2)
    ax2.axhline(ACCURACY_TARGET, color="#ff7f0e", linestyle=":", alpha=0.5,
                label=f"{ACCURACY_TARGET}% target")
    ax2.set_ylabel("Accuracy (%)", color="#ff7f0e")
    ax2.tick_params(axis="y", labelcolor="#ff7f0e")
    ax2.set_ylim(0, 110)

    ax1.set_xticks(range(len(RULE_KEYS_MAIN)))
    ax1.set_xticklabels(lbls, rotation=20, ha="right", fontsize=9)
    ax1.set_title(f"Coverage vs Accuracy — {title}\n"
                  f"(MIN_SUPPORT={MIN_SUPPORT}, K_OF_N={K_OF_N}/{N_NEIGHBOURS})")
    ax1.grid(True, ls="--", alpha=0.35)

    lines1, l1 = ax1.get_legend_handles_labels()
    lines2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, l1 + l2,
               loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
    plt.tight_layout()
    plt.show()


def plot_support_sweep(sweep_results, total_ann, title):
    """
    Two-panel plot:
      Left  — Coverage % vs MIN_SUPPORT threshold for R1
      Right — Accuracy % vs MIN_SUPPORT threshold for R1
    """
    thresholds = sorted(sweep_results.keys())
    covs = [sweep_results[t]["cov"] for t in thresholds]
    accs = [sweep_results[t]["acc"] for t in thresholds]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(thresholds, covs, marker="o", color="#1f77b4", linewidth=2.2)
    ax1.set_xlabel("MIN_SUPPORT threshold")
    ax1.set_ylabel("Coverage (%)")
    ax1.set_title(f"R1 Coverage vs Support threshold\n{title}")
    ax1.set_xticks(thresholds)
    ax1.grid(True, ls="--", alpha=0.35)
    for x, y in zip(thresholds, covs):
        ax1.annotate(f"{y:.1f}%", (x, y), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=8)

    ax2.plot(thresholds, accs, marker="s", color="#ff7f0e", linewidth=2.2)
    ax2.axhline(ACCURACY_TARGET, color="#ff7f0e", linestyle=":", alpha=0.5,
                label=f"{ACCURACY_TARGET}% target")
    ax2.set_xlabel("MIN_SUPPORT threshold")
    ax2.set_ylabel("Accuracy (%)")
    ax2.set_title(f"R1 Accuracy vs Support threshold\n{title}")
    ax2.set_xticks(thresholds)
    ax2.grid(True, ls="--", alpha=0.35)
    ax2.legend()
    for x, y in zip(thresholds, accs):
        ax2.annotate(f"{y:.1f}%", (x, y), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=8)

    plt.tight_layout()
    plt.show()


def plot_tier_breakdown(tier_counts, tier_correct, title):
    """Bar chart showing eligible count and accuracy per confidence tier."""
    tiers  = list(TIER_LABELS.keys())
    counts = [tier_counts.get(t, 0)  for t in tiers]
    accs   = [
        tier_correct.get(t, 0) / tier_counts[t] * 100
        if tier_counts.get(t, 0) else 0
        for t in tiers
    ]
    lbls = [TIER_LABELS[t] for t in tiers]

    fig, ax1 = plt.subplots(figsize=(10, 5))
    x = range(len(tiers))

    bars = ax1.bar(x, counts, color="#4c72b0", alpha=0.75, label="Eligible smORFs")
    ax1.set_ylabel("Eligible smORFs")
    ax1.set_xticks(x)
    ax1.set_xticklabels(lbls, rotation=20, ha="right", fontsize=9)
    ax1.set_title(f"Tiered Confidence — {title}")

    ax2 = ax1.twinx()
    ax2.plot(x, accs, marker="D", color="#dd8452", linewidth=2.2,
             label="Accuracy %")
    ax2.axhline(ACCURACY_TARGET, color="#dd8452", linestyle=":", alpha=0.5)
    ax2.set_ylabel("Accuracy (%)")
    ax2.set_ylim(0, 110)

    lines1, l1 = ax1.get_legend_handles_labels()
    lines2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, l1 + l2, loc="upper right")

    for bar, acc in zip(bars, accs):
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 1, f"{acc:.0f}%",
                 ha="center", va="bottom", fontsize=9, color="#dd8452")

    plt.tight_layout()
    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    for folder in RANGE_DIRS:
        root = BASE / folder

        print(f"\n{'═'*74}")
        print(f"  smORF Rule Evaluator v3  —  {folder}")
        print(f"  MIN_SUPPORT={MIN_SUPPORT}  N_NEIGHBOURS={N_NEIGHBOURS}  K_OF_N={K_OF_N}")
        print(f"{'═'*74}\n")

        # ── collect evidence for every smORF once ─────────────────────────
        all_evidence = []
        for d in root.iterdir():
            if d.is_dir() and d.name.startswith("SHD1_SM.100AA"):
                all_evidence.append(aggregate_smorf_evidence(d))

        total_sm = len(all_evidence)
        print(f"Loaded evidence for {total_sm:,} smORF clusters.\n")

        # ══════════════════════════════════════════════════════════════════
        #  1.  Main rule table (default MIN_SUPPORT)
        # ══════════════════════════════════════════════════════════════════
        counters = {k: {"eligible": 0, "correct": 0} for k in RULE_KEYS_MAIN}
        ann_sm   = 0
        pred_methods: Counter = Counter()
        tier_counts:  Counter = Counter()
        tier_correct: Counter = Counter()

        # unannotated neighbour tracking
        unann_any = unann_up = unann_down = unann_both_same = unann_both_diff = 0
        operon_boundary: list = []   # smORFs where neighbours disagree

        for ev in all_evidence:
            rules, is_ann, tier = evaluate_rules(ev)

            if is_ann:
                ann_sm += 1
                for k in RULE_KEYS_MAIN:
                    if rules[k]["eligible"]:
                        counters[k]["eligible"] += 1
                        counters[k]["correct"]  += int(rules[k]["correct"])

                if tier:
                    tier_counts[tier]  += 1
                    if rules.get("R1_both_minsup", {}).get("correct") or \
                       rules.get("R5_conservative", {}).get("correct"):
                        tier_correct[tier] += 1
                    # re-evaluate correctness per tier
                    tc, uc = ev["target"], ev["upstream"]
                    if tier in ("tier1", "tier2", "tier3") and tc == uc:
                        tier_correct[tier] += 0   # already counted above
            else:
                has_up   = bool(ev["upstream"])
                has_down = bool(ev["downstream"])
                if has_up or has_down:
                    unann_any += 1
                    if has_up:   unann_up   += 1
                    if has_down: unann_down += 1
                    if has_up and has_down:
                        if ev["upstream"] == ev["downstream"]:
                            unann_both_same += 1
                        else:
                            unann_both_diff += 1
                            operon_boundary.append(ev)

                pred = predict_unannotated(ev)
                pred_methods[pred["method"]] += 1

        # ── recompute tier accuracy properly ──────────────────────────────
        tier_counts_r  = Counter()
        tier_correct_r = Counter()
        for ev in all_evidence:
            rules, is_ann, tier = evaluate_rules(ev)
            if not is_ann or tier is None:
                continue
            tc, uc, dc = ev["target"], ev["upstream"], ev["downstream"]
            tier_counts_r[tier] += 1
            if tier in ("tier1", "tier2", "tier3") and tc == uc:
                tier_correct_r[tier] += 1
            elif tier == "tier4" and (tc == uc or tc == dc):
                tier_correct_r[tier] += 1

        # ── print rule table ──────────────────────────────────────────────
        rule_labels_short = {
            "R1_both_minsup":  f"Both neighbours same (support ≥ {MIN_SUPPORT})",
            "R2_down":          "Downstream neighbour matches",
            "R3_up":            "Upstream neighbour matches",
            "R4_either":        "Either neighbour matches",
            "R5_conservative":  f"Conservative: ≥{K_OF_N}/{N_NEIGHBOURS} positions agree",
        }

        print(f"{'Rule':<50} {'Eligible':>9} {'Correct':>9} {'Cov%':>7} {'Acc%':>7}")
        print("─" * 74)
        for k in RULE_KEYS_MAIN:
            el  = counters[k]["eligible"]
            co  = counters[k]["correct"]
            cov = el / ann_sm * 100 if ann_sm else 0
            acc = co / el * 100     if el     else 0
            flag = " ✓" if acc >= ACCURACY_TARGET else ""
            print(f"{rule_labels_short[k]:<50} {el:>9} {co:>9} {cov:>7.1f} {acc:>7.1f}{flag}")

        print(f"\n{'─'*74}")
        print(f"Total smORFs:   {total_sm:,}")
        print(f"Annotated:      {ann_sm:,}  ({ann_sm/total_sm*100:.1f}%)")
        print(f"Unannotated:    {total_sm - ann_sm:,}")

        # ── tiered breakdown ──────────────────────────────────────────────
        print(f"\n{'═'*74}")
        print("  Tiered Confidence Breakdown")
        print(f"{'═'*74}")
        print(f"{'Tier':<55} {'Count':>7} {'Correct':>9} {'Acc%':>7}")
        print("─" * 74)
        for t in TIER_LABELS:
            cnt = tier_counts_r[t]
            cor = tier_correct_r[t]
            acc = cor / cnt * 100 if cnt else 0
            flag = " ✓" if acc >= ACCURACY_TARGET else ""
            print(f"{TIER_LABELS[t]:<55} {cnt:>7} {cor:>9} {acc:>7.1f}{flag}")

        # ── unannotated / operon-boundary summary ─────────────────────────
        print(f"\n{'═'*74}")
        print("  Unannotated smORFs")
        print(f"{'═'*74}")
        print(f"  Has any annotated neighbour:    {unann_any:,}")
        print(f"  Upstream annotated:             {unann_up:,}")
        print(f"  Downstream annotated:           {unann_down:,}")
        print(f"  Both same function:             {unann_both_same:,}")
        print(f"  Both annotated but DIFFERENT    {unann_both_diff:,}  ← operon-boundary candidates")

        print(f"\nFrequency-weighted predictions for unannotated smORFs:")
        for method, cnt in pred_methods.most_common():
            print(f"  • {method:<25}: {cnt:,}")

        # ══════════════════════════════════════════════════════════════════
        #  2.  MIN_SUPPORT sweep for R1
        # ══════════════════════════════════════════════════════════════════
        print(f"\n{'═'*74}")
        print(f"  MIN_SUPPORT sweep for R1 (Both neighbours same + support threshold)")
        print(f"{'Threshold':>12} {'Eligible':>10} {'Correct':>10} {'Cov%':>8} {'Acc%':>8}")
        print("─" * 74)

        sweep_results = {}
        for sup in MIN_SUPPORT_SWEEP:
            el = co = 0
            for ev in all_evidence:
                rules, is_ann, _ = evaluate_rules(ev, support_threshold=sup)
                if is_ann and rules["R1_both_minsup"]["eligible"]:
                    el += 1
                    co += int(rules["R1_both_minsup"]["correct"])
            cov = el / ann_sm * 100 if ann_sm else 0
            acc = co / el * 100     if el     else 0
            flag = " ✓" if acc >= ACCURACY_TARGET else ""
            print(f"{sup:>12} {el:>10} {co:>10} {cov:>8.1f} {acc:>8.1f}{flag}")
            sweep_results[sup] = {"cov": cov, "acc": acc, "eligible": el, "correct": co}

        # ══════════════════════════════════════════════════════════════════
        #  3.  Plots
        # ══════════════════════════════════════════════════════════════════
        plot_main_rules(counters, ann_sm, folder)
        plot_support_sweep(sweep_results, ann_sm, folder)
        plot_tier_breakdown(tier_counts_r, tier_correct_r, folder)


if __name__ == "__main__":
    main()
