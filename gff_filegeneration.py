import pandas as pd
from pathlib import Path

BASE_DIR = "/work/microbiome/users/kruthi/SmORF_neighbourhoods_26_30"

OVERWRITE = True 

missing_neighbors = []
missing_eggnog = []

def create_gff_for_contig(contig_dir):
    contig_path = Path(contig_dir)

    contig_name = contig_path.name.split('_', 1)[1]
    gff_file = contig_path / f"{contig_name}.gff"

    if gff_file.exists() and not OVERWRITE:
        print(f"  ⏩ Skipping {contig_name} (GFF already exists)")
        return

    neighbors_file = contig_path / "neighbors.tsv"
    eggnog_file = contig_path / "EggNOG_annotations.tsv"

    if not neighbors_file.exists():
        missing_neighbors.append(str(contig_path))
        print(f"   Missing neighbors.tsv → SKIPPING {contig_path.name}")
        return

    if not eggnog_file.exists():
        missing_eggnog.append(str(contig_path))
        print(f"   Missing EggNOG_annotations.tsv → SKIPPING {contig_path.name}")
        return
    neighbors = pd.read_csv(neighbors_file, sep="\t")

    target_contig = None
    if "Note" in neighbors.columns:
        target_row = neighbors[neighbors["Note"].str.contains("TARGET", na=False)]
        if not target_row.empty:
            target_contig = target_row.iloc[0]["Contig"]

    with open(eggnog_file, "r") as f:
        header = f.readline().lstrip("#").strip().split("\t")

    eggnog = pd.read_csv(eggnog_file, sep="\t", skiprows=1, names=header)

    annotation_map = {}
    for _, row in eggnog.iterrows():
        cog_name = "Unknown"
        if pd.notna(row.get("eggNOG_OGs")) and "COG" in str(row.get("eggNOG_OGs")):
            cog_raw = str(row["eggNOG_OGs"]).split("@")[0].split(",")[0]
            if cog_raw.startswith("COG"):
                cat = str(row["COG_category"]).strip() if pd.notna(row.get("COG_category")) else ""
                cog_name = f"{cog_raw}-{cat}" if cat and cat != "-" else cog_raw
        annotation_map[row["Contig"]] = cog_name

    # Write GFF (overwrite if exists)
    with open(gff_file, "w") as gff:
        gff.write("##gff-version 3\n")

        for _, row in neighbors.iterrows():
            contig_id = row["Contig"]
            name = annotation_map.get(contig_id, "Unknown")
            score = "1" if name != "Unknown" else "0.0"
            strand = "+" if row["Strand"] == 1 else "-"

            attributes = [f"ID={contig_id}", f"Name={name}"]

            if contig_id == target_contig:
                attributes.extend(["target=1", "Note=TARGET_sMORF", "colour=1"])
            else:
                attributes.append("target=0")

            # IMPORTANT: Use ASCII ';' not '；'
            gff.write(
                f"{contig_id}\tProdigal_v2.6.3\tCDS\t{row['Start']}\t{row['End']}\t"
                f"{score}\t{strand}\t0\t{';'.join(attributes)}\n"
            )

    print(f"  ✔ Wrote {gff_file.name}")


def process_all_smorfs():
    base_path = Path(BASE_DIR)

    print("=== Starting GFF generation ===")
    print(f"OVERWRITE={OVERWRITE}")

    for smorf_dir in sorted(base_path.glob("SHD1_SM.100AA.*")):
        print(f"\nProcessing {smorf_dir.name}...")
        for contig_dir in sorted(smorf_dir.glob("D*_contig_*")):
            create_gff_for_contig(contig_dir)

    print("\n= SUMMARY =")

    print(f"\nMissing neighbors.tsv: {len(missing_neighbors)}")
    for m in missing_neighbors:
        print("  -", m)

    print(f"\nMissing EggNOG_annotations.tsv: {len(missing_eggnog)}")
    for m in missing_eggnog:
        print("  -", m)

    print("\nDone.")

if __name__ == "__main__":
    process_all_smorfs()

