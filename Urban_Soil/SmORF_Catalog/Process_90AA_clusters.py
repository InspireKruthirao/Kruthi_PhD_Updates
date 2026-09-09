import os
import gzip
import re
from collections import defaultdict
import pandas as pd
import importlib.util

# Import helper modules
FASTA_PATH = "/work/microbiome/shanghai_dogs/resource_generation/fasta.py"
LIB_PATH = "/work/microbiome/shanghai_dogs/resource_generation/lib.py"

spec = importlib.util.spec_from_file_location("fasta", FASTA_PATH)
fasta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fasta)

spec = importlib.util.spec_from_file_location("lib", LIB_PATH)
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)

# paths
base_dir = "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil/UrbanSoil_SMORF_resource"
input_fasta = os.path.join(base_dir, "UrbanSoil_100AA_SMORFs_sequences.faa.gz")
origin_tsv = os.path.join(base_dir, "UrbanSoil_100AA_SMORFs_origins.tsv.gz")
output_tsv = os.path.join(base_dir, "UrbanSoil_Clusters.tsv.gz")

cluster_file = os.path.join(base_dir, "UrbanSoil_90AA_SMORFs.clstr")


# Parse CD-HIT cluster file
def parse_cdhit_clusters(cluster_file):
    cluster_mappings = defaultdict(list)
    current_cluster = None

    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">Cluster"):
                current_cluster = line

            else:
                match = re.match(r"\d+\s+\d+aa,\s+>(.*?)\.\.\.\s*(.*)", line)

                if match:
                    seq_id = match.group(1)
                    similarity = match.group(2).strip()

                    if similarity == "*":
                        cluster_mappings[current_cluster].append((seq_id, "*"))

                    else:
                        perc_match = re.search(r"at\s+(\d+\.\d+%)", similarity)

                        if perc_match:
                            percentage = perc_match.group(1)
                            cluster_mappings[current_cluster].append((seq_id, percentage))

    return cluster_mappings


# cluster file
cluster_mappings = parse_cdhit_clusters(cluster_file)

origins = pd.read_csv(origin_tsv, sep="\t")
count100aa = origins["SmORF ID"].value_counts().to_dict()


def smorf_sorting_key(mapping_item):
    # Sort by the number of SmORFs in the cluster and then by the representative
    # 100AA SmORF ID as a tiebreaker
    [rep] = [
        seq_id
        for seq_id, similarity in mapping_item[1]
        if similarity == "*"
    ]

    return (
        -sum(
            count100aa[seq_id]
            for seq_id, _ in mapping_item[1]
        ),
        rep,
    )


clusters = list(cluster_mappings.items())
clusters.sort(key=smorf_sorting_key)

# Map 100AA to 90AA IDs
mappings = []

for idx, (_, seqs) in enumerate(clusters):
    rep_id = lib.pad9("US_SM.90AA", idx)

    for seq_id, similarity in seqs:
        mappings.append((seq_id, rep_id, similarity))

mappings = pd.DataFrame(
    mappings,
    columns=[
        "100AA SmORF ID",
        "90AA SmORF ID",
        "Similarity",
    ],
)

mappings.sort_values(by=["100AA SmORF ID"], inplace=True)

mappings.to_csv(
    output_tsv,
    sep="\t",
    index=False,
    compression="gzip",
)

assert mappings.eval('Similarity == "*"').sum() == len(set(mappings["90AA SmORF ID"]))

mappings = mappings.query('Similarity == "*"')

assert len(set(mappings["90AA SmORF ID"])) == len(set(mappings["100AA SmORF ID"]))

map_90aa_100aa = mappings.set_index(
    "100AA SmORF ID"
)["90AA SmORF ID"].to_dict()

seqs = []

for h, seq in fasta.fasta_iter(input_fasta):
    if h in map_90aa_100aa:
        seqs.append((map_90aa_100aa[h], h, seq))

assert len(seqs) == len(map_90aa_100aa)

seqs.sort()

with gzip.open(
    f"{base_dir}/UrbanSoil_90AA_SMORFs.faa.gz",
    "wt",
) as out_f:

    for h90, h100, seq in seqs:
        out_f.write(f">{h90} {h100}\n{seq}\n")
