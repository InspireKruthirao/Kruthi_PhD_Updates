import sys
import importlib.util
import os
import csv
import gzip
from collections import defaultdict

FASTA_PATH = "/work/microbiome/shanghai_dogs/resource_generation/fasta.py"
LIB_PATH   = "/work/microbiome/shanghai_dogs/resource_generation/lib.py"

spec = importlib.util.spec_from_file_location("fasta", FASTA_PATH)
fasta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fasta)

spec = importlib.util.spec_from_file_location("lib", LIB_PATH)
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)

print("Using:")
print(fasta.__file__)
print(lib.__file__)
print(lib.pad6("US_SM.100AA", 0))

base_dir     = "/work/microbiome/users/kruthi/GMSC_MAPPER_urbansoil"
resource_dir = os.path.join(base_dir, "UrbanSoil_SMORF_resource")
os.makedirs(resource_dir, exist_ok=True)

output_fasta          = os.path.join(resource_dir, "UrbanSoil_100AA_SMORFs_sequences.faa.gz")
sequence_metadata_tsv = os.path.join(resource_dir, "UrbanSoil_100AA_SMORFs_origins.tsv.gz")
habitat_taxonomy_tsv  = os.path.join(resource_dir, "UrbanSoil_100AA_SMORFs_habitat_taxonomy.tsv")
log_file              = os.path.join(resource_dir, "merge_log.txt")

sequence_dict = defaultdict(list)

for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_path):
        continue
    if not sample.endswith("_medaka_polypolish"):
        continue
    fasta_file = os.path.join(sample_path, "mapped.smorfs.faa")
    if os.path.exists(fasta_file):
        for header, seq in fasta.fasta_iter(fasta_file, full_header=True):
            clean_seq = seq.replace(" ", "").replace("\n", "").upper()
            sequence_dict[clean_seq].append(header)

sorted_sequences = sorted(sequence_dict.items(), key=lambda x: (-len(x[1]), x[0]))
seq_to_id = {}

with gzip.open(output_fasta, "wt") as out_f:
    for idx, (seq, headers) in enumerate(sorted_sequences):
        new_id = lib.pad6("US_SM.100AA", idx)
        out_f.write(f">{new_id}\n{seq}\n")
        seq_to_id[seq] = new_id

print(f"FASTA file created: {output_fasta}")

seq_metadata      = defaultdict(list)
seq_counts        = defaultdict(int)
smorf_to_habitat  = defaultdict(set)
smorf_to_taxonomy = defaultdict(str)

for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_path):
        continue
    if not sample.endswith("_medaka_polypolish"):
        continue

    sample_id = sample.replace("_medaka_polypolish", "")

    pred_file = os.path.join(sample_path, "predicted.filtered.smorf.faa")
    if os.path.exists(pred_file):
        for header, seq in fasta.fasta_iter(pred_file, full_header=True):
            clean_seq = seq.replace(" ", "").replace("\n", "").upper()
            seq_counts[clean_seq] += 1
            if clean_seq not in seq_to_id:
                continue
            parts = header.split("#")
            if len(parts) < 5:
                continue
            smorf_id_raw = parts[0].strip()
            contig       = parts[1].strip()
            start        = parts[2].strip()
            end          = parts[3].strip()
            strand       = "+" if parts[4].strip() == "1" else "-"
            coord        = f"{start}-{end}"
            seq_metadata[clean_seq].append({
                "sample_id": sample_id,
                "contig":    contig,
                "coords":    coord,
                "strand":    strand,
                "smorf_id":  smorf_id_raw
            })

    habitat_file = os.path.join(sample_path, "habitat.out.smorfs.tsv")
    if os.path.exists(habitat_file):
        with open(habitat_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                smorf_id_raw = row[0].strip()
                habitats = [h.strip() for h in row[1].split(",") if h.strip()]
                smorf_to_habitat[(smorf_id_raw, sample_id)] = set(habitats)

    taxonomy_file = os.path.join(sample_path, "taxonomy.out.smorfs.tsv")
    if os.path.exists(taxonomy_file):
        with open(taxonomy_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                smorf_id_raw = row[0].strip()
                taxonomy     = row[1].strip()
                smorf_to_taxonomy[(smorf_id_raw, sample_id)] = taxonomy

smorf_order = sorted(seq_to_id.items(), key=lambda x: (-seq_counts[x[0]], x[0]))

with gzip.open(sequence_metadata_tsv, "wt", newline='') as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["SmORF ID", "Sample ID", "Contig", "Coordinates", "Strand"])
    for seq, us_id in smorf_order:
        for entry in seq_metadata[seq]:
            clean_contig = entry["contig"].rsplit("_", maxsplit=1)[0]
            writer.writerow([
                us_id,
                entry["sample_id"],
                clean_contig,
                entry["coords"],
                entry["strand"]
            ])

with open(habitat_taxonomy_tsv, "w", newline='') as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["SmORF ID", "Sample ID", "Habitat", "Taxonomy"])
    for seq, us_id in smorf_order:
        for entry in seq_metadata[seq]:
            sample_id    = entry["sample_id"]
            smorf_id_raw = entry["smorf_id"]
            habitat  = ",".join(sorted(smorf_to_habitat.get((smorf_id_raw, sample_id), set())))
            taxonomy = smorf_to_taxonomy.get((smorf_id_raw, sample_id), "")
            writer.writerow([us_id, sample_id, habitat, taxonomy])

with open(log_file, "w") as log_f:
    log_f.write(f"Total unique sequences: {len(seq_to_id)}\n")
    duplicate_count = sum(1 for count in seq_counts.values() if count > 1)
    log_f.write(f"Number of sequences with duplicates: {duplicate_count}\n")
    log_f.write("Detailed duplicate counts with associated smORFs:\n")
    for seq, us_id in smorf_order:
        count = seq_counts[seq]
        if count > 1:
            smorfs = ", ".join(sorted(set(entry["smorf_id"] for entry in seq_metadata[seq])))
            log_f.write(f"Sequence with ID {us_id} appeared {count} times (smORFs: {smorfs})\n")

print(f"FASTA created:          {output_fasta}")
print(f"Origins TSV:            {sequence_metadata_tsv}")
print(f"Habitat & Taxonomy TSV: {habitat_taxonomy_tsv}")
print(f"Log:                    {log_file}")
