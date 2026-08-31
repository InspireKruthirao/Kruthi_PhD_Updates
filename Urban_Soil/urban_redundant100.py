import importlib.util
from collections import defaultdict


# Paths to helper modules from the Shanghai Dogs workflow
FASTA_PATH = "/work/microbiome/shanghai_dogs/resource_generation/fasta.py"
LIB_PATH = "/work/microbiome/shanghai_dogs/resource_generation/lib.py"


# Load fasta.py
spec = importlib.util.spec_from_file_location("fasta", FASTA_PATH)
fasta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fasta)


# Load lib.py
spec = importlib.util.spec_from_file_location("lib", LIB_PATH)
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)


# Urban Soil Prodigal directory
BASE = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"


# Input ORF catalogue
ORF_FNA = f"{BASE}/UrbanSoil.ORF.fna.xz"

# Output 100NT representative sequences
OUT_100 = f"{BASE}/US.100NT.fna.xz"

# Output redundancy relationships
OUT_MATCHES = f"{BASE}/US.100NT.matches.xz"


# Load all ORF sequences
seqs = list(fasta.fasta_iter(ORF_FNA))


# Sort sequences by length, sequence and header
def k(h_s):
    h, s = h_s
    return (len(s), s, h)


print(f"Loaded {len(seqs)}", flush=True)

seqs.sort(key=k)


# Store unique sequences and redundancy matches
deduped = []
matches = []

prev = ""
h_prev = None


# Remove exact duplicate sequences
for h, seq in seqs:

    if seq == prev:

        # "=" means exact duplicate
        matches.append((h, "=", h_prev))

    else:

        deduped.append((h, seq))
        prev = seq
        h_prev = h


# Remove original sequence list from memory
del seqs


# Length of the shortest unique sequence
min_len = len(deduped[0][1])


# Generate hashes for sequence windows
def rolling_hashes(seq, window_size):

    for i in range(len(seq) - window_size + 1):
        yield hash(seq[i:i + window_size])


# Store sequences sharing matching windows
hash_matches = defaultdict(list)

# Store sequences removed because they are contained
eliminated = set()


print(f"Deduped: {len(deduped)}", flush=True)


# Check whether shorter sequences are contained in longer sequences
for ix, (h, seq) in enumerate(deduped):

    if ix % 100000 == 0:
        print(lib.pad9("Dedupe iteration", ix), flush=True)

    candidates = set()

    # Find possible matching sequences
    for ha in set(rolling_hashes(seq, min_len)):

        candidates.update(hash_matches[ha])
        hash_matches[ha].append(ix)

    # Check candidate sequences
    for ix2 in candidates:

        if ix2 in eliminated:
            continue

        h2, seq2 = deduped[ix2]

        assert ix2 < ix

        # Check if the shorter sequence is inside the current sequence
        if seq2 in seq:

            # "C" means contained sequence
            matches.append((h2, "C", h))

            assert h2 != h

            eliminated.add(ix2)


# Write non-redundant 100NT representative sequences
with lib.xz_out(OUT_100) as f:

    n_ix = 0

    for ix, (h, seq) in enumerate(deduped):

        if ix in eliminated:
            continue

        # Create new 100NT ID
        n_h = lib.pad9("US.100NT", n_ix)

        # Keep original ORF ID in the FASTA header
        f.write(
            f">{n_h} {h}\n{seq}\n".encode("utf-8")
        )

        n_ix += 1


# Write exact duplicate and containment relationships
with lib.xz_out(OUT_MATCHES) as f:

    for h1, relation, h2 in matches:

        f.write(
            f"{h1}\t{relation}\t{h2}\n".encode("utf-8")
        )


# Final summary
print(f"100NT representatives: {n_ix}", flush=True)

print(f"Created: {OUT_100}", flush=True)

print(f"Created: {OUT_MATCHES}", flush=True)
