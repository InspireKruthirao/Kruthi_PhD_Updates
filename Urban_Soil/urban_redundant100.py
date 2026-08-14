import importlib.util
from collections import defaultdict

# Import helper modules
FASTA_PATH = "/work/microbiome/shanghai_dogs/resource_generation/fasta.py"
LIB_PATH = "/work/microbiome/shanghai_dogs/resource_generation/lib.py"

spec = importlib.util.spec_from_file_location("fasta", FASTA_PATH)
fasta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fasta)

spec = importlib.util.spec_from_file_location("lib", LIB_PATH)
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)

BASE = "/work/microbiome/users/kruthi/intermediate_results/urban_soil/Prodigal"

ORF_FNA = f"{BASE}/UrbanSoil.ORF.fna.xz"
OUT_100 = f"{BASE}/US.100NT.fna.xz"
OUT_MATCHES = f"{BASE}/US.100NT.matches.xz"


# Load ORF sequences
seqs = list(fasta.fasta_iter(ORF_FNA))


def k(h_s):
    h, s = h_s
    return (len(s), s, h)


print(f"Loaded {len(seqs)}", flush=True)

seqs.sort(key=k)

deduped = []
matches = []
prev = ""
h_prev = None

for h, seq in seqs:
    if seq == prev:
        matches.append((h, "=", h_prev))
    else:
        deduped.append((h, seq))
        prev = seq
        h_prev = h

del seqs

min_len = len(deduped[0][1])


def rolling_hashes(seq, window_size):
    for i in range(len(seq) - window_size + 1):
        yield hash(seq[i:i + window_size])


hash_matches = defaultdict(list)
eliminated = set()

print(f"Deduped: {len(deduped)}", flush=True)

for ix, (h, seq) in enumerate(deduped):

    if ix % 100000 == 0:
        print(lib.pad9("Dedupe iteration", ix), flush=True)

    candidates = set()

    for ha in set(rolling_hashes(seq, min_len)):
        candidates.update(hash_matches[ha])
        hash_matches[ha].append(ix)

    for ix2 in candidates:

        if ix2 in eliminated:
            continue

        h2, seq2 = deduped[ix2]

        assert ix2 < ix

        if seq2 in seq:
            matches.append((h2, "C", h))
            assert h2 != h
            eliminated.add(ix2)


# Write 100NT representative sequences
with lib.xz_out(OUT_100) as f:

    n_ix = 0

    for ix, (h, seq) in enumerate(deduped):

        if ix in eliminated:
            continue

        n_h = lib.pad9("US.100NT", n_ix)

        f.write(
            f">{n_h} {h}\n{seq}\n".encode("utf-8")
        )

        n_ix += 1


# Write ORF redundancy mappings
with lib.xz_out(OUT_MATCHES) as f:

    for h1, relation, h2 in matches:

        f.write(
            f"{h1}\t{relation}\t{h2}\n".encode("utf-8")
        )


print(f"100NT representatives: {n_ix}", flush=True)
print(f"Created: {OUT_100}", flush=True)
print(f"Created: {OUT_MATCHES}", flush=True)
