Step 1: External MAGs (cohorts)
        ↓
Step 2: Run SKANI vs SHD catalog
        ↓
        Output:
        Skani_Quality_Results/*.tsv
        (contains ANI values)

        ↓
Step 3: Filter using ANI ≥95
        ↓
        Matched (≥95)        Unmatched (<95)
        (same species)       → write to:
                             *_unmatched_paths.txt
                             (e.g., 52 MAGs)

        ↓
Step 4: Run FastANI on unmatched
        ↓
        Input:
        *_unmatched_paths.txt

        ↓
        Output:
        FastANI_results/*.tsv

        ↓
Step 5: Build Excel
        ↓
        Excel shows:
        ✔ MAGs with FastANI hits
        ❌ MAGs with no hits are missing
