#!/bin/bash

# Build the combined Urban Soil ORF catalogue
python3 /work/microbiome/users/kruthi/intermediate_results/build_urban_soil_orf_catalog.py


# Submit with:
# mqsub --mem 20 -t 6 --cpus 12 --bg --segregated-log-file \
#     --name urban_orf_catalog \
#     -- bash /work/microbiome/users/kruthi/intermediate_results/build_urban_soil_orf_catalog.sh
