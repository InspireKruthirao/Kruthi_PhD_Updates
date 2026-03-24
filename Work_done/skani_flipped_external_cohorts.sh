cd /work/microbiome/users/kruthi/MAGs_Onehealth/external_dog_cohorts/Quality_MAGs/Skani

QUALITY_DIR=/work/microbiome/users/kruthi/MAGs_Onehealth/external_dog_cohorts/Quality_MAGs
CATALOG=/work/microbiome/users/kruthi/MAGs_Onehealth/dog_mags_list.txt

for cohort in Coelho_2018_dog Wang_2019_dogs Yarlagadda_2022_global_dog Allaway_2020_dogs Liu_2021_Canidae Xu_2019_dogs Worsley-Tonks_2020_dog; do
    echo "Running HQ: $cohort"
    skani dist \
        --ql ${QUALITY_DIR}/${cohort}_HQ_list.txt \
        --rl ${CATALOG} \
        --min-af 50 -t 40 \
        -o ${cohort}_HQ_ani_flipped.tsv
    echo "Running MQ: $cohort"
    skani dist \
        --ql ${QUALITY_DIR}/${cohort}_MQ_list.txt \
        --rl ${CATALOG} \
        --min-af 50 -t 40 \
        -o ${cohort}_MQ_ani_flipped.tsv
    echo "Done: $cohort"
done
echo "ALL DONE!"
