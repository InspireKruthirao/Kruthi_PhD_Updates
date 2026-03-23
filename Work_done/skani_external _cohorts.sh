cd /work/microbiome/users/kruthi/MAGs_Onehealth/external_dog_cohorts/Quality_MAGs/Skani

QUALITY_DIR=/work/microbiome/users/kruthi/MAGs_Onehealth/external_dog_cohorts/Quality_MAGs
QUERY_LIST=/work/microbiome/users/kruthi/MAGs_Onehealth/dog_mags_list.txt

for cohort in Coelho_2018_dog Wang_2019_dogs Yarlagadda_2022_global_dog Allaway_2020_dogs Liu_2021_Canidae; do
    echo "Running HQ: $cohort"
    skani dist \
        --ql ${QUERY_LIST} \
        --rl ${QUALITY_DIR}/${cohort}_HQ_list.txt \
        --min-af 50 -t 40 \
        -o ${cohort}_HQ_ani.tsv
    echo "Done HQ: $cohort"
done

for cohort in Coelho_2018_dog Wang_2019_dogs Yarlagadda_2022_global_dog Allaway_2020_dogs Liu_2021_Canidae; do
    echo "Running MQ: $cohort"
    skani dist \
        --ql ${QUERY_LIST} \
        --rl ${QUALITY_DIR}/${cohort}_MQ_list.txt \
        --min-af 50 -t 40 \
        -o ${cohort}_MQ_ani.tsv
    echo "Done MQ: $cohort"
done

echo "ALL DONE!"
