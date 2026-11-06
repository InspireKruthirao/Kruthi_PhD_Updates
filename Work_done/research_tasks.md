27/11/2025

- luis: Add cog category alphabets (beside each cog category, so that ranking of smorfs can be done)
- Sam- Redherons presentation- at what level if species get changes, amino acid level changes that difference can be seen in functions
- Luis- evrything to be on git, no docs at any cost-
- Luis- Plot amino acids plot(refer -SEP work) 
- Luis- Get the fraction of the SmORFs that are overlapping -let me divide the catgories and do it
- Luis- What if the contig is having same neighborhoods-how will you decide about it(whether to keep it or remove it / would i just look at 5+_ neighbors or more?)

  28/11/2025
- Luis -Compare the GMSC-Mapper and Paul wilmers's paper data- Run GMSC Mapper on them and their code --see the difference -how many smorfs are being called -Longterm (its not for now)
- Luis - EggNOG Mapper - predicts only  70-80%

  2/12/2025
  - Kruthi - Chcek the over lapping , how many base pairs are being overlapped , which data frames
  - Kruthi - Plot the bar graph based on the how many SmORFs appearing each bin 1-10. 20-30 - check the accuracy and plot based on the relavance
  - KRUTHI - When the SmORF has two annotated neighbors - Think about checking the conserved with BED FORMAT Conversion and then checking the clusters that have repeated conserveness.
  - Kruthi - Check the EggNOG version for the SmORFs that have been annotated by Luis and why is the difference between mostly "L"
  
  4/12/2025
   - Luis - Do you looks at the annotated SmORFs? Kruthi to add the stats of such SmORFs
   - Luis - Add clusters to the duplicated Contigs in the  SmORFs PLots
   - Kruthi - Add Taxonomy (Include the MAG's taxonomy whereever the data is available)
   - Luis - Remove the duplicates of the contigs in the SmORF plots.
   - Luis - Create personal folder /Microbiome/Users/Kruthi
   - Kruthi - do the colors properly and hardcore the colors for the COGs
   - Kruthi- Try plottoing the annotated SmORFs in differnt way

  12/10/2026
   - Kruthi - Added rules C101-C150 to sweep_db.py (conserved neighbourhood, sandwiched, two methods agree, enrichment, conservation-weighted vote, very steep cutoffs); Luis's code and C001-C100 unchanged. Best: sandwiched excluding S, 88.2% accuracy at 13.7% coverage (C130)
   - Kruthi - Fixed the eggNOG annotations of the neighbourhoods: genes matched to eggNOG by sample + Original_ID, one EggNOG_annotations.tsv line per gene, correct Contig column
   - Kruthi - build_smorf_neighbourhoods.py (SmORF_Catalog): one script that writes all neighbourhoods into one database, smorf_neighbourhoods.sqlite, instead of millions of files; stores neighbors.tsv, EggNOG_annotations.tsv and the GFF of every occurrence, plus gene coordinates, strand and distances for distance-based rules
   - Kruthi - Running it on the HPC as 100 PBS jobs at the same time, merged automatically at the end (run_smorf_neighbourhoods.sh)
   - Next - distance and strand rules (C151 onwards) on the new database; add the QUT project code (-P) to PBS jobs

  13/10/2026
   - Kruthi - smorf_neighbourhoods.sqlite built: 156,109 smORFs, 1,609,063 occurrences, 17.3M neighbour genes, 80.9% with eggNOG (eggNOG issue fixed); 100 pieces ran in ~30 min
   - Kruthi - Same sweep results as the old database (C130 88.2%, C027 38.9%), so the new pipeline is correct
   - Kruthi - Added rules C151-C184 using gene strand and distance: operon vote, sandwiched with strand/distance, distance-weighted vote
   - Kruthi - Best: C172 (sandwiched, same strand, within 100 bp, excl. S and R) 93.8% accuracy at 9.7% coverage; C164 (same strand) 91.1% at 13.7%; C175 (vote by bp distance) 42.1% at 99.2%, better than the original method
   - To do - cascade of rules (most accurate rule first, then the next) to predict ~99% of smORFs with accuracy per step
   - To do - check the best rules on held-out bins (avoid picking the best of 184 on the same data)
   - To do - remove the old GFF scripts (Create_UrbanSoil_GFFs.py and its PBS files), now replaced by build_smorf_neighbourhoods.py
   - To do - ask Luis about deleting the old neighbourhood folders; add the QUT project code (-P) to PBS jobs

  14/10/2026
   - Kruthi - Fixed _sql_submit.pbs (Analysis): correct Contig column, Original_ID column added, EggNOG_annotations.tsv now one line per gene ("-" if no eggNOG), real eggNOG match-rate check, second occurrence on a contig no longer overwritten, writes the GFF itself (Create_UrbanSoil_GFFs.py no longer needed)
   - Kruthi - Reran it as 100 PBS jobs: all Neighbourhood_Analysis_SQL_ge6 folders rewritten, files identical to smorf_neighbourhoods.sqlite
   - Kruthi - eggNOG seed_orthologs files have qcov, scov, pident, evalue, bitscore for each hit; wrote eggnog_hit_quality.py (SmORF_Catalog) to measure how well neighbour/target genes match eggNOG
   - To do - run eggnog_hit_quality.py, choose a coverage/identity threshold from the table
   - To do - rerun the 184 rules with weak eggNOG hits removed and compare accuracy
   - To do - cascade of rules + check on held-out bins; remove Create_UrbanSoil_GFFs.py and its PBS files

   20/10/26
   Cluster instead of hand-made rules. Turn each occurrence’s ±5 neighbourhood into a gene string, compare occurrences by similarity (Jaccard), and group similar contexts into clusters.
Make it fast with similarity-search tools like Spotify’s Annoy, or MinHash.
Show the results as a network: each cluster is a node, and overlapping clusters are joined by weaker links.
Don’t collapse genes to COG letters too early, because one letter covers many different genes.
Think about context size: the same smORF can sit in a long gene block in one genome and a short one in another.
Watch for blind spots: fixed rules can miss gene arrangements they weren’t designed for.
Citations matter more than journal impact factor.

What to do:

Prototype the clustering on smorf_neighbourhoods.sqlite: compare each smORF’s occurrences by their neighbours’ COG IDs (Jaccard), then cluster them.
Give each cluster a consensus COG category and a size, i.e. how many occurrences it holds.
Add “biggest cluster’s category” as rule C185 and score it against C172 and the cascade, using the held-out check.
Use COG IDs for matching and COG letters for reporting, which answers the “too coarse” concern.
Note the ±5 window as a limitation, and check whether long and short gene blocks end up in separate clusters.
Later: a network view of the clusters for the interface, and repeating the analysis at 90% and 95% smORF identity.

What’s still pending from before:

Get the held-out and cascade results from Friday’s sweep.
Run the three eggNOG-filtered sweeps and compare them.
Plot accuracy against coverage for all 184 rules.

  22/10/2026
   - Kruthi - Added rules C185-C188 to sweep_db.py: nearest-neighbour vote with Annoy (Spotify), as suggested by Andy. Each smORF is described by the COG IDs within 3 genes (top 512 COG IDs); the 10 or 25 smORFs with the most similar neighbourhoods vote for their known category. smORFs sharing a target gene cannot vote for each other

  26/10/2026
   - Kruthi - plot_rules.py: accuracy against coverage for all rules, with the best trade-offs (Pareto front); run automatically by run_sweep.sh
   - Kruthi - permutation_test.py: the neighbourhood signal is real. Real accuracy 39-69% vs 16-21% when each smORF gets another smORF's neighbourhood, and vs 35-41% when gene order is shuffled within the window; p < 0.01 for all 3 rules (100 randomisations). Sandwich rule: 68.8% real vs 40.6% order-shuffled, so gene position matters
   - Kruthi - distance_decay.py: the signal fades with distance. A neighbour shares the smORF's category 40.5% of the time at rank 1, 19.6% at rank 5, ~11% at random; same strand about twice opposite strand (46.0% vs 22.7% at rank 1); 51.2% for same-strand genes within 100 bp. More neighbours raise coverage (73% -> 94%) but lower accuracy (45.4% -> 36.4%)
   - Kruthi - sweep_db.py shows the results on screen; files only with --save (run_sweep.sh uses --save)
   - Kruthi - Added rules C189-C192: the same nearest-neighbour vote with Voyager (Spotify's successor to Annoy, more accurate: 100% vs 71% of the true 10 nearest found in a test), for scaling up later (SPIRE). Python 3.14 -> 3.13 in the pixi environment for Voyager
   - Kruthi - Results (192 rules): Annoy C185 57.7% accuracy at 87.7% coverage, Voyager C189 58.0% at 80.0%, vs 42.3% at 88.5% for the best rule using only the smORF's own neighbours (C055)
   - Kruthi - Held-out check (cascade.py): best rules equally good on TRAIN and TEST (C172 93.8% / 93.8%), so no overfitting. Cascade: 97.8% of smORFs predicted at 56.1% accuracy; first 3 steps 17.4% at 78.3%

  27/10/2026
   - Kruthi - sweep_db.py faster, same results (checked on 4 test databases): database read once instead of 3 times, Annoy/Voyager vote with numpy, search shared by the rules that differ only in excluded categories, conserved-neighbourhood counts made once, Annoy searches on all CPUs; merged duplicate functions (sandwiched/sandwiched_by, annoy_vote/voyager_vote). Full sweep 28.9 -> 12.1 min
   - Kruthi - Annoy/Voyager: a neighbour with nothing in common (similarity 0) no longer votes, so no prediction without evidence (removed the RuntimeWarning)
   - Kruthi - pixi.toml: [workspace] instead of the old [project] (no more warning); deleted run_signal_tests.sh (the tests are run directly when needed)
   - Kruthi - C172 in numbers: 9,451 smORFs predicted, 7,353 checked at 93.8% accuracy, 1,251 smORFs without a known function get a high-confidence prediction
   - Kruthi - KEGG: eggNOG already gives KEGG annotations. Of genes with eggNOG: COG category 94.5%, KO 52.5%, Pathway 32.5%, Module 21.4%. kegg_analysis.py (Analysis) counts how many smORFs could be evaluated with each label (running)
   - Kruthi - Read Reed, Yuan & de Crecy-Lagard 2026 (preprint, tRNA modification enzymes, COG vs KO): broad COG families mix functionally distinct paralogs (51% COG-KO agreement vs 96% for one-to-one families), so COG labels are imperfect ground truth and KO can be a stricter check
   - To do - Annoy on all CPUs was slower on the HPC (C185 143 -> 164 s): check how many CPUs Python sees in the PBS job
   - To do - dereplication: split TRAIN/TEST by 90% smORF cluster, and stop near-identical smORFs voting for each other in Annoy/Voyager
   - To do - run the 192 rules with KEGG Pathway/Module labels as a second check; run the eggNOG-filtered sweep
   - To do - remove the old GFF scripts in Analysis; README for the smORF neighbourhood pipeline; add the QUT project code (-P) to PBS jobs before 2 November

  28/10/2026 - Supervisor meeting (Luis, Ben)
   - Luis - Framework delivered as planned for this semester; Aim 1 can be packaged as two methods: C172 strict (high quality) and C175 broad (high coverage). Keep exploring until early 2027, then decide papers/thesis chapters and aim higher than the minimum
   - Ben - Decide the goal: best predictor (then use machine learning, e.g. random forest on one-hot encoded neighbours) or understanding (then change one factor at a time, e.g. same strand vs not; C172 vs C166 changes two things)
   - Ben/Luis - Annoy/Voyager not needed at 150k smORFs: do exact nearest neighbours (one matrix multiplication); use approximate search only at SPIRE scale and measure the loss
   - Ben - Vector idea: one-hot COG letters per neighbour, closer genes weighted more; Luis: this can learn that some categories are more informative
   - Luis - Baseline should be "always predict the most common category" (C097 26.1%, C098 25.6%), not random (C099)
   - Luis - Write precisely how occurrences are combined and conflicts handled (e.g. 5 occurrences: 2 say K, 3 say L)
   - Ben - Contigs from the same species are not independent evidence (strain-level copies); matters at SPIRE scale
   - Ben - Requiring >= 6 samples throws away most data; cluster smORFs at 90% (CD-HIT clusters exist) for more coverage; scale up with COG rather than new function databases; stick with COG first, KEGG later (KEGG is more specific, not broader)
   - Luis - Say contigs/contexts, not genomes, when that is what is meant
   - Luis/Ben - Benchmark the rules on all proteins, not only smORFs
   - Luis - Report each prediction with its estimated accuracy; open question: does it transfer to short-read assemblies?
   - Hongyu wants to run it on cow microbiome data
   - Next meeting around January 2027 (2-year mark 21 Jan 2027); second-year seminar and report; Spain visit: ask the consulate to amend the dates, aim Feb-Apr 2027
   - To do - write the aggregation/conflict description; use C097/C098 as baseline; one-factor-at-a-time rule comparisons; exact nearest neighbours vs Annoy/Voyager; ML classifier on one-hot neighbours; rerun with smORFs clustered at 90%; rules on all proteins; make the pipeline easy for Hongyu

  30/10/2026
   - Kruthi - Added rules C193-C228 to sweep_db.py: coupled/overlapping genes, upstream-only and downstream-only votes, same strand at any distance, operon + distance weighting, margin cutoff, tie-breaking (228 rules)
   - Kruthi - Results (228 rules): none of the new rules entered the cascade (still 56.1% at 97.8%); downstream genes are more informative than upstream genes

  04/11/2026 - Feedback from the 28 October meeting
   - Kruthi - Added rules C229-C236: exact nearest neighbours (Ben/Luis; 51 s for all 151k smORFs) and similar COG letter profiles (Ben's vector idea: one-hot letters of 5 genes per side, closer genes weighted more)
   - Kruthi - Exact vs approximate: Annoy finds 54% and Voyager 56% of the true 10 nearest neighbours; exact search is 4-6 points more accurate (C229 62.8% vs Voyager C189 58.2% on TEST)
   - Kruthi - Letter profiles: C233 61.2% at 94.6% coverage (highest coverage of any single rule)
   - Kruthi - Random forest (Ben) on one-hot COG letters per neighbour position + strand + distance: 60.0% at 100% coverage on TEST; at confidence >= 0.5 93.8% at 26.8% (C172: 93.8% at 9.7%). Split TRAIN/TEST by 90% cluster (52,413 clusters): 57.0%, so near-identical smORFs add only ~3 points
   - Kruthi - Baseline (Luis): always the most common category 26.1% (C097)
   - Kruthi - One factor at a time (Ben), 592 rule pairs differing in one setting (median change): same strand +3.8 accuracy / -1.1 coverage; downstream vs upstream +4.4; operon gap 0 -> 200 bp -7.1 / +49.3; neighbours 1 -> 5 -6.4 / +17.6; no change beyond 5 neighbours; dividing votes by how common a category is -8.1; tie-breaking and minimum occurrences ~0
   - Kruthi - Cascade with the 236 rules: 60.8% at 99.9% (was 56.1% at 97.8%); exact search C229 replaced Voyager
   - Kruthi - PBS jobs need -P DFAZCB7230 (lab RPID, from mqsub); GitHub push over SSH from the HPC (key id_ed25519_github)

  05/11/2026
   - Kruthi - Random forest added to sweep_db.py as rules C237-C244 (all / excluding S and R; confidence cutoffs none, 40, 50, 60%). Each half of the annotated smORFs is predicted by a model trained on the other half (the cascade's TRAIN/TEST halves), so accuracies are not inflated
   - Kruthi - Normal proteins (Luis/Ben: "rules on all proteins"; also for data without annotated smORFs): protein_neighbourhoods_sqlite.py built a NEW database of 215,985 random proteins > 100 aa with a COG category (Protein_Neighbourhoods/); protein_tests.py
   - Kruthi - Trained on proteins only: 44% correct on smORFs (vs 60% trained on smORFs, 39% original method, 26% baseline); at confidence >= 0.8 still 97% (5% of smORFs). Usable when a dataset has no annotated smORFs
   - Kruthi - Proteins + smORFs together: 58% (no better than smORFs alone). Neighbours predict normal proteins less well (44%) than smORFs (60%): smORF neighbourhoods are more informative
   - Kruthi - Rules on proteins vs smORFs: accuracies correlate (r = 0.66, 206 rules); sandwiched C166 87% on proteins vs 92% on smORFs; similar-neighbourhood rules much worse on proteins (C229 31% vs 63%)
   - Kruthi - Found: 97% of smORFs overlap a Prodigal gene of about the same length (their known category is their own eggNOG annotation); 3% are intergenic and take the next gene's category
   - To do - leave the 3% intergenic smORFs out of the scoring

  06/11/2026
   - Kruthi - Results (244 rules): random forest C243 (confidence >= 60%) 98.2% at 18.6% coverage, the most accurate rule (C172 93.8% at 9.7%); C241 (>= 50%) 94.0% at 26.6%; C237 60.0% at 100%
   - Kruthi - Cascade (244 rules): 63.2% at 100% coverage; first 5 steps are random forest rules, 49.2% of smORFs at 81.4%
   - Kruthi - Renamed sweep_db.py to predict_function.py; feedback_tests.py merged into it as --report (random forest by confidence, one factor at a time). Dropped the one-off checks (exact vs Annoy/Voyager, 90% cluster split): answered on 04/11
   - Kruthi - Results on GitHub: Pixi_Env/combinations_results.csv (244 rules), feedback_results/, protein_results/
   - To do - check whether the random forest's near-100% cases are trivial (e.g. a smORF inside a bigger gene with the same annotation); write the aggregation/conflict description (Luis); rerun with smORFs clustered at 90% (Ben); simple steps for Hongyu (cow microbiome)
