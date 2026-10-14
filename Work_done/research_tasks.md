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
