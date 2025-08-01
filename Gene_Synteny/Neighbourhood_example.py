import pandas as pd
import gzip
import lzma

smorf_origins = '/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/SHD_SMORF_resource/100AA_SmORFs_origins.tsv.gz'
gene_catalog = '/work/microbiome/shanghai_dogs/data/ShanghaiDogs_OtherResources/GeneCatalog/SHD.ORF.orig.tsv.xz'
final_outfile = '/home/n12228516/31july/smorf_genes_20_neighbours.tsv'
target_smorfs = ['SHD1_SM.100AA.000_000', 'SHD1_SM.100AA.000_001']

def map_smorfs_neighbors(infile1, infile2, outfile, target_smorfs):
    """Map smORFs to contigs and find upstream/downstream genes in a 20-gene window."""
    #gene catalog
    with lzma.open(infile1, 'rt') as f:
        genes_df = pd.read_csv(f, sep='\t', header=0)

    #ontig dictionary
    contig_dict = {}
    for _, row in genes_df.iterrows():
        contig = f"{row['Sample']}_{row['Original_ID']}"
        if contig not in contig_dict:
            contig_dict[contig] = []
        contig_dict[contig].append({
            'ORF': row['ORF'],
            'Start': int(row['Start']),
            'End': int(row['End']),
            'Strand': int(row['Strand']),
            'Info': f"{row['ORF']}:{row['Start']}-{row['End']}#{row['Strand']}"
        })

    #smORF origins
    with gzip.open(infile2, 'rt') as f:
        smorfs_df = pd.read_csv(f, sep='\t', header=None, names=['sORF', 'Sample', 'Contig', 'Position', 'Strand'])

    #target smORFs
    smorfs_df = smorfs_df[smorfs_df['sORF'].isin(target_smorfs)]

    #output
    with open(outfile, 'wt') as out:
        out.write("sORF\tSample\tContig\tsORF_Start\tsORF_End\tsORF_Strand\tUpstream_Genes\tDownstream_Genes\tNotes\n")
        for _, row in smorfs_df.iterrows():
            smorf = row['sORF']
            sample = row['Sample']
            contig_id = row['Contig']
            contig = f"{sample}_{contig_id}"
            position = row['Position']
            strand = row['Strand']
            start, end = map(int, position.split('-'))

            if contig not in contig_dict:
                out.write(f"{smorf}\t{sample}\t{contig_id}\t{start}\t{end}\t{strand}\tNA\tNA\tContig not found\n")
                continue

            genes = contig_dict[contig]
            genes_sorted = sorted(genes, key=lambda x: x['Start'])

            # Find smORF position
            smorf_idx = None
            closest_gene = None
            min_distance = float('inf')
            for i, gene in enumerate(genes_sorted):
                start_diff = abs(gene['Start'] - start)
                end_diff = abs(gene['End'] - end)
                distance = min(start_diff, end_diff)
                if distance < min_distance:
                    min_distance = distance
                    closest_gene = gene
                if (start_diff <= 50 and 
                    end_diff <= 50 and 
                    gene['Strand'] == (1 if strand == '+' else -1)):
                    smorf_idx = i
                    break

            if smorf_idx is None:
                closest_note = f"Closest gene: {closest_gene['ORF']}:{closest_gene['Start']}-{closest_gene['End']} (distance={min_distance} bp)" if closest_gene else "No genes in contig"
                out.write(f"{smorf}\t{sample}\t{contig_id}\t{start}\t{end}\t{strand}\tNA\tNA\tsmORF not found; {closest_note}\n")
                continue

            #10 upstream and 10 downstream genes
            upstream = []
            downstream = []
            for i, gene in enumerate(genes_sorted):
                if gene['Strand'] == (1 if strand == '+' else -1):
                    if i < smorf_idx:
                        upstream.append(gene['Info'])
                    elif i > smorf_idx:
                        downstream.append(gene['Info'])

            upstream = upstream[-10:]
            downstream = downstream[:10]
            upstream_str = ','.join(upstream) if upstream else 'NA'
            downstream_str = ','.join(downstream) if downstream else 'NA'
            notes = f"Contig has {len(genes_sorted)} genes"

            out.write(f"{smorf}\t{sample}\t{contig_id}\t{start}\t{end}\t{strand}\t{upstream_str}\t{downstream_str}\t{notes}\n")

map_smorfs_neighbors(gene_catalog, smorf_origins, final_outfile, target_smorfs)

