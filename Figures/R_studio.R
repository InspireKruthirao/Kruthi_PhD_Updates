library(gggenomes)
library(dplyr)

base_path <- "Z:/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods_25_30/SHD1_SM.100AA.003_915"

read_gff <- function(gff_file, contig_name) {
  df <- read.table(gff_file, sep="\t", header=FALSE, comment.char="#",
                   col.names=c("seq_id", "source", "type", "start", "end", 
                              "score", "strand", "phase", "attributes"),
                   colClasses=c("character", "character", "character", 
                               "integer", "integer", "character", 
                               "character", "character", "character"))
  
  df$feat_id <- sub(".*ID=([^;；]+).*", "\\1", df$attributes)
  df$name <- sub(".*Name=([^;；]+).*", "\\1", df$attributes)
  df$is_target <- grepl("TARGET_sMORF", df$attributes)
  df$seq_id <- contig_name
  
  df %>% select(seq_id, start, end, strand, feat_id, type, name, is_target)
}

contig_dirs <- list.dirs(base_path, recursive=FALSE, full.names=TRUE)
contig_dirs <- contig_dirs[grepl("D[0-9]+_contig_", contig_dirs)]

all_genes <- data.frame()
for (contig_dir in contig_dirs) {
  gff_files <- list.files(contig_dir, pattern="\\.gff$", full.names=TRUE)
  contig_name <- sub("^D[0-9]+_", "", basename(contig_dir))
  genes <- read_gff(gff_files[1], contig_name)
  all_genes <- rbind(all_genes, genes)
}

all_genes <- all_genes %>% mutate(strand = "+")

seqs_data <- all_genes %>%
  group_by(seq_id) %>%
  summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
  mutate(length = end - start)

unique_cogs <- unique(all_genes$name[!all_genes$is_target & all_genes$name != "Unknown"])

vibrant_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")

cog_colors <- setNames(vibrant_colors[1:length(unique_cogs)], unique_cogs)
cog_colors["Unknown"] <- "#D3D3D3"

all_genes <- all_genes %>%
  mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))

final_colors <- c("TARGET_sMORF" = "#000000", cog_colors)

legend_breaks <- c("TARGET_sMORF", sort(unique_cogs), "Unknown")

p <- gggenomes(genes = all_genes, seqs = seqs_data)

plot_result <- p +
  geom_seq() +
  geom_bin_label(size = 3) +
  geom_gene(aes(fill = gene_category)) +
  scale_fill_manual(values = final_colors, name = "Gene Annotation", breaks = legend_breaks) +
  theme_minimal() +
  theme(legend.position = "right", legend.text = element_text(size = 9),
        legend.title = element_text(size = 11, face = "bold"),
        plot.margin = margin(10, 10, 10, 10), axis.text = element_text(size = 9),
        panel.grid.minor = element_blank())

print(plot_result)

plot_width <- 14
plot_height <- max(6, nrow(seqs_data) * 0.5)

ggsave(file.path(base_path, "all_contigs_plot.jpg"), 
       plot = plot_result, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)

cat("\n=== Summary ===\n")
cat("Total contigs:", nrow(seqs_data), "\n")
cat("Total genes:", nrow(all_genes), "\n")
cat("Target sMORFs:", sum(all_genes$is_target), "\n")
