library(gggenomes)
library(dplyr)

base_path <- "Z:/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods_25_30/SHD1_SM.100AA.003_485"

read_gff <- function(gff_file, contig_name) {
  lines <- readLines(gff_file)
  lines <- lines[!grepl("^#", lines)]
  df_list <- strsplit(lines, "\t")
  max_cols <- max(sapply(df_list, length))
  df_mat <- do.call(rbind, lapply(df_list, function(x) c(x, rep(NA, max_cols - length(x)))))
  df <- as.data.frame(df_mat, stringsAsFactors = FALSE)
  colnames(df)[1:9] <- c("orig_seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  df$feat_id <- sub(".*ID=([^;]+).*", "\\1", df$attributes)
  df$name <- sub(".*Name=([^;]+).*", "\\1", df$attributes)
  df$is_target <- grepl("TARGET_sMORF", df$attributes)
  df$seq_id <- contig_name
  df$strand <- ifelse(df$strand == "-", -1, 1)
  df %>% select(seq_id, start, end, strand, feat_id, type, name, is_target)
}

contig_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
contig_dirs <- contig_dirs[grepl("^D[0-9]+_contig_", basename(contig_dirs))]

all_genes <- NULL
for (d in contig_dirs) {
  gff <- list.files(d, pattern = "\\.gff$", full.names = TRUE)
  contig_name <- basename(d)
  genes <- read_gff(gff[1], contig_name)
  all_genes <- rbind(all_genes, genes)
}

contig_fingerprints <- all_genes %>%
  arrange(seq_id, start) %>%
  group_by(seq_id) %>%
  summarise(fingerprint = paste(name, collapse = "|"), .groups = "drop")

clustered <- contig_fingerprints %>%
  group_by(fingerprint) %>%
  summarise(representative = first(seq_id), count = n(), .groups = "drop")

contig_mapping <- clustered %>%
  select(representative, count) %>%
  mutate(display_name = paste0(representative, " (n=", count, ")"))

filtered_genes <- all_genes %>%
  filter(seq_id %in% clustered$representative) %>%
  left_join(contig_mapping, by = c("seq_id" = "representative")) %>%
  mutate(seq_id = display_name) %>%
  select(-display_name, -count)

seqs_data <- filtered_genes %>%
  group_by(seq_id) %>%
  summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
  mutate(length = end - start)

cog_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f",
                 "#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5","#c49c94",
                 "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")

all_cogs <- unique(filtered_genes$name[filtered_genes$name != "Unknown" & !filtered_genes$is_target])
cog_colors <- setNames(rep(cog_palette, length.out = length(all_cogs)), all_cogs)

filtered_genes <- filtered_genes %>%
  mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))

cog_colors["TARGET_sMORF"] <- "#000000"
cog_colors["Unknown"] <- "#FAFAFA"

p <- gggenomes(genes = filtered_genes, seqs = seqs_data) %>%
  flip(1:nrow(seqs_data))

p <- p +
  geom_seq() +
  geom_bin_label(size = 2) +
  geom_gene(aes(fill = gene_category)) +
  scale_fill_manual(values = cog_colors, name = "Annotation") +
  ggtitle(paste0(basename(base_path), " - Clustered Contigs")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 5, 10, 20),
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.text.x = element_text(size = 7)
  ) +
  guides(fill = guide_legend(ncol = 1))

cat("\nFound", nrow(contig_fingerprints), "total contigs\n")
cat("Clustered into", nrow(clustered), "unique groups\n")
cat("Removed", nrow(contig_fingerprints) - nrow(clustered), "duplicates\n\n")
print(clustered %>% select(representative, count) %>% arrange(desc(count)))

print(p)

plot_width <- 14
plot_height <- 10
plot_dpi <- 300

ggsave(file.path(base_path, "all_contigs_plot_clustered.jpg"), 
       plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, limitsize = FALSE)

cat("\nPlot saved\n")
