library(gggenomes)
library(dplyr)
library(ggplot2)

base_path <- "Z:/microbiome/users/kruthi/SmORF_neighbourhoods_25_30/SHD1_SM.100AA.003_485"

read_gff <- function(gff_file, contig_name) {
  lines <- readLines(gff_file)
  lines <- gsub("；", ";", lines)
  lines <- lines[!grepl("^#", lines)]
  if(length(lines) == 0) return(NULL)
  
  df_list <- strsplit(lines, "\t")
  max_cols <- max(sapply(df_list, length))
  df_mat <- do.call(rbind, lapply(df_list, function(x) c(x, rep(NA, max_cols - length(x)))))
  df <- as.data.frame(df_mat, stringsAsFactors = FALSE)
  colnames(df)[1:9] <- c("orig_seqid","source","type","start","end",
                         "score","strand","phase","attributes")
  
  df$start     <- as.integer(df$start)
  df$end       <- as.integer(df$end)
  df$feat_id   <- sub(".*ID=([^;]+).*", "\\1", df$attributes)
  df$name      <- sub(".*Name=([^;]+).*", "\\1", df$attributes)
  df$is_target <- grepl("TARGET_sMORF", df$attributes)
  df$seq_id    <- contig_name
  df$strand    <- 1
  df %>% select(seq_id, start, end, strand, feat_id, type, name, is_target)
}

all_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
contig_dirs <- all_dirs[grepl("_contig_", basename(all_dirs))]
if(length(contig_dirs) == 0) stop("No contig directories found!")

all_genes <- NULL
for(d in contig_dirs) {
  gff <- list.files(d, pattern = "\\.gff$", full.names = TRUE)
  if(length(gff) == 0) next
  genes <- read_gff(gff[1], basename(d))
  if(!is.null(genes)) all_genes <- rbind(all_genes, genes)
}
if(is.null(all_genes)) stop("No genes were loaded.")

contig_fingerprints <- all_genes %>%
  arrange(seq_id, start) %>%
  group_by(seq_id) %>%
  summarise(fingerprint = paste(name, collapse = "|"), .groups = "drop")

clustered <- contig_fingerprints %>%
  group_by(fingerprint) %>%
  summarise(representative = first(seq_id), count = n(), .groups = "drop")

contig_mapping <- clustered %>%
  select(representative, count) %>%
  mutate(display_name = paste0(representative))

filtered_genes <- all_genes %>%
  filter(seq_id %in% clustered$representative) %>%
  left_join(contig_mapping, by = c("seq_id" = "representative")) %>%
  mutate(seq_id = display_name) %>%
  select(-count, -display_name)

seqs_data <- filtered_genes %>%
  group_by(seq_id) %>%
  summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
  mutate(length = end - start)

cog_palette <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5"
)
all_cogs <- unique(filtered_genes$name[!filtered_genes$is_target & filtered_genes$name != "Unknown"])
cog_colors <- setNames(rep(cog_palette, length.out = length(all_cogs)), all_cogs)
cog_colors["TARGET_sMORF"] <- "#000000"
cog_colors["Unknown"]      <- "#D3D3D3"
filtered_genes <- filtered_genes %>% mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))

p <- gggenomes(genes = filtered_genes, seqs = seqs_data) +
  geom_seq() +
  geom_bin_label(size = 2) +
  geom_gene(aes(fill = gene_category)) +
  scale_fill_manual(values = cog_colors, name = "Annotation") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )

print(p)

ggsave(file.path(base_path, "all_contigs_plot_clustered.jpg"),
       plot = p, width = 24, height = 12, dpi = 300, limitsize = FALSE)

cat("\nPlot saved.\n")
