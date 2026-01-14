library(gggenomes)
library(dplyr)

parent_path <- "Z:/microbiome/users/kruthi/SmORF_neighbourhoods_26_30"
all_plots_dir <- file.path(parent_path, "ALL_PLOTS")

# FIX: create ALL_PLOTS upfront so ggsave doesn't prompt/stop
dir.create(all_plots_dir, showWarnings = FALSE, recursive = TRUE)

smorf_dirs <- list.dirs(parent_path, recursive = FALSE, full.names = TRUE)
smorf_dirs <- smorf_dirs[grepl("SHD1_SM", basename(smorf_dirs))]

cat("Found", length(smorf_dirs), "smORF directories to process\n")

if (length(smorf_dirs) == 0) {
  cat("\nERROR: No smORF directories found!\n")
  cat("Checking what's in the parent directory...\n")
  all_dirs <- list.dirs(parent_path, recursive = FALSE, full.names = FALSE)
  cat("First 10 directories:\n")
  print(head(all_dirs, 10))
  stop("Please check the directory structure")
}

read_gff <- function(gff_file, contig_name) {
  lines <- readLines(gff_file, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (length(lines) == 0) return(NULL)

  df_list <- strsplit(lines, "\t")
  max_cols <- max(sapply(df_list, length))
  df_mat <- do.call(rbind, lapply(df_list, function(x) c(x, rep(NA, max_cols - length(x)))))
  df <- as.data.frame(df_mat, stringsAsFactors = FALSE)

  colnames(df)[1:9] <- c("orig_seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")

  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  df$feat_id <- sub(".*ID=([^;]+).*", "\\1", df$attributes)
  df$name <- sub(".*Name=([^;]+).*", "\\1", df$attributes)
  df$is_target <- grepl("TARGET_sMORF", df$attributes)
  df$seq_id <- contig_name

  df$strand <- 1

  df %>% select(seq_id, start, end, strand, feat_id, type, name, is_target)
}

# color palette
cog_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f",
                 "#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5","#c49c94",
                 "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")

# Loop through each smORF directory
for (i in seq_along(smorf_dirs)) {
  base_path <- smorf_dirs[i]
  smorf_name <- basename(base_path)
  cat("\nProcessing [", i, "/", length(smorf_dirs), "]:", smorf_name, "\n")

  # Get contig directories
  contig_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
  contig_dirs <- contig_dirs[grepl("^D[0-9]+_contig_", basename(contig_dirs))]

  if (length(contig_dirs) == 0) {
    cat("  ⚠ No contig directories found, skipping...\n")
    next
  }

  all_genes <- NULL
  for (d in contig_dirs) {
    gff <- list.files(d, pattern = "\\.gff$", full.names = TRUE)
    if (length(gff) == 0) next
    contig_name <- basename(d)
    genes <- read_gff(gff[1], contig_name)
    if (!is.null(genes)) all_genes <- rbind(all_genes, genes)
  }

  if (is.null(all_genes) || nrow(all_genes) == 0) {
    cat("  No genes found, skipping...\n")
    next
  }

  seqs_data <- all_genes %>%
    group_by(seq_id) %>%
    summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
    mutate(length = end - start)

  all_cogs <- unique(all_genes$name)
  cog_colors <- setNames(rep(cog_palette, length.out = length(all_cogs)), all_cogs)
  cog_colors["Unknown"] <- "#D3D3D3"
  cog_colors["TARGET_sMORF"] <- "#000000"

  all_genes <- all_genes %>%
    mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))

  p <- gggenomes(genes = all_genes, seqs = seqs_data) +
    geom_seq() +
    geom_bin_label(size = 2) +
    geom_gene(aes(fill = gene_category)) +
    scale_fill_manual(values = cog_colors, name = "Annotation") +
    ggtitle(basename(base_path)) +
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

  plot_path_local <- file.path(base_path, "all_contigs_plot.jpg")
  ggsave(plot_path_local, plot = p, width = 14, height = 10, dpi = 300, limitsize = FALSE)
  cat("  ✓ Saved to:", plot_path_local, "\n")

  plot_path_all <- file.path(all_plots_dir, paste0(smorf_name, "_plot.jpg"))
  ggsave(plot_path_all, plot = p, width = 14, height = 10, dpi = 300, limitsize = FALSE)
  cat("  ✓ Saved to:", plot_path_all, "\n")
}

cat("\n==\n")
cat("✓ All plots generated successfully!\n")
cat("Total processed:", length(smorf_dirs), "smORFs\n")
cat("=\n")
