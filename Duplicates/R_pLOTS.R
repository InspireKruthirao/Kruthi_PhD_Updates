#!/usr/bin/env Rscript
library(gggenomes)
library(dplyr)

# Base directory containing all smORF folders
base_dir <- "/work/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods_25_30"

# Create ALL_PLOTS directory
all_plots_dir <- file.path(base_dir, "ALL_PLOTS")
dir.create(all_plots_dir, showWarnings = FALSE, recursive = TRUE)

# Function to read GFF file
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

# Function to generate plot for one smORF
generate_smorf_plot <- function(smorf_path, smorf_name) {
  cat(sprintf("Processing %s...\n", smorf_name))
  
  # Get all contig directories
  contig_dirs <- list.dirs(smorf_path, recursive=FALSE, full.names=TRUE)
  contig_dirs <- contig_dirs[grepl("D[0-9]+_contig_", contig_dirs)]
  
  if (length(contig_dirs) == 0) {
    cat(sprintf("  No contig directories found in %s\n", smorf_name))
    return(NULL)
  }
  
  # Read all genes from all contigs
  all_genes <- data.frame()
  for (contig_dir in contig_dirs) {
    gff_files <- list.files(contig_dir, pattern="\\.gff$", full.names=TRUE)
    
    if (length(gff_files) == 0) {
      next
    }
    
    contig_name <- sub("^D[0-9]+_", "", basename(contig_dir))
    genes <- read_gff(gff_files[1], contig_name)
    all_genes <- rbind(all_genes, genes)
  }
  
  if (nrow(all_genes) == 0) {
    cat(sprintf("  No genes found for %s\n", smorf_name))
    return(NULL)
  }
  
  # Force all strands to positive for visualization
  all_genes <- all_genes %>% mutate(strand = "+")
  
  # Create sequence data
  seqs_data <- all_genes %>%
    group_by(seq_id) %>%
    summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
    mutate(length = end - start)
  
  # Define colors
  unique_cogs <- unique(all_genes$name[!all_genes$is_target & all_genes$name != "Unknown"])
  vibrant_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                      "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
  
  cog_colors <- setNames(vibrant_colors[1:min(length(unique_cogs), length(vibrant_colors))], 
                         unique_cogs[1:min(length(unique_cogs), length(vibrant_colors))])
  cog_colors["Unknown"] <- "#D3D3D3"
  
  all_genes <- all_genes %>%
    mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))
  
  final_colors <- c("TARGET_sMORF" = "#000000", cog_colors)
  legend_breaks <- c("TARGET_sMORF", sort(unique_cogs), "Unknown")
  
  # Create plot
  p <- gggenomes(genes = all_genes, seqs = seqs_data)
  
  plot_result <- p +
    geom_seq() +
    geom_bin_label(size = 3) +
    geom_gene(aes(fill = gene_category)) +
    scale_fill_manual(values = final_colors, name = "Gene Annotation", breaks = legend_breaks) +
    ggtitle(smorf_name) +
    theme_minimal() +
    theme(legend.position = "right", legend.text = element_text(size = 9),
          legend.title = element_text(size = 11, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10), axis.text = element_text(size = 9),
          panel.grid.minor = element_blank())
  
  # Calculate plot dimensions
  plot_width <- 14
  plot_height <- max(6, nrow(seqs_data) * 0.5)
  
  # Save in smORF folder
  ggsave(file.path(smorf_path, "all_contigs_plot.jpg"), 
         plot = plot_result, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
  
  # Save in ALL_PLOTS folder with smORF name
  ggsave(file.path(all_plots_dir, paste0(smorf_name, "_plot.jpg")), 
         plot = plot_result, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
  
  cat(sprintf("  Created plots for %s (%d contigs, %d genes)\n", 
              smorf_name, nrow(seqs_data), nrow(all_genes)))
  
  return(list(contigs = nrow(seqs_data), genes = nrow(all_genes), targets = sum(all_genes$is_target)))
}

# Main execution
cat("=== Starting smORF Plot Generation ===\n\n")

# Get all smORF directories
smorf_dirs <- list.dirs(base_dir, recursive=FALSE, full.names=TRUE)
smorf_dirs <- smorf_dirs[grepl("SHD1_SM\\.100AA\\.", basename(smorf_dirs))]

cat(sprintf("Found %d smORF directories\n\n", length(smorf_dirs)))

# Process each smORF
total_contigs <- 0
total_genes <- 0
total_targets <- 0
successful_plots <- 0

for (i in seq_along(smorf_dirs)) {
  smorf_path <- smorf_dirs[i]
  smorf_name <- basename(smorf_path)
  
  cat(sprintf("[%d/%d] ", i, length(smorf_dirs)))
  
  result <- tryCatch({
    generate_smorf_plot(smorf_path, smorf_name)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    total_contigs <- total_contigs + result$contigs
    total_genes <- total_genes + result$genes
    total_targets <- total_targets + result$targets
    successful_plots <- successful_plots + 1
  }
  
  cat("\n")
}

# Summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("Processed: %d/%d smORFs\n", successful_plots, length(smorf_dirs)))
cat(sprintf("Total contigs plotted: %d\n", total_contigs))
cat(sprintf("Total genes plotted: %d\n", total_genes))
cat(sprintf("Total target sMORFs: %d\n", total_targets))
cat(sprintf("\nPlots saved in:\n"))
cat(sprintf("  - Individual folders: <smORF_folder>/all_contigs_plot.jpg\n"))
cat(sprintf("  - Collection folder: %s/\n", all_plots_dir))
cat("\nDone!\n")
