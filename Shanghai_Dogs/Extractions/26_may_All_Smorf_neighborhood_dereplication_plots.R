suppressPackageStartupMessages({
  library(gggenomes)
  library(dplyr)
  library(grid)
})

# =========================================================
# DEREPLICATED NEIGHBORHOOD PLOTS
# RUN FOR ALL smORFs IN DIRECTORY
# =========================================================

parent_path <- getwd()

all_plots_dir <- file.path(parent_path, "ALL_PLOTS")

dir.create(
  all_plots_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Parent:", parent_path, "\n")
cat("ALL_PLOTS:", all_plots_dir, "\n")

# =========================================================
# READ GFF
# =========================================================

read_gff <- function(gff_file, contig_name) {

  lines <- readLines(gff_file, warn = FALSE)

  lines <- lines[!grepl("^#", lines)]

  if (length(lines) == 0) return(NULL)

  df_list <- strsplit(lines, "\t", fixed = TRUE)

  max_cols <- max(vapply(df_list, length, integer(1)))

  df_mat <- do.call(
    rbind,
    lapply(df_list, function(x) {
      c(x, rep(NA, max_cols - length(x)))
    })
  )

  df <- as.data.frame(df_mat, stringsAsFactors = FALSE)

  colnames(df)[1:9] <- c(
    "orig_seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand_chr",
    "phase",
    "attributes"
  )

  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)

  df$feat_id <- ifelse(
    grepl("ID=", df$attributes),
    sub(".*ID=([^;]+).*", "\\1", df$attributes),
    "unknown_id"
  )

  df$name <- ifelse(
    grepl("Name=", df$attributes),
    sub(".*Name=([^;]+).*", "\\1", df$attributes),
    "Unknown"
  )

  df$is_target <- grepl("TARGET_sMORF", df$attributes)

  df$strand <- ifelse(df$strand_chr == "+", "+", "-")

  df$seq_id <- contig_name

  df %>%
    select(
      seq_id,
      start,
      end,
      strand,
      feat_id,
      type,
      name,
      is_target
    )
}

# =========================================================
# COLOR PALETTE
# =========================================================

cog_palette <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728",
  "#9467bd","#8c564b","#e377c2","#7f7f7f",
  "#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94",
  "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5"
)

# =========================================================
# FIND ALL smORF DIRECTORIES
# =========================================================

smorf_dirs <- list.dirs(
  parent_path,
  recursive = FALSE,
  full.names = TRUE
)

smorf_dirs <- smorf_dirs[
  grepl("^SHD1_SM\\.100AA\\.", basename(smorf_dirs))
]

cat("Found", length(smorf_dirs), "smORF directories\n")

if (length(smorf_dirs) == 0) {
  stop("No SHD1_SM.100AA.* directories found")
}

# =========================================================
# SUMMARY COUNTERS
# =========================================================

n_done <- 0
n_skipped_no_contigs <- 0
n_skipped_no_gff <- 0
n_skipped_no_genes <- 0

skipped <- list()

# =========================================================
# LOOP THROUGH ALL smORFs
# =========================================================

for (i in seq_along(smorf_dirs)) {

  base_path <- smorf_dirs[i]

  smorf_name <- basename(base_path)

  cat(
    "\n[",
    i,
    "/",
    length(smorf_dirs),
    "] Processing: ",
    smorf_name,
    "\n",
    sep = ""
  )

  # =========================================================
  # FIND CONTIG DIRECTORIES
  # =========================================================

  contig_dirs <- list.dirs(
    base_path,
    recursive = FALSE,
    full.names = TRUE
  )

  contig_dirs <- contig_dirs[
    grepl("^D[0-9]+_contig_", basename(contig_dirs))
  ]

  if (length(contig_dirs) == 0) {

    cat("  ⚠ No contig directories found\n")

    n_skipped_no_contigs <- n_skipped_no_contigs + 1

    skipped[[smorf_name]] <- "no_contigs"

    next
  }

  # =========================================================
  # READ ALL GFFS
  # =========================================================

  all_genes <- NULL

  gff_found_any <- FALSE

  for (d in contig_dirs) {

    gff <- list.files(
      d,
      pattern = "\\.gff$",
      full.names = TRUE
    )

    if (length(gff) == 0) next

    gff_found_any <- TRUE

    contig_name <- basename(d)

    genes <- read_gff(gff[1], contig_name)

    if (!is.null(genes) && nrow(genes) > 0) {

      all_genes <- rbind(all_genes, genes)
    }
  }

  if (!gff_found_any) {

    cat("  ⚠ No GFF files found\n")

    n_skipped_no_gff <- n_skipped_no_gff + 1

    skipped[[smorf_name]] <- "no_gff"

    next
  }

  if (is.null(all_genes) || nrow(all_genes) == 0) {

    cat("  ⚠ No genes parsed\n")

    n_skipped_no_genes <- n_skipped_no_genes + 1

    skipped[[smorf_name]] <- "no_genes"

    next
  }

  # =========================================================
  # CREATE SIGNATURES
  # =========================================================

  signature_df <- all_genes %>%
    arrange(seq_id, start) %>%
    mutate(
      gene_label = ifelse(
        is_target,
        paste0("TARGET_", name, "_", strand),
        paste0(name, "_", strand)
      )
    ) %>%
    group_by(seq_id) %>%
    summarise(
      signature = paste(gene_label, collapse = "|"),
      .groups = "drop"
    )

  # =========================================================
  # CLUSTER IDENTICAL NEIGHBORHOODS
  # =========================================================

  clusters <- signature_df %>%
    group_by(signature) %>%
    summarise(
      representative = first(seq_id),
      cluster_size = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(cluster_size))

  cat(
    "  Unique neighborhood clusters:",
    nrow(clusters),
    "\n"
  )

  # =========================================================
  # KEEP REPRESENTATIVES ONLY
  # =========================================================

  representatives <- clusters$representative

  plot_genes <- all_genes %>%
    filter(seq_id %in% representatives)

  # =========================================================
  # BUILD SEQ DATA
  # =========================================================

  seqs_data <- plot_genes %>%
    group_by(seq_id) %>%
    summarise(
      start = min(start) - 500,
      end   = max(end) + 500,
      .groups = "drop"
    ) %>%
    mutate(
      length = end - start
    )

  # =========================================================
  # ADD CLUSTER COUNTS TO LABELS
  # =========================================================

  cluster_sizes <- clusters %>%
    select(representative, cluster_size)

  seqs_data <- seqs_data %>%
    left_join(
      cluster_sizes,
      by = c("seq_id" = "representative")
    ) %>%
    mutate(
      seq_id = paste0(
        seq_id,
        " (n=",
        cluster_size,
        ")"
      )
    )

  plot_genes <- plot_genes %>%
    left_join(
      cluster_sizes,
      by = c("seq_id" = "representative")
    ) %>%
    mutate(
      seq_id = paste0(
        seq_id,
        " (n=",
        cluster_size,
        ")"
      )
    )

  # =========================================================
  # COLORS
  # =========================================================

  plot_genes <- plot_genes %>%
    mutate(
      gene_category = ifelse(
        is_target,
        paste0("TARGET (", name, ")"),
        name
      )
    )

  all_cats <- unique(plot_genes$gene_category)

  base_cats <- all_cats[
    !grepl("^TARGET", all_cats)
  ]

  cog_colors <- setNames(
    rep(cog_palette, length.out = length(base_cats)),
    base_cats
  )

  if ("Unknown" %in% names(cog_colors)) {
    cog_colors["Unknown"] <- "#D3D3D3"
  }

  target_cats <- all_cats[
    grepl("^TARGET", all_cats)
  ]

  if (length(target_cats) > 0) {
    cog_colors[target_cats] <- "#000000"
  }

  # =========================================================
  # PLOT
  # =========================================================

  p <- gggenomes(
    genes = plot_genes,
    seqs = seqs_data
  ) +
    geom_seq() +
    geom_bin_label(size = 2.5) +
    geom_gene(aes(fill = gene_category)) +
    scale_fill_manual(
      values = cog_colors,
      name = "Annotation"
    ) +
    ggtitle(smorf_name) +
    theme_minimal() +
    theme(
      legend.position = "right",

      legend.title = element_text(
        size = 10,
        face = "bold"
      ),

      legend.text = element_text(
        size = 8
      ),

      legend.key.size = unit(0.4, "cm"),

      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),

      axis.text.y = element_text(
        size = 7
      ),

      axis.text.x = element_text(
        size = 7
      )
    ) +
    guides(fill = guide_legend(ncol = 1))

  # =========================================================
  # SAVE
  # =========================================================

  local_plot_path <- file.path(
    base_path,
    paste0(smorf_name, ".jpg")
  )

  global_plot_path <- file.path(
    all_plots_dir,
    paste0(smorf_name, ".jpg")
  )

  ggsave(
    filename = local_plot_path,
    plot = p,
    width = 14,
    height = 10,
    dpi = 300,
    limitsize = FALSE
  )

  ggsave(
    filename = global_plot_path,
    plot = p,
    width = 14,
    height = 10,
    dpi = 300,
    limitsize = FALSE
  )

  cat("  ✓ Saved locally:", local_plot_path, "\n")
  cat("  ✓ Saved globally:", global_plot_path, "\n")

  n_done <- n_done + 1
}

# =========================================================
# SUMMARY
# =========================================================

cat("\n================ SUMMARY ================\n")

cat("Total smORF dirs found:", length(smorf_dirs), "\n")

cat("Successfully plotted:", n_done, "\n")

cat("Skipped (no contigs):", n_skipped_no_contigs, "\n")

cat("Skipped (no GFF):", n_skipped_no_gff, "\n")

cat("Skipped (no genes):", n_skipped_no_genes, "\n")

if (length(skipped) > 0) {

  cat("\nSkipped list:\n")

  for (nm in names(skipped)) {

    cat(
      "  -",
      nm,
      ":",
      skipped[[nm]],
      "\n"
    )
  }
}

cat("\nDone.\n")
