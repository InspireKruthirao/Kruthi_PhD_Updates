suppressPackageStartupMessages({
  library(gggenomes)
  library(dplyr)
  library(grid)
})

# =========================================================
# DEREPLICATED NEIGHBORHOOD PLOT
# ONE smORF ONLY
# =========================================================

parent_path <- getwd()

smorf_name <- "SHD1_SM.100AA.003_480"

base_path <- file.path(parent_path, smorf_name)

# =========================================================
# OUTPUT PATHS
# =========================================================

all_plots_dir <- file.path(parent_path, "ALL_PLOTS")

dir.create(
  all_plots_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

local_plot_path <- file.path(
  base_path,
  paste0(smorf_name, ".jpg")
)

global_plot_path <- file.path(
  all_plots_dir,
  paste0(smorf_name, ".jpg")
)

cat("Processing:", smorf_name, "\n")

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

cat("Found", length(contig_dirs), "contigs\n")

# =========================================================
# READ ALL GFFS
# =========================================================

all_genes <- NULL

for (d in contig_dirs) {

  gff <- list.files(
    d,
    pattern = "\\.gff$",
    full.names = TRUE
  )

  if (length(gff) == 0) next

  contig_name <- basename(d)

  genes <- read_gff(gff[1], contig_name)

  if (!is.null(genes) && nrow(genes) > 0) {
    all_genes <- rbind(all_genes, genes)
  }
}

if (is.null(all_genes)) {
  stop("No genes found")
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

cat("Unique neighborhood clusters:",
    nrow(clusters), "\n")

# =========================================================
# KEEP ONLY REPRESENTATIVES
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

cog_palette <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728",
  "#9467bd","#8c564b","#e377c2","#7f7f7f",
  "#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94",
  "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5"
)

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
  )

# =========================================================
# SAVE TO BOTH LOCATIONS
# =========================================================

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

cat("\nSaved locally:\n")
cat(local_plot_path, "\n")

cat("\nSaved in ALL_PLOTS:\n")
cat(global_plot_path, "\n")

cat("\nDone.\n")
