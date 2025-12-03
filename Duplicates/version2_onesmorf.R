library(gggenomes)

library(dplyr)

# --------------------------------------------------------------

# Your folder

# --------------------------------------------------------------

base_path <- "Z:/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods_25_30/SHD1_SM.100AA.003_615"

# --------------------------------------------------------------

# Robust GFF reader

# --------------------------------------------------------------

read_gff <- function(gff_file, contig_name) {

  lines <- readLines(gff_file)

  lines <- lines[!grepl("^#", lines)]

  if (length(lines) == 0) return(NULL)

  

  df_list <- strsplit(lines, "\t")

  max_cols <- max(sapply(df_list, length))

  df_mat  <- do.call(rbind, lapply(df_list, function(x) c(x, rep(NA, max_cols - length(x)))))

  df <- as.data.frame(df_mat, stringsAsFactors = FALSE)

  

  colnames(df)[1:9] <- c("orig_seqid", "source", "type", "start", "end",

                         "score", "strand", "phase", "attributes")

  

  df$start <- as.integer(df$start)

  df$end   <- as.integer(df$end)

  

  df$feat_id   <- sub(".*ID=([^;]+).*", "\\1", df$attributes)

  df$name      <- sub(".*Name=([^;]+).*", "\\1", df$attributes)

  df$is_target <- grepl("TARGET_sMORF", df$attributes)

  df$seq_id    <- contig_name

  df$strand    <- 1  # ALL POINT RIGHT

  

  df %>% select(seq_id, start, end, strand, feat_id, type, name, is_target)

}

# --------------------------------------------------------------

# Read all contigs

# --------------------------------------------------------------

contig_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)

contig_dirs <- contig_dirs[grepl("^D[0-9]+_contig_", basename(contig_dirs))]

all_genes <- NULL

for (d in contig_dirs) {

  gff <- list.files(d, pattern = "\\.gff$", full.names = TRUE)

  if (length(gff) == 0) next

  contig_name <- sub("^D[0-9]+_", "", basename(d))

  genes <- read_gff(gff[1], contig_name)

  if (!is.null(genes)) all_genes <- rbind(all_genes, genes)

}

cat(sprintf("Read %d genes from %d contigs\n", nrow(all_genes), length(unique(all_genes$seq_id))))

# --------------------------------------------------------------

# Sequence info

# --------------------------------------------------------------

seqs_data <- all_genes %>%

  group_by(seq_id) %>%

  summarise(start = min(start) - 500,

            end   = max(end)   + 500, .groups = "drop") %>%

  mutate(length = end - start)

# --------------------------------------------------------------

# Color palette

# --------------------------------------------------------------

cog_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f",

                 "#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5","#c49c94",

                 "#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")

all_cogs <- unique(all_genes$name)

cog_colors <- setNames(rep(cog_palette, length.out = length(all_cogs)), all_cogs)

cog_colors["Unknown"] <- "#D3D3D3"

all_genes <- all_genes %>%

  mutate(gene_category = ifelse(is_target, "TARGET_sMORF", name))

cog_colors["TARGET_sMORF"] <- "#000000"

# --------------------------------------------------------------

# PLOT WITH CUSTOMIZABLE SIZE & SMALLER LEGEND

# --------------------------------------------------------------

p <- gggenomes(genes = all_genes, seqs = seqs_data) +

  geom_seq() +

  geom_bin_label(size = 2) +  # Smaller bin labels

  geom_gene(aes(fill = gene_category)) +

  scale_fill_manual(values = cog_colors, name = "Annotation") +

  ggtitle(basename(base_path)) +

  theme_minimal() +

  theme(

    legend.position = "right",

    legend.title = element_text(size = 10, face = "bold"),      # Smaller legend title

    legend.text = element_text(size = 8),                      # Smaller legend text

    legend.key.size = unit(0.4, "cm"),                         # Smaller legend keys

    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),

    plot.margin = margin(10, 10, 10, 40),                      # Margins: top, right, bottom, left

    axis.text.y = element_text(size = 7),                      # Smaller y-axis labels

    axis.text.x = element_text(size = 7)                       # Smaller x-axis labels

  ) +

  guides(fill = guide_legend(ncol = 1))  # Single column legend

# Show

print(p)

# --------------------------------------------------------------

# SAVE WITH CUSTOM DIMENSIONS

# --------------------------------------------------------------

# Adjust these values to change plot size:

plot_width <- 12         # Width in inches (change this)

plot_height <- 10        # Height in inches (change this)

plot_dpi <- 300          # Resolution (300 is publication quality)

ggsave(

  file.path(base_path, "all_contigs_plot_CUSTOM.jpg"),

  plot = p, 

  width = plot_width, 

  height = plot_height, 

  dpi = plot_dpi, 

  limitsize = FALSE

)

cat(sprintf("DONE! Plot saved as %d x %d inches with smaller legend\n", plot_width, plot_height))
