# Load required libraries
library(gggenomes)
library(dplyr)
library(RColorBrewer)

# Set base path
base_path <- "Z:/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods/SHD1_SM.100AA.006_627"

# Function to read GFF file
read_gff <- function(gff_file, contig_name) {
  df <- read.table(gff_file, sep="\t", header=FALSE, comment.char="#",
                   col.names=c("seq_id", "source", "type", "start", "end", 
                              "score", "strand", "phase", "attributes"),
                   colClasses=c("character", "character", "character", 
                               "integer", "integer", "character", 
                               "character", "character", "character"))
  
  df$feat_id <- sub(".*ID=([^;]+).*", "\\1", df$attributes)
  df$name <- sub(".*Name=([^;]+).*", "\\1", df$attributes)
  df$seq_id <- contig_name
  
  df %>% select(seq_id, start, end, strand, feat_id, type, name)
}

# Find all contig directories
contig_dirs <- list.dirs(base_path, recursive=FALSE, full.names=TRUE)
contig_dirs <- contig_dirs[grepl("D[0-9]+_contig_", contig_dirs)]

# Read all GFF files
all_genes <- data.frame()

for (contig_dir in contig_dirs) {
  gff_files <- list.files(contig_dir, pattern="\\.gff$", full.names=TRUE)
  
  if (length(gff_files) > 0) {
    contig_name <- sub("^D[0-9]+_", "", basename(contig_dir))
    genes <- read_gff(gff_files[1], contig_name)
    all_genes <- rbind(all_genes, genes)
  }
}

# Create seqs track
seqs_data <- all_genes %>%
  group_by(seq_id) %>%
  summarise(
    start = min(start) - 500,
    end = max(end) + 500,
    .groups = "drop"
  ) %>%
  mutate(length = end - start)

# Define colors for COG categories
cog_colors <- c(
  "COG0237" = "#E41A1C", "COG0266" = "#377EB8", "COG0587" = "#4DAF4A",
  "COG0749" = "#984EA3", "COG1477" = "#FF7F00", "COG1501" = "#FFFF33",
  "COG2878" = "#A65628", "Unknown" = "#999999"
)

# Add extra colors if needed
unique_cogs <- unique(all_genes$name)
unique_cogs <- unique_cogs[unique_cogs != "Unknown"]

if (length(unique_cogs) > length(cog_colors) - 1) {
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_cogs))
  names(extra_colors) <- unique_cogs
  cog_colors <- c(cog_colors, extra_colors[!names(extra_colors) %in% names(cog_colors)])
}

# Create plot
p <- gggenomes(genes = all_genes, seqs = seqs_data)

plot_result <- p +
  geom_seq() +
  geom_bin_label(size = 3) +
  geom_gene(aes(fill = name)) +
  scale_fill_manual(values = cog_colors, name = "Gene Annotation") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.3, "cm"),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_blank()
  )

print(plot_result)

# Save plot
plot_width <- 12
plot_height <- max(6, nrow(seqs_data) * 0.4)

ggsave(file.path(base_path, "all_contigs_plot.pdf"), 
       plot = plot_result, width = plot_width, height = plot_height, limitsize = FALSE)

ggsave(file.path(base_path, "all_contigs_plot.png"), 
       plot = plot_result, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
