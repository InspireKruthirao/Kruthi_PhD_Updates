library(gggenomes)
library(dplyr)
library(RColorBrewer)

base_path <- "Z:/microbiome/shanghai_dogs/intermediate-outputs/SmORF_neighbourhoods/SHD1_SM.100AA.003_637"

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

contig_dirs <- list.dirs(base_path, recursive=FALSE, full.names=TRUE)
contig_dirs <- contig_dirs[grepl("D[0-9]+_contig_", contig_dirs)]

all_genes <- data.frame()
for (contig_dir in contig_dirs) {
  gff_files <- list.files(contig_dir, pattern="\\.gff$", full.names=TRUE)
  if (length(gff_files) > 0) {
    contig_name <- sub("^D[0-9]+_", "", basename(contig_dir))
    genes <- read_gff(gff_files[1], contig_name)
    all_genes <- rbind(all_genes, genes)
  }
}

# Force all genes to point right (positive strand)
all_genes <- all_genes %>%
  mutate(strand = "+")

seqs_data <- all_genes %>%
  group_by(seq_id) %>%
  summarise(start = min(start) - 500, end = max(end) + 500, .groups = "drop") %>%
  mutate(length = end - start)

unique_cogs <- unique(all_genes$name)
unique_cogs <- unique_cogs[unique_cogs != "Unknown"]

dark2_colors <- brewer.pal(8, "Dark2")
cog_colors <- setNames(dark2_colors[1:min(length(unique_cogs), 8)], 
                       unique_cogs[1:min(length(unique_cogs), 8)])
cog_colors["Unknown"] <- "#999999"

if (length(unique_cogs) > 8) {
  extra_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_cogs) - 8)
  names(extra_colors) <- unique_cogs[9:length(unique_cogs)]
  cog_colors <- c(cog_colors, extra_colors)
}

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

plot_width <- 12
plot_height <- max(6, nrow(seqs_data) * 0.4)

ggsave(file.path(base_path, "all_contigs_plot.jpg"), 
       plot = plot_result, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
