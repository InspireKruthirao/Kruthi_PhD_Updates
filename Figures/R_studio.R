# Load required libraries
library(gggenomes)
library(dplyr)

# Hardcoded data from your GFF
genes_data <- data.frame(
  seq_id = rep("contig_614", 9),
  start = c(130061, 130661, 131893, 133227, 133444, 134603, 134900, 135272, 138332),
  end = c(130369, 131560, 133038, 133385, 134529, 134845, 135211, 138319, 140941),
  strand = rep("+", 9),
  feat_id = paste0("contig_614_polypolish_", 157:165),
  type = rep("CDS", 9),
  name = c("COG2878", "Unknown", "COG1501", "Unknown", "COG1477", 
           "Unknown", "Unknown", "COG0587", "COG0749"),
  stringsAsFactors = FALSE
)

# Create a seq track (one row per sequence/contig)
# Using the actual range of your genes for better visualization
seqs_data <- data.frame(
  seq_id = "contig_614",
  start = 130000,  # Start just before first gene
  end = 141000,    # End just after last gene
  length = 11000,
  stringsAsFactors = FALSE
)

# Create the gggenomes plot
p <- gggenomes(
  genes = genes_data,
  seqs = seqs_data
)

# Inspect the tracks
print(p %>% track_info())

# Create the plot with better visibility
p +
  geom_seq() +  # Draw the sequence track
  geom_gene(aes(fill = name), size = 8) +  # Draw genes larger and colored by annotation
  scale_fill_manual(
    values = c("COG2878" = "steelblue", 
               "COG1501" = "coral", 
               "COG1477" = "darkgreen",
               "COG0587" = "purple",
               "COG0749" = "orange",
               "Unknown" = "grey70"),
    name = "Gene Annotation"
  ) +
  theme_minimal() +
  labs(title = "Contig 614 Gene Features",
       subtitle = "Region: 130,000 - 141,000 bp") +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )
