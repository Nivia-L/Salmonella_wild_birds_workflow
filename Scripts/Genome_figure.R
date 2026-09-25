# ============================================================
# FIGURE 3 — Vertical phylogenetic tree
# Salmonella Infantis wild birds
# ============================================================

# ------------------------------------------------------------
# 1. PROJECT SETUP
# ------------------------------------------------------------

setup_candidates <- c(
  file.path("R", "00_project_setup.R"),
  file.path("..", "R", "00_project_setup.R")
)

setup_file <- setup_candidates[file.exists(setup_candidates)][1]

if (is.na(setup_file)) {
  stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
}

source(setup_file)


# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------

library(ggtree)
library(ggplot2)
library(ape)
library(treeio)


# ------------------------------------------------------------
# 3. READ TREE AND METADATA
# ------------------------------------------------------------

tree <- read.tree(
  data_file("snp_tree_2023-09-04.nwk")
)

metadata <- read.csv(
  data_file("Supplementary_Table_S1_metadata.csv")
)

metadata$label <- metadata$Isolate_ID


# ------------------------------------------------------------
# 4. PREPARE TREE
# ------------------------------------------------------------

# Remove the outlier that distorts the tree scale
tree_final <- drop.tip(tree, "b")

# Order branches for clearer visualization
tree_final <- ladderize(tree_final)

# Prepare metadata for joining with the tree
metadata_clean <- metadata

names(metadata_clean)[10] <- "Label_Original"

metadata_clean <- metadata_clean[
  !is.na(metadata_clean$Source),
]


# ------------------------------------------------------------
# 5. TREE SCALE
# ------------------------------------------------------------

max_x <- max(
  fortify(tree_final)$x,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 6. SOURCE COLORS
# ------------------------------------------------------------

source_palette <- c(
  "Broilers" = "#1f78b4",
  "Human" = "#e31a1c",
  "Layers" = "#33a02c",
  "Pigs" = "#ff7f00",
  "Supermarket environment" = "#fdbf6f",
  "Wild Bird" = "green"
)


# ------------------------------------------------------------
# 7. FIGURE 3 — VERTICAL TREE
# ------------------------------------------------------------

p <- ggtree(
  tree_final,
  layout = "rectangular",
  linewidth = 0.6,
  color = "grey40"
) %<+% metadata_clean +

  geom_tippoint(
    aes(color = Source),
    size = 3,
    alpha = 0.8
  ) +

  geom_tiplab(
    size = 3
  ) +

  scale_color_manual(
    values = source_palette,
    na.translate = FALSE
  ) +

  theme_tree2() +

  theme(
    legend.position = "right",

    legend.title = element_text(
      face = "bold",
      size = 12
    ),

    legend.text = element_text(
      family = "serif",
      size = 12
    ),

    axis.title.x = element_text(
      face = "bold",
      size = 12
    ),

    plot.margin = margin(
      10, 10, 10, 10
    )
  ) +

  labs(
    x = "Genetic distance"
  ) +

  xlim(
    0,
    max_x * 1.55
  )


# ------------------------------------------------------------
# 8. DISPLAY FIGURE
# ------------------------------------------------------------

print(p)


# ------------------------------------------------------------
# 9. EXPORT FIGURE 3
# ------------------------------------------------------------

# ------------------------------------------------------------
# 9. EXPORT FIGURE 3
# ------------------------------------------------------------

# PNG
ggsave(
  figure_file("Salmonella_fig3.png"),
  plot = p,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# PDF
ggsave(
  figure_file("Salmonella_fig3.pdf"),
  plot = p,
  width = 14,
  height = 10,
  units = "in",
  device = "pdf",
  bg = "white"
)

# JPG
ggsave(
  figure_file("Salmonella_fig3.jpg"),
  plot = p,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
