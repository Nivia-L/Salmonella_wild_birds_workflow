# ============================================================
# TABLE AND SNP DISTANCE MATRIX
# Salmonella wild birds
#
# Pairwise SNP distances calculated DIRECTLY from Data_PCA2.snp
# ============================================================

# ------------------------------------------------------------
# PROJECT SETUP
# ------------------------------------------------------------

setup_candidates <- c(
  file.path("R", "00_project_setup.R"),
  file.path("..", "R", "00_project_setup.R")
)

setup_file <- setup_candidates[file.exists(setup_candidates)][1]

if (is.na(setup_file)) {
  stop(
    "Cannot find R/00_project_setup.R. ",
    "Open Salmonella_wild_birds.Rproj before running."
  )
}

source(setup_file)


# ------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------

library(adegenet)
library(dplyr)
library(ggplot2)
library(reshape2)


# ------------------------------------------------------------
# 1. READ SNP MATRIX AND OFFICIAL METADATA
# ------------------------------------------------------------

# IMPORTANT:
# Distances are calculated directly from the 75-genotype,
# 442-locus biallelic SNP matrix.
#
# The previous version of this script used:
#   snp_tree_2023-09-04.nwk
#   cophenetic.phylo()
#   genome_size <- 4800000
#
# Those steps have been removed. Branch lengths are NOT SNP
# distances and are therefore not converted to SNP counts.

snp_data <- read.snp(
  data_file("Data_PCA2.snp")
)

G <- as.matrix(snp_data)
ids <- indNames(snp_data)

metadata <- read.csv(
  data_file("Supplementary_Table_S1_metadata.csv"),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. VERIFY SNP DATA
# ------------------------------------------------------------

cat("\n========================================\n")
cat("SNP DATA VERIFICATION\n")
cat("========================================\n")

cat("Number of genotypes:", nrow(G), "\n")
cat("Number of SNP loci:", ncol(G), "\n")
cat("Missing genotype values:", sum(is.na(G)), "\n")

observed_values <- sort(
  unique(
    as.vector(G[!is.na(G)])
  )
)

cat(
  "Observed genotype states:",
  paste(observed_values, collapse = ", "),
  "\n"
)

if (!all(observed_values %in% c(0, 1))) {
  stop(
    "Unexpected genotype coding. Expected only 0 and 1."
  )
}


# ------------------------------------------------------------
# 3. ALIGN METADATA TO SNP GENOTYPES
# ------------------------------------------------------------

source_vec <- metadata$Source[
  match(
    ids,
    metadata$Isolate_ID
  )
]

if (any(is.na(source_vec))) {

  cat("\nIDs without Source:\n")

  print(
    ids[is.na(source_vec)]
  )

  stop(
    "Some SNP isolate IDs could not be matched to Supplementary Table S1."
  )
}


# ------------------------------------------------------------
# 4. OFFICIAL SOURCE CATEGORIES
# ------------------------------------------------------------

source_levels <- c(
  "Broilers",
  "Human",
  "Layers",
  "Pigs",
  "Supermarket environment",
  "Wild Bird"
)

unexpected_sources <- setdiff(
  unique(source_vec),
  source_levels
)

if (length(unexpected_sources) > 0) {

  cat("\nUnexpected Source categories:\n")
  print(unexpected_sources)

  stop(
    "The metadata contains Source categories not included in the final analysis."
  )
}

cat("\n========================================\n")
cat("SOURCE COUNTS\n")
cat("========================================\n")

source_counts <- table(
  factor(
    source_vec,
    levels = source_levels
  )
)

print(source_counts)

expected_counts <- c(
  Broilers = 30,
  Human = 4,
  Layers = 17,
  Pigs = 8,
  `Supermarket environment` = 5,
  `Wild Bird` = 11
)

if (!all(
  as.integer(source_counts) ==
    as.integer(expected_counts[names(source_counts)])
)) {
  stop(
    "Source counts do not match the expected 75-genotype dataset."
  )
}


# ------------------------------------------------------------
# 5. VERIFY STUDY S. INFANTIS ISOLATES
# ------------------------------------------------------------

study_ids <- c(
  "CN-Y-11",
  "FC-Y-11",
  "U154s",
  "U168s"
)

study_source <- source_vec[
  match(
    study_ids,
    ids
  )
]

cat("\n========================================\n")
cat("STUDY S. INFANTIS ISOLATES\n")
cat("========================================\n")

print(
  data.frame(
    Isolate_ID = study_ids,
    Source = study_source,
    stringsAsFactors = FALSE
  )
)

if (any(is.na(study_source))) {
  stop(
    "One or more study S. Infantis isolate IDs are absent from Data_PCA2.snp."
  )
}

if (any(study_source != "Wild Bird")) {
  stop(
    "The four study S. Infantis isolates are not all classified as Wild Bird."
  )
}


# ------------------------------------------------------------
# 6. PAIRWISE SNP DISTANCE FUNCTION
# ------------------------------------------------------------
#
# For each isolate pair:
#
#   SNP distance = number of loci with different states
#   (0 versus 1).
#
# Loci that are missing (NA) in either isolate are excluded
# from that pairwise comparison.
#
# This is a direct Hamming-type distance across the biallelic
# SNP matrix and is reported as a number of differing SNP loci.
# ------------------------------------------------------------

pairwise_snp_distance <- function(a, b) {

  keep <- !is.na(a) & !is.na(b)

  if (!any(keep)) {
    return(
      list(
        distance = NA_real_,
        n_shared_loci = 0L
      )
    )
  }

  list(
    distance = sum(
      a[keep] != b[keep]
    ),
    n_shared_loci = sum(keep)
  )
}


# ------------------------------------------------------------
# 7. CALCULATE PAIRWISE SNP DISTANCE MATRIX
# ------------------------------------------------------------

D <- matrix(
  NA_real_,
  nrow = nrow(G),
  ncol = nrow(G),
  dimnames = list(ids, ids)
)

N_shared <- matrix(
  NA_integer_,
  nrow = nrow(G),
  ncol = nrow(G),
  dimnames = list(ids, ids)
)

for (i in seq_len(nrow(G))) {

  for (j in i:nrow(G)) {

    d <- pairwise_snp_distance(
      G[i, ],
      G[j, ]
    )

    D[i, j] <- d$distance
    D[j, i] <- d$distance

    N_shared[i, j] <- d$n_shared_loci
    N_shared[j, i] <- d$n_shared_loci
  }
}


# ------------------------------------------------------------
# 8. VERIFY PAIRWISE DISTANCES
# ------------------------------------------------------------

cat("\n========================================\n")
cat("PAIRWISE SNP DISTANCE SUMMARY\n")
cat("========================================\n")

print(
  summary(
    D[upper.tri(D)]
  )
)

cat("\nShared loci per pair:\n")

print(
  summary(
    N_shared[upper.tri(N_shared)]
  )
)


# ------------------------------------------------------------
# 9. SAVE INDIVIDUAL PAIRWISE DISTANCES
# ------------------------------------------------------------

pair_idx <- which(
  upper.tri(D),
  arr.ind = TRUE
)

pairwise_df <- data.frame(
  Isolate_1 = rownames(D)[pair_idx[, 1]],
  Isolate_2 = colnames(D)[pair_idx[, 2]],
  Source_1 = source_vec[pair_idx[, 1]],
  Source_2 = source_vec[pair_idx[, 2]],
  SNP_distance = D[pair_idx],
  Shared_loci = N_shared[pair_idx],
  stringsAsFactors = FALSE
)

write.csv(
  pairwise_df,
  table_file(
    "Pairwise_SNP_distances_direct.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. MEAN DISTANCE BETWEEN SOURCE GROUPS
# ------------------------------------------------------------

result <- matrix(
  NA_real_,
  nrow = length(source_levels),
  ncol = length(source_levels),
  dimnames = list(
    source_levels,
    source_levels
  )
)

n_comparisons <- matrix(
  0L,
  nrow = length(source_levels),
  ncol = length(source_levels),
  dimnames = list(
    source_levels,
    source_levels
  )
)

for (i in seq_along(source_levels)) {

  for (j in seq_along(source_levels)) {

    ids_i <- which(
      source_vec == source_levels[i]
    )

    ids_j <- which(
      source_vec == source_levels[j]
    )

    if (i == j) {

      # Within-source mean:
      # unique isolate pairs only; diagonal excluded.
      sub <- D[
        ids_i,
        ids_i,
        drop = FALSE
      ]

      vals <- sub[
        upper.tri(sub)
      ]

    } else {

      # Between-source mean:
      # all pairwise comparisons between the two groups.
      vals <- as.vector(
        D[
          ids_i,
          ids_j,
          drop = FALSE
        ]
      )
    }

    vals <- vals[
      !is.na(vals)
    ]

    if (length(vals) > 0) {

      result[i, j] <- mean(
        vals
      )

      n_comparisons[i, j] <- length(vals)
    }
  }
}


# ------------------------------------------------------------
# 11. FINAL MEAN SNP DISTANCE MATRIX
# ------------------------------------------------------------

result_rounded <- round(
  result,
  2
)

cat("\n========================================\n")
cat("MEAN PAIRWISE SNP DISTANCE MATRIX\n")
cat("========================================\n")

print(
  as.data.frame(result_rounded)
)

cat(
  "\nDiagonal values represent mean distances among unique\n",
  "pairs within each Source group.\n",
  sep = ""
)

cat(
  "Off-diagonal values represent mean distances across all\n",
  "pairwise comparisons between the two Source groups.\n"
)


# ------------------------------------------------------------
# 12. SAVE SQUARE SNP MATRIX
# ------------------------------------------------------------

result_snp_df <- as.data.frame(
  result_rounded
)

result_snp_df <- cbind(
  Source = rownames(result_snp_df),
  result_snp_df
)

write.csv(
  result_snp_df,
  table_file(
    "Table_SNP_distances_matrix.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. LONG TABLE OF UNIQUE SOURCE PAIRS
# ------------------------------------------------------------

dist_long <- data.frame(
  Source_1 = character(),
  Source_2 = character(),
  Mean_SNP_distance = numeric(),
  n_comparisons = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(source_levels)) {

  for (j in i:length(source_levels)) {

    if (is.na(result[i, j])) {
      next
    }

    dist_long <- rbind(
      dist_long,
      data.frame(
        Source_1 = source_levels[i],
        Source_2 = source_levels[j],
        Mean_SNP_distance = round(
          result[i, j],
          2
        ),
        n_comparisons = n_comparisons[i, j],
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n========================================\n")
cat("UNIQUE SOURCE-PAIR TABLE\n")
cat("========================================\n")

print(
  dist_long,
  row.names = FALSE
)

write.csv(
  dist_long,
  table_file(
    "Table_SNP_distances.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. HEATMAP
# ------------------------------------------------------------

heatmap_df <- melt(
  result_rounded,
  varnames = c(
    "Source_1",
    "Source_2"
  ),
  value.name = "Mean_SNP_distance"
)

heatmap_df$Source_1 <- factor(
  heatmap_df$Source_1,
  levels = rev(source_levels)
)

heatmap_df$Source_2 <- factor(
  heatmap_df$Source_2,
  levels = source_levels
)

p_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = Source_2,
    y = Source_1,
    fill = Mean_SNP_distance
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +

  geom_text(
    aes(
      label = ifelse(
        is.na(Mean_SNP_distance),
        "\u2014",
        sprintf(
          "%.2f",
          Mean_SNP_distance
        )
      )
    ),
    size = 4,
    color = "white"
  ) +

  scale_fill_viridis_c(
    name = "Mean SNP\ndistance",
    option = "viridis",
    limits = c(0, 40),
    breaks = c(0, 10, 20, 30, 40),
    oob = scales::squish,
    na.value = "grey90"
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),

    axis.text.y = element_text(
      size = 10
    ),

    axis.ticks = element_blank(),

    legend.title = element_text(
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ------------------------------------------------------------
# 15. EXPORT HEATMAP
# ------------------------------------------------------------

pdf(
  figure_file(
    "Fig_SNP_distances_heatmap.pdf"
  ),
  width = 18 / 2.54,
  height = 14 / 2.54
)

print(
  p_heatmap
)

dev.off()


jpeg(
  figure_file(
    "Fig_SNP_distances_heatmap.jpg"
  ),
  width = 18,
  height = 14,
  units = "cm",
  res = 600,
  bg = "white"
)

print(
  p_heatmap
)

dev.off()


# ------------------------------------------------------------
# 16. FINAL VERIFICATION
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FINAL VERIFICATION\n")
cat("========================================\n")

cat(
  "Genotypes analyzed: ",
  nrow(G),
  "\n",
  sep = ""
)

cat(
  "SNP loci analyzed: ",
  ncol(G),
  "\n",
  sep = ""
)

cat(
  "Source groups: ",
  length(source_levels),
  "\n",
  sep = ""
)

cat(
  "Source categories: ",
  paste(
    source_levels,
    collapse = " | "
  ),
  "\n",
  sep = ""
)

cat(
  "Wild Duck present: NO\n"
)

cat(
  "Distance source: Data_PCA2.snp\n"
)

cat(
  "Distance definition: number of differing 0/1 SNP states\n"
)

cat(
  "Missing loci: excluded separately for each isolate pair\n"
)

cat(
  "Tree branch lengths used: NO\n"
)

cat(
  "Genome-size scaling used: NO\n"
)

cat(
  "Expected matrix dimension: 6 x 6\n"
)

cat(
  "Actual matrix dimension: ",
  nrow(result),
  " x ",
  ncol(result),
  "\n",
  sep = ""
)

cat("\nExpected mean-distance matrix:\n")

print(
  round(
    result,
    2
  )
)

cat(
  "\nAnalysis completed successfully.\n"
)
