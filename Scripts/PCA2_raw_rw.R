# Project paths only; no analytical changes.
# ============================================================
# Project paths only; no analytical changes.
# ============================================================

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


# ============================================================
# FIGURA 2 – PCA
# PCA de los 75 genotipos de S. Infantis
# Los cuatro aislamientos propios se identifican como:
# "Wild Birds (This study)"
# ============================================================

library(adegenet)
library(ggplot2)
library(grid)


# ============================================================
# 1. CARGAR DATOS
# ============================================================

Data_PCA2 <- read.snp(
  data_file("Data_PCA2.snp")
)


# ============================================================
# 2. PCA
# ============================================================

pca_res <- glPca(
  Data_PCA2,
  nf = 3
)

var <- pca_res$eig /
  sum(pca_res$eig) * 100


# ============================================================
# 3. DATA FRAME DE LA PCA
#    Source procede de la S1 corregida
# ============================================================

pca_df <- as.data.frame(pca_res$scores)

pca_df$ID <- indNames(Data_PCA2)


# ============================================================
# 4. CARGAR METADATA OFICIAL DE S1
# ============================================================

metadata <- read.csv(
  data_file("Supplementary_Table_S1_metadata.csv"),
  stringsAsFactors = FALSE
)


# ============================================================
# 5. ASIGNAR SOURCE DESDE S1
# ============================================================

pca_df$Source <- metadata$Source[
  match(
    pca_df$ID,
    metadata$Isolate_ID
  )
]


# ============================================================
# 6. IDENTIFICAR LOS CUATRO AISLAMIENTOS PROPIOS
# ============================================================

this_study_ids <- c(
  "CN-Y-11",
  "FC-Y-11",
  "U154s",
  "U168s"
)


# ============================================================
# 7. VARIABLE STUDY
# ============================================================

pca_df$Study <- ifelse(
  pca_df$ID %in% this_study_ids,
  "Wild Birds (This study)",
  "Other isolates"
)


# ============================================================
# 8. VERIFICACIONES
# ============================================================

cat("\n========================================\n")
cat("PCA – VERIFICACIÓN DE METADATA\n")
cat("========================================\n")

cat("\nNúmero de genotipos:\n")
print(nrow(pca_df))

cat("\nConteo por Source:\n")
print(
  table(
    pca_df$Source,
    useNA = "ifany"
  )
)

cat("\nConteo por Study:\n")
print(
  table(
    pca_df$Study,
    useNA = "ifany"
  )
)

cat("\nAislamientos propios:\n")

print(
  pca_df[
    pca_df$ID %in% this_study_ids,
    c(
      "ID",
      "Source",
      "Study"
    )
  ]
)

cat("\nIDs sin Source:\n")

print(
  pca_df$ID[
    is.na(pca_df$Source)
  ]
)


# ============================================================
# 9. PALETA DE SOURCE
# ============================================================

source_palette_pca <- c(

  "Broilers" =
    "#1f78b4",

  "Human" =
    "#e31a1c",

  "Layers" =
    "#33a02c",

  "Pigs" =
    "#ff7f00",

  "Supermarket environment" =
    "#fdbf6f",

  "Wild Birds" =
    "green"

)


# ============================================================
# 10. GRAFICO PCA
# ============================================================

# Paleta para Source + aislamientos propios
pca_palette <- c(
  "Broilers" = "#1f78b4",
  "Human" = "#e31a1c",
  "Layers" = "#33a02c",
  "Pigs" = "#ff7f00",
  "Supermarket environment" = "#fdbf6f",
  "Wild Bird" = "grey50",
  "Wild Birds (This study)" = "green"
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +

  geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = 0,
    color = "grey70",
    linewidth = 0.4
  ) +

  # ----------------------------------------------------------
  # AISLAMIENTOS EXTERNOS
  # Color según Source
  # ----------------------------------------------------------

  geom_point(
    data = subset(
      pca_df,
      Study != "Wild Birds (This study)"
    ),
    aes(color = Source),
    size = 2.5,
    alpha = 0.85
  ) +

  # ----------------------------------------------------------
  # AISLAMIENTOS PROPIOS
  # Color y etiqueta propia
  # ----------------------------------------------------------

  geom_point(
    data = subset(
      pca_df,
      Study == "Wild Birds (This study)"
    ),
    aes(color = Study),
    shape = 17,
    size = 4
  ) +

  scale_color_manual(
    name = "Source",
    values = pca_palette,
    breaks = c(
      "Broilers",
      "Human",
      "Layers",
      "Pigs",
      "Supermarket environment",
      "Wild Bird",
      "Wild Birds (This study)"
    ),
    drop = FALSE
  ) +

  coord_equal() +

  theme_classic(
    base_size = 12
  ) +

  labs(
    x = paste0(
      "PC1 (",
      round(var[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var[2], 1),
      "%)"
    )
  ) +

  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),

    legend.position = "bottom",

    legend.title = element_text(
      face = "bold",
      size = 11
    ),

    legend.text = element_text(size = 10),

    legend.key.size = unit(
      0.45,
      "cm"
    ),

    legend.spacing.x = unit(
      0.3,
      "cm"
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

# ============================================================
# 11. MOSTRAR FIGURA
# ============================================================

print(p_pca)


# ============================================================

# 12. EXPORTAR FIGURA
# ============================================================
pdf(
  file = figure_file("Figure2_PCA.pdf"),
  width = 8,
  height = 6.5
)
print(p_pca)
dev.off()

png(
  filename = figure_file("Figure2_PCA.png"),
  width = 8,
  height = 6.5,
  units = "in",
  res = 600
)
print(p_pca)
dev.off()

jpeg(
  filename = figure_file("Figure2_PCA.jpg"),
  width = 8,
  height = 6.5,
  units = "in",
  res = 600,
  quality = 100
)
print(p_pca)
dev.off()
# ============================================================
# FIN DEL SCRIPT
# ============================================================
