# =============================================================================
# 07_comparacion_cohortes.R – Integración multi-cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# Objetivo general del script
# ---------------------------
# En este script hacemos un análisis completo a nivel multi-cohorte. La idea es
# sintetizar y comparar los patrones de respuesta a distintos niveles:
#
#   PARTE A – Vías Hallmark (nivel funcional)
#     A1) I-SPY2 (ssGSEA Hallmark + limma) vs GSE241876 (fgsea Hallmark)
#     A2) I-SPY2 (ssGSEA) vs HTG_BrighTNess (ssGSEA)
#     A3) GSE241876 (fgsea) vs HTG_BrighTNess (ssGSEA)
#       Generamos tablas de vías comunes, anotamos significación y calculamos
#       correlaciones de efectos + scatterplots comparativos.
#
#   PARTE B – Comparación a nivel de gen (limma)
#     B1) I-SPY2 vs GSE241876 (formalizamos la comparación de logFC a nivel
#         de gen, que ya habíamos empezado a explorar de forma más dispersa).
#
#   PARTE B2 – UpSet de genes diferencialmente expresados (DEGs)
#       Construimos un diagrama de tipo UpSet para visualizar el solapamiento
#       de genes DE entre cohortes:
#            ISPY2, HTG_MEDI4736, HTG_BrighTNess, GSE241876.
#
#   PARTE C – Vista integrativa de la firma IFN-γ entre cohortes
#     C1) Resumen del rendimiento de IFN-γ en I-SPY2 (AUC ROC, salida del script 06).
#     C2) Resumen del comportamiento de IFN-γ en HTG_BrighTNess (logFC de vía).
#     C3) Resumen del nivel basal de IFN-γ en HTG_SCANB (distribución de scores).
#       Construimos una tabla sintética que resume cómo se comporta la “firma
#       IFN-γ” en los distintos contextos.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    install.packages("UpSetR")
  }
  library(UpSetR)})

# =============================================================================
# Rutas y funciones auxiliares
# =============================================================================

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

write_tsv_safe <- function(df, filename) {
  # Helper para guardar tablas en la carpeta de resultados de forma consistente:
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

save_plot <- function(p, filename, w = 7, h = 6, dpi = 300) {
  # Helper para guardar figuras con el mismo formato en todo el proyecto
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}

# Scatter genérico para comparación entre cohortes (lo reutilizamos varias veces):
plot_scatter_comparison <- function(df, x_var, y_var, col_var,
                                    x_lab, y_lab, title) {
  ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var), color = !!sym(col_var))) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = title,
      x     = x_lab,
      y     = y_lab,
      color = "Categoría"
    ) +
    theme_bw(base_size = 11)}

# =============================================================================
# BLOQUE 0 – Resumen global: nº de genes DE por cohorte
# (a partir de los resúmenes generados en 04 y 05)
# =============================================================================

resumen_files <- c(
  GSE173839_ISPY2 = "H4_ISPY2_resumen.tsv",
  HTG_MEDI4736    = "H4_HTG_MEDI4736_resumen.tsv",
  HTG_BrighTNess  = "H4_HTG_BrighTNess_resumen.tsv",
  GSE241876       = "H5_GSE241876_limma_resumen.tsv")

resumen_list <- lapply(names(resumen_files), function(coh) {
  f <- file.path(path_tabs, resumen_files[[coh]])
  if (!file.exists(f)) {
    warning("No se encontró el fichero de resumen para ", coh, ": ", f)
    return(NULL)}
  df <- readr::read_tsv(f, show_col_types = FALSE)
  
  # Nos aseguramos de que haya una columna 'cohort' uniforme
  if (!"cohort" %in% names(df)) {
    df$cohort <- coh}
  
  df})

# Filtramos las entradas NULL (cohortes para las que no hay resumen disponible):
resumen_list <- Filter(Negate(is.null), resumen_list)

if (length(resumen_list) > 0) {
  resumen_global_deg <- dplyr::bind_rows(resumen_list) %>%
    dplyr::select(cohort, dplyr::everything())
  
  write_tsv_safe(resumen_global_deg,
                 "H7_resumen_genes_DE_por_cohorte.tsv")
  message("Resumen global de genes DE por cohorte:")
  print(resumen_global_deg)
} else {
  message("No se encontró ningún resumen de genes DE para combinar.")}

# =============================================================================
# PARTE A – Comparación a nivel de vías Hallmark entre cohortes
# =============================================================================

message("=== PARTE A: Comparación de vías Hallmark entre cohortes ===")

# -----------------------------
# A1) I-SPY2 vs GSE241876
# -----------------------------
# Aquí comparamos los efectos a nivel de vía Hallmark entre I-SPY2 (ssGSEA+limma)
# y GSE241876 (fgsea). Para cada vía miramos logFC en I-SPY2 (score ssGSEA R vs NR)
# y NES en GSE241876 (fgsea R vs NR); y anotamos qué vías son significativas en
# cada cohorte:

# I-SPY2: H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv:
file_ispy_path <- file.path(path_tabs,
                            "H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv")
stopifnot(file.exists(file_ispy_path))

ispy2_path <- readr::read_tsv(file_ispy_path, show_col_types = FALSE)
if (!("pathway" %in% names(ispy2_path))) {
  stop("El archivo H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv no tiene columna 'pathway'.")}

ispy2_clean <- ispy2_path %>%
  dplyr::select(pathway, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val)

# GSE241876: H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv:
file_gse_path <- file.path(path_tabs,
                           "H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv")
stopifnot(file.exists(file_gse_path))

gse_path <- readr::read_tsv(file_gse_path, show_col_types = FALSE)
if (!("pathway" %in% names(gse_path))) {
  stop("El archivo H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv no tiene columna 'pathway'.")}

gse_clean <- gse_path %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(
    NES_GSE = NES,
    FDR_GSE = padj)

# Unimos por "pathway" para quedarnos sólo con las vías comunes:
path_ISPY2_vs_GSE <- dplyr::inner_join(
  ispy2_clean, gse_clean, by = "pathway"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_ambas = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_ambas              ~ "Significativa en ambas",
      sig_ISPY2 & !sig_GSE   ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE   ~ "Solo GSE241876",
      TRUE                   ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo I-SPY2",
        "Solo GSE241876",
        "No significativa")))

cor_paths_all <- suppressWarnings(
  cor(path_ISPY2_vs_GSE$logFC_ISPY2,
      path_ISPY2_vs_GSE$NES_GSE,
      method = "spearman",
      use    = "complete.obs"))

summary_paths_ISPY2_GSE <- tibble::tibble(
  comparacion               = "ISPY2_vs_GSE241876",
  n_pathways_comunes        = nrow(path_ISPY2_vs_GSE),
  n_sig_ISPY2               = sum(path_ISPY2_vs_GSE$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE                 = sum(path_ISPY2_vs_GSE$sig_GSE,   na.rm = TRUE),
  n_sig_ambas               = sum(path_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE),
  cor_spearman_logFC_vs_NES = cor_paths_all)

write_tsv_safe(path_ISPY2_vs_GSE,
               "H7_pathways_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_paths_ISPY2_GSE,
               "H7_pathways_ISPY2_vs_GSE241876_resumen.tsv")

p_paths_ISPY2_GSE <- plot_scatter_comparison(
  df      = path_ISPY2_vs_GSE,
  x_var   = "logFC_ISPY2",
  y_var   = "NES_GSE",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR) en ssGSEA (Hallmark)",
  y_lab   = "GSE241876 – NES (fgsea, Hallmark)",
  title   = "Vías Hallmark: I-SPY2 vs GSE241876")

save_plot(p_paths_ISPY2_GSE,
          "H7_ISPY2_vs_GSE241876_Hallmark_scatter.png")

# -----------------------------
# A2) I-SPY2 vs HTG_BrighTNess
# -----------------------------
# Repetimos el mismo esquema de comparación pero ahora entre I-SPY2 y
# HTG_BrighTNess, usando logFC de vías en ambos casos (ssGSEA+limma):

file_brtn_path <- file.path(path_tabs,
                            "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv")
stopifnot(file.exists(file_brtn_path))

brtn_path <- readr::read_tsv(file_brtn_path, show_col_types = FALSE)
if (!("pathway" %in% names(brtn_path))) {
  stop("El archivo H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv no tiene columna 'pathway'.")}

brtn_clean <- brtn_path %>%
  dplyr::select(pathway, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_BRTN = logFC,
    FDR_BRTN   = adj.P.Val)

path_ISPY2_vs_BRTN <- dplyr::inner_join(
  ispy2_clean, brtn_clean, by = "pathway"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_BRTN  = FDR_BRTN  < 0.05,
    sig_ambas = sig_ISPY2 & sig_BRTN,
    categoria = dplyr::case_when(
      sig_ambas              ~ "Significativa en ambas",
      sig_ISPY2 & !sig_BRTN  ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_BRTN  ~ "Solo BrighTNess",
      TRUE                   ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo I-SPY2",
        "Solo BrighTNess",
        "No significativa")))

cor_paths_ISPY2_BRTN <- suppressWarnings(
  cor(path_ISPY2_vs_BRTN$logFC_ISPY2,
      path_ISPY2_vs_BRTN$logFC_BRTN,
      method = "spearman",
      use    = "complete.obs"))

summary_paths_ISPY2_BRTN <- tibble::tibble(
  comparacion                  = "ISPY2_vs_HTG_BrighTNess",
  n_pathways_comunes           = nrow(path_ISPY2_vs_BRTN),
  n_sig_ISPY2                  = sum(path_ISPY2_vs_BRTN$sig_ISPY2, na.rm = TRUE),
  n_sig_BRTN                   = sum(path_ISPY2_vs_BRTN$sig_BRTN,   na.rm = TRUE),
  n_sig_ambas                  = sum(path_ISPY2_vs_BRTN$sig_ambas,  na.rm = TRUE),
  cor_spearman_logFC_vs_logFC  = cor_paths_ISPY2_BRTN)

write_tsv_safe(path_ISPY2_vs_BRTN,
               "H7_pathways_ISPY2_vs_HTG_BrighTNess_merged.tsv")
write_tsv_safe(summary_paths_ISPY2_BRTN,
               "H7_pathways_ISPY2_vs_HTG_BrighTNess_resumen.tsv")

p_paths_ISPY2_BRTN <- plot_scatter_comparison(
  df      = path_ISPY2_vs_BRTN,
  x_var   = "logFC_ISPY2",
  y_var   = "logFC_BRTN",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR) en vías Hallmark",
  y_lab   = "HTG_BrighTNess – log2(FC) (R vs NR) en vías Hallmark",
  title   = "Vías Hallmark: I-SPY2 vs HTG_BrighTNess")

save_plot(p_paths_ISPY2_BRTN,
          "H7_ISPY2_vs_HTG_BrighTNess_Hallmark_scatter.png")

# -----------------------------
# A3) GSE241876 vs HTG_BrighTNess
# -----------------------------
# Por último, cruzamos GSE241876 (NES de fgsea) con HTG_BrighTNess (logFC).
# De nuevo anotamos la significación en cada cohorte y calculamos la
# correlación entre NES y logFC:

path_GSE_vs_BRTN <- dplyr::inner_join(
  gse_clean, brtn_clean, by = "pathway"
) %>%
  dplyr::mutate(
    sig_GSE   = FDR_GSE   < 0.05,
    sig_BRTN  = FDR_BRTN  < 0.05,
    sig_ambas = sig_GSE & sig_BRTN,
    categoria = dplyr::case_when(
      sig_ambas              ~ "Significativa en ambas",
      sig_GSE & !sig_BRTN    ~ "Solo GSE241876",
      !sig_GSE & sig_BRTN    ~ "Solo BrighTNess",
      TRUE                   ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo GSE241876",
        "Solo BrighTNess",
        "No significativa")))

cor_paths_GSE_BRTN <- suppressWarnings(
  cor(path_GSE_vs_BRTN$NES_GSE,
      path_GSE_vs_BRTN$logFC_BRTN,
      method = "spearman",
      use    = "complete.obs"))

summary_paths_GSE_BRTN <- tibble::tibble(
  comparacion               = "GSE241876_vs_HTG_BrighTNess",
  n_pathways_comunes        = nrow(path_GSE_vs_BRTN),
  n_sig_GSE                 = sum(path_GSE_vs_BRTN$sig_GSE,  na.rm = TRUE),
  n_sig_BRTN                = sum(path_GSE_vs_BRTN$sig_BRTN, na.rm = TRUE),
  n_sig_ambas               = sum(path_GSE_vs_BRTN$sig_ambas, na.rm = TRUE),
  cor_spearman_NES_vs_logFC = cor_paths_GSE_BRTN)

write_tsv_safe(path_GSE_vs_BRTN,
               "H7_pathways_GSE241876_vs_HTG_BrighTNess_merged.tsv")
write_tsv_safe(summary_paths_GSE_BRTN,
               "H7_pathways_GSE241876_vs_HTG_BrighTNess_resumen.tsv")

p_paths_GSE_BRTN <- plot_scatter_comparison(
  df      = path_GSE_vs_BRTN,
  x_var   = "NES_GSE",
  y_var   = "logFC_BRTN",
  col_var = "categoria",
  x_lab   = "GSE241876 – NES (fgsea, Hallmark)",
  y_lab   = "HTG_BrighTNess – log2(FC) (R vs NR) en vías Hallmark",
  title   = "Vías Hallmark: GSE241876 vs HTG_BrighTNess")

save_plot(p_paths_GSE_BRTN,
          "H7_GSE241876_vs_HTG_BrighTNess_Hallmark_scatter.png")

# =============================================================================
# PARTE B – Comparación a nivel de gen (limma) entre I-SPY2 y GSE241876
# =============================================================================

message("=== PARTE B: Comparación a nivel de gen (I-SPY2 vs GSE241876) ===")

# -----------------------------
# B1) I-SPY2 vs GSE241876 (genes)
# -----------------------------
# En este bloque cruzamos directamente los resultados limma a nivel de gen
# entre I-SPY2 y GSE241876, y calculamos el número de genes DE por cohorte, el
# número de genes significativos en ambas, y la correlación de logFC y restringida
# a genes significativos en ambas:


# I-SPY2: limma R vs NR:
file_ispy_genes <- file.path(path_tabs, "H4_ISPY2_limma_R_vs_NR.tsv")
if (!file.exists(file_ispy_genes)) {
  # fallback por si el nombre es H3_...
  file_ispy_genes_alt <- file.path(path_tabs, "H3_ISPY2_limma_R_vs_NR.tsv")
  if (file.exists(file_ispy_genes_alt)) {
    file_ispy_genes <- file_ispy_genes_alt
  } else {
    stop("No se encontró ni H4_ISPY2_limma_R_vs_NR.tsv ni H3_ISPY2_limma_R_vs_NR.tsv.")}}

ispy2_genes <- readr::read_tsv(file_ispy_genes, show_col_types = FALSE)
if (!("Gene" %in% names(ispy2_genes))) {
  ispy2_genes <- ispy2_genes %>% tibble::rownames_to_column("Gene")}

ispy2_genes_clean <- ispy2_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val)

# GSE241876: limma-voom R vs NR:
file_gse_genes <- file.path(path_tabs, "H5_GSE241876_limma_R_vs_NR.tsv")
stopifnot(file.exists(file_gse_genes))

gse_genes <- readr::read_tsv(file_gse_genes, show_col_types = FALSE)
if (!("Gene" %in% names(gse_genes))) {
  gse_genes <- gse_genes %>% tibble::rownames_to_column("Gene")}

gse_genes_clean <- gse_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_GSE = logFC,
    FDR_GSE   = adj.P.Val)

genes_ISPY2_vs_GSE <- dplyr::inner_join(
  ispy2_genes_clean, gse_genes_clean, by = "Gene"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_ambas = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_ambas              ~ "Significativo en ambas",
      sig_ISPY2 & !sig_GSE   ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE   ~ "Solo GSE241876",
      TRUE                   ~ "No significativo"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativo en ambas",
        "Solo I-SPY2",
        "Solo GSE241876",
        "No significativo")))

cor_genes_all <- suppressWarnings(
  cor(genes_ISPY2_vs_GSE$logFC_ISPY2,
      genes_ISPY2_vs_GSE$logFC_GSE,
      method = "spearman",
      use    = "complete.obs"))

if (any(genes_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE)) {
  df_sig_both <- genes_ISPY2_vs_GSE %>%
    dplyr::filter(sig_ambas,
                  is.finite(logFC_ISPY2),
                  is.finite(logFC_GSE))
  if (nrow(df_sig_both) > 1) {
    cor_genes_sig_both <- suppressWarnings(
      cor(df_sig_both$logFC_ISPY2,
          df_sig_both$logFC_GSE,
          method = "spearman",
          use    = "complete.obs"))
  } else {
    cor_genes_sig_both <- NA_real_}
} else {
  cor_genes_sig_both <- NA_real_}

summary_genes_ISPY2_GSE <- tibble::tibble(
  comparacion                 = "ISPY2_vs_GSE241876",
  n_genes_comunes             = nrow(genes_ISPY2_vs_GSE),
  n_sig_ISPY2                 = sum(genes_ISPY2_vs_GSE$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE                   = sum(genes_ISPY2_vs_GSE$sig_GSE,   na.rm = TRUE),
  n_sig_ambas                 = sum(genes_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE),
  cor_spearman_logFC_global   = cor_genes_all,
  cor_spearman_logFC_sig_both = cor_genes_sig_both)

write_tsv_safe(genes_ISPY2_vs_GSE,
               "H7_genes_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_genes_ISPY2_GSE,
               "H7_genes_ISPY2_vs_GSE241876_resumen.tsv")

genes_plot <- genes_ISPY2_vs_GSE %>%
  dplyr::filter(is.finite(logFC_ISPY2), is.finite(logFC_GSE))

p_genes_ISPY2_GSE <- plot_scatter_comparison(
  df      = genes_plot,
  x_var   = "logFC_ISPY2",
  y_var   = "logFC_GSE",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR)",
  y_lab   = "GSE241876 – log2(FC) (R vs NR)",
  title   = "Efectos a nivel de gen: I-SPY2 vs GSE241876")

save_plot(p_genes_ISPY2_GSE,
          "H7_ISPY2_vs_GSE241876_genes_logFC_scatter.png")

# =============================================================================
# PARTE B2 – UpSet de solapamiento de genes diferencialmente expresados
# =============================================================================

message("=== PARTE B2: UpSet de genes diferencialmente expresados (DEGs) ===")

# Definimos los umbrales de significación que vamos a usar como "DE":
deg_threshold_fdr   <- 0.05
deg_threshold_logfc <- log2(1.5)

deg_files <- c(
  ISPY2          = "H4_ISPY2_limma_R_vs_NR.tsv",
  HTG_MEDI4736   = "H4_HTG_MEDI4736_limma_R_vs_NR.tsv",
  HTG_BrighTNess = "H4_HTG_BrighTNess_limma_R_vs_NR.tsv",
  GSE241876      = "H5_GSE241876_limma_R_vs_NR.tsv")

deg_sets <- lapply(names(deg_files), function(coh) {
  f <- file.path(path_tabs, deg_files[[coh]])
  if (!file.exists(f)) {
    warning("No se encontró el fichero limma para ", coh, ": ", f)
    return(NULL)}
  
  tt <- readr::read_tsv(f, show_col_types = FALSE)
  
  # Aseguramos que exista columna 'Gene':
  if (!("Gene" %in% colnames(tt))) {
    tt <- tt %>% tibble::rownames_to_column("Gene")}
  
  # Filtramos genes DE según FDR y |logFC|:
  sig_genes <- tt %>%
    dplyr::filter(
      !is.na(adj.P.Val),
      adj.P.Val < deg_threshold_fdr,
      !is.na(logFC),
      abs(logFC) >= deg_threshold_logfc
    ) %>%
    dplyr::pull(Gene) %>%
    unique()
  
  message("Cohorte ", coh, ": ", length(sig_genes),
          " genes DE (FDR < ", deg_threshold_fdr,
          ", |logFC| >= ", round(deg_threshold_logfc, 3), ").")
  
  sig_genes})

names(deg_sets) <- names(deg_files)

# Eliminamos cohortes sin genes DE o sin fichero:
deg_sets <- Filter(function(x) !is.null(x) && length(x) > 0, deg_sets)

if (length(deg_sets) >= 2) {
  # Construimos una tabla booleana de pertenencia (Gene x cohorte)
  all_genes <- sort(unique(unlist(deg_sets)))
  membership_mat <- sapply(deg_sets, function(s) all_genes %in% s)
  membership_tbl <- membership_mat %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "Gene")
  
  # Guardamos la tabla de pertenencia por cohorte:
  write_tsv_safe(membership_tbl,
                 "H7_DEGs_membership_por_cohorte.tsv")
  
  # UpSet plot basado en los conjuntos de genes DE de cada cohorte:
  outfile <- file.path(path_figs, "H7_DEGs_UpSet_genes.png")
  png(outfile, width = 2200, height = 1400, res = 300)
  grid::grid.newpage()
  
  UpSetR::upset(
    upset_input,
    nsets   = length(deg_sets),
    order.by = "freq")
  
  dev.off()
  message("Figura UpSet de genes DE guardada en: ", outfile)
  
} else {
  message("No hay suficientes cohortes con genes DE significativos para construir el UpSet.")}

# =============================================================================
# PARTE C – Vista integrativa de la firma IFN-γ entre cohortes
# =============================================================================

message("=== PARTE C: Vista integrativa de IFN-γ (Hallmark) entre cohortes ===")

# I-SPY2: AUC IFN-γ (salida del script 06_validacion_sensibilidad.R)
file_ifng_ispy_auc <- file.path(path_tabs, "H7_ISPY2_IFNG_AUC.tsv")
if (file.exists(file_ifng_ispy_auc)) {
  ifng_ispy_auc <- readr::read_tsv(file_ifng_ispy_auc, show_col_types = FALSE)
} else {
  warning("No se encontró H7_ISPY2_IFNG_AUC.tsv; se dejará NA para AUC IFN-γ en I-SPY2.")
  ifng_ispy_auc <- tibble::tibble(
    marker = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    AUC    = NA_real_,
    cohort = "GSE173839_ISPY2",
    n_R    = NA_real_,
    n_NR   = NA_real_)}

# HTG_BrighTNess: efecto IFN-γ a nivel de vía (logFC R vs NR en Hallmark IFN-γ):
path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
ifng_brtn_row <- brtn_clean %>%
  dplyr::filter(pathway == path_ifng)

if (nrow(ifng_brtn_row) == 0) {
  warning("En HTG_BrighTNess no se encontró la vía IFN-γ en la tabla de vías.")
  ifng_brtn_summary <- tibble::tibble(
    cohort  = "HTG_BrighTNess",
    pathway = path_ifng,
    logFC   = NA_real_,
    FDR     = NA_real_)
} else {
  ifng_brtn_summary <- tibble::tibble(
    cohort  = "HTG_BrighTNess",
    pathway = path_ifng,
    logFC   = ifng_brtn_row$logFC_BRTN,
    FDR     = ifng_brtn_row$FDR_BRTN)}

# HTG_SCANB: resumen descriptivo de scores IFN-γ (nivel basal en una cohorte grande):
file_ifng_scanb <- file.path(path_tabs, "H5_HTG_SCANB_IFNG_score_resumen.tsv")
if (file.exists(file_ifng_scanb)) {
  ifng_scanb_summary <- readr::read_tsv(file_ifng_scanb, show_col_types = FALSE)
} else {
  warning("No se encontró H5_HTG_SCANB_IFNG_score_resumen.tsv; se dejarán NA para SCAN-B.")
  ifng_scanb_summary <- tibble::tibble(
    cohort     = "HTG_SCANB",
    pathway    = path_ifng,
    n_samples  = NA_real_,
    mean_score = NA_real_,
    sd_score   = NA_real_,
    q25        = NA_real_,
    q50        = NA_real_,
    q75        = NA_real_)}

# Construimos una tabla integrativa que resume el "comportamiento IFN-γ"
# en las distintas cohortes, con el tipo de medida y algo de contexto extra:
ifng_integrativo <- tibble::tibble(
  cohorte        = c("GSE173839_ISPY2", "HTG_BrighTNess", "HTG_SCANB"),
  tipo_medida    = c(
    "AUC ROC (score vía IFN-γ)",
    "log2(FC) R vs NR (vía IFN-γ)",
    "Score basal ssGSEA IFN-γ"
  ),
  valor_principal = c(
    ifng_ispy_auc$AUC[1],
    ifng_brtn_summary$logFC[1],
    ifng_scanb_summary$mean_score[1]
  ),
  info_adicional = c(
    paste0("n_R = ", ifng_ispy_auc$n_R[1],
           ", n_NR = ", ifng_ispy_auc$n_NR[1]),
    paste0(
      "FDR = ",
      if (!is.na(ifng_brtn_summary$FDR[1]))
        signif(ifng_brtn_summary$FDR[1], 3) else "NA"
    ),
    paste0(
      "n = ", ifng_scanb_summary$n_samples[1],
      "; q25/q50/q75 = ",
      if (!is.na(ifng_scanb_summary$q25[1]))
        paste0(
          signif(ifng_scanb_summary$q25[1], 3), "/",
          signif(ifng_scanb_summary$q50[1], 3), "/",
          signif(ifng_scanb_summary$q75[1], 3)
        ) else "NA")))

write_tsv_safe(ifng_integrativo,
               "H7_IFNG_integrativo_cohortes.tsv")
print(ifng_integrativo)

message("=== Script 07 (comparación de cohortes: vías + genes + IFN-γ + UpSet) COMPLETADO ===")
