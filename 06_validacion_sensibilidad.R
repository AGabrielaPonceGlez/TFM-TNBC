# =============================================================================
# 06_validacion_sensibilidad.R – Validación y análisis de sensibilidad en I-SPY2
# Proyecto:
#   Análisis transcriptómico para la identificación de patrones de
#   respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
#   En este script nos centramos en la cohorte principal (GSE173839_ISPY2) y
#   evaluamos hasta qué punto la firma de IFN-γ (y otras vías inmunes) es:
#     - discriminativa (capacidad para separar R vs NR),
#     - estable frente al desequilibrio R:NR (sensibilidad),
#     - biológicamente coherente a nivel génico y de vías.
#
#   *Por coherencia y continuidad, se hizo un ajuste entre script 6/7 con respecto
#   a la planificación inicial. Este Script está etiquetado como 06, pero realmente
#   se corresponde con el bloque de validación / sensibilidad, es decir, hito 7:
#
# Contenido:
#   BLOQUE 1 – ROC para IFN-γ (Hallmark) y resumen de AUC.
#   BLOQUE 2 – Screen de AUC para todas las vías Hallmark (mini-screening).
#   BLOQUE 3 – Modelo "multivía" combinando firmas inmunes (IFN-γ, IFN-α, TNF).
#   BLOQUE 4 – Submuestreos balanceados R/NR para evaluar la estabilidad del AUC.
#   BLOQUE 5 – Proxies génicos de IFN-γ en I-SPY2 + validación en GSE241876.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(pROC)
  library(msigdbr)})

# =============================================================================
# Rutas de trabajo y helpers generales
# =============================================================================

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

write_tsv_safe <- function(df, filename) {
  # Función de ayuda para guardar tablas en la carpeta estándar de resultados.
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

save_plot <- function(p, filename, w = 7, h = 6, dpi = 300) {
  # Función de ayuda para guardar figuras con un tamaño y resolución coherentes
  # en todo el proyecto.
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}

# Scatter genérico para comparaciones de efectos (lo reutilizamos en varios bloques)
plot_scatter_comparison <- function(df,
                                    x_var, y_var,
                                    col_var,
                                    x_lab, y_lab,
                                    title) {
  ggplot(df, aes(x = !!sym(x_var),
                 y = !!sym(y_var),
                 color = !!sym(col_var))) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = "Categoría"
    ) +
    theme_bw(base_size = 11)}

# Función para obtener la expresión de un gen (columna GeneSymbol)
# para un conjunto concreto de muestras. Si un gen aparece varias veces,
# promediamos sus filas para quedarnos con un único vector por gen:
get_gene_expr <- function(expr_tbl, gene_symbol, sample_ids) {
  sub <- expr_tbl %>%
    dplyr::filter(GeneSymbol == gene_symbol)
  
  if (nrow(sub) == 0) {
    # Si el gen no está en el panel, devolvemos NA para todas las muestras.
    return(rep(NA_real_, length(sample_ids)))}
  
  mat <- as.matrix(sub[, sample_ids, drop = FALSE])
  mode(mat) <- "numeric"
  
  # Si hay filas duplicadas de un mismo gen, usamos la media por muestra.
  if (nrow(mat) > 1) {
    vals <- colMeans(mat, na.rm = TRUE)
  } else {
    vals <- as.numeric(mat[1, ])}
  
  vals[match(sample_ids, names(vals))]}

# -----------------------------------------------------------------------------
# 0) Cargar metadatos y scores ssGSEA para I-SPY2
# -----------------------------------------------------------------------------
# En este primer bloque preparamos los metadatos consolidados (meta_master)
# y la matriz de scores ssGSEA Hallmark que calculamos antes en el script 05.
# A partir de aquí trabajaremos siempre en la subcohorte R/NR:

meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character))

# Scores ssGSEA Hallmark (salida del script 05 para I-SPY2)
ssgsea_ispy <- readr::read_tsv(
  file.path(path_tabs, "H5_ISPY2_ssGSEA_scores.tsv"),
  show_col_types = FALSE)

# Matriz de scores (vías x muestras): convertimos la tabla a matrix
stopifnot("pathway" %in% names(ssgsea_ispy))
scores_mat <- as.matrix(ssgsea_ispy[, -1, drop = FALSE])
rownames(scores_mat) <- ssgsea_ispy[["pathway"]]
samples_scores <- colnames(scores_mat)

# Fenotipo I-SPY2 (R/NR) restringido a las muestras presentes en los scores
pheno_ispy <- meta_master %>%
  dplyr::filter(
    cohort   == "GSE173839_ISPY2",
    sample_id %in% samples_scores,
    response %in% c("R", "NR")
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R")))

# Reordenamos la matriz de scores para que el orden de columnas coincida
# exactamente con el orden de muestras en pheno_ispy:
common_samples <- pheno_ispy$sample_id
scores_mat     <- scores_mat[, common_samples, drop = FALSE]
stopifnot(identical(colnames(scores_mat), pheno_ispy$sample_id))

message("I-SPY2 – muestras R/NR para validación/sensibilidad: ",
        ncol(scores_mat))

# =============================================================================
# BLOQUE 1 – ROC para IFN-γ (Hallmark): capacidad discriminativa básica
# =============================================================================

message("=== BLOQUE 1 – ROC para IFN-γ (Hallmark) ===")

path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
if (!(path_ifng %in% rownames(scores_mat))) {
  stop("No se encontró la vía HALLMARK_INTERFERON_GAMMA_RESPONSE en los scores ssGSEA.")}

# Extraemos el score de la vía IFN-γ y lo usamos como marcador continuo:
score_ifng <- as.numeric(scores_mat[path_ifng, ])
names(score_ifng) <- colnames(scores_mat)

roc_df_base <- tibble::tibble(
  sample_id  = pheno_ispy$sample_id,
  response   = pheno_ispy$response,
  score_ifng = score_ifng)

# Curva ROC con pROC (NR como referencia, R como evento):
roc_ifng <- pROC::roc(
  response ~ score_ifng,
  data   = roc_df_base,
  levels = c("NR", "R"),
  direction = ">")
auc_ifng <- as.numeric(pROC::auc(roc_ifng))

roc_df <- tibble::tibble(
  sens = roc_ifng$sensitivities,
  spec = roc_ifng$specificities)

p_roc_ifng <- ggplot(roc_df, aes(x = 1 - spec, y = sens)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = paste0("I-SPY2 – ROC IFN-γ (AUC = ",
                   round(auc_ifng, 3), ")"),
    x = "1 - Especificidad",
    y = "Sensibilidad"
  ) +
  theme_bw(base_size = 11)

save_plot(p_roc_ifng, "H6_ISPY2_ROC_IFNG.png")

summary_ifng <- tibble::tibble(
  marker   = path_ifng,
  AUC      = auc_ifng,
  cohort   = "GSE173839_ISPY2",
  n_R      = sum(pheno_ispy$response == "R"),
  n_NR     = sum(pheno_ispy$response == "NR"))

# Guardamos el resumen:
write_tsv_safe(summary_ifng, "H6_ISPY2_IFNG_AUC.tsv")
message("AUC IFN-γ (Hallmark) en I-SPY2: ", round(auc_ifng, 3))

# =============================================================================
# BLOQUE 2 – AUC por vía Hallmark (mini-screening de capacidad discriminativa)
# =============================================================================

message("=== BLOQUE 2 – AUC por vía Hallmark (screening) ===")

# Función auxiliar: dado el vector de scores de una vía y el fenotipo R/NR,
# calculamos el AUC asociando esa vía a la respuesta clínica:
calc_auc_pathway <- function(scores_vec, pheno_df) {
  df <- tibble::tibble(
    sample_id = pheno_df$sample_id,
    response  = pheno_df$response,
    score     = as.numeric(scores_vec[pheno_df$sample_id]))
  
  # Si por cualquier motivo sólo queda una clase, no podemos definir un AUC.
  if (length(unique(df$response[!is.na(df$score)])) < 2) {
    return(NA_real_)}
  
  roc_obj <- pROC::roc(
    response ~ score,
    data   = df,
    levels = c("NR", "R"),
    direction = ">")
  as.numeric(pROC::auc(roc_obj))}

all_pathways <- rownames(scores_mat)

# Recorremos todas las vías Hallmark y obtenemos un AUC por vía.
auc_vec <- sapply(all_pathways, function(pw) {
  calc_auc_pathway(scores_mat[pw, ], pheno_ispy)})

auc_tbl <- tibble::tibble(
  pathway = all_pathways,
  AUC     = as.numeric(auc_vec)
) %>%
  dplyr::arrange(dplyr::desc(AUC))

write_tsv_safe(auc_tbl, "H6_ISPY2_Hallmark_AUCs.tsv")

# Representamos las 15 vías con mayor AUC
auc_top15 <- auc_tbl %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(
    pathway = factor(pathway, levels = rev(pathway)))

p_auc_top15 <- ggplot(auc_top15,
                      aes(x = pathway, y = AUC)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    title = "I-SPY2 – AUC por vía Hallmark (top 15)",
    x = "Vía Hallmark",
    y = "AUC (R vs NR)"
  ) +
  theme_bw(base_size = 11)

save_plot(p_auc_top15, "H6_ISPY2_Hallmark_AUCs_top15.png")
message("BLOQUE 2 COMPLETADO – AUC calculado para todas las vías Hallmark.")

# =============================================================================
# BLOQUE 3 – Modelo multivía (combinación de firmas inmunes)
# =============================================================================

message("=== BLOQUE 3 – Modelo multivía (firmas inmunes combinadas) ===")

# Seleccionamos algunas vías inmunes "clásicas" como candidatas a combinar:
candidate_paths <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

present_paths <- candidate_paths[candidate_paths %in% rownames(scores_mat)]
if (length(present_paths) < 2) {
  warning("Hay menos de 2 vías inmunes candidatas presentes; ",
          "el modelo multivía será muy limitado.")}

# Construimos un data.frame con la respuesta y los scores de las vías candidatas:
multi_df <- pheno_ispy %>%
  dplyr::select(sample_id, response)
for (pw in present_paths) {
  multi_df[[pw]] <- as.numeric(scores_mat[pw, multi_df$sample_id])}

# Ajustamos un modelo logístico con respuesta R vs NR en función de los scores
form_str <- paste("response ~", paste(present_paths, collapse = " + "))
model_multi <- glm(
  formula = as.formula(form_str),
  data    = multi_df,
  family  = binomial)

multi_df$prob_R_multi <- stats::predict(model_multi, type = "response")

# Comparamos el AUC del modelo multivía con el de IFN-γ sola:
roc_multi <- pROC::roc(
  response ~ prob_R_multi,
  data   = multi_df,
  levels = c("NR", "R"),
  direction = ">")
auc_multi <- as.numeric(pROC::auc(roc_multi))

message("AUC modelo multivía: ", round(auc_multi, 3),
        " | AUC IFN-γ: ", round(auc_ifng, 3))

roc_df_multi <- tibble::tibble(
  sens   = roc_multi$sensitivities,
  spec   = roc_multi$specificities,
  marker = "Multivía")
roc_df_ifng  <- tibble::tibble(
  sens   = roc_ifng$sensitivities,
  spec   = roc_ifng$specificities,
  marker = "IFN-γ")
roc_both <- dplyr::bind_rows(roc_df_multi, roc_df_ifng)

p_roc_multi <- ggplot(roc_both,
                      aes(x = 1 - spec, y = sens, color = marker)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = paste0(
      "I-SPY2 – ROC comparativa (AUC IFN-γ = ",
      round(auc_ifng, 3),
      "; AUC multivía = ",
      round(auc_multi, 3), ")"
    ),
    x = "1 - Especificidad",
    y = "Sensibilidad",
    color = "Modelo"
  ) +
  theme_bw(base_size = 11)

save_plot(p_roc_multi, "H6_ISPY2_ROC_IFNG_vs_multivia.png")

summary_auc_multi <- tibble::tibble(
  modelo = c("IFN-γ", "Multivía"),
  AUC    = c(auc_ifng, auc_multi))
write_tsv_safe(summary_auc_multi, "H6_ISPY2_AUC_IFNG_vs_multivia.tsv")

message("BLOQUE 3 COMPLETADO – modelo multivía ajustado y comparado.")

# =============================================================================
# BLOQUE 4 – Submuestreos balanceados R/NR (estabilidad / sensibilidad)
# =============================================================================

message("=== BLOQUE 4 – Submuestreos balanceados R/NR ===")

# Aquí estudiamos hasta qué punto el AUC de IFN-γ es estable frente al
# desequilibrio R:NR. Para ello generamos muchos subconjuntos balanceados
# y recalculamos el AUC en cada uno:

set.seed(123)
labels <- pheno_ispy$response
scores <- score_ifng
idx_R  <- which(labels == "R")
idx_NR <- which(labels == "NR")
n_R    <- length(idx_R)
n_NR   <- length(idx_NR)
n_min  <- min(n_R, n_NR)

message("I-SPY2 – R: ", n_R, " | NR: ", n_NR,
        " | tamaño de subconjuntos balanceados: ", n_min)

n_iter <- 200
auc_boot <- numeric(n_iter)

for (b in seq_len(n_iter)) {
  sel_R  <- sample(idx_R, n_min, replace = FALSE)
  sel_NR <- sample(idx_NR, n_min, replace = FALSE)
  idx_b  <- c(sel_R, sel_NR)
  
  df_b <- tibble::tibble(
    response = labels[idx_b],
    score    = scores[idx_b])
  
  roc_b <- pROC::roc(
    response ~ score,
    data   = df_b,
    levels = c("NR", "R"),
    direction = ">")
  
  auc_boot[b] <- as.numeric(pROC::auc(roc_b))}

boot_tbl <- tibble::tibble(
  iter = seq_len(n_iter),
  AUC  = auc_boot)

boot_summary <- boot_tbl %>%
  dplyr::summarise(
    AUC_media  = mean(AUC, na.rm = TRUE),
    AUC_sd     = sd(AUC, na.rm = TRUE),
    AUC_q25    = quantile(AUC, 0.25, na.rm = TRUE),
    AUC_q50    = quantile(AUC, 0.50, na.rm = TRUE),
    AUC_q75    = quantile(AUC, 0.75, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    AUC_original = auc_ifng,
    n_iter       = n_iter,
    n_R          = n_R,
    n_NR         = n_NR,
    n_bal        = n_min)

write_tsv_safe(boot_tbl,     "H6_ISPY2_IFNG_AUC_bootstrap_balanceado.tsv")
write_tsv_safe(boot_summary, "H6_ISPY2_IFNG_AUC_bootstrap_resumen.tsv")

p_boot <- ggplot(boot_tbl, aes(x = AUC)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = auc_ifng, linetype = "dashed") +
  labs(
    title = "I-SPY2 – Distribución de AUC IFN-γ\nsubmuestreos balanceados R/NR",
    x = "AUC (IFN-γ, subconjuntos balanceados)",
    y = "Frecuencia"
  ) +
  theme_bw(base_size = 11)

save_plot(p_boot, "H6_ISPY2_IFNG_AUC_bootstrap_balanceado_hist.png")
message("BLOQUE 4 COMPLETADO – análisis de estabilidad AUC con submuestreos.")

# =============================================================================
# BLOQUE 5 – Proxies de IFN-γ en I-SPY2 + validación génica en MEDI
# =============================================================================

message("=== BLOQUE 5 – Proxies IFN-γ + validación génica en MEDI ===")

# ------------------------------------------------------------------
# 5.1) Proxies de la firma IFN-γ en I-SPY2 (genes Hallmark presentes)
# ------------------------------------------------------------------

# Obtenemos la lista de genes de la vía Hallmark IFN-γ desde MSigDB:
msig_h <- msigdbr::msigdbr(
  species    = "Homo sapiens",
  collection = "H")

genes_ifng_hallmark <- msig_h %>%
  dplyr::filter(gs_name == path_ifng) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

# Cargamos la matriz de expresión Agilent de I-SPY2:
ispy_expr_tbl <- data.table::fread(
  file.path("data", "GSE173839",
            "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"),
  sep = "\t", header = TRUE
) %>% tibble::as_tibble()

# Nos aseguramos de que la primera columna se llame GeneSymbol:
colnames(ispy_expr_tbl)[1] <- "GeneSymbol"

# Nos quedamos sólo con las muestras que estamos usando en pheno_ispy:
expr_samples_ord <- pheno_ispy$sample_id[
  pheno_ispy$sample_id %in% colnames(ispy_expr_tbl)]

ispy_expr_tbl <- ispy_expr_tbl %>%
  dplyr::select(GeneSymbol, dplyr::all_of(expr_samples_ord))

# Reordenamos los scores IFN-γ para que coincidan con estas muestras:
score_ifng_ord <- score_ifng[expr_samples_ord]

# Elegimos los genes de la vía IFN-γ que realmente están en el panel Agilent:
genes_proxies <- intersect(genes_ifng_hallmark, ispy_expr_tbl$GeneSymbol)

message("Nº de genes Hallmark IFN-γ presentes en I-SPY2 (proxies): ",
        length(genes_proxies))

# Para cada gen proxy calculamos la correlación Spearman entre su expresión y
# y el score IFN-γ; y el AUC individual del gen como marcador R vs NR:
proxy_list <- lapply(genes_proxies, function(g) {
  expr_g <- get_gene_expr(ispy_expr_tbl, g, expr_samples_ord)
  
  if (all(is.na(expr_g))) {
    tibble::tibble(
      Gene             = g,
      cor_IFNG_score   = NA_real_,
      pval_cor_IFNG    = NA_real_,
      AUC_gene         = NA_real_)
  } else {
    # Correlación entre score IFN-γ y expresión del gen
    ct <- suppressWarnings(
      cor.test(score_ifng_ord, expr_g, method = "spearman"))
    
    # AUC del gen como marcador univariante
    df_gene <- tibble::tibble(
      response = pheno_ispy$response[match(expr_samples_ord, pheno_ispy$sample_id)],
      expr     = expr_g)
    
    if (length(unique(df_gene$response[!is.na(df_gene$expr)])) < 2) {
      auc_g <- NA_real_
    } else {
      roc_g <- pROC::roc(
        response ~ expr,
        data   = df_gene,
        levels = c("NR", "R"),
        direction = ">")
      auc_g <- as.numeric(pROC::auc(roc_g))}
    
    tibble::tibble(
      Gene             = g,
      cor_IFNG_score   = unname(ct$estimate),
      pval_cor_IFNG    = ct$p.value,
      AUC_gene         = auc_g)}})

proxies_ifng_tbl <- dplyr::bind_rows(proxy_list) %>%
  dplyr::arrange(dplyr::desc(abs(cor_IFNG_score)))

write_tsv_safe(proxies_ifng_tbl,
               "H6_ISPY2_proxies_IFNG_Hallmark_genes.tsv")

message("Proxies IFN-γ en I-SPY2 calculados (correlación con score IFN-γ y AUC individual).")

# ------------------------------------------------------------------
# 5.2) Validación de genes cardinales de IFN-γ en GSE241876 (MEDI)
# ------------------------------------------------------------------

# Definimos el panel de genes "cardinales" de la vía IFN-γ que queremos
# seguir entre cohortes:
genes_cardinal <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "CD274", "PDCD1")

# Resultados limma a nivel de gen para I-SPY2
# Aquí usamos el fichero H4_ISPY2_limma_R_vs_NR.tsv generado en el script 04:

file_ispy_genes <- file.path(path_tabs, "H4_ISPY2_limma_R_vs_NR.tsv")
stopifnot(file.exists(file_ispy_genes))

ispy2_genes <- readr::read_tsv(
  file_ispy_genes,
  show_col_types = FALSE)

# Nos aseguramos de disponer de una columna explícita 'Gene'. Si no está,
# la creamos a partir de los rownames.
if (!("Gene" %in% names(ispy2_genes))) {
  ispy2_genes <- ispy2_genes %>%
    tibble::rownames_to_column("Gene")}

# Nos quedamos sólo con la información relevante para la comparación:
# logFC y FDR para I-SPY2, renombrando las columnas para que sean explícitas.
ispy2_genes_clean <- ispy2_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val)


# Resultados limma para GSE241876 (MEDI), generados en el script 05:
gse_genes <- readr::read_tsv(
  file.path(path_tabs, "H5_GSE241876_limma_R_vs_NR.tsv"),
  show_col_types = FALSE)

if (!("Gene" %in% names(gse_genes))) {
  gse_genes <- gse_genes %>%
    tibble::rownames_to_column("Gene")}

gse_genes_clean <- gse_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_GSE = logFC,
    FDR_GSE   = adj.P.Val)

# Construimos una tabla conjunta para estos genes cardinales en ambas cohortes:
cardinal_tbl <- tibble::tibble(Gene = genes_cardinal) %>%
  dplyr::left_join(ispy2_genes_clean, by = "Gene") %>%
  dplyr::left_join(gse_genes_clean,   by = "Gene") %>%
  dplyr::mutate(
    sig_ISPY2 = !is.na(FDR_ISPY2) & FDR_ISPY2 < 0.05,
    sig_GSE   = !is.na(FDR_GSE)   & FDR_GSE   < 0.05,
    same_dir  = dplyr::case_when(
      is.na(logFC_ISPY2) | is.na(logFC_GSE) ~ NA,
      TRUE ~ sign(logFC_ISPY2) == sign(logFC_GSE)))

write_tsv_safe(cardinal_tbl,
               "H6_IFNG_cardinal_genes_I-SPY2_vs_GSE241876.tsv")

message("Validación génica en MEDI (GSE241876) completada para genes cardinales IFN-γ.")

message("=== Script 06 (validación + sensibilidad + proxies/validación génica) COMPLETADO ===")
