# =============================================================================
# 05_analisis_funcional.R – Análisis funcional por cohorte
# Proyecto:
#   Análisis transcriptómico para la identificación de patrones de
#   respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este script realizamos el análisis funcional en cuatro niveles:
#
#   A) – I-SPY2:
#     - Cálculo de actividades de vías Hallmark mediante ssGSEA (GSVA, API nueva).
#     - Modelo limma a nivel de vía (R vs NR) ajustando por brazo de tratamiento.
#
#   B) – GSE241876 (RNA-seq crudo, cohorte MEDI)
#     - Análisis diferencial génico (limma-voom, R vs NR).
#     - GSEA clásica (fgsea) usando las estadísticas de gen (t).
#
#   C) – HTG_BrighTNess
#     - ssGSEA Hallmark sobre matriz HTG (OBP).
#     - Modelo limma a nivel de vía (R vs NR).
#
#   D) – HTG_SCANB
#     - ssGSEA Hallmark descriptivo (sin R/NR).
#     - Resumen de la vía de respuesta a IFN-γ.
#
# Flujo general:
#   1) Cargar meta_master y las matrices de expresión necesarias.
#   2) Cargar colecciones Hallmark de MSigDB (msigdbr).
#   3) PARTE A: ssGSEA + limma sobre scores (I-SPY2).
#   4) PARTE B: limma-voom + GSEA (GSE241876).
#   5) PARTE C: ssGSEA + limma (HTG_BrighTNess).
#   6) PARTE D: ssGSEA descriptivo (HTG_SCANB).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(GSVA)
  library(msigdbr)
  library(fgsea)})


#### 1) Rutas y funciones auxiliares ####

# En primer lugar, definimos las rutas de trabajo. Si falta algo,
# detenemos la ejecución con un mensaje claro:

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")

rutas_necesarias <- c(path_tabs, path_figs)
faltan_rutas     <- rutas_necesarias[!dir.exists(rutas_necesarias)]

if (length(faltan_rutas) > 0) {
  stop(
    "Las siguientes carpetas necesarias no existen:\n",
    paste(" -", faltan_rutas, collapse = "\n"),
    "\nPor favor, crea la estructura de carpetas (p.ej. ejecutando 01_setup.R) antes de continuar."
  )
} else {
  message("Estructura de carpetas verificada: ",
          paste(rutas_necesarias, collapse = ", "))}

# Helper para armonizar identificadores (muestras, pacientes, etc.):
to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)}

# Guardado de tablas TSV con mensajito:
write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

# Guardado homogéneo de figuras:
save_plot <- function(p, filename, w = 8, h = 6, dpi = 300) {
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}

# A partir de una tabla de expresión (genes en primera columna, resto muestras)
# devolvemos una lista con la matriz de expresión y el vector de genes:
prepare_expr <- function(tbl, gene_col = 1) {
  stopifnot(is.data.frame(tbl))
  if (ncol(tbl) < 2) stop("La tabla de expresión no tiene columnas suficientes.")
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)}

# Emparejamos matriz de expresión y meta_master para una cohorte concreta:
match_expr_meta <- function(expr_list, meta_df, cohort_tag) {
  X  <- expr_list$expr
  cn <- colnames(X)
  
  meta_ids <- meta_df %>%
    dplyr::filter(cohort == cohort_tag) %>%
    dplyr::pull(sample_id) %>%
    unique()
  
  cn_up <- to_upper_trim(cn)
  keep  <- cn_up %in% meta_ids
  
  if (!any(keep)) {
    # Intentamos limpiar sufijos típicos de ficheros si no hay match directo
    cn2  <- gsub("\\.CEL$|\\.GZ$|\\.TXT$|\\.FASTQ$|\\.BAM$", "", cn_up)
    keep <- cn2 %in% meta_ids
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn2[keep]
  } else {
    X <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]}
  
  samples <- colnames(X)
  
  pheno <- meta_df %>%
    dplyr::filter(cohort == cohort_tag, sample_id %in% samples) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  pheno <- pheno[match(samples, pheno$sample_id), , drop = FALSE]
  stopifnot(identical(pheno$sample_id, samples))
  
  list(expr = X, genes = expr_list$genes, pheno = pheno)}

# Colapsamos genes duplicados por símbolo usando la media:
collapse_by_gene <- function(X, genes) {
  ord   <- order(genes)
  genes <- genes[ord]
  X     <- X[ord, , drop = FALSE]
  X2    <- rowsum(X, group = genes) / as.vector(table(genes))
  list(expr = X2, genes = rownames(X2))
}

# Gráfico de barras de vías (espera columnas logFC y adj.P.Val):
plot_pathway_bar <- function(tt, title, top_n = 15) {
  tt2 <- tt %>%
    dplyr::mutate(
      pathway = rownames(.),
      sig     = adj.P.Val < 0.05
    ) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ggplot(tt2, aes(x = pathway, y = logFC, fill = sig)) +
    geom_col() +
    coord_flip() +
    labs(
      title = title,
      x     = "Vía Hallmark",
      y     = "log2(FC) (R vs NR)"
    ) +
    theme_bw(base_size = 11)}

# A partir de un objeto HTG devolvemos una matriz genes x muestras con
# rownames = símbolo de gen (o ID disponible):
get_htg_matrix <- function(obj) {
  if (is.matrix(obj)) {
    if (is.null(rownames(obj))) {
      stop("La matriz HTG no tiene rownames; no podemos identificar los genes.")}
    return(obj)}
  if (is.data.frame(obj)) {
    if (ncol(obj) < 2) {
      stop("El data.frame HTG tiene menos de 2 columnas.")}
    genes <- obj[[1]]
    mat   <- as.matrix(obj[, -1, drop = FALSE])
    rownames(mat) <- as.character(genes)
    suppressWarnings(mode(mat) <- "numeric")
    return(mat)}
  stop("Objeto HTG no soportado (ni matrix ni data.frame).")}


#### 2) Cargar meta_master y matrices de expresión ####

# Cargamos meta_master desde el TSV generado en 02_curacion_datos.R:
meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character))

message("meta_master cargado. Dim: ", paste(dim(meta_master), collapse = " x "))

# Matriz de expresión I-SPY2 (microarrays Agilent):
ispy_tbl <- data.table::fread(
  file.path(
    "data", "GSE173839",
    "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
  ),
  sep = "\t", header = TRUE
) %>% tibble::as_tibble()

message("I-SPY2 (Agilent) cargado. Dim: ", paste(dim(ispy_tbl), collapse = " x "))

# Matriz raw de GSE241876 (RNA-seq crudo):
gse_tbl <- data.table::fread(
  file.path("data", "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
) %>% tibble::as_tibble()

message("GSE241876 raw cargado. Dim: ", paste(dim(gse_tbl), collapse = " x "))

# Entorno HTG con BrighTNess, MEDI y SCAN-B:
htg_env <- new.env()
path_htg_rdata <- file.path("data", "G9_HTG", "valid_datasets_HTG.RData")
if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
} else {
  stop("No se encontró valid_datasets_HTG.RData; revisar ruta o ejecutar scripts previos.")}


#### 3) Colecciones Hallmark de MSigDB ####

msig_h <- msigdbr::msigdbr(
  species    = "Homo sapiens",
  collection = "H")

hallmark_list <- msig_h %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  { setNames(.$genes, .$gs_name) }

message("Vías Hallmark cargadas: ", length(hallmark_list))

# Nos guardamos el nombre de la vía de IFN-γ para usarlo luego en SCAN-B:
path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
if (!(path_ifng %in% names(hallmark_list))) {
  warning("No se encontró la vía HALLMARK_INTERFERON_GAMMA_RESPONSE en Hallmark.")}


# =============================================================================
# PARTE A – I-SPY2: ssGSEA (Hallmark) + limma a nivel de vía
# =============================================================================

message("=== PARTE A – I-SPY2: ssGSEA Hallmark + limma vías (R vs NR) ===")

#### A1) Preparación de expresión I-SPY2 ####

prep_ispy <- prepare_expr(ispy_tbl, gene_col = 1)
X_all     <- prep_ispy$expr
genes_all <- prep_ispy$genes

# Filtramos genes con varianza > 0 para evitar problemas:
v_all <- apply(X_all, 1, stats::var)
keep  <- is.finite(v_all) & v_all > 0
X     <- X_all[keep, , drop = FALSE]
genes <- genes_all[keep]

# Emparejamos con meta_master para la cohorte I-SPY2:
matched_ispy <- match_expr_meta(
  list(expr = X, genes = genes),
  meta_master,
  "GSE173839_ISPY2")

X          <- matched_ispy$expr
genes      <- matched_ispy$genes
pheno_ispy <- matched_ispy$pheno

# Nos quedamos solo con R/NR y preparamos 'response' y 'arm':
pheno_ispy <- pheno_ispy %>%
  dplyr::filter(response %in% c("R", "NR")) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R")),   # coef 'responseR' = R > NR
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),       # agrupamos brazos muy minoritarios
    arm      = stats::relevel(factor(arm), ref = "CONTROL"))

samples_keep_ispy <- pheno_ispy$sample_id
X <- X[, samples_keep_ispy, drop = FALSE]

pheno_ispy <- pheno_ispy[match(colnames(X), pheno_ispy$sample_id), , drop = FALSE]
stopifnot(identical(pheno_ispy$sample_id, colnames(X)))

# Colapsamos posibles genes duplicados por símbolo:
collapsed_ispy <- collapse_by_gene(X, genes)
X_ssg           <- collapsed_ispy$expr
rownames(X_ssg) <- collapsed_ispy$genes

message("I-SPY2: genes tras colapso: ", nrow(X_ssg))
message("I-SPY2: muestras (R/NR) tras emparejar: ", ncol(X_ssg))


#### A2) ssGSEA con GSVA (API nueva) ####

ssgsea_param_ispy <- GSVA::ssgseaParam(
  exprData = X_ssg,
  geneSets = hallmark_list)

scores_ispy <- GSVA::gsva(
  ssgsea_param_ispy,
  verbose = FALSE)

message("Dim scores ssGSEA I-SPY2 (vías x muestras): ",
        paste(dim(scores_ispy), collapse = " x "))

scores_ispy_df <- as.data.frame(scores_ispy) %>%
  tibble::rownames_to_column("pathway")

write_tsv_safe(scores_ispy_df, "H5_ISPY2_ssGSEA_scores.tsv")


#### A3) Modelo limma a nivel de vía (R vs NR, ajustando por 'arm') ####

design_ispy <- stats::model.matrix(~ response + arm, data = pheno_ispy)

fit_ispy <- limma::lmFit(scores_ispy, design_ispy)
fit_ispy <- limma::eBayes(fit_ispy, trend = TRUE, robust = TRUE)

coef_name <- "responseR"
if (!(coef_name %in% colnames(coef(fit_ispy)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para I-SPY2 (vías).")}

tt_path_ispy <- limma::topTable(
  fit_ispy,
  coef    = coef_name,
  number  = Inf,
  sort.by = "P")

tt_path_ispy_out <- tt_path_ispy %>%
  tibble::rownames_to_column("pathway")

write_tsv_safe(tt_path_ispy_out, "H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv")
write_tsv_safe(
  tt_path_ispy_out %>% dplyr::slice_head(n = 30),
  "H5_ISPY2_ssGSEA_R_vs_NR_pathways_top30.tsv")

summary_ispy_path <- tibble::tibble(
  cohort                 = "GSE173839_ISPY2",
  n_pathways             = nrow(tt_path_ispy),
  n_sig_fdr_05           = sum(tt_path_ispy$adj.P.Val < 0.05),
  n_sig_absLFC_0.2_FDR05 = sum(
    tt_path_ispy$adj.P.Val < 0.05 &
      abs(tt_path_ispy$logFC) >= 0.2))

write_tsv_safe(summary_ispy_path, "H5_ISPY2_ssGSEA_resumen.tsv")
print(summary_ispy_path)

p_bar_ispy <- plot_pathway_bar(
  tt_path_ispy,
  "I-SPY2 – ssGSEA Hallmark: vías asociadas a R vs NR (ajustado por brazo)",
  top_n = 15)

save_plot(p_bar_ispy, "H5_ISPY2_ssGSEA_pathways_barplot_top15.png")

message("PARTE A (I-SPY2 ssGSEA + limma vías) COMPLETADA.")


# =============================================================================
# PARTE B – GSE241876: limma-voom + GSEA clásica (fgsea)
# =============================================================================

message("=== PARTE B – GSE241876: limma-voom + GSEA Hallmark ===")

#### B1) Preparar matriz genes x muestras ####

# Nos quedamos con símbolo de gen + columnas de muestras (R020, R022, etc.)
gse_expr_tbl <- gse_tbl %>%
  dplyr::select(
    GeneSymbol,
    dplyr::starts_with("R"))

message("GSE241876 – expresión (GeneSymbol + muestras). Dim: ",
        paste(dim(gse_expr_tbl), collapse = " x "))

genes_g_all <- gse_expr_tbl$GeneSymbol
X_gse_all   <- as.matrix(gse_expr_tbl[, -1, drop = FALSE])
mode(X_gse_all) <- "numeric"
rownames(X_gse_all) <- genes_g_all

message("GSE241876 – matriz inicial: ",
        paste(dim(X_gse_all), collapse = " x "))

# Filtro por varianza > 0:
v_gse    <- apply(X_gse_all, 1, stats::var)
keep_var <- is.finite(v_gse) & v_gse > 0
X_gse    <- X_gse_all[keep_var, , drop = FALSE]
genes_g  <- genes_g_all[keep_var]

message("GSE241876 – genes tras varianza > 0: ", nrow(X_gse))

# Eliminamos genes sin símbolo:
keep_not_na <- !is.na(genes_g)
X_gse       <- X_gse[keep_not_na, , drop = FALSE]
genes_g     <- genes_g[keep_not_na]

message("GSE241876 – genes tras eliminar NA: ", nrow(X_gse))

# Colapsamos genes duplicados por símbolo:
collapsed_gse <- collapse_by_gene(X_gse, genes_g)
X_gse_coll    <- collapsed_gse$expr
genes_g_coll  <- collapsed_gse$genes
rownames(X_gse_coll) <- genes_g_coll

message("GSE241876 – genes tras colapso por símbolo: ", nrow(X_gse_coll))
message("GSE241876 – muestras (todas): ", ncol(X_gse_coll))


#### B2) Emparejar con meta_master y filtrar R/NR ####

samples_expr_gse <- colnames(X_gse_coll)

pheno_gse <- meta_master %>%
  dplyr::filter(
    cohort    == "GSE241876",
    sample_id %in% samples_expr_gse,
    response  %in% c("R", "NR")
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R")))

samples_keep_gse <- pheno_gse$sample_id
X_gse_coll       <- X_gse_coll[, samples_keep_gse, drop = FALSE]

pheno_gse <- pheno_gse[match(colnames(X_gse_coll), pheno_gse$sample_id), , drop = FALSE]
stopifnot(identical(pheno_gse$sample_id, colnames(X_gse_coll)))

message("GSE241876 – muestras con R/NR: ", ncol(X_gse_coll))
print(table(pheno_gse$response))


#### B3) limma-voom (análisis diferencial génico) ####

dge_gse <- edgeR::DGEList(
  counts = X_gse_coll,
  genes  = data.frame(Gene = rownames(X_gse_coll)))

keep_edge <- edgeR::filterByExpr(dge_gse, group = pheno_gse$response)
dge_gse   <- dge_gse[keep_edge, , keep.lib.sizes = FALSE]
dge_gse   <- edgeR::calcNormFactors(dge_gse, method = "TMM")

message("GSE241876 – genes tras filterByExpr: ", nrow(dge_gse))

design_gse <- stats::model.matrix(~ response, data = pheno_gse)

v_gse_voom <- limma::voom(dge_gse, design_gse, plot = FALSE)
fit_gse    <- limma::lmFit(v_gse_voom, design_gse)
fit_gse    <- limma::eBayes(fit_gse, trend = TRUE, robust = TRUE)

coef_name_gse <- "responseR"
if (!(coef_name_gse %in% colnames(coef(fit_gse)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para GSE241876.")}

tt_gse <- limma::topTable(
  fit_gse,
  coef    = coef_name_gse,
  number  = Inf,
  sort.by = "P")

if (!("Gene" %in% colnames(tt_gse))) {
  tt_gse <- tt_gse %>% tibble::rownames_to_column("Gene")}

write_tsv_safe(tt_gse, "H5_GSE241876_limma_R_vs_NR.tsv")
write_tsv_safe(
  tt_gse %>% dplyr::slice_head(n = 50),
  "H5_GSE241876_limma_R_vs_NR_top50.tsv")

summary_gse <- tibble::tibble(
  cohort                 = "GSE241876",
  n_genes                = nrow(tt_gse),
  n_sig_fdr_05           = sum(tt_gse$adj.P.Val < 0.05, na.rm = TRUE),
  n_sig_absLFC_1.5_FDR05 = sum(
    tt_gse$adj.P.Val < 0.05 &
      abs(tt_gse$logFC) >= log2(1.5),
    na.rm = TRUE))

write_tsv_safe(summary_gse, "H5_GSE241876_limma_resumen.tsv")
print(summary_gse)


#### B4) GSEA clásica (fgsea) con estadísticas de gen ####

stats_gse <- tt_gse$t
names(stats_gse) <- tt_gse$Gene

stats_gse <- stats_gse[!is.na(names(stats_gse))]
stats_gse <- sort(stats_gse, decreasing = TRUE)

fgsea_res <- fgsea::fgsea(
  pathways = hallmark_list,
  stats    = stats_gse,
  minSize  = 10,
  maxSize  = 500,
  nperm    = 10000)

fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
fgsea_res_tbl <- fgsea_res %>% as_tibble()

write_tsv_safe(fgsea_res_tbl, "H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv")

summary_gsea <- tibble::tibble(
  cohort                 = "GSE241876",
  n_pathways             = nrow(fgsea_res_tbl),
  n_sig_fdr_05           = sum(fgsea_res_tbl$padj < 0.05, na.rm = TRUE),
  n_sig_absNES_1.5_FDR05 = sum(
    fgsea_res_tbl$padj < 0.05 &
      abs(fgsea_res_tbl$NES) >= 1.5,
    na.rm = TRUE))

write_tsv_safe(summary_gsea, "H5_GSE241876_fgsea_resumen.tsv")
print(summary_gsea)

# Para reutilizar el barplot, adaptamos NES/padj a columnas logFC/adj.P.Val:
tt_gsea_path <- fgsea_res_tbl %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(logFC = NES, adj.P.Val = padj) %>%
  tibble::column_to_rownames("pathway")

p_bar_gse <- plot_pathway_bar(
  tt_gsea_path,
  "GSE241876 – GSEA Hallmark: vías asociadas a R vs NR",
  top_n = 15)

save_plot(p_bar_gse, "H5_GSE241876_fgsea_Hallmark_barplot_top15.png")

message("PARTE B (GSE241876 limma-voom + GSEA) COMPLETADA.")


# =============================================================================
# PARTE C – HTG_BrighTNess: ssGSEA Hallmark + limma a nivel de vía (R vs NR)
# =============================================================================

message("=== PARTE C – HTG_BrighTNess: ssGSEA Hallmark + limma vías ===")

if (!exists("brtn.htg.Rseq", envir = htg_env)) {
  warning("No se encontró 'brtn.htg.Rseq' en htg_env; se omite análisis HTG_BrighTNess.")
} else {
  # C1) Matriz de expresión HTG BrighTNess
  brtn_htg   <- get("brtn.htg.Rseq", envir = htg_env)
  X_brtn_raw <- get_htg_matrix(brtn_htg)   # genes x muestras
  
  colnames(X_brtn_raw) <- to_upper_trim(colnames(X_brtn_raw))
  colnames(X_brtn_raw) <- sub("^X", "", colnames(X_brtn_raw))
  
  genes_brtn <- rownames(X_brtn_raw)
  
  collapsed_brtn <- collapse_by_gene(X_brtn_raw, genes_brtn)
  X_brtn         <- collapsed_brtn$expr
  genes_brtn     <- collapsed_brtn$genes
  
  message("HTG_BrighTNess – genes tras colapso: ", nrow(X_brtn))
  message("HTG_BrighTNess – muestras totales: ", ncol(X_brtn))
  
  # C2) Fenodata R/NR
  pheno_brtn <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_BrighTNess",
      response %in% c("R", "NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR", "R"))
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  common_ids <- intersect(colnames(X_brtn), pheno_brtn$sample_id)
  message("HTG_BrighTNess – muestras con fenotipo R/NR: ", length(common_ids))
  
  if (length(common_ids) < 4 ||
      length(unique(pheno_brtn$response[pheno_brtn$sample_id %in% common_ids])) < 2) {
    warning("HTG_BrighTNess: muy pocas muestras R/NR tras emparejar; no se realiza limma a nivel de vía.")
  } else {
    X_brtn     <- X_brtn[, common_ids, drop = FALSE]
    pheno_brtn <- pheno_brtn[match(common_ids, pheno_brtn$sample_id), , drop = FALSE]
    
    stopifnot(identical(colnames(X_brtn), pheno_brtn$sample_id))
    
    # C3) ssGSEA Hallmark
    message("Calculando ssGSEA (Hallmark) para HTG_BrighTNess...")
    
    keep_genes <- !is.na(rownames(X_brtn)) & rownames(X_brtn) != ""
    X_brtn_ssg <- X_brtn[keep_genes, , drop = FALSE]
    
    ssgsea_param_brtn <- GSVA::ssgseaParam(
      exprData = X_brtn_ssg,
      geneSets = hallmark_list)
    
    scores_brtn <- GSVA::gsva(
      ssgsea_param_brtn,
      verbose = FALSE)
    
    message("HTG_BrighTNess – ssGSEA dim (vías x muestras): ",
            paste(dim(scores_brtn), collapse = " x "))
    
    scores_brtn_df <- as.data.frame(scores_brtn) %>%
      tibble::rownames_to_column("pathway")
    
    write_tsv_safe(scores_brtn_df,
                   "H5_HTG_BrighTNess_ssGSEA_scores.tsv")
    
    # C4) limma a nivel de vía (~ response)
    design_brtn <- stats::model.matrix(~ response, data = pheno_brtn)
    
    fit_brtn <- limma::lmFit(scores_brtn, design_brtn)
    fit_brtn <- limma::eBayes(fit_brtn, trend = TRUE, robust = TRUE)
    
    coef_name_brtn <- "responseR"
    if (!(coef_name_brtn %in% colnames(coef(fit_brtn)))) {
      stop("HTG_BrighTNess: no se encontró el coeficiente 'responseR' en el modelo de vías.")}
    
    tt_path_brtn <- limma::topTable(
      fit_brtn,
      coef    = coef_name_brtn,
      number  = Inf,
      sort.by = "P")
    
    tt_path_brtn_out <- tt_path_brtn %>%
      tibble::rownames_to_column("pathway")
    
    write_tsv_safe(tt_path_brtn_out,
                   "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv")
    write_tsv_safe(
      tt_path_brtn_out %>% dplyr::slice_head(n = 30),
      "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways_top30.tsv")
    
    summary_brtn_path <- tibble::tibble(
      cohort                 = "HTG_BrighTNess",
      n_pathways             = nrow(tt_path_brtn),
      n_sig_fdr_05           = sum(tt_path_brtn$adj.P.Val < 0.05, na.rm = TRUE),
      n_sig_absLFC_0.2_FDR05 = sum(
        tt_path_brtn$adj.P.Val < 0.05 &
          abs(tt_path_brtn$logFC) >= 0.2,
        na.rm = TRUE))
    
    write_tsv_safe(summary_brtn_path,
                   "H5_HTG_BrighTNess_ssGSEA_resumen.tsv")
    print(summary_brtn_path)
    
    p_bar_brtn <- plot_pathway_bar(
      tt_path_brtn,
      "HTG_BrighTNess – ssGSEA Hallmark: vías asociadas a R vs NR",
      top_n = 15)
    
    save_plot(p_bar_brtn,
              "H5_HTG_BrighTNess_ssGSEA_pathways_barplot_top15.png")
    
    message("PARTE C (HTG_BrighTNess ssGSEA + limma vías) COMPLETADA.")}}


# =============================================================================
# PARTE D – HTG_SCANB: ssGSEA Hallmark descriptivo (sin R/NR)
# =============================================================================

message("=== PARTE D – HTG_SCANB: ssGSEA Hallmark descriptivo ===")

if (!exists("scanb.htg.Rseq", envir = htg_env)) {
  warning("No se encontró 'scanb.htg.Rseq' en htg_env; se omite análisis HTG_SCANB.")
} else {
  # D1) Matriz de expresión HTG SCAN-B
  scanb_htg   <- get("scanb.htg.Rseq", envir = htg_env)
  X_scanb_raw <- get_htg_matrix(scanb_htg)  # genes x muestras
  
  colnames(X_scanb_raw) <- to_upper_trim(colnames(X_scanb_raw))
  colnames(X_scanb_raw) <- sub("^X", "", colnames(X_scanb_raw))
  
  genes_scanb <- rownames(X_scanb_raw)
  
  collapsed_scanb <- collapse_by_gene(X_scanb_raw, genes_scanb)
  X_scanb         <- collapsed_scanb$expr
  genes_scanb     <- collapsed_scanb$genes
  
  message("HTG_SCANB – genes tras colapso: ", nrow(X_scanb))
  message("HTG_SCANB – muestras totales: ", ncol(X_scanb))
  
  # D2) ssGSEA Hallmark descriptivo
  keep_genes  <- !is.na(rownames(X_scanb)) & rownames(X_scanb) != ""
  X_scanb_ssg <- X_scanb[keep_genes, , drop = FALSE]
  
  message("Calculando ssGSEA (Hallmark) para HTG_SCANB...")
  
  ssgsea_param_scanb <- GSVA::ssgseaParam(
    exprData = X_scanb_ssg,
    geneSets = hallmark_list)
  
  scores_scanb <- GSVA::gsva(
    ssgsea_param_scanb,
    verbose = FALSE)
  
  message("HTG_SCANB – ssGSEA dim (vías x muestras): ",
          paste(dim(scores_scanb), collapse = " x "))
  
  scores_scanb_df <- as.data.frame(scores_scanb) %>%
    tibble::rownames_to_column("pathway")
  
  write_tsv_safe(scores_scanb_df,
                 "H5_HTG_SCANB_ssGSEA_scores.tsv")
  
  # D3) Resumen descriptivo de la vía IFN-γ (si está disponible)
  if (path_ifng %in% rownames(scores_scanb)) {
    ifng_scanb <- as.numeric(scores_scanb[path_ifng, ])
    
    summary_ifng_scanb <- tibble::tibble(
      cohort    = "HTG_SCANB",
      pathway   = path_ifng,
      n_samples = length(ifng_scanb),
      mean_score = mean(ifng_scanb, na.rm = TRUE),
      sd_score   = sd(ifng_scanb, na.rm = TRUE),
      q25        = quantile(ifng_scanb, 0.25, na.rm = TRUE),
      q50        = quantile(ifng_scanb, 0.50, na.rm = TRUE),
      q75        = quantile(ifng_scanb, 0.75, na.rm = TRUE))
    
    write_tsv_safe(summary_ifng_scanb,
                   "H5_HTG_SCANB_IFNG_score_resumen.tsv")
    print(summary_ifng_scanb)
    
    df_ifng_scanb <- tibble::tibble(
      sample_id  = colnames(scores_scanb),
      IFNG_score = ifng_scanb)
    
    p_ifng_scanb <- ggplot(df_ifng_scanb, aes(x = IFNG_score)) +
      geom_histogram(bins = 30, fill = "grey70", color = "black") +
      labs(
        title = "HTG_SCANB – Distribución de scores IFN-γ (Hallmark, ssGSEA)",
        x     = "Score ssGSEA IFN-γ",
        y     = "Frecuencia"
      ) +
      theme_bw(base_size = 11)
    
    save_plot(p_ifng_scanb,
              "H5_HTG_SCANB_IFNG_ssGSEA_hist.png")
    
  } else {
    warning("En HTG_SCANB no se encontró la vía HALLMARK_INTERFERON_GAMMA_RESPONSE en los scores ssGSEA.")
    }
  
  message("PARTE D (HTG_SCANB ssGSEA descriptivo) COMPLETADA.")
}

message("=== Script 05 (análisis funcional multi-cohorte) COMPLETADO ===")
