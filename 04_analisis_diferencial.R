# =============================================================================
# 04_analisis_diferencial.R – R vs NR por cohorte
# Proyecto:
#   Análisis transcriptómico para la identificación de patrones de
#   respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este script vamos a desarrollar el análisis diferencial de expresión entre
# respondedores (R) y no respondedores (NR) para tres cohortes concretas del
# proyecto. Por un lado, trabajamos con la cohorte I-SPY2 (GSE173839), en la
# que disponemos de datos de microarrays Agilent; en este caso aplicamos limma
# clásico sobre las intensidades, ajustando además por el brazo de tratamiento
# para tener en cuenta el efecto del brazo control frente a durvalumab/olaparib.
# Por otro lado, analizamos las cohortes HTG MEDI4736 y HTG BrighTNess, ambas
# basadas en la tecnología HTG EdgeSeq OBP y con matrices de expresión ya
# normalizadas. En estas dos cohortes aplicamos limma de forma directa sobre
# la matriz normalizada utilizando la respuesta como covariable principal.
#
# Es importante destacar que la cohorte GSE241876 (RNA-seq) no se incluye en
# este script porque solo dispone de dos pacientes no respondedores, lo que
# hace que un análisis limma clásico no sea adecuado. Esta cohorte se utilizará
# más adelante de forma exploratoria en el contexto de ssGSEA/GSVA en otro
# script, donde el foco estará en patrones de enriquecimiento de firmas y no
# en un análisis diferencial gen a gen.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)})

# Empezamos definiendo las rutas donde vamos a guardar tanto las tablas como las
# figuras que se generen en este script:

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

showWarnings = FALSE


#### Funciones de apoyo ####

# Definimos primero una función sencilla para limpiar y homogeneizar
# identificadores pasando todo a carácter, recortando espacios en blanco y
# convirtiendo a mayúsculas. Esto lo usamos tanto en los IDs de muestra como
# en algunos campos de metadatos:
to_upper_trim <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)}

# Definimos también una función auxiliar para guardar tablas en formato TSV en
# la carpeta de resultados, mostrando un mensaje cada vez que escribimos un
# fichero:
write_tsv_safe <- function(df, filename){
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

# También utilizamos esta función para guardar gráficos en PNG de
# manera homogénea, evitando repetir la llamada a ggsave con los mismos
# argumentos en distintos puntos del código:
save_plot <- function(p, filename, w = 8, h = 6, dpi = 300){
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}

# A continuación definimos una función que, a partir de una tabla donde la
# primera columna contiene los identificadores de gen y el resto son columnas
# de expresión, construye una lista con la matriz numérica de expresión y el
# vector de genes. Así podemos estandarizar el punto de partida para
# limma en todas las cohortes donde la matriz viene en formato tabla:
prepare_expr <- function(tbl, gene_col = 1){
  stopifnot(is.data.frame(tbl))
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)}

# Esta función alinea la matriz de expresión con los metadatos de meta_master
# para una cohorte concreta. A partir de la lista creada con prepare_expr
# recuperamos la matriz, homogeneizamos los nombres de sus columnas, buscamos
# qué columnas tienen IDs presentes en meta_master y, si hace falta, probamos
# a eliminar sufijos típicos de ficheros. El resultado es una matriz X
# filtrada a las muestras que realmente están anotadas:
match_expr_meta <- function(expr_list, meta_df, cohort_tag){
  X  <- expr_list$expr
  cn <- colnames(X)
  
  meta_ids <- meta_df %>%
    dplyr::filter(cohort == cohort_tag) %>%
    dplyr::pull(sample_id) %>%
    unique()
  
  cn_up <- to_upper_trim(cn)
  keep  <- cn_up %in% meta_ids
  
  if (!any(keep)) {
    cn2  <- gsub("\\.CEL$|\\.GZ$|\\.TXT$|\\.FASTQ$|\\.BAM$", "", cn_up)
    keep <- cn2 %in% meta_ids
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn2[keep]
  } else {
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]
  }
  X}

# Para visualizar los resultados de limma utilizamos una función que genera
# un volcano plot. A partir de la tabla de resultados marcamos como
# "significativos" aquellos genes que superan simultáneamente un umbral de
# FDR y de tamaño de efecto (log2FC) y dibujamos las líneas de corte
# correspondientes en el gráfico:
volcano_plot <- function(tt, title, lfc_thr = log2(1.5), p_thr = 0.05){
  tt2 <- tt %>%
    dplyr::mutate(
      sig = ifelse(
        adj.P.Val < p_thr & abs(logFC) >= lfc_thr,
        "Significativo", "No significativo"))
  
  ggplot(tt2, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(shape = sig), alpha = 0.85) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
    geom_hline(yintercept = -log10(p_thr), linetype = 2) +
    labs(title = title, x = "log2(FC)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)}

# Además del volcano plot, generamos también un MA plot clásico, donde
# representamos el logFC frente a la expresión media, marcando la línea
# horizontal en cero para identificar genes sobreexpresados o infraexpresados
# en R frente a NR:
ma_plot <- function(tt, title){
  ggplot(tt, aes(x = AveExpr, y = logFC)) +
    geom_point(alpha = .6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = title,
         x = "Media log2 expresión",
         y = "log2(FC) R vs NR") +
    theme_bw(base_size = 12)}

# Por último, definimos una pequeña función para extraer matrices de expresión
# a partir de los objetos que tenemos en htg_env:
get_htg_matrix <- function(obj){
  if (is.matrix(obj)) {
    return(obj)}
  if (is.data.frame(obj)) {
    genes <- obj[[1]]
    mat   <- as.matrix(obj[, -1, drop = FALSE])
    rownames(mat) <- as.character(genes)
    return(mat)}
  stop("Objeto HTG no soportado (ni matrix ni data.frame).")}


#### Cargamos metadatos armonizados ####

# En este punto cargamos la tabla maestra de metadatos generada en el script 02.
# Sobre ella construiremos los diseños de limma para cada cohorte:
meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character))

# También cargamos el archivo RData que contiene las matrices y metadatos
# relacionados con las cohortes HTG (MEDI4736, BrighTNess, SCAN-B) en un
# entorno aislado. A partir de aquí extraeremos las matrices necesarias:
htg_env <- new.env()
path_htg_rdata <- file.path("data", "G9_HTG", "valid_datasets_HTG.RData")
if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
} else {
  warning("No se encontró valid_datasets_HTG.RData; no se podrán analizar MEDI/BrighTNess.")}


# ---------------------------------------------------------------------------
# 4.a) Cohorte I-SPY2 (GSE173839) – microarrays Agilent
# ---------------------------------------------------------------------------

# Empezamos con la cohorte I-SPY2. Cargamos la tabla de expresión desde el fichero
# original y la convertimos a tibble para manejarla mejor:
ispy_tbl <- data.table::fread(
  file.path(
    "data", "GSE173839",
    "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
  ),
  sep = "\t", header = TRUE
) %>% tibble::as_tibble()

message("I-SPY2 – preparación del diseño (R vs NR con covariable 'arm').")

# A partir de la tabla de expresión construimos la matriz numérica y el vector
# de genes, filtramos los genes con varianza estrictamente positiva para evitar
# problemas numéricos y nos quedamos con esa matriz reducida para el análisis:
prep   <- prepare_expr(ispy_tbl, gene_col = 1)
X_all  <- prep$expr
genes  <- prep$genes
v_all  <- apply(X_all, 1, stats::var)
keep   <- is.finite(v_all) & v_all > 0
X      <- X_all[keep, , drop = FALSE]
genes  <- genes[keep]

# A continuación emparejamos las columnas de la matriz de expresión con los
# metadatos de meta_master para la cohorte I-SPY2. Con esto nos aseguramos de
# trabajar únicamente con muestras que tienen anotación clínica disponible:
X <- match_expr_meta(list(expr = X, genes = genes), meta_master, "GSE173839_ISPY2")
samples <- colnames(X)

# Construimos ahora el fenodata para I-SPY2 a partir de meta_master. Nos
# quedamos solo con las muestras clasificadas como R o NR, homogenizamos el
# brazo de tratamiento, agrupamos brazos muy raros y definimos la variable
# response como factor con NR como categoría de referencia. Así el coeficiente
# que extraemos de limma corresponde a R frente a NR:
pheno <- meta_master %>%
  dplyr::filter(
    cohort    == "GSE173839_ISPY2",
    sample_id %in% samples,
    response  %in% c("R","NR")
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR","R")),
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),
    arm      = stats::relevel(factor(arm), ref = "CONTROL"))

# Reordenamos las filas de pheno para que sigan exactamente el mismo orden que
# las columnas de la matriz de expresión, y comprobamos que los IDs coinciden:
pheno <- pheno[match(samples, pheno$sample_id), , drop = FALSE]
stopifnot(identical(pheno$sample_id, samples))

# Con los datos ya alineados, construimos la matriz de diseño de limma
# incluyendo la respuesta como efecto principal y el brazo de tratamiento como
# covariable. Ajustamos el modelo lineal para cada gen y aplicamos eBayes con
# las opciones de tendencia y robustez activadas:
design <- stats::model.matrix(~ response + arm, data = pheno)
fit    <- limma::lmFit(X, design)
fit    <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# Extraemos ahora específicamente el contraste asociado a R frente a NR,
# comprobando primero que el coeficiente existe en el modelo:
coef_name <- "responseR"
if (!(coef_name %in% colnames(coef(fit)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para I-SPY2")
}

tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("rowid") %>%
  dplyr::mutate(Gene = genes[as.integer(rowid)]) %>%
  dplyr::select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Guardamos la tabla completa de resultados, una versión resumida con los 50
# genes más significativos, y generamos el volcano plot y el MA plot para esta
# cohorte:
write_tsv_safe(tt, "H4_ISPY2_limma_R_vs_NR.tsv")
write_tsv_safe(
  tt %>% dplyr::slice_head(n = 50),
  "H4_ISPY2_limma_R_vs_NR_top50.tsv")

save_plot(
  volcano_plot(tt, "I-SPY2 – R vs NR (limma, ~ response + arm)"),
  "H4_ISPY2_volcano.png")
save_plot(
  ma_plot(tt, "I-SPY2 – MA plot (limma)"),
  "H4_ISPY2_MA.png")

# Por último, construimos un pequeño resumen con el número total de genes
# analizados, cuántos superan FDR < 0.05 y cuántos, además, presentan un
# cambio de expresión absoluto superior a 1.5 en escala lineal (log2FC
# >= log2(1.5)):
summary_ispy <- tibble::tibble(
  cohort                 = "GSE173839_ISPY2",
  n_genes                = nrow(tt),
  n_sig_fdr_05           = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
  n_sig_absLFC_1.5_FDR05 = sum(
    tt$adj.P.Val < 0.05 &
      abs(tt$logFC) >= log2(1.5),
    na.rm = TRUE))
write_tsv_safe(summary_ispy, "H4_ISPY2_resumen.tsv")
print(summary_ispy)

message("Análisis diferencial para I-SPY2 completado.")


# ---------------------------------------------------------------------------
# 4.b) Cohorte HTG MEDI4736 – limma directo (datos HTG normalizados)
# ---------------------------------------------------------------------------

# A continuación analizamos la cohorte HTG MEDI4736, que es la cohorte
# principal del trabajo basada en HTG EdgeSeq OBP. La matriz de expresión ya
# está normalizada en el objeto medi.htg.Rseq dentro de htg_env, por lo que
# podemos aplicar limma directamente sin pasar por edgeR:

if (exists("medi.htg.Rseq", envir = htg_env)) {
  
  medi_htg   <- get("medi.htg.Rseq", envir = htg_env)
  X_medi_htg <- get_htg_matrix(medi_htg)
  
  # Homogeneizamos los nombres de las columnas:
  colnames(X_medi_htg) <- to_upper_trim(colnames(X_medi_htg))
  colnames(X_medi_htg) <- sub("^X", "", colnames(X_medi_htg))
  
  # Construimos el fenodata a partir de meta_master, quedándonos solo con
  # aquellas muestras de la cohorte HTG_MEDI4736 con anotación clara de
  # respuesta R o NR. Definimos de nuevo NR como categoría de referencia:
  pheno_medi <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_MEDI4736",
      response %in% c("R","NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR","R"))
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Emparejamos ahora por la intersección de IDs entre la matriz de expresión
  # y el fenodata para quedarnos únicamente con las muestras comunes:
  common_ids_medi <- intersect(colnames(X_medi_htg), pheno_medi$sample_id)
  
  if (length(common_ids_medi) == 0) {
    warning("HTG MEDI4736: no hay IDs comunes entre expresión y meta_master.")
  } else {
    
    X_medi_htg  <- X_medi_htg[, common_ids_medi, drop = FALSE]
    pheno_medi  <- pheno_medi[match(common_ids_medi, pheno_medi$sample_id), , drop = FALSE]
    
    message("HTG MEDI4736 – muestras R/NR para limma: ", ncol(X_medi_htg))
    
    # Solo tiene sentido ajustar limma si tenemos un número razonable de
    # muestras y ambas categorías de respuesta están presentes. Si se cumple,
    # construimos un diseño simple con la respuesta como único efecto:
    if (ncol(X_medi_htg) >= 4 &&
        length(unique(pheno_medi$response)) == 2) {
      
      design_medi <- stats::model.matrix(~ response, data = pheno_medi)
      
      fit_medi <- limma::lmFit(X_medi_htg, design_medi)
      fit_medi <- limma::eBayes(fit_medi, trend = TRUE, robust = TRUE)
      
      coef_name_medi <- "responseR"
      if (!(coef_name_medi %in% colnames(coef(fit_medi)))) {
        stop("No se encontró el coeficiente 'responseR' en el modelo para HTG MEDI4736")}
      
      tt_medi <- limma::topTable(
        fit_medi, coef = coef_name_medi,
        number = Inf, sort.by = "P"
      ) %>%
        tibble::rownames_to_column("Gene")
      
      write_tsv_safe(tt_medi, "H4_HTG_MEDI4736_limma_R_vs_NR.tsv")
      write_tsv_safe(
        tt_medi %>% dplyr::slice_head(n = 50),
        "H4_HTG_MEDI4736_limma_R_vs_NR_top50.tsv")
      
      save_plot(
        volcano_plot(tt_medi, "HTG MEDI4736 – R vs NR (limma)"),
        "H4_HTG_MEDI4736_volcano.png")
      save_plot(
        ma_plot(tt_medi, "HTG MEDI4736 – MA plot (limma)"),
        "H4_HTG_MEDI4736_MA.png")
      
      summary_medi <- tibble::tibble(
        cohort                 = "HTG_MEDI4736",
        n_genes                = nrow(tt_medi),
        n_sig_fdr_05           = sum(tt_medi$adj.P.Val < 0.05, na.rm = TRUE),
        n_sig_absLFC_1.5_FDR05 = sum(
          tt_medi$adj.P.Val < 0.05 & abs(tt_medi$logFC) >= log2(1.5),
          na.rm = TRUE))
      write_tsv_safe(summary_medi, "H4_HTG_MEDI4736_resumen.tsv")
      print(summary_medi)
      
    } else {
      warning("HTG MEDI4736 no tiene suficientes muestras R/NR para limma (tras emparejar IDs).")}}
  
} else {
  message("Objeto 'medi.htg.Rseq' no encontrado en htg_env; se omite análisis HTG MEDI4736.")}


# ---------------------------------------------------------------------------
# 4.c) Cohorte HTG BrighTNess – limma directo (datos HTG normalizados)
# ---------------------------------------------------------------------------

# Repetimos ahora lo mismo para la cohorte HTG BrighTNess, que utilizamos como
# validación adicional de respuesta pero en un contexto sin inmunoterapia:

if (exists("brtn.htg.Rseq", envir = htg_env)) {
  
  brtn_htg   <- get("brtn.htg.Rseq", envir = htg_env)
  X_brtn_htg <- get_htg_matrix(brtn_htg)
  
  # Homogeneizamos los nombres de las columnas eliminando prefijos y
  # pasando a mayúsculas, igual que hicimos con MEDI4736.
  colnames(X_brtn_htg) <- to_upper_trim(colnames(X_brtn_htg))
  colnames(X_brtn_htg) <- sub("^X", "", colnames(X_brtn_htg))
  
  # Construimos el fenodata a partir de meta_master para la cohorte
  # HTG_BrighTNess, quedándonos solo con muestras R y NR.
  pheno_brtn <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_BrighTNess",
      response %in% c("R","NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR","R"))
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Alineamos por la intersección de IDs entre la matriz de expresión y los
  # metadatos de BrighTNess, y reordenamos para que coincidan:
  common_ids <- intersect(colnames(X_brtn_htg), pheno_brtn$sample_id)
  
  if (length(common_ids) == 0) {
    warning("HTG BrighTNess: no hay IDs comunes entre la matriz de expresión y meta_master.")
  } else {
    X_brtn_htg <- X_brtn_htg[, common_ids, drop = FALSE]
    pheno_brtn <- pheno_brtn[match(common_ids, pheno_brtn$sample_id), , drop = FALSE]
    
    message("HTG BrighTNess – muestras R/NR para limma (IDs comunes): ", ncol(X_brtn_htg))
    
    # Como antes, comprobamos que tenemos suficientes muestras y ambas
    # categorías de respuesta antes de ajustar el modelo:
    if (ncol(X_brtn_htg) >= 4 &&
        length(unique(pheno_brtn$response)) == 2) {
      
      design_brtn <- stats::model.matrix(~ response, data = pheno_brtn)
      
      fit_b <- limma::lmFit(X_brtn_htg, design_brtn)
      fit_b <- limma::eBayes(fit_b, trend = TRUE, robust = TRUE)
      
      coef_name_b <- "responseR"
      if (!(coef_name_b %in% colnames(coef(fit_b)))) {
        stop("No se encontró el coeficiente 'responseR' en el modelo para HTG BrighTNess")}
      
      tt_brtn <- limma::topTable(
        fit_b, coef = coef_name_b,
        number = Inf, sort.by = "P"
      ) %>%
        tibble::rownames_to_column("Gene")
      
      write_tsv_safe(tt_brtn, "H4_HTG_BrighTNess_limma_R_vs_NR.tsv")
      write_tsv_safe(
        tt_brtn %>% dplyr::slice_head(n = 50),
        "H4_HTG_BrighTNess_limma_R_vs_NR_top50.tsv")
      
      save_plot(
        volcano_plot(tt_brtn, "HTG BrighTNess – R vs NR (limma)"),
        "H4_HTG_BrighTNess_volcano.png")
      save_plot(
        ma_plot(tt_brtn, "HTG BrighTNess – MA plot (limma)"),
        "H4_HTG_BrighTNess_MA.png")
      
      summary_brtn <- tibble::tibble(
        cohort                 = "HTG_BrighTNess",
        n_genes                = nrow(tt_brtn),
        n_sig_fdr_05           = sum(tt_brtn$adj.P.Val < 0.05, na.rm = TRUE),
        n_sig_absLFC_1.5_FDR05 = sum(
          tt_brtn$adj.P.Val < 0.05 & abs(tt_brtn$logFC) >= log2(1.5),
          na.rm = TRUE))
      write_tsv_safe(summary_brtn, "H4_HTG_BrighTNess_resumen.tsv")
      print(summary_brtn)
      
    } else {
      warning("HTG BrighTNess no tiene suficientes muestras R/NR para limma (tras emparejar IDs).")}}
  
} else {
  message("Objeto 'brtn.htg.Rseq' no encontrado en htg_env; se omite análisis HTG BrighTNess.")}

message("Análisis diferencial COMPLETADO para I-SPY2, HTG MEDI4736 y HTG BrighTNess (limma / limma directo HTG).")
