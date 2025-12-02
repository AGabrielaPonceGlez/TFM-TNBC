# =============================================================================
# 03_exploracion_inicial.R – QA (NA/duplicados), boxplots y PCA por cohorte
# Proyecto:
#   Análisis transcriptómico para la identificación de patrones de
#   respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este tercer script nos centramos en una exploración inicial de calidad de
# las distintas cohortes incluidas en el trabajo. Trabajamos con la cohorte
# MEDI4736 (HTG) como cohorte principal, con I-SPY2 (GSE173839, microarrays
# Agilent) como validación cruzada, con BrighTNess (HTG) como una validación
# adicional de respuesta en un contexto sin inmunoterapia, con GSE241876
# (RNA-seq) como validación exploratoria con un número limitado de pacientes
# evaluables y con SCAN-B (HTG) como cohorte destinada a análisis de
# supervivencia.
#
# Lo que hacemos de forma general es, en primer lugar, realizar un control
# básico de calidad sobre los metadatos combinados (meta_master), buscando
# posibles duplicados por cohorte y sample_id y cuantificando los valores
# perdidos (NA) por campo. A continuación, exploramos la cohorte I-SPY2
# mediante boxplots de distribución de expresión por muestra y un PCA
# coloreado por respuesta y estratificado por brazo de tratamiento. Después
# trabajamos con la cohorte GSE241876 aplicando una transformación log1p
# sobre la matriz de expresión y generando de nuevo boxplots y PCA de forma
# exploratoria para las pacientes evaluables R/NR. Finalmente, exploramos las
# matrices HTG disponibles (BrighTNess, MEDI4736, SCAN-B) de manera más
# genérica a través de boxplots y PCA, con el objetivo de detectar posibles
# artefactos globales o patrones de variación llamativos.


suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)})

# -----------------------------------------------------------------------------
# Rutas de salida para resultados y figuras
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Helpers generales
# -----------------------------------------------------------------------------
# A continuación definimos algunas funciones auxiliares que vamos a utilizar de
# forma repetida a lo largo del script. Con esto buscamos centralizar en un solo sitio
# todo lo que tenga que ver con el guardado de figuras y tablas, la preparación
# de matrices de expresión y el emparejamiento entre expresión y metadatos:

# Función para guardar figuras de ggplot2 de forma homogénea:
save_plot <- function(p, filename, w = 9, h = 6, dpi = 300) {
  ggsave(
    filename = file.path(path_figs, filename),
    plot     = p,
    width    = w,
    height   = h,
    dpi      = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}

# Función para guardar tablas en formato TSV con un mensaje de confirmación:
write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

# Función para armonizar identificadores (muestras, pacientes, etc.) a un
# formato consistente. Convertimos a carácter, recortamos espacios y pasamos
# todo a mayúsculas:
to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)}

# A partir de una tabla de expresión donde los genes se encuentran en la primera
# columna y el resto de columnas corresponden a muestras, construimos una lista
# con la matriz de expresión numérica, el vector de genes y el vector de IDs de
# muestra. De esta manera trabajamos con una estructura estándar en el resto del
# script:
prepare_expr <- function(tbl, gene_col = 1) {
  stopifnot(is.data.frame(tbl))
  if (ncol(tbl) < 2) stop("La tabla de expresión no tiene columnas suficientes.")
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes, samples = colnames(mat))}

# Esta función empareja una matriz de expresión con meta_master a partir de
# los sample_id. Partimos de la lista devuelta por prepare_expr, buscamos qué
# columnas de la matriz tienen un identificador presente en meta_master para
# la cohorte de interés y, si no hay coincidencia directa, intentamos también
# limpiar sufijos típicos de ficheros como .CEL, .GZ, etc. Devolvemos una
# versión filtrada de la matriz de expresión con solo las muestras emparejadas:
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
    cn2  <- gsub("\\.CEL$|\\.GZ$|\\.TXT$|\\.FASTQ$|\\.BAM$", "", cn_up)
    keep <- cn2 %in% meta_ids
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn2[keep]
  } else {
    X <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]}
  
  list(expr = X, genes = expr_list$genes, samples = colnames(X))}

# Esta función genera un boxplot de la distribución de expresión por muestra:
plot_box <- function(X, title) {
  df <- as_tibble(X) %>%
    dplyr::mutate(.gene = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      -.gene,
      names_to  = "sample",
      values_to = "value")
  
  ggplot(df, aes(x = sample, y = value)) +
    geom_boxplot(outlier.size = 0.25) +
    labs(title = title, x = "Muestra", y = "Expresión") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank())}

# Esta función realiza un PCA por cohorte y lo representa en el plano de las
# dos primeras componentes principales (PC1 y PC2). Filtramos genes con
# varianza positiva, calculamos el PCA sobre las muestras y unimos el resultado
# con meta_master para colorear por respuesta (R/NR/NE) y utilizar la forma del
# punto para distinguir los brazos de tratamiento:
plot_pca <- function(X, meta_df, cohort_tag, title) {
  v    <- apply(X, 1, stats::var)
  keep <- is.finite(v) & v > 0
  X2   <- X[keep, , drop = FALSE]
  if (ncol(X2) < 3) stop("Menos de 3 muestras tras emparejamiento; no se puede hacer PCA.")
  
  pcs <- prcomp(t(X2), center = TRUE, scale. = TRUE)
  
  pdat <- as_tibble(pcs$x[, 1:2]) %>%
    dplyr::mutate(sample_id = rownames(pcs$x)) %>%
    dplyr::left_join(
      meta_df %>% dplyr::filter(cohort == cohort_tag),
      by = dplyr::join_by(sample_id)
    )
  
  ggplot(pdat, aes(x = PC1, y = PC2, color = response, shape = arm)) +
    geom_point(size = 2, alpha = 0.95) +
    scale_color_manual(
      values = c(
        "R"  = "#1B9E77",
        "NR" = "#D95F02",
        "NE" = "#7570B3"
      ),
      na.value = "grey70"
    ) +
    labs(title = title, color = "Respuesta", shape = "Brazo") +
    theme_bw(base_size = 12)}

# -----------------------------------------------------------------------------
# Cargar objetos del Script 2 (o reconstruir desde disco)
# -----------------------------------------------------------------------------
need    <- c("meta_master", "ispy_expr", "medi_raw", "medi_norm", "htg_env", "htg_meta_h")
missing <- need[!need %in% ls(envir = .GlobalEnv)]
if (length(missing)) {
  message("Faltan en memoria: ", paste(missing, collapse = ", "))}

# Si meta_master no está en memoria, lo cargamos desde el TSV generado en el
# script 02. Si tampoco existe el fichero, detenemos la ejecución (habría que ejecutar
# primero el script 02):
if (!exists("meta_master", envir = .GlobalEnv)) {
  meta_master_path <- file.path(path_tabs, "meta_master.tsv")
  if (file.exists(meta_master_path)) {
    meta_master <- readr::read_tsv(meta_master_path, show_col_types = FALSE) %>%
      dplyr::mutate(across(
        c(cohort, technology, sample_id, arm, response, timepoint),
        ~ ifelse(is.na(.x), NA, as.character(.x))
      ))
    message("meta_master cargado desde: ", meta_master_path)
  } else {
    stop("No se encontró 'meta_master.tsv'. Ejecutar antes 02_curacion_datos.R.")}}

# Reconstruimos I-SPY2 desde disco si no lo tenemos ya en memoria. De nuevo,
# simplemente leemos la matriz de expresión de Agilent y dejamos que el resto
# del script se encargue de emparejarla con meta_master:
if (!exists("ispy_expr", envir = .GlobalEnv)) {
  path_ispy_expr <- file.path(
    "data", "GSE173839",
    "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz")
  if (file.exists(path_ispy_expr)) {
    ispy_expr <- data.table::fread(path_ispy_expr, sep = "\t", header = TRUE) |> as_tibble()
    message("I-SPY2 (Agilent) cargado desde disco. Dim: ",
            paste(dim(ispy_expr), collapse = " x "))
  } else {
    warning("No se encontró la matriz de I-SPY2 en disco.")}}

# Para GSE241876 intentamos cargarmos la matriz normalizada, o si falta, la raw:
base_medi      <- file.path("data", "GSE241876")
path_medi_norm <- file.path(base_medi, "GSE241876_DeseqNormalizedCount.csv.gz")
path_medi_raw  <- file.path(base_medi, "GSE241876_raw_count_matrix.csv.gz")

medi_norm <- medi_raw <- NULL

if (file.exists(path_medi_norm)) {
  medi_norm <- data.table::fread(path_medi_norm) |> as_tibble()
  message("GSE241876 – matriz DESeq2 cargada. Dim: ",
          paste(dim(medi_norm), collapse = " x "))
} else if (file.exists(path_medi_raw)) {
  medi_raw <- data.table::fread(path_medi_raw) |> as_tibble()
  message("GSE241876 – matriz raw cargada. Dim: ",
          paste(dim(medi_raw), collapse = " x "))
} else {
  stop("No se encontró ninguna matriz de GSE241876 en disco.")}


# -----------------------------------------------------------------------------
# QA de metadatos: duplicados y NA
# -----------------------------------------------------------------------------
# Antes de pasar a las matrices de expresión, realizamos un control básico de
# calidad sobre meta_master:

meta_master <- tibble::as_tibble(meta_master)

dup <- meta_master %>%
  dplyr::group_by(cohort, sample_id) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1)

write_tsv_safe(dup, "H3_duplicados_meta_master.tsv")
print(dup, n = 20)

na_tbl <- meta_master %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to  = "campo",
    values_to = "n_na"
  ) %>%
  dplyr::arrange(dplyr::desc(n_na))

write_tsv_safe(na_tbl, "H3_NA_por_campo_meta_master.tsv")
print(na_tbl, n = nrow(na_tbl))

# ============================================================================
# 1) Exploración I-SPY2 (Agilent): boxplot y PCA
# ============================================================================
# En esta sección exploramos la cohorte I-SPY2. Primero preparamos la matriz de
# expresión, la emparejamos con meta_master para asegurarnos de que solo
# trabajamos con muestras anotadas, evaluamos la escala de expresión para
# decidir si aplicamos o no una transformación log1p y finalmente generamos
# tanto los boxplots por muestra como un PCA anotado por respuesta y brazo de
# tratamiento:

if (exists("ispy_expr") && is.data.frame(ispy_expr)) {
  ispy_prep    <- prepare_expr(ispy_expr, gene_col = 1)
  ispy_matched <- match_expr_meta(ispy_prep, meta_master, "GSE173839_ISPY2")
  Xispy        <- ispy_matched$expr
  
  message("I-SPY2 – muestras tras emparejar con meta: ", ncol(Xispy))
  
  med_all <- median(Xispy, na.rm = TRUE)
  if (is.finite(med_all) && med_all > 50) {
    Xispy_t         <- log1p(Xispy)
    trans_used_ispy <- "log1p"
  } else {
    Xispy_t         <- Xispy
    trans_used_ispy <- "as-is"}
  
  p1 <- plot_box(Xispy_t, "I-SPY2 (Agilent) – Distribución por muestra")
  save_plot(p1, "H3_ISPY2_boxplot.png")
  
  p2 <- plot_pca(
    Xispy_t,
    meta_master,
    "GSE173839_ISPY2",
    "I-SPY2 (Agilent) – PCA")
  save_plot(p2, "H3_ISPY2_PCA.png")
  
  tibble(
    cohort    = "GSE173839_ISPY2",
    n_genes   = nrow(Xispy_t),
    n_samples = ncol(Xispy_t),
    transform = trans_used_ispy
  ) |> write_tsv_safe("H3_ISPY2_resumen.tsv")
  
} else {
  message("I-SPY2 no disponible en memoria/disco para exploración.")}

# ============================================================================
# 2) Exploración GSE241876 (RNA-seq, 15 pacientes evaluables)
# ============================================================================
# Ahora exploramos la cohorte GSE241876, centrándonos únicamente en las
# 15 pacientes evaluables que ya hemos definido en el script 02. El objetivo es
# comprobar, de manera exploratoria, si existe algún patrón claro entre R y NR
# a nivel transcriptómico, siempre teniendo presente que el tamaño muestral es
# limitado:

medi_tbl <- if (exists("medi_norm") && is.data.frame(medi_norm)) {
  medi_norm
} else if (exists("medi_raw") && is.data.frame(medi_raw)) {
  medi_raw
} else {
  NULL}

if (!is.null(medi_tbl)) {
  
  drop_annot_cols_local <- function(tbl) {
    bad  <- c("ENSEMBLEID","ENTREZID","GENESYMBOL","SYMBOL",
              "GENE","DESCRIPTION","COLUMN1")
    keep <- !(toupper(names(tbl)) %in% bad)
    tbl[, keep]
  }
  medi_tbl <- drop_annot_cols_local(medi_tbl)
  
  medi_prep    <- prepare_expr(medi_tbl, gene_col = 1)
  medi_matched <- match_expr_meta(medi_prep, meta_master, "GSE241876")
  Xmedi        <- medi_matched$expr
  samples_medi <- colnames(Xmedi)
  
  message("GSE241876 – muestras tras emparejar con meta: ", ncol(Xmedi))
  
  pheno_medi <- meta_master %>%
    dplyr::filter(
      cohort   == "GSE241876",
      sample_id %in% samples_medi
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  pheno_medi <- pheno_medi[match(samples_medi, pheno_medi$sample_id),
                           , drop = FALSE]
  stopifnot(identical(pheno_medi$sample_id, samples_medi))
  
  pheno_medi %>%
    dplyr::count(timepoint, response, name = "n_muestras") %>%
    write_tsv_safe("H3_GSE241876_conteo_response.tsv")
  
  Xmedi_t <- log1p(Xmedi)
  
  p3 <- plot_box(Xmedi_t, "GSE241876 (RNA-seq) – Distribución por muestra")
  save_plot(p3, "H3_GSE241876_boxplot.png")
  
  p4 <- plot_pca(
    Xmedi_t,
    pheno_medi,
    "GSE241876",
    "GSE241876 (RNA-seq) – PCA (pacientes evaluables)"
  )
  save_plot(p4, "H3_GSE241876_PCA.png")
  
  tibble(
    cohort    = "GSE241876",
    n_genes   = nrow(Xmedi_t),
    n_samples = ncol(Xmedi_t),
    transform = "log1p"
  ) |> write_tsv_safe("H3_GSE241876_resumen.tsv")
  
} else {
  message("GSE241876 no disponible en memoria/disco para exploración.")}

# ============================================================================
# 3) Exploración HTG (BrighTNess, MEDI4736, SCAN-B)
# ============================================================================
# Para terminar, hacemos una exploración global de las matrices de expresión
# HTG disponibles en htg_env. No distinguimos aquí todavía entre cohortes
# específicas, sino que buscamos objetos que parezcan matrices de expresión y
# generamos boxplots y PCA genéricos. Esto nos ayuda a detectar posibles
# outliers o efectos globales de lote antes de entrar en análisis más
# específicos.:

if (exists("htg_env")) {
  objs            <- ls(htg_env)
  expr_candidates <- objs[grepl("Rseq|expr|mat|counts", objs, ignore.case = TRUE)]
  did_any         <- FALSE
  
  for (o in expr_candidates) {
    obj <- get(o, envir = htg_env)
    
    if (is.matrix(obj) || is.data.frame(obj)) {
      
      if (is.matrix(obj)) {
        X    <- obj
        Xnum <- suppressWarnings(apply(X, 2, as.numeric))
        if (!is.matrix(Xnum)) Xnum <- as.matrix(X)
        colnames(Xnum) <- colnames(X)
        X_t <- log1p(abs(Xnum))
      } else {
        prep <- prepare_expr(as_tibble(obj), gene_col = 1)
        X_t  <- log1p(abs(prep$expr))}
      
      tag <- gsub("[^A-Za-z0-9_]+", "_", o)
      
      p5 <- plot_box(X_t, paste0("HTG EdgeSeq – Distribución (", o, ")"))
      save_plot(p5, paste0("H3_HTG_", tag, "_boxplot.png"))
      
      v    <- apply(X_t, 1, stats::var)
      keep <- is.finite(v) & v > 0
      if (sum(keep) > 10 && ncol(X_t) >= 3) {
        pcs <- prcomp(t(X_t[keep, , drop = FALSE]),
                      center = TRUE, scale. = TRUE)
        pdat <- as_tibble(pcs$x[, 1:2]) %>%
          dplyr::mutate(sample = rownames(pcs$x))
        
        p6 <- ggplot(pdat, aes(PC1, PC2)) +
          geom_point(size = 2) +
          labs(title = paste0("HTG EdgeSeq – PCA (", o, ")")) +
          theme_bw(base_size = 12)
        
        save_plot(p6, paste0("H3_HTG_", tag, "_PCA.png"))}
      
      tibble(
        cohort    = paste0("HTG_", toupper(tag)),
        n_genes   = nrow(X_t),
        n_samples = ncol(X_t),
        transform = "log1p"
      ) |> write_tsv_safe(paste0("H3_HTG_", tag, "_resumen.tsv"))
      
      did_any <- TRUE}}
  
  if (!did_any) {
    message("No se identificaron matrices de expresión utilizables dentro de htg_env.")
  }
} else {
  message("htg_env no está disponible.")}

message("Exploración inicial (QA) completada.")
