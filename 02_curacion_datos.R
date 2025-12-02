# =======================================================================
# 02_curacion_datos.R – Carga, curación y diccionario de metadatos
# Proyecto:
#  Análisis transcriptómico para la identificación de patrones de
#  respuesta a inmunoterapia en cáncer de mama triple negativo.
# =======================================================================

# En este segundo script procedemos a la curación de la información clínica y de
# anotación de las distintas cohortes que vamos a utilizar en el TFM. Nuestro
# objetivo es cargar los metadatos clínicos y de diseño de cada cohorte, armonizar
# las columnas para que todas compartan el mismo formato, unir todos los metadatos
# en una tabla maestra y construir tanto un diccionario de metadatos como una
# tabla resumen por cohorte, tecnología y respuesta.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(readxl)
  library(SummarizedExperiment)
  library(POMA)})

# -----------------------------
# Rutas principales y helpers
# -----------------------------

# En primer lugar, definimos las rutas del proyecto y nos aseguramos de que exista
# la carpeta donde vamos a guardar las tablas de salida de este script. En esta
# carpeta iremos almacenando los ficheros intermedios y finales relacionados con
# los metadatos:

path_data    <- "data"
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
dir.create(path_tabs, showWarnings = FALSE, recursive = TRUE)

# Definimos una pequeña función auxiliar para guardar ficheros TSV de forma
# segura. De esta manera, centralizamos el manejo de posibles errores de escritura
# y además mostramos un mensaje informativo cada vez que guardamos una tabla:

write_tsv_safe <- function(df, path) {
  tryCatch({
    readr::write_tsv(df, path)
    message("Guardado: ", path)
  }, error = function(e) {
    warning("No se pudo guardar ", path, " -> ", e$message)})}

# También definimos una función auxiliar para homogeneizar identificadores. La
# idea es convertir lo que recibimos a carácter, eliminar espacios en blanco al
# inicio y al final y pasar todo a mayúsculas:

to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)}

# =======================================================================
# 1) COHORTES HTG (BrighTNess, MEDI4736, SCAN-B)
# =======================================================================

# En esta primera sección trabajamos con las cohortes HTG incluidas en el archivo
# 'valid_datasets_HTG.RData'. Lo que queremos aquí es extraer, para cada cohorte,
# los metadatos clínicos relevantes, traducir la respuesta a un formato binario R / NR
#y construir una tabla armonizada 'htg_meta_h' con columnas estándar que luego
# podamos combinar con el resto de estudios:

path_htg_rdata <- file.path(path_data, "G9_HTG", "valid_datasets_HTG.RData")

# Para evitar llenar el entorno global con objetos intermedios, cargamos el
# .RData en un entorno aislado. A partir de ahí, vamos extrayendo lo que nos
# interesa cohorte por cohorte:

htg_env    <- new.env()
htg_meta_h <- NULL

if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
  
  # Utilizamos la lista 'combo' para ir almacenando los metadatos de cada
  # cohorte HTG por separado y, al final, combinarlos en una sola tabla:
  combo <- list()
  
  # ---------------------------
  # HTG BrighTNess
  # ---------------------------
  # En la cohorte BrighTNess disponemos de la variable
  # 'pathologic_complete_response', que indica si el paciente alcanzó pCR.
  # La transformamos a R/NR para unificarla con el resto del proyecto:
  
  if (exists("brtn.htg.cli.dat", envir = htg_env)) {
    df <- get("brtn.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG BrighTNess: ", nrow(df), " muestras")
    
    combo[["HTG_BrighTNess"]] <- df %>%
      mutate(
        resp_raw = toupper(trimws(pathologic_complete_response))
      ) %>%
      transmute(
        cohort     = "HTG_BrighTNess",
        technology = "HTG_EdgeSeq_OBP",
        # Armonizamos el ID de muestra a mayúsculas y sin espacios.
        sample_id  = to_upper_trim(Sample),
        # 'planned_arm' indica el brazo de tratamiento planeado.
        arm        = to_upper_trim(planned_arm),
        # Reetiquetamos la respuesta a R/NR.
        response   = case_when(
          resp_raw == "PCR" ~ "R",
          resp_raw %in% c("RD", "NON-PCR", "NON PCR", "NO PCR") ~ "NR",
          TRUE ~ NA_character_
        ),
        # En esta cohorte no diferenciamos explícitamente momentos PRE/POST.
        timepoint  = NA_character_)}
  
  # ---------------------------
  # HTG MEDI4736
  # ---------------------------
  # En la cohorte MEDI4736 la respuesta patológica completa aparece codificada
  # con distintos formatos (YES/NO, 0/1, PCR/R, etc.). Aquí la traducimos a un
  # formato homogéneo R/NR para poder agrupar los resultados en el análisis:
  
  if (exists("medi.htg.cli.dat", envir = htg_env)) {
    df <- get("medi.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG MEDI4736: ", nrow(df), " muestras")
    
    combo[["HTG_MEDI4736"]] <- df %>%
      mutate(
        resp_raw = toupper(trimws(pathologic.complete.response))
      ) %>%
      transmute(
        cohort     = "HTG_MEDI4736",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(sample_name),
        arm        = to_upper_trim(treatment),
        response   = case_when(
          resp_raw %in% c("YES","Y","1","PCR","R","RESPONDER") ~ "R",
          resp_raw %in% c("NO","N","0","NR","NONRESPONDER")   ~ "NR",
          TRUE ~ NA_character_
        ),
        timepoint  = NA_character_)}
  
  # ---------------------------
  # HTG SCAN-B (sin R/NR; para supervivencia)
  # ---------------------------
  # En el caso de SCAN-B, la cohorte se ha diseñado principalmente para análisis
  # de supervivencia y no disponemos de una clasificación clara en R/NR. Por
  # tanto dejamos 'response' como NA y mantenemos únicamente la información de
  # cohorte, tecnología e ID de muestra:
  
  if (exists("scanb.htg.cli.dat", envir = htg_env)) {
    df <- get("scanb.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG SCAN-B: ", nrow(df), " muestras")
    
    combo[["HTG_SCANB"]] <- df %>%
      transmute(
        cohort     = "HTG_SCANB",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(SAMPLE),
        arm        = NA_character_,
        response   = NA_character_,
        timepoint  = NA_character_)}
  
  # Una vez construidas las tablas de metadatos para cada cohorte HTG, las
  # combinamos en un solo objeto, filtrando muestras sin identificador válido
  # y mostrando un pequeño resumen del número de muestras por cohorte y respuesta:
  
  combo_df <- Filter(\(x) is.data.frame(x) || tibble::is_tibble(x), combo)
  
  if (length(combo_df) > 0) {
    htg_meta_h <- dplyr::bind_rows(combo_df) %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(sample_id) & sample_id != "")
    
    message("Metadatos combinados HTG: ", nrow(htg_meta_h), " muestras totales.")
    htg_meta_h %>%
      dplyr::count(cohort, response) %>%
      print()
  } else {
    message("No se identificaron cohortes HTG con metadatos válidos.")}
  
} else {
  warning("No se encontró archivo HTG: ", path_htg_rdata)}

# =======================================================================
# 2) GSE173839 (I-SPY2, Agilent) – expresión + metadatos clínicos
# =======================================================================

# En esta sección trabajamos con el estudio GSE173839 (I-SPY2). Cargamos tanto
# la matriz de expresión de Agilent como el fichero de biomarcadores, y a partir
# de este último construimos una tabla armonizada de metadatos 'ispy_meta_h' con
# las columnas clave en el formato del trabajo:

path_ispy_expr <- file.path(
  path_data, "GSE173839",
  "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
)
path_ispy_meta <- file.path(
  path_data, "GSE173839",
  "GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz")

ispy_expr <- ispy_meta <- NULL

# --- Matriz de expresión I-SPY2 ---
# Cargamos la matriz de expresión de Agilent desde el archivo de texto y
# mostramos un pequeño resumen de dimensiones para comprobar que la lectura
# se ha realizado correctamente:

if (file.exists(path_ispy_expr)) {
  ispy_expr <- data.table::fread(path_ispy_expr, sep = "\t", header = TRUE) %>%
    tibble::as_tibble()
  message("I-SPY2 matriz (Agilent) cargada. Dim: ", paste(dim(ispy_expr), collapse = " x "))
} else {
  warning("No se encontró matriz I-SPY2: ", path_ispy_expr)}

# --- Metadatos clínicos I-SPY2 ---
# El fichero de biomarcadores contiene las variables clínicas, de tratamiento
# y de respuesta que necesitamos armonizar. Lo cargamos y revisamos sus
# dimensiones para asegurarnos de que todo está en orden:

if (file.exists(path_ispy_meta)) {
  ispy_meta <- data.table::fread(path_ispy_meta) %>% tibble::as_tibble()
  message("I-SPY2 metadatos cargados. Dim: ", paste(dim(ispy_meta), collapse = " x "))
} else {
  warning("No se encontró metadatos I-SPY2: ", path_ispy_meta)
}

# --- Estandarización de columnas clave para meta_master ---
# A continuación identificamos las columnas candidatas a contener el ID de
# muestra, el brazo de tratamiento, la respuesta y el momento de la biopsia.
# Seleccionamos la primera coincidencia para cada tipo de variable y la
# mapeamos a la estructura común del trabajo:

ispy_meta_h <- NULL

if (!is.null(ispy_meta)) {
  cand_id   <- intersect(
    names(ispy_meta),
    c("sample_id","Sample_ID","ResearchID"))
  cand_arm  <- intersect(
    names(ispy_meta),
    c("arm","Arm","Treatment"))
  cand_resp <- intersect(
    names(ispy_meta),
    c("response","Response","pCR.status"))
  cand_tp   <- intersect(
    names(ispy_meta),
    c("timepoint","Timepoint","biopsy_time","Visit"))
  
  ispy_meta_h <- ispy_meta %>%
    dplyr::mutate(
      sample_id = if (length(cand_id))   .[[cand_id[1]]]   else NA,
      arm       = if (length(cand_arm))  .[[cand_arm[1]]]  else NA,
      response  = if (length(cand_resp)) .[[cand_resp[1]]] else NA,
      timepoint = if (length(cand_tp))   .[[cand_tp[1]]]   else NA
    ) %>%
    dplyr::transmute(
      cohort     = "GSE173839_ISPY2",
      technology = "Agilent_microarray",
      sample_id  = to_upper_trim(sample_id),
      arm        = dplyr::if_else(is.na(arm), NA_character_, to_upper_trim(arm)),
      response   = dplyr::if_else(
        is.na(response),
        NA_character_,
        toupper(as.character(response))
      ),
      timepoint  = dplyr::if_else(
        is.na(timepoint),
        NA_character_,
        as.character(timepoint)
      )
    ) %>%
    # Finalmente, traducimos las distintas codificaciones de respuesta
    # (1/0, R/NR, YES/NO, etc.) al formato binario R/NR que utilizaremos
    # en todo el TFM:
    dplyr::mutate(
      response = dplyr::case_when(
        response %in% c("1","PCR","R","RESPONDER","YES","SI","SÍ") ~ "R",
        response %in% c("0","NO_PCR","NONRESPONDER","NR","NO","NO RESPONDER","-1") ~ "NR",
        TRUE ~ response))}

# =======================================================================
# 3) GSE241876 (RNA-seq) – expresión + clínico manual (15 pacientes)
# =======================================================================

# En esta sección trabajamos con la cohorte GSE241876. Vamos a construir una tabla
# de metadatos 'medi_meta_h' armonizada con el resto de cohortes, pero limitándonos
# únicamente a los pacientes evaluables para los que ya hemos curado manualmente la
# información clínica y de respuesta:

path_medi_raw    <- file.path(path_data, "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
path_medi_norm   <- file.path(path_data, "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz")
path_medi_readme <- file.path(path_data, "GSE241876", "GSE241876_readme.txt")
path_medi_clin   <- file.path(path_data, "GSE241876", "GSE241876_clinical_manual.csv")

medi_raw  <- medi_norm <- NULL
medi_clin <- NULL

# --- Cargar matrices de expresión ---
# Cargamos tanto la matriz de recuentos brutos como la normalizada con DESeq2
# si los ficheros están disponibles, simplemente para tenerlos accesibles en el
# entorno cuando los necesitemos en scripts posteriores:

if (file.exists(path_medi_raw)) {
  medi_raw <- data.table::fread(path_medi_raw) |> as_tibble()
  message("GSE241876 raw cargado. Dim: ", paste(dim(medi_raw), collapse = " x "))}
if (file.exists(path_medi_norm)) {
  medi_norm <- data.table::fread(path_medi_norm) |> as_tibble()
  message("GSE241876 DESeq2 cargado. Dim: ", paste(dim(medi_norm), collapse = " x "))}
if (file.exists(path_medi_readme)) {
  message("Readme GSE241876: ", path_medi_readme)}

# --- Cargar CSV clínico (YA LIMPIO: solo 15 pacientes evaluables R/NR) ---
# El fichero 'GSE241876_clinical_manual.csv' es el resultado de una curación
# previa donde ya hemos definido qué pacientes son evaluables y cómo se
# codifica la respuesta (R/NR) en base a la información contenida en el estudio.
# Aquí simplemente nos encargamos de armonizar el formato de las variables a
# mayúsculas y sin espacios:

if (file.exists(path_medi_clin)) {
  medi_clin <- readr::read_delim(
    file = path_medi_clin,
    delim = ";",
    show_col_types = FALSE,
    trim_ws = TRUE
  ) |>
    dplyr::mutate(
      patient_id  = to_upper_trim(patient_id),
      response    = toupper(trimws(response)),   # R / NR
      recist      = toupper(trimws(recist)),
      PDL1_status = toupper(trimws(PDL1_status)),
      alive       = toupper(trimws(alive)))
  
  message("Metadatos clínicos manuales (GSE241876) cargados. n = ", nrow(medi_clin))
  print(medi_clin)
} else {
  warning("No se encontró el fichero clínico manual: ", path_medi_clin)}

# --- Limpieza de columnas de anotación ---
# Definimos una función que elimine las columnas de anotación génica para
# quedarnos únicamente con las columnas que representan muestras. Esto facilita
# el manejo posterior de las matrices de expresión:

drop_annot_cols <- function(tbl) {
  if (is.null(tbl)) return(tbl)
  bad  <- c("ENSEMBLEID","ENTREZID","GENESYMBOL","SYMBOL","GENE","DESCRIPTION","COLUMN1")
  keep <- !(toupper(names(tbl)) %in% bad)
  tbl[, keep]}

medi_raw  <- drop_annot_cols(medi_raw)
medi_norm <- drop_annot_cols(medi_norm)

# --- Extraer IDs de muestra ---
# A partir de las matrices de expresión extraemos los identificadores de
# muestra tomando los nombres de columna, asumiendo que la primera columna
# corresponde a genes y el resto a muestras:

extract_sample_ids <- function(expr_tbl) {
  if (is.null(expr_tbl)) return(NULL)
  if (ncol(expr_tbl) <= 1) return(NULL)
  colnames(expr_tbl)[-1]
}

medi_ids <- extract_sample_ids(medi_norm)
if (is.null(medi_ids)) medi_ids <- extract_sample_ids(medi_raw)

# --- Crear metadatos armonizados ---
# A partir de los IDs de muestra y de la tabla clínica manual, construimos la
# tabla 'medi_meta_h' con la misma estructura que el resto de cohortes
# (cohort, technology, sample_id, arm, response, timepoint), quedándonos
# únicamente con pacientes evaluables:

medi_meta_h <- NULL

if (!is.null(medi_ids) && !is.null(medi_clin)) {
  
  medi_meta_h <- tibble(
    cohort     = "GSE241876",
    technology = "RNAseq",
    sample_id  = to_upper_trim(medi_ids)
  ) |>
    mutate(
      # Derivamos el ID de paciente a partir de la muestra, quitando el sufijo
      # de tiempo cuando existe:
      patient_id = gsub("_(PRE|POST)$", "", sample_id),
      timepoint  = case_when(
        grepl("_PRE$",  sample_id) ~ "PRE",
        grepl("_POST$", sample_id) ~ "POST",
        TRUE ~ NA_character_
      ),
      # Todos los pacientes de esta cohorte reciben la misma combinación de
      # tratamiento neoadyuvante según el diseño del estudio:
      arm        = "CNP_PEMBROLIZUMAB+CARBOPLATINO+NAB_PACLITAXEL"
    ) |>
    # Unimos con la tabla clínica manual (15 pacientes evaluables R/NR).
    left_join(medi_clin, by = "patient_id") |>
    # Nos quedamos solo con pacientes presentes en el CSV clínico y con
    # respuesta conocida:
    filter(
      patient_id %in% medi_clin$patient_id,
      !is.na(response)
    ) |>
    # Priorizamos PRE > POST > NA para quedarnos con una única muestra por
    # paciente, siempre que sea posible:
    arrange(
      patient_id,
      match(timepoint, c("PRE","POST"))
    ) |>
    distinct(patient_id, .keep_all = TRUE) |>
    # Normalizamos 'response' explícitamente a R/NR por seguridad:
    mutate(
      response = case_when(
        response %in% c("R","NR") ~ response,
        TRUE ~ NA_character_
      )
    )
  
  message("Metadatos armonizados GSE241876 (solo pacientes evaluables): ",
          nrow(medi_meta_h), " filas.")
  medi_meta_h %>%
    dplyr::count(response, name = "n_pacientes") %>%
    print()
  
} else {
  message("No se pudieron construir metadatos armonizados para GSE241876 (faltan IDs o clínico).")}

# =======================================================================
# 4) TABLA MAESTRA + DICCIONARIO + RESUMEN
# =======================================================================

# Para terminar, combinamos los metadatos armonizados de todas las cohortes
# disponibles (HTG, I-SPY2 y GSE241876) en una tabla maestra 'meta_master'.
# A partir de ella construimos un diccionario de metadatos que describe el
# significado de cada campo y generamos una tabla resumen por cohorte,
# tecnología y respuesta para tener una visión global de cuántas muestras R/NR
# tenemos en cada caso:

meta_list <- Filter(
  Negate(is.null),
  list(ispy_meta_h, medi_meta_h, htg_meta_h)
)

if (length(meta_list) == 0) {
  meta_master <- tibble::tibble()
} else {
  meta_master <- dplyr::bind_rows(meta_list) %>%
    # Nos aseguramos de que cada combinación cohort/technology/sample_id sea
    # única, por si hubiese duplicados indeseados.
    dplyr::distinct(cohort, technology, sample_id, .keep_all = TRUE) %>%
    dplyr::arrange(cohort, technology, sample_id) %>%
    tibble::as_tibble()}

message("Dim meta_master: ", paste(dim(meta_master), collapse = " x "))
print(utils::head(meta_master, 10))

# Guardamos la tabla maestra:
write_tsv_safe(meta_master, file.path(path_tabs, "meta_master.tsv"))

# --- Diccionario de metadatos ---
# Creamos una pequeña tabla que describa el significado de cada campo de
# 'meta_master' para facilitar la interpretación de los resultados:

dict_meta <- tibble::tibble(
  campo = c("cohort","technology","sample_id","arm","response","timepoint"),
  descripcion = c(
    "Cohorte o estudio de origen",
    "Tecnología (RNAseq / Agilent_microarray / HTG_EdgeSeq_OBP)",
    "Identificador de muestra (armonizado en mayúsculas y sin espacios)",
    "Brazo de tratamiento (si aplica)",
    "Respuesta: R/NR/NA (pCR, RECIST u otros criterios, según la cohorte)",
    "Momento de biopsia o visita (si aplica; por ejemplo, PRE/POST tratamiento)"))

write_tsv_safe(dict_meta, file.path(path_tabs, "diccionario_metadatos.tsv"))

# Por último, generamos una tabla resumen para tener un recuento de cuántas
# muestras R/NR (o NA) tenemos por cohorte y tecnología:

if (nrow(meta_master) > 0) {
  resumen1 <- meta_master %>%
    dplyr::count(cohort, technology, response, name = "n_muestras") %>%
    dplyr::arrange(cohort, technology, dplyr::desc(n_muestras))
  
  write_tsv_safe(resumen1, file.path(path_tabs, "resumen_cohorte_tecnologia_response.tsv"))
  print(resumen1)}

message("Carga y curación inicial COMPLETADA.")

