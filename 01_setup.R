# ============================================================
# 01_setup.R – Configuración del entorno reproducible del TFM
# Proyecto:
#   Análisis transcriptómico para la identificación de patrones de
#   respuesta a inmunoterapia en cáncer de mama triple negativo.
# ============================================================


# En este primer script vamos a preparar el entorno de trabajo de todo el análisis.
# Para ello, vamos a emplear "renv", configuramos los repositorios de Bioconductor,
# instalamos y cargamos los paquetes y creamos la estructura de carpetas para el
# almacenamiento de datos/resultados. De esta manera nos aseguramos que el análisis
# es completamente reproducible:

# Primero, vamos a instalar "renv" e inicializar el entorno:
install.packages("renv")

if (!dir.exists("renv")) {
  renv::init(bare = TRUE)
} else {
  renv::activate()}

# A continuación, configuramos BiocManager y los repositorios para que todas las
# instalaciones vengan de las fuentes correctas, evitando errores. En Windows además
# indicamos que preferimos instalar binarios, lo cual acelera mucho la instalación
# de los paquetes:

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(repos = BiocManager::repositories())

if (tolower(Sys.info()[["sysname"]]) == "windows") {
  options(pkgType = "binary")}

# Ahora, verificamos si los paquetes de CRAN que vamos a usar en el análisis están
# instalados. De lo contrario, los instalamos y cargamos:

paquetes_cran <- c(
  "tidyverse", "data.table", "readr", "readxl",
  "ggplot2", "ROCR", "msigdbr", "pheatmap")

for (p in paquetes_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))}

# Hacemos lo mismo para los paquetes de Bioconductor. En primer lugar tenemos
# los paquetes base, luego el SE y POMA (para matrices de expresión y metadatos),
# y por último, GSVA, para ejecutar ssGSEA/GSVA en las matrices de expresión:

bioc_base <- c("S4Vectors", "IRanges", "GenomeInfoDb")
bioc_obj  <- c("SummarizedExperiment", "POMA")
bioc_func <- c("GSVA")

for (pkg in c(bioc_base, bioc_obj, bioc_func)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)}
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))}

# Vamos a definir la estructura de carpetas para todo el análisis:

path_data    <- file.path("data")
path_scripts <- file.path("scripts")
path_results <- file.path("results")
path_figs    <- file.path("results", "figuras")
path_tabs    <- file.path("results", "tablas")
path_docs    <- file.path("docs")

dirs <- c(path_data, path_scripts, path_results, path_figs, path_tabs, path_docs)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Ahora, verificamos que los archivos de datos que esperamos utilizar durante
# el análisis estén disponibles en la carpeta correspondiente. Definimos una
# lista explícita con todos los ficheros necesarios, tanto los de HTG como los
# de I-SPY2 y los de GSE241876, y comprobamos uno por uno si existen. Así evitamos
# que se ejecute el resto del pipeline sin algún archivo. Si falta algo, mostramos
# un aviso indicando exactamente qué archivo no está disponible:

esperados <- c(
  file.path(path_data, "G9_HTG", "valid_datasets_HTG.RData"),
  file.path(path_data, "G9_HTG", "G9_combined_results_all_genes_24JUL2023.xlsx"),
  file.path(path_data, "GSE173839",
            "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"),
  file.path(path_data, "GSE173839",
            "GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_raw_count_matrix.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_readme.txt"),
  file.path(path_data, "GSE241876", "GSE241876_clinical_manual.csv"))

faltan <- esperados[!file.exists(esperados)]

if (length(faltan) > 0) {
  warning("Los siguientes archivos esperados no se encontraron en 'data':\n",
          paste(" -", faltan, collapse = "\n"))
} else {
  message("Todos los archivos esperados están disponibles en 'data'.")}

# Finalmente, ejecutamos un snapshot de renv para dejar registrado el estado exacto
# de las versiones de todos los paquetes usados. Envolvemos esta llamada en un "try"
# para evitar errores que interrumpan la ejecución del script:

try(renv::snapshot(prompt = FALSE, force = TRUE), silent = TRUE)

message("Setup COMPLETADO. Archivo 'renv.lock' actualizado y entorno listo.")
