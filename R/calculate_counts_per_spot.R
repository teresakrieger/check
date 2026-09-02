#' Obtain counts per Visium spot
#' 
#' This function is intended as a diagnostic step for checking whether high 
#' counts per spot match the tissue position.
#'
#' @param visium_dir Character scalar. Base directory containing the Visium
#'   sample directories.
#' @param visium_sample Character scalar. Name or identifier of the Visium
#'   sample directory.
#' @param matrix_name Character scaler. Name of the matrix used for calculating
#' the counts per spot. Either `raw` or `filtered`.  
#' @param exclude_deprecated_genes Logical. Whether features marked as 
#' `DEPRECATED` should be excluded from the calculation. Defaults to `TRUE`.
#' @param remove_background_spots Logical. Whether background spots should 
#' be removed from the returned object. Defaults to `FALSE`, so background spots
#' are retained and can be compared to tissue positions. Selecting `FALSE` requires
#' the raw matrix to be available, while `TRUE` enables using the filtered matrix,
#' only.
#'   
#' @export
get_cps <- function(visium_dir, visium_sample, matrix_name=c("filtered", "raw"), 
                    exclude_deprecated_genes=TRUE, remove_background_spots=FALSE){
  matrix_name <- match.arg(matrix_name)
  
  # Load raw matrix counts
  file_dir <- file.path(visium_dir, visium_sample, "outs")
  raw_folder <- file.path(file_dir,"raw_feature_bc_matrix")
  fm_features <- file.path(file_dir,"filtered_feature_bc_matrix", "features.tsv.gz")
  h5_raw <- file.path(file_dir,"raw_feature_bc_matrix.h5")
  
  if ((dir.exists(raw_folder) && matrix_name == "raw") || 
      (dir.exists(raw_folder) && matrix_name == "filtered" && file.exists(fm_features) )){
    message("Loading ", raw_folder)
    # Load the raw matrix
    mtx <- readMM(file = file.path(raw_folder, "matrix.mtx.gz"))
    features <- read.delim(file.path(raw_folder, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
    colnames(mtx) <- read.delim(file.path(raw_folder, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)$V1
    rownames(mtx) <- features$V1
    
    # subset to filtered matrix
    if (matrix_name == "filtered"){
      features <- read.delim(fm_features, header = FALSE, stringsAsFactors = FALSE)$V1
      mtx <- mtx[rownames(mtx) %in% features, ]}
    
  }else if (file.exists(h5_raw)){
    message("Loading ", h5_raw)
    h5 <- rhdf5::h5read(h5_raw, name = "matrix")
    # convert to mtx matrix
    mtx <- methods::new(
      "dgCMatrix",
      x = as.numeric(h5$data),
      i = as.integer(h5$indices),
      p = as.integer(h5$indptr),
      Dim = as.integer(h5$shape),
      Dimnames = list(h5$features$id, h5$barcodes))
    
    if (matrix_name == "filtered"){
      # subset to filtered matrix
      target_idx <- h5$features$target_sets$`Visium Human Transcriptome Probe Set` + 1
      mtx <- mtx[as.numeric(target_idx), , drop = FALSE]}
    
  }else{
    stop("Can't find ", h5_raw, ". ")}
  
  # exclude as deprecated marked genes
  mtx <- mtx[grepl("^DEPRECATED", rownames(mtx)), , drop=FALSE]
  cps <- Matrix::colSums(mtx)
  names(cps) <- colnames(mtx)
  
  return(cps)}
