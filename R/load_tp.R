#' Load Visium tissue positions
#'
#' Load and optionally validate a Visium tissue-position CSV file.
#'
#' @param visium_dir Character scalar. Base directory containing the Visium
#'   sample directories.
#' @param visium_sample Character scalar. Name or identifier of the Visium
#'   sample directory.
#' @param tp_file_name Character scalar. Name of the tissue-position file in
#'   the sample's \file{outs/spatial} directory. Defaults to
#'   \code{"tissue_positions.csv"}.
#' @param L_capture_area Positive numeric scalar defining the side lengths
#'  of the squared capture-area. Supported values are \code{6500} and 
#'  \code{11000}. Defaults to \code{NULL}.
#' @param check_values Logical scalar. Whether required values and coordinate
#'   types should be validated and \code{NA} values should be excluded. 
#'   Defaults to \code{TRUE}.
#'
#' @return Tissue positions as dataframe.
#'
#' @examples
#' \dontrun{
#' tp <- load_tp(
#'   visium_dir = "/path/to/visium",
#'   visium_sample = "sample_1",
#' )
#' }
#'
#' @export
load_tp <- function(visium_dir, visium_sample, tp_file_name="tissue_positions.csv", L_capture_area=NULL, check_values=TRUE){
 
  # check logical arguments
  if (!is.logical(check_values) || length(check_values) != 1L || is.na(check_values)) {
    stop("'check_values' must be one non-missing logical value.")}
  
  # check directory & file
  file_dir <- file.path(visium_dir, visium_sample, "outs", "spatial")
  tp_file <- file.path(file_dir, tp_file_name)
  
  if (!dir.exists(file_dir) || !file.exists(tp_file)){
    stop("Can't find directory ", file_dir, " with ", tp_file_name,".")}
  
  # check for tp header
  tp_cols <-  c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  first_line <- readLines(tp_file, n = 1L)
  got_header <- identical(first_line, paste(tp_cols, collapse = ","))
  
  # load tp  
  tp <- utils::read.csv(file = tp_file, stringsAsFactors = FALSE, header = got_header) 
  
  n_cols <- length(tp_cols)
  if (ncol(tp) != n_cols){
    stop("Expected ",n_cols," columns in the tissue position file, got ", ncol(tp))}
    
    # got expected tp column names
    if (!got_header){
      if (is.character(tp[,1])){
        colnames(tp)[1] <- tp_cols[1]}
      if ( all(is.numeric(unlist(tp[,2:6])))){
        if (all(tp$V2 %in% c(0,1))){
          colnames(tp)[2] <- tp_cols[2]}
        if (max(tp$V3) < max(tp$V4) && max(c(tp$V3, tp$V4)) < abs(max(c(tp$V5, tp$V6)))){
          colnames(tp)[3] <- tp_cols[3]
          colnames(tp)[4] <- tp_cols[4]
          colnames(tp)[5] <- tp_cols[5]
          colnames(tp)[6] <- tp_cols[6]}}}
    
    if (check_values){
      tp <- check_tp(tp = tp, L_capture_area = L_capture_area)}
  
  return(tp)}

