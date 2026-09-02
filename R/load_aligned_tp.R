#' Load aligned Visium tissue positions
#'
#' Load and optionally validate the aligned Visium tissue positions and 
#' orientation-adjustment metadata.
#'
#' @param file_dir Character scalar. Directory from where the file(s) are loaded.
#' @param file_name Character scalar. Base filename without extension.
#'   Defaults to \code{"aligned_tissue_positions"}.
#' @param visium_sample Optional character scalar. Sample identifier appended
#'   to the base filename. Defaults to \code{NULL}.
#' @param attr_name Character scalar. Name of the attribute containing
#'   orientation-adjustment information. Defaults to
#'   \code{"orientation_adjustments"}.
#' @param file_type Character scalar specifying the file format:
#'   Either \code{"rds"} (default) or \code{"csv"}.   
#' @param L_capture_area Positive numeric scalar defining the side lengths
#'  of the squared capture-area. Supported values are \code{6500} and 
#'  \code{11000}. Defaults to \code{NULL}.
#' @param check_values Logical scalar. Whether required values and coordinate
#'   types should be validated and \code{NA} values should be excluded. 
#'   Defaults to \code{TRUE}.
#'
#' @return Data frame containing aligned Visium tissue positions and orientation
#' adjustments as attribute information.
#' 
#' @examples
#' \dontrun{
#' tp <- load_aligned_tp(
#'   file_dir = "/path/to/file",
#'   visium_sample = "Sample_1",
#' )
#' }
#'
#' @export
load_aligned_tp <- function(file_dir, file_name="aligned_tissue_positions", visium_sample = NULL, 
                            attr_name = "orientation_adjustments", file_type="rds", L_capture_area= NULL,
                            check_values=TRUE){
  
  # Character inputs
  check_string(file_dir, "file_dir")
  check_string(file_name, "file_name")
  check_string(attr_name, "attr_name")
  check_string(visium_sample, "visium_sample", allow_null = TRUE)
  
  file_type <- match.arg(file_type, choices = c("rds", "csv"))
  
  # check logical arguments
  if (!is.logical(check_values) || length(check_values) != 1L || is.na(check_values)) {
    stop("'check_values' must be one non-missing logical value.")}
  
  # check directory & file
  if (!is.null(visium_sample)){
    visium_sample <- paste0("_", visium_sample)}
  
  tp_file <- file.path(file_dir, paste0(file_name, visium_sample,".", file_type))
  if (!dir.exists(file_dir)){
    stop("Can't find the file directory ", file_dir, ".")}
  
  if (!file.exists(tp_file)){
    stop("Can't find aligned tissue position file ", tp_file, ".")}
  
  # --------------------------------------------------------------------------#
  
  # Load tp
  # expected adjustments
  expected_attr_names <- c("adjust_array", "adjust_pixels")
  # attribute information
  valid_adjustments <- c("swap", "x",  "y")
  
  message("Loading ", tp_file)
  
  # RDS
  if (file_type == "rds"){
    tp <- readRDS(file = tp_file)
    attr_info <- attr(tp, attr_name)
    if (is.null(attr_info)){
      stop("Found no valid attribute information for ", tp_file, ".")}}
  
  # CSV
  if (file_type == "csv"){
    tp <- utils::read.csv(file = tp_file, stringsAsFactors = FALSE)
    
    # attribute information
    tp_attr <- file.path(file_dir, paste0(file_name, visium_sample, "_",attr_name,".", file_type))
    if (!file.exists(tp_attr)){
      stop("Can't find attribute information ",tp_attr,".")}
    
    message("Loading ", tp_attr)
    tp_attr <- utils::read.csv(file = tp_attr, stringsAsFactors = FALSE)

    # check availability of expected columns
    expected_cols <- c("coordinate_type", "orientation_adjustments")
    missing_cols <- setdiff(expected_cols, colnames(tp_attr))
    if (length(missing_cols) > 0L) {
      stop("Attribute file is missing required columns: ",
           paste(missing_cols, collapse = ", "), ".")}
    
    # only accept expected_attr_names
    unknown_types <- setdiff(unique(tp_attr$coordinate_type), expected_attr_names)
    if (length(unknown_types) > 0L) {
      stop("Attribute file contains unknown coordinate type(s): ",
        paste(unknown_types, collapse = ", "),".")}
    
    # reconstruct attribute list
    attr_info <- setNames(
      lapply(expected_attr_names, function(x) {
          values <- tp_attr$orientation_adjustments[
            tp_attr$coordinate_type == x]
          if (length(values) == 0L) {
            return(NULL)}
          values}),
      expected_attr_names)}
  
  # --------------------------------------------------------------------------#
  
  # sanity checks
  if (!is.list(attr_info) || !setequal(names(attr_info), expected_attr_names)) {
    stop("Attribute '", attr_name,"' must contain 'adjust_array' and 'adjust_pixels'.")}
  
  for (x in expected_attr_names) {
    adjustments <- attr_info[[x]]
    if (is.null(adjustments)) {
      next}
    
    if (!is.character(adjustments) || anyNA(adjustments) ||
        any(!adjustments %in% valid_adjustments) ||
        anyDuplicated(adjustments)) {
      stop("'", attr_name, "$", x, "' must contain unique values selected from: ",
           paste(valid_adjustments, collapse = ", "),".")}}
  
  # tissue positions
  if (isTRUE(check_values)){
    tp <- check_tp(tp = tp, L_capture_area = L_capture_area, 
            array_columns= c("array_row", "array_col","array_x", "array_y"), 
            pixel_columns = c("pxl_row_in_fullres", "pxl_col_in_fullres", "pixels_x","pixels_y"))}

  # ensure availability of attribute information
  attr(tp, attr_name) <- attr_info
  
  return(tp)}

