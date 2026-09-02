#' Save aligned Visium tissue positions
#'
#' Save aligned Visium tissue positions and orientation-adjustment metadata
#' as an RDS file, or as separate tissue-position and metadata CSV files.
#'
#' @param tp Data frame containing aligned Visium tissue positions.
#' @param out_dir Character scalar. Directory where the file(s) are saved.
#' @param file_name Character scalar. Base filename without extension.
#'   Defaults to \code{"aligned_tissue_positions"}.
#' @param visium_sample Optional character scalar. Sample identifier appended
#'   to the base filename. Defaults to \code{NULL}.
#' @param attr_name Character scalar. Name of the attribute containing
#'   orientation-adjustment information. Defaults to
#'   \code{"orientation_adjustments"}.
#' @param file_type Character vector specifying one or both output formats:
#'   \code{"rds"} (default) and \code{"csv"}.
#'
#' @return Invisibly returns the paths of files written.
#' 
#' @examples
#' \dontrun{
#' tp <- save_aligned_tp(
#'   tp = tp_df,
#'   out_dir = "/path/to/output",
#'   visium_sample = "Sample_1",
#' )
#' }
#'
#' @export
save_aligned_tp <- function(tp, out_dir,
    file_name = "aligned_tissue_positions",
    visium_sample = NULL,
    attr_name = "orientation_adjustments",
    file_type = "rds") {
  
  if (!is.data.frame(tp)) {
    stop("'tp' must be a data.frame object.")}
  
  # Character inputs
  check_string(out_dir, "out_dir")
  check_string(file_name, "file_name")
  check_string(attr_name, "attr_name")
  check_string(visium_sample, "visium_sample", allow_null = TRUE)
  
  file_type <- match.arg(file_type, choices = c("rds", "csv"), several.ok = TRUE)
  
  # Check orientation metadata
  attr_info <- attr(tp, attr_name)
  
  if (is.null(attr_info)) {
    stop("'tp' does not contain the '", attr_name, "' attribute. ",
      "Use utils::write.csv() if only the data frame should be saved.")}
  
  expected_attr_names <- c("adjust_array", "adjust_pixels")
  
  if (!is.list(attr_info) || !setequal(names(attr_info), expected_attr_names)) {
    stop("Attribute '", attr_name,"' must contain 'adjust_array' and 'adjust_pixels'.")}
  
  valid_adjustments <- c("swap", "x",  "y")
  
  for (x in expected_attr_names) {
    adjustments <- attr_info[[x]]
    if (is.null(adjustments)) {
      next}
    
    if (!is.character(adjustments) || anyNA(adjustments) ||
      any(!adjustments %in% valid_adjustments) ||
      anyDuplicated(adjustments)) {
      stop("'", attr_name, "$", x, "' must contain unique values selected from: ",
        paste(valid_adjustments, collapse = ", "),".")}}
  
  # Output directory
  if (!dir.exists(out_dir)) {
    message("Creating output directory ", out_dir, ".")
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)}
  
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  
  if (!is.null(visium_sample)){
    visium_sample <- paste0("_", visium_sample)}
  
  saved_files <- c()
  # RDS
  if ("rds" %in% file_type){
    rds_file <- file.path(out_dir,paste0(file_name, visium_sample,".rds"))
    saveRDS(tp, file = rds_file)
    saved_files <- rds_file}
  
  # CSV
  if ("csv" %in% file_type) {
    csv_file <- file.path(out_dir,paste0(file_name, visium_sample,".csv"))
    
    # save data frame
    utils::write.csv(tp, file = csv_file, row.names = FALSE)
    
    attr_df <- do.call(rbind, lapply(expected_attr_names, function(x) {
          adjustments <- attr_info[[x]]
          if (is.null(adjustments) ||length(adjustments) == 0L) {
            return(NULL)}
          data.frame(coordinate_type = x, orientation_adjustments = adjustments, stringsAsFactors = FALSE)}))
    
    if (is.null(attr_df)) {
      attr_df <- data.frame(
        coordinate_type = character(),
        orientation_adjustments = character())}
    
    # save attributes
    attr_file <- file.path(out_dir,paste0(file_name, visium_sample, "_", attr_name,".csv"))
    utils::write.csv(attr_df, file = attr_file, row.names = FALSE)
  
    saved_files <- c(saved_files, csv_file, attr_file)}
  
  invisible(saved_files)}

