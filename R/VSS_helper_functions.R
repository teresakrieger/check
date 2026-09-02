
#' capitalize the first letter
#' @keywords internal
to_capital <- function(x){
   paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))}

#####################################################################

#' @keywords internal
calc_abs_diff <- function(tp, col_name){
  abs(diff(c(min(tp[[col_name]]), max(tp[[col_name]]))))}

#####################################################################

# load image file

#' @keywords internal
load_img <- function(img_file, reduce_img_size=TRUE, img_height = 600, use_terra=TRUE, resize_method="bilinear"){
  message("Loading ", img_file)
  file_type <- tolower(tools::file_ext(img_file))
  if (file_type %in% c("jpg", "jpeg")){
  if (!requireNamespace("jpeg", quietly = TRUE)){
    stop(
      "Package 'jpeg' is required to load image file.",
      "It can be installed with install.packages(\"jpeg\").")}
  p1 <- jpeg::readJPEG(img_file)
}else if (identical(file_type,"png")){
  if (!requireNamespace("png", quietly = TRUE)){
    stop("Package 'png' is required to load image file.",
         "It can be installed with install.packages(\"png\").")}
  p1 <- png::readPNG(img_file)
}else{
  stop("Unsupported image extension: '.", file_type,
       "'. Supported formats are .jpg/.jpeg or .png")}
  
  # reduce image size
  if (reduce_img_size){
    p1 <- adjust_img_size(img=p1, h = img_height, use_terra = use_terra, method = resize_method)}
  
  # convert into normalized array
  if (inherits(p1, "SpatRaster")) {
    img_arr <- terra::as.array(p1)
  }else{
    img_arr <- as.array(p1)}
  max_val <- max(img_arr, na.rm = TRUE)
  if (!is.finite(max_val)|| !max_val > 0){
    stop("The detected image contains no valid positive entries.")}
  if (max_val > 1){
    img_arr <- img_arr / max_val}
  img_arr[img_arr < 0] <- 0
  img_arr[img_arr > 1] <- 1
  
  # create plot
  p1 <-   ggplot2::ggplot() +
    ggplot2::annotation_raster(
      raster = img_arr,
      xmin = 0,
      xmax = dim(img_arr)[2],
      ymin = 0,
      ymax = dim(img_arr)[1],
      interpolate = FALSE) +
    ggplot2::coord_fixed(
      xlim = c(0,  dim(img_arr)[2]),
      ylim = c(0,  dim(img_arr)[1]),
      expand = FALSE)
  
return(p1)} 

#####################################################################

#' Check numeric input
#'
#' @param x Object to validate.
#' @param name Character scalar naming the argument.
#' @param n Optional positive integer specifying the required length.
#' @param value_range Permitted numeric range: \code{"positive"},
#'   \code{"negative"}, or \code{"any"}.
#' @param value_type Required numeric type: \code{"numeric"} or
#'   \code{"integer"}.
#' @param allow_zero Logical scalar indicating whether zero is permitted for
#'   positive or negative ranges. Defaults to \code{FALSE}.
#'
#' @return The validated numeric vector.
#'
#' @keywords internal
check_numeric <- function(x, name, n=NULL, value_range = c("positive", "negative", "any"), value_type = c("numeric", "integer"), allow_zero=FALSE) {
  value_range <- match.arg(value_range)
  value_type <- match.arg(value_type)
  
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("'name' must be one non-empty character value.")}
  
  if (is.null(x)){
    stop("No benefit of checking NULL as numeric. Please provide a proper value for 'x'.")}
  
  if (value_range != "any"){
      if (!is.logical(allow_zero) || length(allow_zero) != 1L || is.na(allow_zero)) {
    stop("'allow_zero' must be one non-missing logical value.")}}

  if (is.null(n)) {
    n <- length(x)
  } else {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) || n < 1 || n != as.integer(n)) {
      stop("'n' must be NULL or one positive integer-valued number.")}
    n <- as.integer(n)}
  
  if (is.null(n) || !is.integer(n) || length(n) != 1L || n < 1){
    warning("No valid n provided. Using length of ", name, ".")
    n <- length(x)}
  
  x <- suppressWarnings(as.numeric(x))
  if (value_type == "integer"){
    y <- as.integer(x)
  }else{
    y <- x}
  
  range_valid <- switch(
    value_range,
    positive = if (allow_zero) x >= 0 else x > 0,
    negative = if (allow_zero) x <= 0 else x < 0,
    any = rep(TRUE, length(x)))
  
  if (length(x) != n || any(is.na(x)) || any(!is.finite(x)) || any(y != x) || any(!range_valid)) {
    range_description <- switch(
      value_range,
      positive = if (allow_zero) "non-negative " else "positive ",
      negative = if (allow_zero) "non-positive " else "negative ",
      any = "")
    
    stop(name," must contain exactly ", n, " valid ",range_description, value_type," value",
         if (n == 1L) "" else "s",".")}
  return(y)}

#####################################################################

#' @keywords internal
check_string <- function(x, name, allow_null = FALSE) {
  
  if (allow_null && is.null(x)) {
    return(invisible(NULL))}
  
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("'", name, "' must be one non-empty character value.")}}

#####################################################################

#' Check Visium tissue positions
#'
#' @keywords internal
check_tp <- function(tp, L_capture_area=NULL, 
                     tp_columns=c("barcode", "in_tissue"), 
                     array_columns =c("array_row","array_col"),
                     pixel_columns = c("pxl_row_in_fullres","pxl_col_in_fullres")){
  
  if (!is.data.frame(tp)) {
    stop("'tp' must be a data frame.")}
  
  all_columns <- c(tp_columns, array_columns, pixel_columns)
  missing_columns <- setdiff(all_columns, colnames(tp))
  if (length(missing_columns) > 0L) {
    stop("'tp' is missing required columns: ", paste(missing_columns, collapse = ", "), ".")}
  
  n_cols <- length(all_columns)
  if (ncol(tp) != n_cols){
    stop("Expected ",n_cols," columns in the tissue position file, got ", ncol(tp))}
  
  if (nrow(tp) == 0L){
    stop("'tp' doesn't contain any entry.")}
  
  if (anyNA(tp[all_columns])) {
    na_entry <- paste(all_columns[unlist(lapply(all_columns, function(x) any(is.na(tp[x]))))], collapse = ", ")
    n_before <- nrow(tp)
    complete_rows <- stats::complete.cases(tp[, all_columns, drop = FALSE])
    tp <- tp[complete_rows, ,drop = FALSE]
    if (nrow(tp) == 0L){
      stop("After excluding NA entries, no tissue positions remain.")}
    warning("Detected NA values in ",na_entry," column. 
               \n", sum(!complete_rows), " incomplete tissue-position row(s) were excluded; ",
            nrow(tp), " of ", n_before," rows remain.")}
  
  if (anyDuplicated(tp[c("array_row", "array_col")]) || anyDuplicated(tp$barcode)) {
    stop("Each barcode and each combination of 'array_row' and 'array_col' must be unique.")}
  
  # Valid Visium positions have matching row/column parity.
  invalid_parity <- (tp$array_row %% 2L != tp$array_col %% 2L)
  if (any(invalid_parity)) {
    warning("'tp' contains ", sum(invalid_parity)," row(s) with incompatible  'array_row' and 'array_col' parity.")}
  
  # value type
  # check columns
  if (!is.character(tp$barcode) || anyNA(tp$barcode) || any(!nzchar(tp$barcode))) {
    stop("'barcode' must contain non-empty character values.")}
  
  for (x_name in pixel_columns){
    tp[,x_name] <-  check_numeric(x=tp[,x_name], name=x_name, n= nrow(tp), value_range="any", value_type="numeric")}
  
  for (x_name in array_columns){
    tp[,x_name] <-  check_numeric(x=tp[,x_name], name=x_name,  n= nrow(tp), value_range="positive", value_type="integer", allow_zero = TRUE)}
  
  tp$in_tissue <- check_numeric(x = as.character(tp$in_tissue), name = "in_tissue", n = nrow(tp),
                                value_range = "positive", value_type = "integer", allow_zero = TRUE)
  
  if (!all(tp$in_tissue %in% c(0L, 1L))) {
    stop("'in_tissue' must contain only 0 and 1.")}
  
  if (!is.null(L_capture_area)){
    L_capture_area <- check_numeric(x= L_capture_area, name="L_capture_area", n=1L, value_range="positive", value_type="numeric") 
    if(L_capture_area %in% c(6500, 11000)){
      if (L_capture_area == 6500){
        rows <- 0:77
        columns <- 0:127}
      if (L_capture_area == 11000){
        rows <- 0:127
        columns <- 0:223}
      
      n_rows <- length(rows)
      n_cols <- length(columns)
      
      # Ensure compatibility
      expected_spots <- n_rows * n_cols / 2L
      if (nrow(tp) != expected_spots) {
        warning("'tp' contains ", nrow(tp), " positions, whereas ",expected_spots," spots are
        expected for the full (",L_capture_area,"µm)²-array.")}
      
      unknown_rows <- c(setdiff(unique(tp$array_row), rows))
      unknown_cols <- c(setdiff(unique(tp$array_col), columns))
      if (length(unknown_rows) > 0L) {
        stop("'tp$array_row' contains values not represented in 'rows': ",
             paste(unknown_rows, collapse = ", "), ".")}
      
      if (length(unknown_cols) > 0L) {
        stop("'tp$array_col' contains values not represented in 'columns': ",
             paste(unknown_cols, collapse = ", "),".")}
        
      }else{
        warning("Expected a capture area of 6500 or 11000µm.")}}
  
  return(tp)}

#####################################################################

#' @keywords internal
fill_df <- function(df, row_index, img){
  df[row_index,"xmin_img"] <- terra::xmin(img)
  df[row_index,"xmax_img"] <- terra::xmax(img)
  df[row_index,"ymin_img"] <- terra::ymin(img)
  df[row_index,"ymax_img"] <- terra::ymax(img)
  return(df)}

#####################################################################

# adjust rgb image to be used in ggplot
#' @keywords internal
raster_to_ggdf <- function(r) {
  stopifnot(terra::nlyr(r) == 3)
  
  df <- as.data.frame(r, xy = TRUE)
  names(df) <- c("x", "y", "R", "G", "B")
  
  # Normalize to 0-1
  max_val <- max(df$R, df$G, df$B, na.rm = TRUE)
  df <- df %>%
    dplyr::mutate(
      Rn = R / max_val,
      Gn = G / max_val,
      Bn = B / max_val,
      hex = rgb(Rn, Gn, Bn))
  return(df)}

#####################################################################

#' @keywords internal
flip_image_df <- function(df_img, flip_x=FALSE, flip_y=FALSE) {
  if (flip_x) {
    df_img$x <- sum(range(df_img$x, na.rm = TRUE)) - df_img$x}
  if (flip_y) {
    df_img$y <- sum(range(df_img$y, na.rm = TRUE)) - df_img$y}
  return(df_img)}

#####################################################################

#' @keywords internal
validate_dimension <- function(x, argument_name) {
  if (is.null(x)) {
    return(invisible(NULL))}
  
  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < 1 ||
      x != round(x)) {
    stop("'", argument_name,
         "' must be NULL or one positive integer.")}
  
  return(invisible(NULL))}

#####################################################################
