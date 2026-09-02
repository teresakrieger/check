# helper functions for spatial alignment

#######################################################################################

#--------------------------------------------------------------#
# check_transformation_matrix_json()
#--------------------------------------------------------------#


#' Check the transformation matrix in manual alignment JSON file from 10X Loupe Browser
#'
#' @param json_dir Character Scalar. Directory, with the JSON file from 
#' 10X Loupe Manual Alignment. 
#' @param slide_id Character Scalar. Name of the JSON manual alignment file.
#' 
#' @return A named list containing:
#' \describe{
#'   \item{\code{json_file}}{The manual alignment JSON file information from 
#'    10X Loupe Browser comprising a named list with `transform`, `oligo`, and
#'    `fiducial` slots..}
#'   \item{\code{slide2image}}{A named list with a transformation matrix `A` 
#'    for converting slide (`x`, `y`) to image coordinates (`imageX`, `imageY`),
#'    along with the converted `oligo` and `fiducial` data.frame objects and 
#'    their maximal errors `error_oligo` and `error_fiducial` for each axis.}
#'   \item{\code{image2slide}}{A named list with the inverted transformation 
#'    matrix `A_inv` for converting image (`imageX`, `imageY`) to slide 
#'    coordinates (`x`, `y`), along with the converted `oligo_inv` and
#'    `fiducial_inv` data.frame objects and their maximal errors `error_oligo` 
#'    and `error_fiducial` for each axis.}
#' }
#' 
#' @details
#' This function reads the JSON file exported by 10X Loupe Browser for manual
#' fiducial alignment. This JSON structure is not a documented stable Space
#' Ranger output interface and may change between 10X Loupe Browser versions.
#' The function therefore validates the expected structure before extracting
#' transformation information.
#'
#' @examples
#' \dontrun{
#' check_transformation_matrix_json(
#'   json_dir = "path/to/json/file",
#'   slide_id="V11F11-111-A1"
#' )
#' }
#'
#' @keywords internal
check_transformation_matrix_json <- function(json_dir, slide_id){
  
  # check slide id
  check_string(slide_id, "slide_id", allow_null = FALSE)
  
  # load json
  if (!dir.exists(json_dir)) {
    stop("Can't find directory ", json_dir, ".")}
  
  if(!grepl("\\.json$", slide_id)){
    file_name <- paste0(slide_id, ".json")
  }else{
    file_name <- slide_id}
  
  json_file <- file.path(json_dir, file_name)
  if (!file.exists(json_file)) {
    stop("Can't find JSON file ", file_name, " in directory ", json_dir,".")}  
  
  json_file <- jsonlite::fromJSON(json_file)
  
  # check required json slots
  required_json_slots <- c("transform", "oligo", "fiducial")
  
  missing_slots <- setdiff(required_json_slots,names(json_file))
  
  if (length(missing_slots) > 0L) {
    stop("The JSON file is missing required slots: ",
         paste(missing_slots, collapse = ", "), 
         ". The Loupe manual-alignment JSON format may have changed.")}
  
  A <- json_file$transform
  oligo <- as.data.frame(json_file$oligo)
  fid <- as.data.frame(json_file$fiducial)
  
  # check transformation matrix
  if (!is.matrix(A) || !is.numeric(A) || !identical(dim(A), c(3L, 3L)) || anyNA(A) || any(!is.finite(A))) {
    stop("'transform' must be a finite numeric 3 x 3 matrix.")}
  if (!isTRUE(all.equal(unname(A[3, ]), c(0, 0, 1), tolerance = 1e-8))) {
    warning("Expected c(0, 0, 1) in transformation-matrix row 3 for verified affine transformation.")}
  # ensure transformation can be performed
  det_A <- det(A)
  if (!is.finite(det_A) || abs(det_A) < sqrt(.Machine$double.eps)) {
    stop("The transformation matrix is singular or numerically close to singular.")}
  
  # check oligo & fiducial cols
  if(!inherits(oligo, "data.frame")){
    stop("'oligo' should be a data.frame object.")}
  if(!inherits(fid, "data.frame")){
    stop("'fiducial' should be a data.frame object.")}
  
  # check required cols
  required_cols <- c("x", "y", "row", "col","imageX", "imageY")
  
  missing_cols <- setdiff(required_cols, colnames(oligo))
  if (length(missing_cols) > 0L) {
    stop("'oligo' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  missing_cols <- setdiff(required_cols, colnames(fid))
  if (length(missing_cols) > 0L) {
    stop("'fiducial' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  if (diff(range(oligo$col)) < diff(range(oligo$row))){
    warning("Observed a smaller 'col' than 'row' index range,
            which may indicate inverted axes.")}
  
  if (anyDuplicated(paste(oligo$row, oligo$col, sep = "_")) > 0L){
    stop("Observed duplicated spot positions.")}
  
  # check column values
  numeric_columns <- c("x", "y", "imageX", "imageY")
  array_columns <- c("row", "col")
  
  for (x_name in numeric_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="any", value_type="numeric")
    fid[,x_name] <-  check_numeric(x=fid[,x_name], name= paste0("fid$",x_name), n= nrow(fid), value_range="any", value_type="numeric")}
  
  for (x_name in array_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="positive", value_type="integer", allow_zero=TRUE)
    fid[,x_name] <-  check_numeric(x=fid[,x_name], name= paste0("fid$",x_name), n= nrow(fid), value_range="any", value_type="integer")}
  
  #--------------------------------------------------------------#
  
  # Forward transformation: slide -> image
  # Convert slide coordinates to 3x3 coordinates
  v_oligo <- cbind(oligo$x, oligo$y, 1)
  v_fid <- cbind(fid$x, fid$y, 1)
  
  # Transformation of spatial coordinates into pixel coordinates
  # Perform matrix-vector multiplication & re-combine with original output
  # oligo
  t_oligo <- v_oligo %*% t(A)
  pred_oligo <- oligo
  pred_oligo$pred_imageX <- t_oligo[,1L]
  pred_oligo$pred_imageY <- t_oligo[,2L] 
  # fiducials
  t_fid <- v_fid %*% t(A)
  pred_fid <- fid
  pred_fid$pred_imageX <- t_fid[,1L]
  pred_fid$pred_imageY <- t_fid[,2L] 
  
  #--------------------------------------------------------------#
  
  # linear model
  fit_imageX <- stats::lm(imageX ~ x+y, data = oligo)
  fit_imageY <- stats::lm(imageY ~ x+y, data = oligo) 
  
  A_imageX <- stats::coef(fit_imageX)
  A_imageY <- stats::coef(fit_imageY)
  
  A_lm <- matrix(c(A_imageX["x"], A_imageX["y"], A_imageX["(Intercept)"],
                       A_imageY["x"], A_imageY["y"], A_imageY["(Intercept)"],
                       0, 0, 1), nrow = 3, byrow = TRUE)
  
  A_identical <- all.equal(A, A_lm)
  if (!isTRUE(A_identical)){
    message(sample_id, ": Got distinct transformation matrices A from 'json_file$transform' and by linear model.")}
  
  #--------------------------------------------------------------#
  
  # save to list
  S2I <- list("A"=A, 
              "lm_imageX" = fit_imageX,
              "lm_imageY" = fit_imageY,
              "pred_oligo"= pred_oligo,
              "A_identical" =A_identical,
              "error_oligo"=c(max_abs_imageX = max(abs(pred_oligo$pred_imageX - pred_oligo$imageX)),
                              max_abs_imageY = max(abs(pred_oligo$pred_imageY - pred_oligo$imageY))),
              "pred_fiducial"=pred_fid,
              "error_fiducial"= c(max_abs_imageX = max(abs(pred_fid$pred_imageX - pred_fid$imageX)),
                                  max_abs_imageY = max(abs(pred_fid$pred_imageY - pred_fid$imageY))))
  
  #--------------------------------------------------------------#
  
  # Backward transformation: image -> slide
  # calculate inverse transformation matrix
  A_inv <- solve(A)
  
  # Convert pixel coordinates to 3x3 coordinates
  v_oligo <- cbind(oligo$imageX, oligo$imageY, 1)
  v_fid <- cbind(fid$imageX, fid$imageY, 1)
  
  # Transformation of spatial coordinates into pixel coordinates
  # Perform matrix-vector multiplication
  # oligo
  t_oligo <- v_oligo %*% t(A_inv)
  pred_oligo <- oligo
  pred_oligo$pred_x <- t_oligo[,1L]
  pred_oligo$pred_y <- t_oligo[,2L]
  # fiducials
  t_fid <- v_fid %*% t(A_inv)
  pred_fid <- fid
  pred_fid$pred_x <- t_fid[,1L]
  pred_fid$pred_y <- t_fid[,2L] 
  
  #--------------------------------------------------------------#
  
  # linear model
  fit_x <- stats::lm(x ~ imageX + imageY, data = oligo)
  fit_y <- stats::lm(y ~ imageX + imageY, data = oligo)  

  A_x <- stats::coef(fit_x)
  A_y <- stats::coef(fit_y)

  A_lm_inv <- matrix(c(A_x["imageX"], A_x["imageY"], A_x["(Intercept)"],
                       A_y["imageX"], A_y["imageY"], A_y["(Intercept)"],
                       0, 0, 1), nrow = 3, byrow = TRUE)
  
  A_inv_identical <- all.equal(A_inv, A_lm_inv)
  if (!isTRUE(A_inv_identical)){
    message(sample_id, ": Got distinct inverse transformation matrices A_inv from linear model and 'solve(A)'.")}
  
  #--------------------------------------------------------------#
  
  # save to list
  I2S <- list("A_inv"=A_inv,
              "lm_x" = fit_x,
              "lm_y" = fit_y,
              "oligo_inv"= pred_oligo,
              "A_inv_identical" =A_inv_identical,
              "error_oligo_inv"=  c(max_abs_x = max(abs(pred_oligo$pred_x - pred_oligo$x)),
                                    max_abs_y = max(abs(pred_oligo$pred_y - pred_oligo$y))),
              "fiducial_inv"=pred_fid,
              "error_fiducial_inv"= c(max_abs_x = max(abs(pred_fid$pred_x - pred_fid$x)),
                                      max_abs_y = max(abs(pred_fid$pred_y - pred_fid$y))))
  
  #--------------------------------------------------------------#
  
  # save data
  json_transformation <- list(json_file = json_file,
                              slide2image = S2I,
                              image2slide = I2S)
  
  return(json_transformation)}


#######################################################################################

#--------------------------------------------------------------#
# fiducial_symbols_json()
#--------------------------------------------------------------#

#' Obtain corner localization for fiducial symbols in 10X Loupe Browser manual alignment JSON file
#'
#' @param json_file List with fiducial slot, including a data.frame object with 
#' `row`, `col`, and`fidName` columns.
#' @param sample_id Character scalar, optional. Will enable printing 
#' sample-specific warnings and error messages.
#' 
#' @return A named list containing:
#' \describe{
#'   \item{\code{symbol_corners}}{A list of fiducial symbols named by their
#'   corner position.}
#'   \item{\code{symbol_df}}{The `fiducial` data.frame subset for fiducial symbols.}
#' }
#'
#' @keywords internal
fiducial_symbols_json <- function(json_file, sample_id=NULL){
  
  # check sample id
  check_string(sample_id, "sample_id", allow_null = TRUE)
  
  if (!is.null(sample_id)){
    sample_name <- paste0(sample_id, ": ")
  }else{
    sample_name <- ""}
  
  # check required json slots
  if (length(setdiff("fiducial",names(json_file))) > 0L) {
    stop(sample_name, "The JSON file is missing the required fiducial slot")}
  
  fid <- as.data.frame(json_file$fiducial)
  
  if(!inherits(fid, "data.frame")){
    stop(sample_name, "'fiducial' should be a data.frame object.")}
  
  # check required cols
  required_cols <- c("row", "col", "fidName")
  
  missing_cols <- setdiff(required_cols, colnames(fid))
  if (length(missing_cols) > 0L) {
    stop(sample_name, "'fiducial' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  # check column values
  array_columns <- c("row", "col")
  
  for (x_name in array_columns){
    fid[,x_name] <-  check_numeric(x=fid[,x_name], name= paste0("fid$",x_name), n= nrow(fid), value_range="any", value_type="integer")}
  
  
  # select fiducial markers
  fm <- fid[!is.na(fid$fidName) & nzchar(fid$fidName), ,drop = FALSE]
    if (nrow(fm) == 0L) {
    stop(sample_name, "No named fiducial symbols were found.")}

  # array coordinates
  # low row -> upper
  # high row -> lower
  # low col -> left
  # high col -> right
  row_sorted <- which.max(diff(sort(fm$row, decreasing = FALSE)))
  col_sorted <- which.max(diff(sort(fm$col, decreasing = FALSE)))
  low_row <- fm$fidName[order(fm$row, decreasing = FALSE)][seq_len(row_sorted)]
  low_col <- fm$fidName[order(fm$col, decreasing = FALSE)][seq_len(col_sorted)]
  high_row <- fm$fidName[order(fm$row, decreasing = FALSE)][(row_sorted+1):length(fm$fidName)]
  high_col <- fm$fidName[order(fm$col, decreasing = FALSE)][(col_sorted+1):length(fm$fidName)]
  
  symbol_corners <- list(Upper_left = intersect(low_row, low_col),
                         Lower_left  = intersect(high_row, low_col),
                         Lower_right = intersect(high_row, high_col),
                         Upper_right = intersect(low_row, high_col))
  symbol_df <- fm
  symbol_df$array_coordinates <- NA_character_
  for (marker in names(symbol_corners)){
    symbol_df$array_coordinates[symbol_df$fidName == symbol_corners[[marker]]] <- marker}
  return(list(symbol_corners = symbol_corners,
              symbol_df = symbol_df))}

#######################################################################################

#--------------------------------------------------------------#
# `derive_slide_coordinates_json()`
#--------------------------------------------------------------#

#' Obtain the slide coordinate information provided by 10X Loupe Browser's manual alignment JSON file
#'
#' @param json_file List object. Data from 10X Loupe Browser's manual alignment 
#' JSON file, with named list slots `transform`, `oligo`, and `fiducial`. The 
#' `transform` slot is not required by this function.
#' 
#' @return A named list with
#' \decribe{
#' \item{\code{oligo}}{Linear models (`lm_x` `lm_y`) fitting the slide coordinates 
#'  (`x`, `y`) and the array coordinates (`row`, `col`), axis-specific increments
#'  (`dx`, `dy`), and `median_ratio` with `dx`/`dy` for the oligo data.frame.}
#' \item{\code{fiducials}}{Linear models (`lm_x` `lm_y`) fitting the slide coordinates 
#'  (`x`, `y`) and the array coordinates (`row`, `col`), axis-specific increments
#'  (`dx`, `dy`), and  `median_ratio` with `dx`/`dy`for the fiducial data.frame.}
#' \item{\code{coordinate_stats}}{Relation between slide and array coordinates,
#' and axis determination.}
#' }
#'
#' @examples
#' \dontrun{
#' derive_slide_coordinates_json(
#'   json_file = json_data
#' )
#' 
#' # or if returned object from 'check_transformation_matrix_json()' is used:
#' derive_slide_coordinates_json(
#'   json_file = json_data$json_file
#' )
#' } 
#'
#' @keywords internal
derive_slide_coordinates_json <- function(json_file){
  
  # check required json slots
  required_json_slots <- c("oligo", "fiducial")
  
  missing_slots <- setdiff(required_json_slots,names(json_file))
  
  if (length(missing_slots) > 0L) {
    stop("The JSON file is missing required slots: ",
         paste(missing_slots, collapse = ", "), ".")}
  
  oligo <- as.data.frame(json_file$oligo)
  fid <- as.data.frame(json_file$fiducial)
  
  # check oligo & fiducial cols
  if(!inherits(oligo, "data.frame")){
    stop("'oligo' should be a data.frame object.")}
  if(!inherits(fid, "data.frame")){
    stop("'fiducial' should be a data.frame object.")}
  
  # check required cols
  required_cols <- c("x", "y", "row", "col")
  
  missing_cols <- setdiff(required_cols, colnames(oligo))
  if (length(missing_cols) > 0L) {
    stop("'oligo' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  missing_cols <- setdiff(required_cols, colnames(fid))
  if (length(missing_cols) > 0L) {
    stop("'fiducial' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  if (diff(range(oligo$col)) < diff(range(oligo$row))){
    warning("Observed a smaller 'col' than 'row' index range, which may indicate inverted axes.")}
  
  if (anyDuplicated(paste(oligo$row, oligo$col, sep = "_")) > 0L){
    stop("Observed duplicated spot positions.")}
  
  # check column values
  numeric_columns <- c("x", "y")
  array_columns <- c("row", "col")
  
  for (x_name in numeric_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="any", value_type="numeric")
    fid[,x_name] <-  check_numeric(x=fid[,x_name], name= paste0("fid$",x_name), n= nrow(fid), value_range="any", value_type="numeric")}
  
  for (x_name in array_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="positive", value_type="integer", allow_zero=TRUE)
    fid[,x_name] <-  check_numeric(x=fid[,x_name], name= paste0("fid$",x_name), n= nrow(fid), value_range="any", value_type="integer")}
  
  #--------------------------------------------------------------#
  
  # coordinate system
  x_stats <- oligo %>%
    dplyr::group_by(col) %>%
    dplyr::reframe(
      range_x = diff(range(x)),
      range_y = diff(range(y)),
      expected = diff(range(y)) > diff(range(x)),
      inverted = diff(range(x)) > diff(range(y)))
  
  y_stats <- oligo %>%
    dplyr::group_by(row) %>%
    dplyr::reframe(
      range_x = diff(range(x)),
      range_y = diff(range(y)),
      expected = diff(range(x)) > diff(range(y)),
      inverted = diff(range(y)) > diff(range(x)))
  
  if (median(x_stats$range_y) > median(x_stats$range_x) &&
      median(y_stats$range_x) > median(y_stats$range_y)){
    # expected axes
    x_axis <- "col"
    y_axis <- "row"
  }else if (median(x_stats$range_x) > median(x_stats$range_y) &&
            median(y_stats$range_y) > median(y_stats$range_x)) {
    # inverted axes
    x_axis <- "row"
    y_axis <- "col"
  }else{
    stop("Unexected array grid.")}
  
  oligo$arrayX <- oligo[[x_axis]]
  fid$arrayX <- fid[[x_axis]]
  oligo$arrayY <- oligo[[y_axis]]
  fid$arrayY <- fid[[y_axis]]
  
  if (sum(x_stats$inverted, y_stats$inverted)>0L){
    warning("Observed ", sum(x_stats$inverted, y_stats$inverted), " inverted slide distances across coordinate array axes.")}
  
  # spatial models
  oligo_x <- stats::lm(x ~ arrayX, data = oligo)
  oligo_y <- stats::lm(y ~ arrayY, data = oligo)
  fid_x <- stats::lm(x ~ arrayX, data = fid)
  fid_y <- stats::lm(y ~ arrayY, data = fid)
  
  slide_coordinates <- list(oligo= list(lm_x = oligo_x,
                                        dx = abs(unname(coef(oligo_x)["arrayX"])),
                                        median_abs_rx = median(abs(stats::residuals(oligo_x))),
                                        max_abs_rx = max(abs(stats::residuals(oligo_x))),
                                        
                                        lm_y = oligo_y,
                                        dy = abs(unname(coef(oligo_y)["arrayY"])),
                                        median_abs_ry = median(abs(stats::residuals(oligo_y))),
                                        max_abs_ry = max(abs(stats::residuals(oligo_y))),
                                        xy_ratio = abs(unname(coef(oligo_x)["arrayX"])/
                                                         unname(coef(oligo_y)["arrayY"]))),
                            fiducials= list(lm_x = fid_x,
                                            dx = abs(unname(coef(fid_x)["arrayX"])),
                                            lm_y = fid_y,
                                            dy = abs(unname(coef(fid_y)["arrayY"])),
                                            xy_ratio = abs(unname(coef(fid_x)["arrayX"])/
                                                             unname(coef(fid_y)["arrayY"]))),
                            coordinate_stats = list(
                              x_stats = x_stats,
                              y_stats = y_stats,
                              standard_mapping = c(x = "col", y = "row"),
                              inferred_mapping = c(x = x_axis, y = y_axis)))
  
  return(slide_coordinates)}


#######################################################################################

#--------------------------------------------------------------#
# `combine_tp_json()`
#--------------------------------------------------------------#

#' Combine the 10X Loupe Browser's manual alignment JSON file and the tissue positions file
#'
#' @param json_file List object. Data from JSON manual alignment file, with named
#'  list slots `transform`, `oligo`, and `fiducial`. Only `oligo` is required for this function.
#' @param base_dir Base directory with a folder `sample_id`containing the tissue position and 
#' `scalefactors_json.json` file.
#' @param sample_id Name of the folder in `base_dir` containing the tissue position and 
#' `scalefactors_json.json` file.
#' 
#' @return A named list with
#' \decribe{
#' \item{\code{tp_json}}{Tissue positions combined with `oligo` information from json file.}
#' \item{\code{scalefactors}}{Data from `scalefactors_json.json` file.}
#' }
#'
#' @examples
#' \dontrun{
#' derive_slide_coordinates_json(
#'   json_file = json_data
#' )
#' 
#' # or if returned object from 'check_transformation_matrix_json()' is used:
#' derive_slide_coordinates_json(
#'   json_file = json_data$json_file
#' )
#' } 
#'
#' @keywords internal
combine_tp_json <- function(json_file, base_dir, sample_id){
  
  # check sample id
  check_string(sample_id, "sample_id", allow_null = TRUE)
  
  # check json_file
  if (!"oligo" %in% names(json_file)) {
    stop(sample_id,": The JSON file is missing the required oligo slot.")}
  
  required_oligo <- c("row", "col")
  oligo <- as.data.frame(json_file$oligo)
  if (!inherits( oligo, "data.frame")|| length(setdiff(required_oligo, colnames(oligo)))> 0L){
    stop(sample_id,":`oligo` should be a data frame with a `row` and `col` column.")}
  
  # check column values
  numeric_columns <- c("x", "y", "imageX", "imageY")
  array_columns <- c("row", "col")
  
  for (x_name in numeric_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="any", value_type="numeric")}
  
  for (x_name in array_columns){
    oligo[,x_name] <-  check_numeric(x=oligo[,x_name], name= paste0("oligo$",x_name), n= nrow(oligo), value_range="positive", value_type="integer", allow_zero=TRUE)}
  
  # check directory & file
  file_dir <- file.path(base_dir, sample_id)
  if (!dir.exists(file_dir)){
    stop(sample_id,": Sample-specific file directory ", file_dir, " can't be found.")}
  # tp file
  tp_file <- grep("tissue_positions", list.files(file_dir), value = TRUE)
  if (length(tp_file) != 1L){
    stop(sample_id,": Can't find exactly one `tissue_positions` file in sample-specific directory ", file_dir, ".")}
  # sf file
  sf_file <- file.path(file_dir, "scalefactors_json.json")
  if (!file.exists(sf_file)){
    stop(sample_id,": Can't find `scalefactors_json.json` in sample-specific directory ", file_dir, ".")}
  
  # load tp & check for tp header
  tp_cols <-  c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  first_line <- readLines(file.path(file_dir, tp_file), n = 1L)
  got_header <- identical(first_line, paste(tp_cols, collapse = ","))
  # load tp  
  tp <- utils::read.csv(file = file.path(file_dir, tp_file), stringsAsFactors = FALSE, header = got_header) 
  if (ncol(tp) != length(tp_cols)){
    stop(sample_id, ": Expected ", length(tp_cols), " got ", ncol(tp))}
  colnames(tp) <- tp_cols
  
  # check column values
  numeric_columns <- c("pxl_row_in_fullres", "pxl_col_in_fullres")
  array_columns <- c("array_row", "array_col")
  
  for (x_name in numeric_columns){
    tp[,x_name] <-  check_numeric(x=tp[,x_name], name= paste0("tp$",x_name), n= nrow(tp), value_range="any", value_type="numeric")}
  
  for (x_name in array_columns){
    tp[,x_name] <-  check_numeric(x=tp[,x_name], name= paste0("tp$",x_name), n= nrow(tp), value_range="positive", value_type="integer", allow_zero=TRUE)}
  
  # load sf & check required slots
  sf <-  jsonlite::fromJSON(sf_file)
  sf_cols <-  c("tissue_hires_scalef", "tissue_lowres_scalef")
  if (!is.list(sf) || length(setdiff(sf_cols, names(sf)))>0L){
    stop(sample_id,": Expected a list with ",paste(sf_cols, collapse= " and "), " slots for scalefactors_json.json. 
         \nGot a ", class(sf), " with ", paste(setdiff(sf_cols, names(sf)), collapse=" and "), ".")} 
  optional_cols <- c("fiducial_diameter_fullres", "spot_diameter_fullres")
  missing_cols <- setdiff(optional_cols, names(sf))
  if (length(missing_cols)>0L){
    warning(sample_id,": scalefactors_json.json is missing expected slots ", 
            paste(missing_cols, collapse=" and "), ".")} 

  # convert full resolution pixels into image resolution pixels
  tp$pxl_col_in_hires <- tp$pxl_col_in_fullres*sf$tissue_hires_scalef
  tp$pxl_row_in_hires <- tp$pxl_row_in_fullres*sf$tissue_hires_scalef
  tp$pxl_col_in_lowres <- tp$pxl_col_in_fullres*sf$tissue_lowres_scalef
  tp$pxl_row_in_lowres <- tp$pxl_row_in_fullres*sf$tissue_lowres_scalef
  
  # combine with json information
  # ensure unique spot positions for reliable merging
  tp_duplicated <- anyDuplicated(tp[, c("array_row", "array_col")]) > 0L
  oligo_duplicated <- anyDuplicated(oligo[,required_oligo]) > 0L
  if (tp_duplicated || oligo_duplicated){
    stop(sample_id, ": Detected duplicated spot positions in ", 
         paste(c("tp"[tp_duplicated], "oligo"[oligo_duplicated]), collapse = " and "), ".")}
  if (nrow(tp) != nrow(oligo)){
    stop(sample_id, ": Observed distinct number of entries in tp (", nrow(tp),") and oligo (", nrow(oligo),").")
  }else{
    n_entries <- nrow(tp)}
  tp  <- merge(tp, oligo, by.x = c("array_row", "array_col"), by.y= c("row", "col"))
  if (nrow(tp) != n_entries){
    stop(sample_id, ": Got ", n_entries, " before merging tp and oligo, but only ", nrow(tp), " could be matched.")}
  return(list("tp_json"= tp,
              "n_spots" = list("tp"= n_entries,
                               "json_oligo"=n_entries,
                               "tp_json"= nrow(tp)),
              "scalefactors"= sf))}


#######################################################################################

#--------------------------------------------------------------#
# `pixel_conversion_json_tp()`
#--------------------------------------------------------------#
#'
#' Calculate the transformation matrices for converting JSON image coordinates
#' into tissue position pixel coordinates (B) and backwards (B_inv)
#'
#' @param tp_json Data.frame with `array_row`, `array_col`,`imageX`,`imageY`,
#'  `pxl_row_in_fullres`, `pxl_col_in_fullres` columns.
#' @param sample_id Character scalar, optional. Will enable printing 
#' sample-specific warnings and error messages.
#' 
#' @return A named list containing:
#' \describe{
#'   \item{\code{B}}{Transformation matrix B for converting json image 
#'   coordinates into tp pixel coordinates.}
#'   \item{\code{B_inv}}{Inverse Transformation matrix B_inv for converting 
#'   tp pixel coordinates into json image coordinates.}
#' }
#'
#' @keywords internal
pixel_conversion_json_tp <- function(tp_json, sample_id=NULL){
  
  # check sample id
  check_string(sample_id, "sample_id", allow_null = TRUE)
  
  if (!is.null(sample_id)){
    sample_name <- paste0(sample_id, ": ")
  }else{
    sample_name <- ""}
  
  # check tp_json
  if(!inherits(tp_json, "data.frame")){
    stop("'tp_json' should be a data.frame object.")}
  
  # check required cols
  required_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres", "imageX", "imageY")
  
  missing_cols <- setdiff(required_cols, colnames(tp_json))
  if (length(missing_cols) > 0L) {
    stop(sample_name, "'tp_json' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  # check column values
  numeric_columns <- c("pxl_row_in_fullres", "pxl_col_in_fullres", "imageX", "imageY")
  
  for (x_name in numeric_columns){
    tp_json[,x_name] <-  check_numeric(x=tp_json[,x_name], name= paste0("tp_json$",x_name), n= nrow(tp_json), value_range="any", value_type="numeric")}
  
  #--------------------------------------------------------------#
  # Transformation of json image coordinates into tp pixel coordinates
  # fit linear model 
  fit_pxl_row <- stats::lm(pxl_row_in_fullres ~ imageX + imageY, data = tp_json)
  fit_pxl_col <- stats::lm(pxl_col_in_fullres ~ imageX + imageY, data = tp_json)
  
  B_pxl_row <- stats::coef(fit_pxl_row)
  B_pxl_col <- stats::coef(fit_pxl_col)
  
  # build transformation matrix B
  B <- matrix(c(B_pxl_row["imageX"], B_pxl_row["imageY"], B_pxl_row["(Intercept)"],
                B_pxl_col["imageX"], B_pxl_col["imageY"], B_pxl_col["(Intercept)"],
                0, 0, 1), nrow = 3, byrow = TRUE)
  
  #--------------------------------------------------------------#
  
  # get max prediction error
  # Convert slide coordinates to 3x3 coordinates
  v_image <- cbind(tp_json$imageX, tp_json$imageY, 1)
  t_pixels <- v_image %*% t(B)
  pred_pixels <- tp_json
  pred_pixels$pred_pxl_row <- t_pixels[,1L]
  pred_pixels$pred_pxl_col <- t_pixels[,2L] 
  
  # save to list
  JSON2TP <- list("B"=B,
                  "lm_pxl_row" = fit_pxl_row,
                  "lm_pxl_col" = fit_pxl_col,
                  "pred_pixels" = pred_pixels,
                  "error_pixels"=  c(max_abs_pxl_row = max(abs(pred_pixels$pred_pxl_row - tp_json$pxl_row_in_fullres)),
                                     max_abs_pxl_col = max(abs(pred_pixels$pred_pxl_col - tp_json$pxl_col_in_fullres))))
  
  #--------------------------------------------------------------#
  
  # build inverse transformation matrix B_inv
  det_B <- det(B)    
  if (!is.finite(det_B) || abs(det_B) < sqrt(.Machine$double.eps)) {
    stop(sample_name, "The fitted image-to-pixel transformation matrix is singular or numerically close to singular.")}
  
  B_inv <- solve(B)
  
  #--------------------------------------------------------------#
  
  # Transformation of  tp pixel coordinates into json image coordinates
  # fit linear model
  fit_imageX <- stats::lm(imageX ~ pxl_row_in_fullres + pxl_col_in_fullres, data = tp_json)
  fit_imageY <- stats::lm(imageY ~ pxl_row_in_fullres + pxl_col_in_fullres, data = tp_json)
  
  B_imageX <- stats::coef(fit_imageX)
  B_imageY <- stats::coef(fit_imageY)
  
  # build inverse transformation matrix B_lm_inv from linear model
  B_lm_inv <- matrix(c(B_imageX["pxl_row_in_fullres"], B_imageX["pxl_col_in_fullres"], B_imageX["(Intercept)"],
                       B_imageY["pxl_row_in_fullres"], B_imageY["pxl_col_in_fullres"], B_imageY["(Intercept)"],
                       0, 0, 1), nrow = 3, byrow = TRUE)
  
  # check identical results
  B_inv_identical <- all.equal(B_inv, B_lm_inv)
  if (!isTRUE(B_inv_identical)){
    warning(sample_name, "The independently fitted reverse transformation differs from solve(B).")}
  
  #--------------------------------------------------------------#
  
  # get max prediction error 
  pixels <- cbind(tp_json$pxl_row_in_fullres, tp_json$pxl_col_in_fullres, 1)
  t_image <- pixels %*% t(B_inv)
  pred_image <- tp_json
  pred_image$pred_imageX <- t_image[,1L]
  pred_image$pred_imageY <- t_image[,2L] 
  
  # save to list
  TP2JSON <- list("B_inv"=B_inv, 
                  "lm_imageX" = fit_imageX,
                  "lm_imageY" = fit_imageY,
                  "B_inv_identical" = B_inv_identical,
                  "pred_image_inv"= pred_image,
                  "error_image_inv"=c(max_abs_imageX = max(abs(pred_image$pred_imageX - tp_json$imageX)),
                                  max_abs_imageY = max(abs(pred_image$pred_imageY - tp_json$imageY))))
                  
                  return(list(
                    JSON2TP = JSON2TP,
                    TP2JSON = TP2JSON))}


#######################################################################################

#--------------------------------------------------------------#
# get ranges from 'symbol_df'
#--------------------------------------------------------------#

get_ranges <- function(symbol_df,col){
  range_col <- strsplit(symbol_df[,col], split = "-")
  names(range_col) <- rownames(symbol_df)
  lapply(range_col, function(x){
    diff(as.numeric(x))})}

