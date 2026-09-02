#' Obtain axis-specific pixel deviations for the tissue positions file
#'
#' @param tp_json Data.frame with `array_row`, `array_col`,`x`,`y`,
#'  `pxl_row_in_fullres`, `pxl_col_in_fullres` columns.
#' 
#' @return A named list with pixel deviation information:
#' \decribe{
#' \item{\code{lm_array}}{List object with linear models fitting the 
#' `pxl_row_in_fullres` and the `pixel_col_in_fullres` coordinates to the 
#' array grid.}
#' \item{\code{lm_slide}}{List object with linear models fitting the 
#' `pxl_row_in_fullres` and the `pixel_col_in_fullres` coordinates to the 
#' slide coordinates.}
#' \item{\code{residual_df}}{`tp_json` data.frame extended by 4 columns, namely
#'  `dev_pxl_row`, `dev_pxl_col`, `dev_pxl_x`, and `dev_pxl_y`, comprising the 
#'  residuals for each linear model.}
#' \item{\code{residual_summary}}{Summary of the `residual_df`.}
#' }
#'
#' @export
pixel_deviations_tp <- function(tp){
  
  # check tp
  tp <- check_tp(tp)
  
  required_cols <- c("array_row", "array_col","x","y", "pxl_row_in_fullres", "pxl_col_in_fullres")
  
  missing_cols <- setdiff(required_cols, colnames(tp_json))
  if (length(missing_cols) > 0L) {
    stop("'tp_json' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  #--------------------------------------------------------------#
  
  # pixel deviations
  lm_row_array <- stats::lm(pxl_row_in_fullres ~ array_row + array_col, data = tp_json)
  lm_row_slide <- stats::lm(pxl_row_in_fullres ~ x + y, data = tp_json)
  lm_col_array <- stats::lm(pxl_col_in_fullres ~ array_row + array_col, data = tp_json)
  lm_col_slide <- stats::lm(pxl_col_in_fullres ~ x + y, data = tp_json)
  
  dev_df <- tp_json
  dev_df$dev_pxl_row_array <- stats::residuals(lm_row_array)
  dev_df$dev_pxl_col_array <- stats::residuals(lm_col_array)
  dev_df$dev_pxl_row_slide <- stats::residuals(lm_row_slide)
  dev_df$dev_pxl_col_slide <- stats::residuals(lm_col_slide)
  
  residual_summary <- list("array"=list(pxl_row_in_fullres = summarize_residuals(dev_df$dev_pxl_row),
                                        pxl_col_in_fullres = summarize_residuals(dev_df$dev_pxl_col)),
                           "slide"=list(pxl_row_in_fullres = summarize_residuals(dev_df$dev_pxl_x),
                                        pxl_col_in_fullres = summarize_residuals(dev_df$dev_pxl_y)))
  
  pixel_deviations <- list(
    lm_array = list(pxl_row_in_fullres = lm_row,
                    pxl_col_in_fullres= lm_col),
    lm_slide = list(pxl_row_in_fullres = lm_x,
                    pxl_col_in_fullres= lm_y),
    residual_df = dev_df,
    residual_summary = residual_summary)

  return(pixel_deviations)}


#######################################################################################

#--------------------------------------------------------------#
# `pixel_deviations()`
#--------------------------------------------------------------#


#' Obtain axis-specific pixel deviations for the tissue positions file
#'
#' @param tp_json Data.frame with `array_row`, `array_col`,`x`,`y`,
#'  `pxl_row_in_fullres`, `pxl_col_in_fullres` columns.
#' 
#' @return A named list with pixel deviation information:
#' \decribe{
#' \item{\code{lm_array}}{List object with linear models fitting the 
#' `pxl_row_in_fullres` and the `pixel_col_in_fullres` coordinates to the 
#' array grid.}
#' \item{\code{lm_slide}}{List object with linear models fitting the 
#' `pxl_row_in_fullres` and the `pixel_col_in_fullres` coordinates to the 
#' slide coordinates.}
#' \item{\code{residual_df}}{`tp_json` data.frame extended by 4 columns, namely
#'  `dev_pxl_row`, `dev_pxl_col`, `dev_pxl_x`, and `dev_pxl_y`, comprising the 
#'  residuals for each linear model.}
#' \item{\code{residual_summary}}{Summary of the `residual_df`.}
#' }
#'
#' @keywords internal
pixel_deviations <- function(tp_json){
  
  if(!inherits(tp_json, "data.frame")){
    stop("'tp_json' should be a data.frame object.")}
  
  # check required cols
  required_cols <- c("array_row", "array_col","x","y", "pxl_row_in_fullres", "pxl_col_in_fullres")
  
  missing_cols <- setdiff(required_cols, colnames(tp_json))
  if (length(missing_cols) > 0L) {
    stop("'tp_json' is missing required columns: '", paste(missing_cols, collapse = ", "))}
  
  #--------------------------------------------------------------#
  
  # pixel deviations
  lm_row <- stats::lm(pxl_row_in_fullres ~ array_row + array_col, data = tp_json)
  lm_x <- stats::lm(pxl_row_in_fullres ~ x + y, data = tp_json)
  lm_col <- stats::lm(pxl_col_in_fullres ~ row + col, data = tp_json)
  lm_y <- stats::lm(pxl_col_in_fullres ~ x + y, data = tp_json)
  
  dev_df <- tp_json
  dev_df$dev_pxl_row <- stats::residuals(lm_row)
  dev_df$dev_pxl_col <- stats::residuals(lm_col)
  dev_df$dev_pxl_x <- stats::residuals(lm_x)
  dev_df$dev_pxl_y <- stats::residuals(lm_y)
  
  residual_summary <- list("array"=list(pxl_row_in_fullres = summarize_residuals(dev_df$dev_pxl_row),
                                        pxl_col_in_fullres = summarize_residuals(dev_df$dev_pxl_col)),
                           "slide"=list(pxl_row_in_fullres = summarize_residuals(dev_df$dev_pxl_x),
                                        pxl_col_in_fullres = summarize_residuals(dev_df$dev_pxl_y)))
  
  pixel_deviations <- list(
    lm_array = list(pxl_row_in_fullres = lm_row,
                    pxl_col_in_fullres= lm_col),
    lm_slide = list(pxl_row_in_fullres = lm_x,
                    pxl_col_in_fullres= lm_y),
    residual_df = dev_df,
    residual_summary = residual_summary)
  
  return(pixel_deviations)}

#--------------------------------------------------------------#
# `summarize_residuals()`
#--------------------------------------------------------------#

# summarize pixel residuals
summarize_residuals <- function(x) {
  abs_x <- abs(x)
  
  c(median_abs = stats::median(abs_x),
    mad = stats::mad(x),
    rmse = sqrt(mean(x^2)),
    q95_abs = unname(stats::quantile(abs_x, 0.95)),
    q99_abs = unname(stats::quantile(abs_x, 0.99)),
    max_abs = max(abs_x))}

# median_abs  median absolute error
# mad         residual spread
# rmse        quadratic error
# q95_abs     q95 upper-tail error
# q99_abs     q99 upper-tail error
# max_abs     max observed error



