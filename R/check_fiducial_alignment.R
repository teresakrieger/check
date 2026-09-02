#' Check fiducial frame detection & align tissue image, array, and pixel positions
#'
#'
#' Load and plot the detected tissue image and the tissue_<image_resolution>_image.
#' Visualize the array-grid positions and full-resolution pixel positions from
#' a 10x Genomics Visium tissue-position file and align array-grid (`adjust_array`)
#' and pixel positions (`adjust_pixels`) to match the tissue_<image_resolution>_image. 
#' All plots will be returned as a \code{ggplot2} object, along with the aligned
#' tissue positions as \code{data.frame}.
#'
#' This function is intended as a diagnostic step, ensuring correct fiducial 
#' marker detection, and for creating a \file{csv} file for aligned tissue
#' positions for `spatial_conversion()`.
#'
#' @param visium_dir Character scalar. Base directory containing the Visium
#'   sample directories.
#' @param visium_sample Character scalar. Name or identifier of the Visium
#'   sample directory.
#' @param adjust_array Character vector or \code{NULL}. Optional transformations
#'   applied to array coordinates to match the image tissue orientation. 
#'   Available transformations are \code{"swap"} to exchange the plotted x 
#'   and y coordinates, \code{"x"} to mirror the plotted x axis, and \code{"y"} 
#'   to mirror the plotted y axis. Transformations are applied in the order: 
#'   swap, x, y. Defaults to \code{NULL}.
#' @param adjust_pixels Character vector or \code{NULL}. Optional transformations
#'   applied to pixel coordinates to match the image tissue orientation. 
#'   Available transformations are \code{"swap"}, \code{"x"}, and \code{"y"}, 
#'   with the same definitions as for \code{adjust_array}. Defaults to \code{NULL}.
#' @param tp_file_name Character scalar. Name of the tissue-position file within
#'   the sample's \file{outs/spatial} directory. Defaults to
#'   \code{"tissue_positions.csv"}.
#' @param spot_color Named character vector defining colors for the
#'   \code{in_tissue} values. The vector names must correspond to values in the
#'   \code{in_tissue} column. Defaults to
#'   \code{c("0" = "lightgrey", "1" = "blue")}.
#' @param image_resolution Character scalar. Resolution of the tissue image file
#'   within the sample's \file{outs/spatial} directory. Available options are
#'   \code{"lowres"} or \code{"hires"}. Defaults to \code{"lowres"}.
#' @param detected_tissue_image Character scalar. Name of the detected tissue image
#'   within the sample's \file{outs/spatial} directory. JPEG and PNG files are
#'   supported. Defaults to \code{"detected_tissue_image"}.
#' @param reduce_img_size Logical scalar. Whether to reduce the image
#'   resolution before constructing the plot. This can substantially reduce
#'   memory usage. Defaults to \code{TRUE}.
#' @param use_terra Logical scalar. Whether to use the package \pkg{terra}
#'   for image reduction before constructing the plot. Defaults to \code{TRUE}.
#' @param img_height Positive integer. Target image height in pixels when
#'   \code{reduce_img_size = TRUE}. The width is calculated while preserving
#'   the original aspect ratio. Defaults to 600.
#' @param resize_method Character scalar. Interpolation method passed to
#'   \code{\link{adjust_img_size}} when resizing through \pkg{terra}. Common
#'   choices are \code{"bilinear"} and \code{"near"}.
#'
#' @return A named list of objects containing:
#' \describe{
#'   \item{\code{detected_tissue_image}}{A \code{ggplot} object showing the detected
#'     tissue with the aligned fiducial frame, or \code{NULL} when the image was 
#'     not requested or could not be loaded.}
#'   \item{\code{array_spots}}{A \code{ggplot} object showing the Visium
#'     array-grid positions.}
#'   \item{\code{pixel_positions}}{A \code{ggplot} object showing the spot
#'     coordinates in full-resolution image pixels.}
#'    \item{\code{tp}}{The tissue positions \code{data.frame} with named attribute 
#'    `orientation_adjustments`, including information about adjustments made to the
#'     array positions `adjust_array` and to the pixel positions `adjust_pixels`}
#'  }
#'
#' @details
#' The function expects the standard Space Ranger directory structure:
#'
#' \preformatted{
#' <visium_dir>/<visium_sample>/outs/spatial/
#' }
#'
#' with a csv tissue position file, the detected image and tissue_<image_resolution>_image
#' files as JPEG or PNG image. With JPEG files requiring the suggested package \pkg{jpeg},
#' and PNG images requiring the suggested package \pkg{png}.
#'
#' @seealso
#' \code{\link{adjust_img_size}}
#'
#' @examples
#' \dontrun{
#' check_detection <- check_fiducial_alignment(
#'   visium_dir = "/path/to/visium",
#'   visium_sample = "sample_1"
#' )
#' 
#' # ggplot objects
#' check_detection$tissue_image
#' check_detection$detected_tissue_image
#' check_detection$array_spots
#' check_detection$pixel_positions
#' 
#' # tissue positions with aligned 
#' # array and pixel positions
#' head(check_detection$tp)
#' }
#'
#' @export
check_fiducial_alignment <- function(visium_dir, visium_sample, adjust_array=NULL, adjust_pixels= NULL,
                                     tp_file_name="tissue_positions.csv", spot_color=c("0"="lightgrey", "1"="blue"), 
                                     image_resolution = "lowres", detected_tissue_image="detected_tissue_image", 
                                     reduce_img_size=TRUE, use_terra=TRUE, img_height=600L, resize_method="bilinear"){
  
  image_resolution <- match.arg(image_resolution, c("lowres", "hires"))
  res_info <- ifelse(image_resolution == "lowres", "Low resolution Image", "High resolution Image")
  img_file_name <- paste0("tissue_",image_resolution,"_image")
  
  # check logical arguments
  logical_arguments <- list(
    reduce_img_size = reduce_img_size,
    use_terra = use_terra)
  
  for (argument_name in names(logical_arguments)) {
    argument_value <- logical_arguments[[argument_name]]
    
    if (!is.logical(argument_value) || length(argument_value) != 1L || is.na(argument_value)) {
      stop("'", argument_name, "' must be one non-missing logical value.")}}

  if (!is.null(adjust_array)){  
    adjust_array <- match.arg(adjust_array, choices = c("x", "y", "swap"), several.ok = TRUE)}
  if (!is.null(adjust_pixels)){
    adjust_pixels <- match.arg(adjust_pixels, choices = c("x", "y", "swap"), several.ok = TRUE)}

  # load tp
  tp <- load_tp(visium_dir=visium_dir, visium_sample=visium_sample, tp_file_name=tp_file_name, check_values=TRUE)

  # add aligned spot coordinates
  missing_colors <- setdiff(
    unique(as.character(tp$in_tissue)),
    names(spot_color))
  
  if (length(missing_colors) > 0L) {
    stop("'spot_color' does not define colors for: ",
      paste(missing_colors, collapse = ", "))}
  
  tp$array_x <- tp$array_col
  tp$array_y <- tp$array_row
  tp$pixels_x <- tp$pxl_col_in_fullres
  tp$pixels_y <- tp$pxl_row_in_fullres
  tp$in_tissue <- factor(as.character(tp$in_tissue), levels = names(spot_color))
    if (!is.null(adjust_array)){
    if ("swap" %in% adjust_array){
      new_col <- tp$array_x    
      tp$array_x <- tp$array_y
      tp$array_y <- new_col}
      if ("x" %in% adjust_array){
        tp$array_x <- max(tp$array_x)-tp$array_x}
      if ("y" %in% adjust_array){
        tp$array_y <- max(tp$array_y)-tp$array_y}}
  
      if (!is.null(adjust_pixels)){
        if ("swap" %in% adjust_pixels){
          new_col <- tp$pixels_x   
          tp$pixels_x <- tp$pixels_y
          tp$pixels_y <- new_col}
        if ("x" %in% adjust_pixels){
          tp$pixels_x <- max(tp$pixels_x)-tp$pixels_x}
        if ("y" %in% adjust_pixels){
          tp$pixels_y <- max(tp$pixels_y)-tp$pixels_y}}
  
  # create plots  
  # define base theme
  base_theme <- ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank(), 
                   axis.ticks.y = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = 'white'))
  
  # determine ratio
  array_ratio <- calc_abs_diff(tp=tp, col_name = "array_x")/calc_abs_diff(tp=tp, col_name = "array_y")
  pixel_ratio <- calc_abs_diff(tp=tp, col_name = "pixels_x")/calc_abs_diff(tp=tp, col_name = "pixels_y")
  inv_ratio <- character()
  if(!is.finite(array_ratio) || !array_ratio > 0){
    inv_ratio <- c(inv_ratio, "array")
    array_ratio <- 1}
  if(!is.finite(pixel_ratio) || !pixel_ratio > 0){
    inv_ratio <- c(inv_ratio, "pixel")
    pixel_ratio <- 1}
  
  if (length(inv_ratio)> 0L) {
    warning("Could not determine a valid coordinate ratio for: ",
            paste(inv_ratio, collapse = ", "),". Using 1.")}
  
   # array spot positions
  p3 <- ggplot2::ggplot(tp, ggplot2::aes(x = array_x,
                                         y = array_y,
                                         color = in_tissue))+
    ggplot2::geom_point(size = 1, pch=19)+ ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=spot_color, breaks = names(spot_color)) + 
    ggplot2::coord_fixed(ratio = array_ratio) + 
    ggplot2::labs(title ="Array Spot Positions") +  
    base_theme
  
  # pixel positions
  p4 <- ggplot2::ggplot(tp, ggplot2::aes(x= pixels_x, 
                                         y= pixels_y, 
                                         color = in_tissue))+
    ggplot2::geom_point(size = 1, pch=19)+ ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=spot_color, breaks = names(spot_color)) + 
    ggplot2::coord_fixed(ratio = pixel_ratio) +
    ggplot2::labs(title = "Pixel Positions") +
    base_theme
  
  # Load image
  file_dir <- file.path(visium_dir, visium_sample, "outs", "spatial")
  available_files <- list.files(file_dir)
  img_file <- file.path(file_dir, grep(img_file_name, available_files, value = TRUE))
  if (length(img_file) != 1L){
    stop("Expected exactly one image file.")}
  detected_tissue_file <- file.path(file_dir, grep(detected_tissue_image, available_files, value = TRUE))
  if (length(detected_tissue_file) != 1L){
    stop("Expected exactly one detected tissue file.")}
      p1 <- tryCatch({
        p1 <- load_img(img_file,  reduce_img_size=reduce_img_size, img_height=img_height, use_terra = use_terra, resize_method = resize_method)
        p1 <- p1 + ggplot2::labs(title = res_info) + base_theme}, 
           error = function(e){
             warning("Could not load ",res_info,": ", conditionMessage(e), ".")})

      p2 <- tryCatch({
        p2 <- load_img(detected_tissue_file,  reduce_img_size=reduce_img_size, img_height=img_height, use_terra = use_terra, resize_method = resize_method)
        p2 <- p2 + ggplot2::labs(title = "Detected Tissue") + base_theme},
           error = function(e){
             warning("Could not load Detected Tissue: ", conditionMessage(e), ".")})

  
  attr(tp, "orientation_adjustments") <- list(
    adjust_array = adjust_array,
    adjust_pixels = adjust_pixels)
  
    return(list("tissue_image"=p1,
                "detected_tissue_image"=p2,
                "array_spots"=p3,
                "pixel_positions"=p4,
                "tp"=tp))}

