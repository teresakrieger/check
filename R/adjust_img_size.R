#' Resize an image array
#'
#' Reduce or resize a two-dimensional image matrix or a three-dimensional
#' image array. When only one target dimension is supplied, the other is
#' calculated from the original image aspect ratio.
#'
#' Resizing can be performed using spatial interpolation implemented by
#' \pkg{terra}, or by selecting approximately evenly spaced source pixels.
#' The latter provides a lightweight nearest-neighbour-style fallback and is
#' intended for downsampling only.
#'
#' @param img Numeric matrix or array containing an image. The first two
#'   dimensions must correspond to image height and width, respectively.
#'   A third dimension may contain image channels.
#' @param h Integer scalar giving the requested output height. May be
#'   \code{NULL} when \code{w} is supplied. Defaults to \code{600L}.
#' @param w Integer scalar giving the requested output width. May be
#'   \code{NULL} when \code{h} is supplied. Defaults to \code{NULL}.
#' @param use_terra Logical scalar. Whether to resize the image using
#'   \code{\link[terra]{resample}}. Defaults to \code{TRUE}.
#' @param method Character scalar specifying the interpolation method passed
#'   to \code{\link[terra]{resample}}. Common choices include
#'   \code{"bilinear"} and \code{"near"}. When \code{NULL}, evenly spaced
#'   source pixels are selected instead. Defaults to \code{"bilinear"}.
#' @param tol Numeric scalar giving the tolerated relative difference between
#'   the original and requested aspect ratios when both \code{h} and
#'   \code{w} are explicitly supplied. Defaults to \code{0.05}.
#'
#' @return A matrix or array containing the resized image. Matrix input
#'   remains two-dimensional; array input retains its channel dimension.
#'
#' @details
#' When only one of \code{h} or \code{w} is supplied, the missing dimension
#' is calculated while preserving the original image aspect ratio, subject
#' to integer rounding.
#'
#' The non-\pkg{terra} fallback selects approximately evenly spaced source
#' pixels and therefore does not interpolate new values. It supports
#' downsampling but not upsampling.
#'
#' @examples
#' \dontrun{
#' img_small <- adjust_img_size(
#'   img = img,
#'   h = 600L
#' )
#'
#' img_small_nearest <- adjust_img_size(
#'   img = img,
#'   h = 600L,
#'   use_terra = FALSE
#' )
#' }
#'
#' @export
adjust_img_size <- function(img, h = 600L, w = NULL, use_terra=TRUE, method = "bilinear", tol = 0.05) {
  
  # sanity checks
  if (is.null(h) && is.null(w)) {
    stop("'h' or 'w' must be specified.")}
  validate_dimension(w, "w")
  validate_dimension(h, "h")
  h <- if (is.null(h)) NULL else as.integer(h)
  w <- if (is.null(w)) NULL else as.integer(w)
  
  img_dim <- dim(img)
  if (is.null(img_dim) || !length(img_dim) %in% c(2L, 3L) || !is.numeric(img)) {
    stop("'img' must be a numeric matrix or three-dimensional array.")}
  
  if (!is.logical(use_terra) || length(use_terra) != 1L || is.na(use_terra)) {
    stop("'use_terra' must be TRUE or FALSE.")}
  
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol < 0) {
    stop("'tol' must be a finite, non-negative numeric value.")}
  
  if (!is.null(method) &&
    (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% c("bilinear", "near"))) {
    stop("'method' must be NULL, 'bilinear', or 'near'.")}
  
  orig_h <- img_dim[1]
  orig_w <- img_dim[2]
  
  expected_ratio <- orig_w / orig_h
  if (!is.null(h) && !is.null(w)){
    requested_ratio <- w / h
    if (abs(requested_ratio / expected_ratio - 1) > tol) {
      warning("Both 'h' and 'w' were specified and alter the original image ratio > ", tol)}}
    
  if (is.null(h)) {
    h <- as.integer(round(orig_h * w / orig_w))
  } else if (is.null(w)) {
    w <- as.integer(round(orig_w * h / orig_h))}
  
  if (h < 1L || w < 1L) {
    stop("'h' and 'w' must both be at least 1.")}
  
  # no changes required
  if (h == orig_h && w == orig_w) {
    return(img)}
  
  # resize using terra
  if (use_terra && !is.null(method)){
    r <- terra::rast(img)
    
    template <- terra::rast(
      nrows = h,
      ncols = w,
      nlyrs = terra::nlyr(r),
      extent = terra::ext(r))
    
    r_small <- terra::resample(
      r, 
      template,
      method = method)
    
    resized <- terra::as.array(r_small)
    
    if (length(img_dim) == 2L) {
      resized <- resized[, , 1L, drop = TRUE]}
    return(resized)
    
  }else{
    
    # resize using source-pixel selection
    if (h > orig_h || w > orig_w) {
      stop("The source-pixel fallback cannot upscale images. Use ",
        "'use_terra = TRUE' with a non-NULL interpolation method.")}
    
    # Select approximately evenly spaced source pixels
    row_idx <- unique(round(seq.int(1L, orig_h, length.out = h)))
    col_idx <- unique(round(seq.int(1L, orig_w, length.out = w)))
    
    if (length(row_idx) != h || length(col_idx) != w){
      warning("Transformed image pixels don't match the original dimensions.")}
    
    if (length(dim(img)) == 2L) {
      img_small <- img[row_idx, col_idx, drop = FALSE]
    } else {
      img_small <- img[row_idx, col_idx, , drop = FALSE]}
    
    return(img_small)}}

