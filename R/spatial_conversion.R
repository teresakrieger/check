#' Convert Visium array and pixel coordinates to physical coordinates
#'
#' Reconstruct physical coordinates for Visium array positions and transform
#' full-resolution image-pixel coordinates into the same coordinate system.
#'
#' @param aligned_tp Data frame containing all Space Ranger tissue-positions 
#' with \code{array_row}, \code{array_col}, \code{pxl_row_in_fullres}, and 
#'  \code{pxl_col_in_fullres} columns, and aligned array (\code{array_x},
#'   \code{array_y}) and pixel positions (\code{pixels_x}, \code{pixels_y}),
#'   as returned by `check_fiducial_alignment()`.
#' @param spex SpatialFeature Object as returned by `load_spex()`.
#' @param L_capture_area Positive numeric scalar defining the side lengths
#'  of the squared capture-area. Supported values are \code{6500} and 
#'  \code{11000}.
#' @param spatial_geometry Character Scalar. Defaults to \code{"hex100"},
#' representing the best supported spatial geometry settings with a nearest
#' neighbor spot distance of 100µm in a hexagonal grid. 
#' Given the limited space of the capture area, exceeded by the default spatial
#' geometry grid `hex100`, we also support \code{k0}, corresponding to zero 
#' distance between the spot's edges and the capture area edges, while 
#' ensuring the containment of all spots.
#' @param image_resolution Character scalar. Image resolution to load, either
#'   \code{"lowres"} or \code{"hires"}.
#'
#' @return A \code{SpatialFeatureExperiment} containing:
#' \describe{
#'   \item{spatial array coordinates}{Spot-center coordinates transformed to the
#'     requested spatial coordinate system.}
#'   \item{\code{spatial_array}}{Buffered \code{sf} spot polygons in the spatial
#'     array coordinate system.}
#'   \item{spatial pixel coordinates}{Pixel coordinates transformed to the
#'     requested spatial coordinate system.}
#'   \item{\code{spatial_pixels}}{Buffered \code{sf} Pixel polygons in the 
#'     spatial pixel coordinate system.}
#'   \item{transformed image}{The selected image with an extent aligned to the
#'     spatial coordinates.}
#'   \item{metadata}{Image and spot bounding boxes, scaling factors, the
#'     original upper image boundary, and the spatial-unit label.}
#' }
#' 
#' @seealso
#' \code{\link[SpatialFeatureExperiment]{toSpatialFeatureExperiment}}
#'   
#' @examples
#' \dontrun{
#' # convert spots spatially aligned tissue position & image data
#' for 6.5mm array (L_capture_area = 6500) into spatially 
#' aligned sfe object   
#' sfe_6500 <- spatial_conversion(
#'   aligned_tp = tp_aligned,
#'   spex = spex,
#'   L_capture_area = 6500,
#'   image_resolution="hires"
#' )
#' }
#'
#'
#' @export
spatial_conversion <- function(aligned_tp, spex, L_capture_area,
                               spatial_geometry = "hex100",
                               image_resolution=c("lowres", "hires")){
  
  # check input 
  tp <- aligned_tp
  spot_diameter <- 55
  spatial_unit <- "µm"
  
  spatial_geometry <- match.arg(spatial_geometry, choices = c("hex100", "k0"))
  
  # -------------------------------------------------------------------------#
  
  # tp
  # check availability of aligned array & pixel positions
  array_aligned <- c("array_x", "array_y")
  pixels_aligned <- c("pixels_x", "pixels_y")
  
  if (!all(array_aligned %in% colnames(tp)) && !!all(pixels_aligned %in% colnames(tp))) {
    stop("Observed 'aligned_tp' column names: ",paste(colnames(tp), collapse = ", "), 
         "\nMissing aligned array ('array_x', 'array_y')a nd pixel coordinates ('pixels_x', 'pixels_y').")}
  
  # get attribute information
  attr_info <- attr(tp, "orientation_adjustments")
  
  # check tp
  tp <- check_tp(tp = tp, L_capture_area = L_capture_area, 
                 array_columns= c("array_row", "array_col",array_aligned), 
                 pixel_columns = c("pxl_row_in_fullres", "pxl_col_in_fullres", pixels_aligned))
  
  # no attribute information -> derive swap information from ranges
  if (is.null(attr_info) || !all(c("adjust_array", "adjust_pixels") %in% names(attr_info))){
    stop("Missing expected attribute information 'orientation_adjustments' with
            \n'adjust_array' and 'adjust_pixels' slots in tp object.")
  }else{
    # attribute information -> swap information
    swapped_array <- "swap" %in%  attr_info$adjust_array
    # swapped_pixels <- "swap" %in%  attr_info$adjust_pixels
    }
  
  # -------------------------------------------------------------------------#
  
  # image resolution
  image_resolution <- match.arg(image_resolution)
  
  # array settings
  spatial_unit <- match.arg(spatial_unit, choices=c("µm", "mm"))
  L_capture_area <- check_numeric(x= L_capture_area, name="L_capture_area", n=1L, value_range="positive", value_type="numeric") 
  spot_diameter <- check_numeric(x= spot_diameter, name="spot_diameter", n=1L, value_range="positive", value_type="numeric") 

  if(!L_capture_area %in% c(6500, 11000) || spot_diameter != 55){
    stop("`spatial_conversion()` is only implemented for a capture area of 6500 or 11000µm and a spot_diameter of 55µm.")}

  if (L_capture_area == 6500){
    rows <- 0:77
    columns <- 0:127}
  if (L_capture_area == 11000){
    rows <- 0:127
    columns <- 0:223}
  
  n_rows <- length(rows)
  n_cols <- length(columns)
  
  # spot radius
  r <- spot_diameter/2
  
  # incorporate spatial information
  if (spatial_geometry == "hex100"){
    d_col <- 50  
    d_row <- 50 * sqrt(3)
    m_row <- L/2 - n_cols*25-2.5
    m_col <- L/2 - (n_rows-1L)*(25*sqrt(3))-r}
  if (spatial_geometry == "k0"){
    d_col <- (L-spot_distance)/(n_col-1)
    d_row <- (L-spot_distance)/(n_row-1)
    m_row <- r
    m_col <- r}
  
  # check containment
  if (m_row < r || m_col < r) {
    message(
      "The selected spatial reconstruction does not keep complete spots ",
      "inside the capture area. Derived margins:
      \nrow = ",round(m_row, 5), ",
      \ncolumn = ", round(m_col, 5), "; 
      \nrequired minimum = ", r, ".")}
  
  # create new df
  df_corners <- data.frame(matrix(nrow=2, ncol=8), row.names = c("spatial_array", "spatial_pixels"))
  colnames(df_corners) <- c("xmin_spots",  "xmax_spots","ymin_spots", "ymax_spots", "xmin_img",  "xmax_img","ymin_img", "ymax_img")
  
  #---------------------------------------------------------------------------#
  
  # match spex & tp
  idx <- match(colnames(spex), tp$barcode)
  
  if (anyNA(idx)) {
    stop("Could not match ", sum(is.na(idx)),
      " SpatialExperiment barcode(s) to 'tp'.")}
  
  tp <- tp[idx, , drop = FALSE]
  
  stopifnot(identical(tp$barcode, colnames(spex)))
  
  rownames(tp) <- tp$barcode

  #---------------------------------------------------------------------------#
  
  # spex
  # get scaling factors
  sf <- SpatialExperiment::scaleFactors(spex)
  # convert to sfe
  sfe <- SpatialFeatureExperiment::toSpatialFeatureExperiment(spex, spotDiameter = spot_diameter, unit =  spatial_unit)
  
  # get image
  # Extract corners for the bounding box of the spatial image
  img <- SpatialExperiment::getImg(sfe, image_id = image_resolution)@image
  
  if (!inherits(img, "SpatRaster")) {
    stop("Expected a terra 'SpatRaster' image. Found ", class(img),".")}
  
  old_ext <- terra::ext(img)
  if (!inherits(old_ext, "SpatExtent")) {
    stop("Expected a terra 'SpatExtent' object for image_corners. Found ", class(old_ext),".")}
  
  image_corners <- list(xmin = terra::xmin(old_ext),
                        xmax = terra::xmax(old_ext),
                        ymin = terra::ymin(old_ext),
                        ymax = terra::ymax(old_ext))
  
  for (corner in names(image_corners)){
      image_corners[corner] <- check_numeric(x=image_corners[corner],name=corner, n=1L, value_range="any", value_type="numeric")  }

  #---------------------------------------------------------------------------#

  # A) array transformation
  if (!swapped_array){
    dx <- d_col
    dy <- d_row
    mx <- m_col
    my <- m_row
    nx <- n_cols
    ny <- n_rows
  }else{
    dx <- d_row
    dy <- d_col
    mx <- m_row
    my <- m_col
    nx <- n_rows
    ny <- n_cols}
  
  tp$spatial_x <- (tp$array_x*dx)+mx
  tp$spatial_y <- (tp$array_y*dy)+my
  
  # ensure valid ranges
  inv_range <- NULL
  if (diff(range(tp$spatial_x)) <= 0 || diff(range(tp$spatial_y)) <= 0){
    inv_range <- c(inv_range, "spatial")}
  
  if (diff(range(tp$pixels_x)) <= 0 || diff(range(tp$pixels_y)) <= 0) {
    inv_range <- c(inv_range, "pixel")}
  
  if (!is.null(inv_range)) {
    stop("The ",paste(inv_range, collapse = " & ")," coordinates must span more than one position on both axes.")}
  
  # add image positions for fullres_image
  tp$fullres_image_x <- tp$pxl_col_in_fullres
  tp$fullres_image_y <- image_corners$ymax - tp$pxl_row_in_fullres
  
  # pixel transformation
  # get scaling factors
  sf_x <- diff(range(tp$spatial_x))/diff(range(tp$fullres_image_x))
  sf_y <- diff(range(tp$spatial_y))/diff(range(tp$fullres_image_y))
  
  # get offsets for positions
  k_x <- min(tp$spatial_x) - min(tp$fullres_image_x) * sf_x
  k_y <- min(tp$spatial_y) - min(tp$fullres_image_y) * sf_y
  
  # Transform aligned image-pixel coordinates into the
  # reconstructed physical coordinate system.
  tp$spatial_pxl_x <- k_x + tp$fullres_image_x * sf_x
  tp$spatial_pxl_y <- k_y + tp$fullres_image_y * sf_y
  
  # add delta
  tp$delta_spatial_x <- tp$spatial_x - tp$spatial_pxl_x
  tp$delta_spatial_y <- tp$spatial_y - tp$spatial_pxl_y
  
  #---------------------------------------------------------------------------#
  
  # adjust raster: spatial transformation of image used to 
  # spatially aligned image with adjusted scaling factors & offsets
  new_ext <- terra::ext(
    k_x + terra::xmin(old_ext)*sf_x,
    k_x + terra::xmax(old_ext)*sf_x,
    k_y + terra::ymin(old_ext)*sf_y,
    k_y + terra::ymax(old_ext)*sf_y)
  
  terra::ext(img) <- new_ext
  
  ##############################################
  
  # consistency
  stopifnot(all(colnames(sfe) %in% rownames(tp)))
  
  # add geometries for spatial coordinates
  spatial_coords <- list(
    "spatial_array" = c("spatial_x", "spatial_y"),
    "spatial_pixels" = c("spatial_pxl_x", "spatial_pxl_y"))
  
  for (coord_id in names(spatial_coords)){
    coord_cols <- spatial_coords[[coord_id]]
    coord_df <- tp[ , coord_cols, drop = FALSE]
    spots_spatial <- sf::st_as_sf(coord_df, coords = coord_cols,
                                  agr = "identity", remove = FALSE)
    # add spot diameter & barcode info
    spots_spatial <- sf::st_buffer(spots_spatial, dist =r)
    rownames(spots_spatial) <- tp$barcode
    
    # add spatial spot ranges to df
    spot_box <- sf::st_bbox(spots_spatial)
    names(spot_box) <- paste0(names(spot_box), "_spots")
    df_corners[coord_id,]  <- spot_box[match(colnames(df_corners), names(spot_box))]
    # add image raster range info to df
    df_corners <- fill_df(df=df_corners,row_index = which(rownames(df_corners)== coord_id), img = img)
    
    # spatially resolved array positions
    sf_matrix <- as.matrix(tp[ , coord_cols])
    SpatialExperiment::spatialCoords(sfe) <- sf_matrix
    SpatialFeatureExperiment::colGeometry(sfe, coord_id) <- spots_spatial}
  
  # add / adjust metadata
  S4Vectors::metadata(sfe)$aligned_tp <- tp
  S4Vectors::metadata(sfe)$corner_information <- df_corners
  S4Vectors::metadata(sfe)$spatial_resolution <- list(
    spatial_unit = spatial_unit,
    L_capture_area = L_capture_area,
    spot_diameter = spot_diameter,
    n_rows = n_rows,
    n_columns = n_cols,
    n_x = nx,
    n_y = ny,
    d_x = dx,
    d_y = dy,
    d_diagonal = sqrt(dx^2 + dy^2),
    margin_center_x = mx,
    margin_center_y = my,
    margin_edge_x = mx-r,
    margin_edge_y = my-r,
    scaling_factor_x = sf_x,
    scaling_factor_y = sf_y,
    offset_x = k_x,
    offset_y = k_y)
  
  return(sfe)}
