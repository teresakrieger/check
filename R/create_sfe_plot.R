#' Plot a spatially aligned Visium image and spot geometry
#'
#' Create a \code{ggplot2} overlay of a spatially transformed Visium image and
#' an \code{sf} spot geometry stored in a
#' \code{SpatialFeatureExperiment}. Spots can be displayed using one fixed
#' outline color or colored according to a geometry metadata column.
#'
#' @param sfe A \code{SpatialFeatureExperiment} containing the requested image
#'   and column geometry.
#' @param image_resolution Character scalar. Image identifier passed to
#'   \code{\link[SpatialExperiment]{getImg}}, for example \code{"lowres"} or
#'   \code{"hires"}.
#' @param geom_id Character scaler. Name of the spatially resolved column
#'   geometry retrieved with \code{\link[SpatialFeatureExperiment]{colGeometry}}.
#'   Available options include \code{"spatial_array"} and \code{"spatial_pixels"}.
#' @param plot_title Character scalar or \code{NULL}. Plot title to be used.
#' @param color_by Character scalar or \code{NULL}. Name of a metadata column
#'   or the selected geometry field used to color spot outlines. If \code{NULL}, 
#'   all spots are drawn using one fixed color.
#' @param colors_use Character vector or \code{NULL}. Colors used for spot
#'   outlines. When \code{color_by = NULL}, this must contain one color. When
#'   \code{color_by} is supplied, it may be either an unnamed vector with one
#'   color per observed value or a named vector whose names match the observed
#'   values. If insufficient colors are supplied, a default qualitative palette
#'   is generated.
#'  @param metadata_slot Optional, Character Scalar. Name of the metadata slot
#'   to look for the `color_by` column, if `metadata` is selected as `color_source`.
#' @param colors_source Character vector or \code{NULL} (default). When 
#'  \code{color_by} is not `NULL`, this must contain the slot in the sfe object, 
#'  where the color information can be retrieved. Available options include `geometry`
#'  and `metadata`.
#'
#' @return A \code{ggplot} object for the selected \code{geom_id}
#'  containing the transformed image raster and the specific spatial spot geometry.
#'
#' Missing values in \code{color_by} are displayed transparently.
#'
#' @seealso
#' \code{\link{load_sfe}},
#' \code{\link{raster_to_ggdf}},
#'
#' @examples
#' \dontrun{
#' p <- create_sfe_plot(
#'   sfe = sfe,
#'   image_resolution = "lowres"
#' )
#'
#' p_selected <- create_sfe_plot(
#'   sfe = sfe,
#'   image_resolution = "lowres",
#'   color_by = "selected",
#'   colors_use = c(
#'     selected = "cornflowerblue",
#'     out = "lightgrey"
#'   ),
#'   
#'   plot_title = "Selected spots"
#' )
#' }
#'
#' @export
create_sfe_plot <- function(sfe, image_resolution, 
                            geom_id="spatial_pixels", 
                            plot_title="Spatially aligned Pixel Positions", 
                            color_by=NULL, colors_use =NULL, 
                            metadata_slot = "spot_metadata",
                            color_source = NULL){
  if (!inherits(sfe, "SpatialFeatureExperiment")){
    stop("Expected sfe to be a SpatialFeatureExperiment. Got ", class(sfe), ".")}
  
  spot_diameter <- 55
  spatial_unit <- "µm"
  
  # check color by
  if (!is.null(color_by)){
    if (is.null(color_source)){
      stop("Please specify a color_source for color_by.")}
  color_source <- match.arg(color_source, choices = c("geometry", "metadata"))}

  # check geometries
  if (is.null(geom_id) || length(geom_id) != 1L || !geom_id %in% names(SpatialFeatureExperiment::colGeometries(sfe))){
    stop("Please specify a geom_id available in colGeometries(sfe).")}
  
  geom_id<- match.arg(geom_id, choices =  c("spatial_array", "spatial_pixels"))
  tp_spatial <- SpatialFeatureExperiment::colGeometry(sfe, geom_id)  
    
  img <- SpatialExperiment::getImg(sfe, image_id = image_resolution)@image
    if (!inherits(img, "SpatRaster")) {
      stop("Expected a terra 'SpatRaster' image. Found ", class(img),".")}
    if (terra::nlyr(img) != 3L) {
      stop("Expected exactly 3 RGB layers. Found ", terra::nlyr(img)," layers.")}
    
    # image
    df_img <- raster_to_ggdf(img)

    if (!is.null(spatial_unit)){
      su <- paste0(" (",spatial_unit, ")")
    }else{
      su <- ""}
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = df_img, ggplot2::aes(x = x, y = y, fill = hex)) +
      ggplot2::scale_fill_identity()
    
    # one color for all spots
    if (is.null(color_by) ){
      if(is.null(colors_use) || length(colors_use) != 1L){
        colors_use <- "cornflowerblue"}
      p <- p +
        ggplot2::geom_sf(data = tp_spatial, fill = NA, color= unlist(unname(colors_use)), linewidth = 0.25, show.legend = "point")}
    
    # color by used
      if (!is.null(color_by)){
        # color_source="geomtry"
        if (color_source =="geometry"){
          if (!color_by %in% colnames(tp_spatial)){
            stop("Can't find ", color_by, " in colGeometry ", geom_id)}}
        
        # color_source="metadata"
        if (color_source == "metadata"){
          sfe_meta <- S4Vectors::metadata(sfe)
          
          metadata_available <- !is.null(metadata_slot) &&
            metadata_slot %in% names(sfe_meta) &&
            is.data.frame(sfe_meta[[metadata_slot]])
          
          meta_has_color <- metadata_available &&
            color_by %in% colnames(sfe_meta[[metadata_slot]])
          
          if (!metadata_available){
            stop("Can't find ", metadata_slot, " in the metadata of the sfe object.")}
          if (!meta_has_color){
            stop("Can't find ", color_by, " in ",metadata_slot, " metadata of the sfe object.")}
          
          sfe_meta <- sfe_meta[[metadata_slot]]

          idx <- match(rownames(tp_spatial), rownames(sfe_meta))
          
          if (anyNA(idx)) {
           warning("Could not match all spots from colGeometry `",
                 geom_id, "` to metadata slot `", metadata_slot,"`.")}
          tp_spatial[[color_by]] <- sfe_meta[[color_by]]}
        
      color_values <- tp_spatial[[color_by]]
      color_levels <- unique(as.character(color_values[!is.na(color_values)]))
      n_colors <- length(color_levels)
      if (n_colors == 0L){
        stop("Column ", color_by, " doesn't contain valid entries.")}
      
      if (is.null(colors_use) || length(colors_use) < n_colors){
        warning("Colors specified doesn't match number of unique entries in ",color_by,".
              \nUsing default colors")
        colors_use <- grDevices::hcl.colors(
          n = n_colors,
          palette = "Pastel 1")
        names(colors_use) <- color_levels
      }else{
        if (is.null(names(colors_use))) {
          colors_use <- colors_use[seq_len(n_colors)]
          names(colors_use) <- color_levels
        }else{
          colors_use <- colors_use[color_levels]}
        # Check that named colors cover every observed value & aren't NA
        missing_colors <- setdiff(color_levels, names(colors_use))
        if (length(missing_colors) > 0L || anyNA(colors_use)) {
          stop("No valid color was supplied for: ", 
               paste(unique(c(missing_colors, names(colors_use)[is.na(colors_use)])), collapse = ", "))}}
      
      # no coord_equal() or related parameter are required, since geom_sf() handles ratios based on values
      p <- p +
        ggplot2::geom_sf(data = tp_spatial, mapping = ggplot2::aes(color= .data[[color_by]]), fill = NA,  linewidth = 0.25, show.legend = "point") +
        ggplot2::scale_color_manual(values = colors_use, breaks = color_levels, na.value = "transparent", drop=FALSE, name=color_by)}
      
    p <- p + ggplot2::coord_sf(expand = FALSE) +
      ggplot2::labs(x = paste0("Spatial X", su), y = paste0("Spatial Y", su),
                    title = plot_title) +
      ggplot2::theme_minimal(base_size = 12)+
      ggplot2::theme(axis.title.y=ggplot2::element_text(margin = ggplot2::margin(r=5)), 
                     axis.title.x=ggplot2::element_text(margin = ggplot2::margin(t=2)))+
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape=1, linewidth=1, stroke=1, size=2)))
  
  return(p)}


