#' Prepare an interactive Visium spot-selection application
#'
#' Construct the user-interface and server components for a Shiny application
#' that displays spatial Visium spots over a tissue image and supports manual
#' lasso or box selection through \pkg{plotly}.
#'
#' Selected spots are reported by barcode and spatial centroid coordinates and
#' can be downloaded as a CSV file.
#'
#' @param sfe A \code{SpatialFeatureExperiment} containing a spatial image and
#'   an \code{sf} column geometry.
#' @param image_resolution Character scalar. Identifier of the image to display, 
#' such as \code{"lowres"} or \code{"hires"}.
#' @param geom_id Character scalar. Name of the column geometry containing the
#'   spatially resolved spot positions. Defaults to \code{"spatial_pixels"}.
#' @param pt_color Character scalar. Marker color used for selectable spots.
#'   Defaults to \code{"blue"}.
#' @param pt_size Positive numeric scalar. Plotly marker size. Defaults to 5.
#' @param pt_opacity Numeric scalar between 0 and 1. Marker opacity. Defaults
#'   to 0.5.
#' @param plot_height Character scalar defining the height of the plotly output,
#'   for example \code{"800px"}. Defaults to \code{"800px"}.
#' @param download_filename Character scalar. Default filename used when
#'   downloading selected spots. Defaults to
#'   \code{"manual_selected_visium_spots.csv"}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{spots_df}}{A data frame containing spot barcodes and centroid
#'     coordinates.}
#'   \item{\code{img_ext}}{The spatial extent of the displayed image.}
#'   \item{\code{user_interface}}{A Shiny user-interface definition.}
#'   \item{\code{server_settings}}{A Shiny server function.}
#' }
#'
#' The returned components can be launched with:
#'
#' \preformatted{
#' app <- prepare_interactive_spot_selection(sfe, "lowres")
#' shiny::shinyApp(app$user_interface, app$server_settings)
#' }
#'
#' @details
#' The image is embedded as a static Plotly background, while only the spot
#' centroids are represented as interactive markers. This avoids converting
#' the complete image raster or spot polygons into interactive Plotly objects.
#'
#' The packages \pkg{png} and \pkg{base64enc} are required at run time.
#'
#' @seealso
#' \code{\link{create_sfe_plot}}
#'
#' @examples
#' \dontrun{
#' app <- prepare_interactive_spot_selection(
#'   sfe = sfe,
#'   image_resolution = "lowres"
#' )
#'
#' shiny::shinyApp(
#'   ui = app$user_interface,
#'   server = app$server_settings
#' )
#' }
#'
#' @export
prepare_interactive_spot_selection <- function(sfe, image_resolution, geom_id="spatial_pixels", barcode_col="barcode", 
                                               pt_color="blue",  # selected color, else transparent
                                               pt_size=5, pt_opacity=0.5, plot_height="800px",
                                               download_filename = "manual_selected_visium_spots.csv"){
  
  if (!inherits(sfe, "SpatialFeatureExperiment")){
    stop("Expected sfe to be a SpatialFeatureExperiment. Got ", class(sfe), ".")}
  
  # Convert sf polygons to centroid coordinates for selection
  spots_spatial <- SpatialFeatureExperiment::colGeometry(sfe, geom_id)
  if (is.null(spots_spatial) || !inherits(spots_spatial, "sf") || nrow(spots_spatial) == 0L) {
    stop(
      "Geometry '", geom_id, "' must be a valid 'sf' object, available in sfe.")}
   spots_df <- spots_spatial %>%
    sf::st_centroid() %>%
    dplyr::mutate(
      barcode = rownames(spots_spatial),
      x = sf::st_coordinates(geometry)[, 1],
      y = sf::st_coordinates(geometry)[, 2])  %>%
    sf::st_drop_geometry()
  
  # transform SpatRaster image into normalized array
  img <- SpatialExperiment::getImg(sfe, image_id = image_resolution)@image
  if (!inherits(img, "SpatRaster")) {
    stop("Expected a terra 'SpatRaster' image. Found ", class(img),".")}
  
  img_arr <- terra::as.array(img)
  if (length(dim(img_arr)) != 3L) {
    stop( "The image could not be converted to a three-dimensional array.")}
  if (!is.finite(max(img_arr, na.rm = TRUE)) ||max(img_arr, na.rm = TRUE)<=0 || !is.finite(max(img_arr, na.rm = TRUE))){
    stop("No valid image values found.")}
  if (max(img_arr, na.rm = TRUE) > 1){
      img_arr <- img_arr / max(img_arr, na.rm = TRUE)}
  img_arr[img_arr < 0] <- 0
  img_arr[img_arr > 1] <- 1
  
  # save temporary image file -> will be saved to working directory
  png_file <- tempfile(fileext = ".png")
  png::writePNG(img_arr, png_file)
  encoded_image <- base64enc::dataURI(file = png_file, mime = "image/png")
  unlink(png_file)
  
  # return image corner information
  img_ext <- terra::ext(img)
  
  ########################################################
  
  # create user interface
  ui <- shiny::fluidPage(
    shiny::titlePanel("Manual Visium spot selection"),
    shiny::fileInput("upload_barcodes", "Upload selected barcodes", accept = c(".csv", "text/csv")),
    plotly::plotlyOutput("spot_plot", height = plot_height),
    shiny::downloadButton("download_selected", "Download selected spots"),
    shiny::tableOutput("selected_table"))
  
  # create plot & server settings
  server <- function(input, output, session) {
    
    # Persistent selection state.
    # The app starts without pre-selected spots.
    selected_barcodes <- shiny::reactiveVal(character(0))
    
    # ----------------------------------------------------
    # Load optional pre-selected barcodes
    # ----------------------------------------------------
    
    shiny::observeEvent(input$upload_barcodes, {
      
      uploaded <- tryCatch(
        utils::read.csv(
          input$upload_barcodes$datapath,
          stringsAsFactors = FALSE,
          check.names = FALSE),
        error = function(e) {
          shiny::showNotification(
            paste("Could not read the barcode file:",
              conditionMessage(e)),
            type = "error")
          return(NULL)})
      
      if (is.null(uploaded)) {
        return()}
      
      if (!barcode_col %in% colnames(uploaded)) {
        shiny::showNotification(paste0("The uploaded CSV must contain a column named '", barcode_col, "'."), type = "error")
        return()
      }else{
          colnames(uploaded)[colnames(uploaded)==barcode_col] <- "barcode"}
      
      uploaded_barcodes <- unique(trimws(as.character(uploaded$barcode)))
      
      uploaded_barcodes <- uploaded_barcodes[!is.na(uploaded_barcodes) & nzchar(uploaded_barcodes)]
      matched_barcodes <- intersect(uploaded_barcodes, spots_df$barcode)
      selected_barcodes(matched_barcodes)
      
      unmatched_barcodes <- setdiff(uploaded_barcodes, spots_df$barcode)      
      if (length(unmatched_barcodes) > 0L) {
        unmatched_info <- paste0(length(unmatched_barcodes),
            " pre-selected barcodes couldn't be matched and were ignored. \n")
      }else{
        unmatched_info <- "All barcodes could be matched: "}
      
      shiny::showNotification(
        paste0(unmatched_info, "Added ", length(matched_barcodes),
          " pre-selected spots."),
        type = "message")})
    
    # ----------------------------------------------------
    # Interactive spot selection
    # ----------------------------------------------------
    
    # add pre-selected spots
    plot_data <- shiny::reactive({
      spots_df %>% dplyr::mutate(
          selected = barcode %in% selected_barcodes(),
          plot_color = ifelse(selected, pt_color, "transparent"),
          plot_opacity = ifelse(selected, 1, pt_opacity))})
    
    # interactive plotly plot
    output$spot_plot <- plotly::renderPlotly({
      current_data <- plot_data()
      plotly::plot_ly(data = current_data, 
                      x = ~x, 
                      y = ~y,
                      key = ~barcode,
                      customdata= ~barcode,
                      text = ~paste0(
                        "barcode: ", barcode,
                        "<br>x: ", round(x, 2),
                        "<br>y: ", round(y, 2),
                        "<br>selected: ", selected),
                      hoverinfo = "text",
                      type = "scatter",
                      mode = "markers",
                      marker = list(
                        size = pt_size,
                        color = current_data$plot_color,
                        opacity = current_data$plot_opacity),
                      source = "visium_spots") %>%
        plotly::layout(
          images = list(
            list(
              source = encoded_image,
              xref = "x",
              yref = "y",
              x = terra::xmin(img_ext),
              y = terra::ymax(img_ext),
              sizex = terra::xmax(img_ext) - terra::xmin(img_ext),
              sizey = terra::ymax(img_ext) - terra::ymin(img_ext),
              sizing = "stretch",
              opacity = 1,
              layer = "below")),
          xaxis = list(
            range = c(terra::xmin(img_ext), terra::xmax(img_ext)),
            title = "Spatial X"),
          yaxis = list(
            range = c(terra::ymin(img_ext), terra::ymax(img_ext)),
            title = "Spatial Y",
            scaleanchor = "x"),
          dragmode = "lasso") %>%
        
        plotly::event_register("plotly_selected")})
    
    
    spot_selection <- shiny::reactive({
      
      plotly::event_data(
        "plotly_selected",
        source = "visium_spots",
        priority = "event")})
    
    shiny::observeEvent(spot_selection(),{
        sel <- spot_selection()
        
        if (is.null(sel) || 
            nrow(sel) == 0L || 
            !"customdata" %in% colnames(sel) || 
            all(is.na(sel$customdata))){
          return()}
        
        sel_barcodes <- unique(as.character(sel$customdata))
        
        current_barcodes <- selected_barcodes()
        
        selected_barcodes(
          union(
            # newly selected
            setdiff(current_barcodes, sel_barcodes),
            # removed
            setdiff(sel_barcodes, current_barcodes)))},
      ignoreInit = TRUE)
    
    # ----------------------------------------------------
    # Return selected spots
    # ----------------------------------------------------

    selected_spots <- shiny::reactive({
      
      spots_df %>%
        dplyr::filter(barcode %in% selected_barcodes())})
    
    output$selected_table <- shiny::renderTable({
      selected_spots()})
    
    output$download_selected <- shiny::downloadHandler(
      filename = function() {
        download_filename},
      content = function(file) {
        readr::write_csv(selected_spots(), file)})}
  
  return(list("spots_df"=spots_df,
              "img_ext" = img_ext,
              "user_interface"=ui, 
              "server_settings"=server))}
