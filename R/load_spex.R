#' Load a `SpatialExperiment` Object for a 10x Genomics Visium sample
#'
#' Load a Visium feature-barcode matrix and image as \code{SpatialExperiment}
#'
#' @param visium_dir Character scalar. Base directory containing the Visium
#'   sample directories.
#' @param visium_sample Character scalar. Name of the Visium sample directory.
#' @param matrix_name Character scalar. Name of Matrix passed to
#'   \code{\link[SpatialExperiment]{read10xVisium}} through its
#'   \code{data} argument. Available options are \code{"filtered"} (default)
#'   or \code{"raw"}.
#' @param matrix_file_format Character scalar. Matrix type used for the \code{type} 
#'  argument in \code{\link[SpatialExperiment]{read10xVisium}}, 
#'  either \code{"h5"} or \code{"sparse"}.
#' @param image_resolution Character scalar. Image resolution of the tissue image
#' file to load, either  \code{"lowres"} or \code{"hires"}.
#'
#' @return A \code{SpatialExperiment} object.
#'
#' @details
#' The function assumes the standard Space Ranger directory structure:
#'
#' \preformatted{
#' <visium_dir>/<visium_sample>/outs/
#' }
#' 
#' @seealso
#' \code{\link[SpatialExperiment]{read10xVisium}}
#'
#' @examples
#' \dontrun{
#' spex <- load_spex(
#'   visium_dir = "/path/to/visium",
#'   visium_sample = "sample_1",
#'   matrix_name = "filtered",
#'   matrix_file_format = "h5",
#'   image_resolution = "lowres"
#' )
#' }
#'
#' @export
load_spex <- function(visium_dir, visium_sample,
                     matrix_name="filtered",matrix_file_format=c("h5", "sparse"),
                     image_resolution=c("lowres", "hires")){

  # sanity check
  matrix_file_format <- match.arg(matrix_file_format, choices = c("h5", "sparse"))
  image_resolution <- match.arg(image_resolution)
  matrix_name <- match.arg(matrix_name, choices = c("filtered", "raw"))
  
  type_use <- switch(matrix_file_format, h5 = "HDF5", sparse = "sparse")

  file_dir <- file.path(visium_dir, visium_sample, "outs")
  
  if (!dir.exists(file_dir)){
    stop("Can't find directory ", file_dir, ".")}
  
  ##############################################
  
  # load data
  message("Loading Visium data from ", file_dir)
  
  # SpatialFeatureExperiment::read10xVisiumSFE() not feasible with 11mm array & filtered matrix 
  # due to internal barcode-location reference for findVisiumGraph() with 4992 original Visium positions
  # furthermore, only "full_res_image_pixel" or microns are accepted as spatial units file.path(visium_dir,visium_sample)
  spex <- SpatialExperiment::read10xVisium(samples = file_dir, sample_id = visium_sample,
                                          type = type_use, data = matrix_name, images = image_resolution, load = FALSE)

    return(spex)}

