#' Download vignette files
#'
#' @param output_dir Parent output directory, where a folder for each visium 
#' sample, matching the expected Space Ranger output structure, will be created.
#' @param visium_sample Character vector. Supported samples are 
#' \describe{
#'   \item{\code{Visium_FFPE_Human_Ovarian_Cancer}}{Ovarian Cancer sample
#'   generated using Space Ranger Version 1.3.0 with an Array Size of 6.5mm.}
#'   \item{\code{CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma}}{Ovarian Cancer
#'   Sample generated using Space Ranger Version 2.0.0 with an Array Size of 11mm.}}
#'   The sample name  \code{visium_sample} has to match the file name from 
#'   the `Batch download` command as specified on the Visium dataset download site with
#'  `https://cf.10xgenomics.com/samples/spatial-exp/<spaceranger_version>/<visium_sample>`.
#'  @param spaceranger_version Optional, named character vector. Spaceranger 
#'  version used to generate the `visium_sample` file. Defaults to NULL, since
#'  the SpaceRanger versions are resolved automatically for the vignette data.
#' @param matrix_file_format Character scalar. Defines the file format of the 
#'  matrix, which will be downloaded. Available options are `"h5"` for 
#'  downloading the `h5` file or `"sparse"` for the `tar.gz` file. The 
#'  compressed `tar.gz` file will be unpacked after the successful download,
#'  and the `tar.gz` file will be removed, thereafter.
#' @param matrix_name Character vector specifying one or both matrices to
#'   download. Available values are \code{"filtered"} (default) and \code{"raw"}.
#' @param download_handler Optional, character scalar. Download handler 
#' used for downloading the requested files. Supported methods are `curl`, 
#'  `libcurl`, `wget`, `auto`. If `NULL` (default) the currently selected 
#'  download handler will be used.
#' @param image_resolution Character vector specifying one or both tissue-image
#'   resolutions to retain after extraction. Available options are `lowres` 
#'   for the low resolution image (default) and `hires` for the 
#'  high resolution image.
#' 
#' @return Output directory where the downloaded files and folders can be found, 
#' matching the Space Ranger output structure `<visium_dir>/<visium_sample>/outs/`
#'
#' @examples
#' \dontrun{
#' # Download Human Ovarian Cancer 6.5mm & Human Ovarian Cancer 11mm datasets
#'  with low & high resolution image
#' download_vignette_files(
#'   output_dir = "path/to/output/directory",
#'   visium_sample = c("Visium_FFPE_Human_Ovarian_Cancer", 
#'   "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma"),
#'   image_resolution=c("lowres","hires")
#' )
#'
#' # Download Human Ovarian Cancer 11mm dataset with low resolution image
#' download_vignette_files(
#'   output_dir = "path/to/output/directory",
#'   visium_sample =  "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma"
#' )
#' }
#'
#' @export
download_vignette_files <- function(output_dir, visium_sample = c("Visium_FFPE_Human_Ovarian_Cancer", 
                                                                   "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma"),
                                    matrix_file_format= "h5", matrix_name = "filtered", 
                                    download_handler=NULL, image_resolution="lowres") {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)

  # sanity checks
  default_samples <- c("Visium_FFPE_Human_Ovarian_Cancer", "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma")
  visium_sample <- match.arg(visium_sample, choices = default_samples, several.ok = TRUE)
  matrix_file_format <- match.arg(matrix_file_format,choices = c("h5", "sparse"), several.ok = TRUE)
  matrix_name <- match.arg(matrix_name, choices = c("filtered", "raw"), several.ok = TRUE)
  image_resolution <- match.arg(image_resolution, choices = c("lowres", "hires"), several.ok = TRUE)  
  
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("'output_dir' must be one non-empty character value.")}
  
  # check download handler
  valid_download_handlers <- c( "curl", "libcurl", "wget", "auto")
  if (is.null(download_handler)) {
    download_handler <- getOption("download.file.method", default = "auto")}
  
  if (!is.character(download_handler) || length(download_handler) != 1L || is.na(download_handler) || !download_handler %in% valid_download_handlers) {
    stop("'download_handler' must be one of ", paste(valid_download_handlers, collapse = ", "), ".")}
  
  if (download_handler %in% c("curl", "wget") && !nzchar(Sys.which(download_handler))) {
    stop("Download method '", download_handler, "' was requested, but the executable was not found.")}
  
  message("Using ", download_handler, " for downloading the requested files.")
  
  # define spatial files
  keep_spatial <- c("detected_tissue_image.jpg", "scalefactors_json.json")
  
  if ("lowres" %in% image_resolution){
    keep_spatial <- c(keep_spatial, "tissue_lowres_image.png")}
  if ("hires" %in% image_resolution){
    keep_spatial <- c(keep_spatial, "tissue_hires_image.png")}
  
    spaceranger_version <- list("Visium_FFPE_Human_Ovarian_Cancer" = "1.3.0",
                                "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma" = "2.0.0")
  # create directory
  for (sample_name in visium_sample){
    output_folder <- file.path(output_dir, sample_name, "outs")
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    output_folder <- normalizePath(output_folder, mustWork = TRUE)
    spaceranger_id <- spaceranger_version[[sample_name]] 
    
    options(timeout = max(300, old_timeout))
    message("Downloading ",sample_name," with
          \n- ",paste(matrix_name, sep=" and ")," matrix in ", paste(matrix_file_format, collapse=" and "), " matrix_file_format," , 
            "\n- detected_tissue_image,",
            "\n- ",paste(image_resolution, collapse=" and ")," image",
            "\n- scalefactors_json.json (Distinct from 10X Loupe Browser Manual Alignment json file!)")
    
    # define download paths
    base_url <- paste0("https://cf.10xgenomics.com/samples/spatial-exp/", spaceranger_id,"/", sample_name)  
    matrix_suffix <- switch(matrix_file_format, h5=".h5", sparse= ".tar.gz")
    online_files <- c(paste0(sample_name, "_",matrix_name,"_feature_bc_matrix",matrix_suffix),
               paste0(sample_name, "_spatial.tar.gz"))
  
    # download & extract tar files
    for (file in online_files){
      local_file <- gsub(paste0("^", sample_name, "_"), "", file)
      dest_file <- file.path(output_folder, local_file)
      status <-  utils::download.file(url = paste0(base_url,"/", file), 
                                      destfile = dest_file, 
                                      method = download_handler, mode = "wb")
      
      if (!identical(status, 0L)) {
        stop("Download failed for '", file, "' with status ", status, ".")}
      if (grepl("\\.tar\\.gz$", local_file)){
        untar_status <- utils::untar(tarfile=dest_file, exdir= output_folder)
        new_folder <- file.path(output_folder, 
                                gsub("\\.tar\\.gz$", "",  local_file))
        
        if (!dir.exists(new_folder) || ! identical(untar_status, 0L)) {
          warning("Extraction of  ",local_file, " didn't work. 
                  \nPlease extract the file manually.")
        }else{
          file.remove(dest_file)}}}
      
      # cleanup not required spatial files
    spatial_dir <- file.path(output_folder, "spatial")
    if (!dir.exists(spatial_dir)){
      stop("Can't find spatial directory: ", spatial_dir)}
      spatial_files <- list.files(spatial_dir)
      remove_spatial <- spatial_files[!spatial_files %in% keep_spatial &
          !grepl("^tissue_positions", spatial_files)]
      file.remove(file.path(spatial_dir, remove_spatial))
      
    message("Downloaded requested vignette files for ",sample_name," to ", output_folder, ".")}
  
  return(invisible(NULL))}


#---------------------------------------------------------------------------#

# matrix_sizes <- list(
#   h5 = c(filtered = 26548409, raw = 40881630),
#   sparse = c(filtered = 82878588, raw = 110985967))
# matrix_extracted <- list(
#   h5 = c(filtered = 26548409, raw = 40881630),
#   sparse = c(filtered = 82690062+155265+19786, raw= 110646809+297011+25674))
# spatial_size <- 16119457
# netto_spatial <- 1817752 + # detected tissue image
#   186470 +  # tissue positions
#   165  # scalefactors_json
# if("lowres" %in% image_resolution){
#   netto_spatial <- netto_spatial + 1045184}  # tissue lowres
# if ("hires" %in% image_resolution){
#   netto_spatial <- netto_spatial + 11409022}
# 
# matrix_sizes <- list(
#   h5 = c(filtered = 37939439, raw = 55712353),
#   sparse = c(filtered = 122922633, raw = 153825900))
# matrix_extracted <- list(
#   h5 = c(filtered = 37939439, raw = 55712353),
#   sparse = c(filtered = 122720911+156435+25627, raw=153429930+297935+76003))
# spatial_size <- 35315564    
# netto_spatial <- 1749037 + # detected tissue image
#   559837 + # tissue positions
#   207  # scalefactors_json
# if("lowres" %in% image_resolution){
#   netto_spatial <- netto_spatial + 318324} # tissue lowres
# if ("hires" %in% image_resolution){
#   netto_spatial <- netto_spatial + 3831187}
# 
# data_size <- sum(matrix_sizes[[matrix_file_format]][matrix_name])
# data_extracted <- sum(matrix_extracted[[matrix_file_format]][matrix_name])
# 
# options(timeout = max(300, old_timeout))
# message("Downloading ",file_name," with
#           \n- Matrix Size: ",format_mb(data_size), 
#         "\n-Spatial Data Size (brutto): ", format_mb(spatial_size), 
#         "\n-Spatial Data Size (netto): ", format_mb(netto_spatial),
#         "\n-total downloaded: ", format_mb(sum(data_size, spatial_size)),
#         "\n-total kept: ", format_mb(sum(data_extracted, netto_spatial)))

