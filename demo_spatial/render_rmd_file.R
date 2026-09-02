# render rmd file to docs
repo_dir <- getwd()
repo_dir <- normalizePath(file.path(repo_dir),mustWork = TRUE)
rmarkdown::render(
  input = file.path(repo_dir,"demo_spatial", "Reconstruction_models.Rmd"),
  output_file = "Reconstruction_models",
  output_dir = file.path(repo_dir, "docs"))

# render rmd file to docs
repo_dir <- getwd()
repo_dir <- normalizePath(file.path(repo_dir),mustWork = TRUE)
rmarkdown::render(
  input = file.path(repo_dir,"vignettes", "Fiducial_Detection.Rmd"),
  output_file = "Fiducial_Detection",
  output_dir = file.path(repo_dir, "docs"))

# render rmd file to docs
repo_dir <- getwd()
repo_dir <- normalizePath(file.path(repo_dir),mustWork = TRUE)
rmarkdown::render(
  input = file.path(repo_dir,"vignettes", "Spatial_Resolution.Rmd"),
  output_file = "Spatial_Resolution",
  output_dir = file.path(repo_dir, "docs"))