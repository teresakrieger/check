
# define directories & slide ids
demo_dir <- "C:/Users/Eli/Desktop/Patho/temp_git/visiumSpatialSolutions/demo_spatial"
base_dir <- file.path("D:", "Visium_spatial")
out_dir <- file.path(base_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
slide_IDs <- list("spatial_Glioblastoma" = "V52Y10-317-B1",
                  "spatial_Human_Ovarian_Cancer" = "V10L13-020-D1",
                  "spatial_Human_Ovarian_Carcinoma" ="V52Y09-003-B1",
                  "spatial_Lung_Cancer" = "V52Y10-286-B1",
                  "spatial_Lung_Squamous_Cell_Carcinoma" = "V42A20-354-D1",
                  "spatial_Skin_Melanoma" = "V42A20-355-A1")

array_size <- list("spatial_Human_Ovarian_Cancer" = "6.5mm",
                   "spatial_Human_Ovarian_Carcinoma" ="11mm",
                   "spatial_Glioblastoma" = "11mm",
                   "spatial_Lung_Cancer" = "11mm",
                   "spatial_Lung_Squamous_Cell_Carcinoma" = "6.5mm",
                   "spatial_Skin_Melanoma" = "6.5mm")

# load functions
source(file.path(demo_dir, "functions_spatial_alignment.R"))

#--------------------------------------------------------------#
# Internal Validation of the 10X Loupe Browser generated Manual Alignment JSON file
#--------------------------------------------------------------#

# Check Transformation Matrices

json_check <- list()
for (sample_id in names(slide_IDs)){
  json_dir <- file.path(base_dir, sample_id)
  json_check[[sample_id]] <- check_transformation_matrix_json(json_dir= json_dir, slide_id = slide_IDs[[sample_id]])}


#--------------------------------------------------------------#

# conversion error
# A: Slide2Image 
# A_inv: Image2Slide
df_A <- data.frame(matrix(nrow=length(slide_IDs), ncol=4), row.names = names(slide_IDs))
df_A_inv <- df_A
colnames(df_A) <- c("Sample_Name", "Array_Size", "max_abs_imageX", "max_abs_imageY")
colnames(df_A_inv) <- c("Sample_Name", "Array_Size", "max_abs_x", "max_abs_y")
for (row in names(slide_IDs)){
  df_A[row,] <- c(row, array_size[[row]], json_check[[row]]$slide2image$error_oligo)
  df_A_inv[row,] <- c(row, array_size[[row]], json_check[[row]]$image2slide$error_oligo_inv)}

utils::write.csv(df_A, file.path(out_dir, "slide2image_transformation_inaccuracy.csv"), row.names = FALSE)
utils::write.csv(df_A_inv, file.path(out_dir, "image2slide_transformation_inaccuracy.csv"), row.names = FALSE)

#--------------------------------------------------------------#

# transformation matrix
df_A <- data.frame(matrix(nrow=length(slide_IDs)*3, ncol=7), 
                   row.names =c(paste0(rep(names(slide_IDs), each=3),c("_row1","_row2", "_row3"))))
colnames(df_A) <- c("Sample_Name", "Array_Size","lm_identical", "row", "col1", "col2", "col3")
df_A_inv <- data.frame(matrix(nrow=length(slide_IDs)*3, ncol=7), 
                       row.names =c(paste0(rep(names(slide_IDs), each=3),c("_row1","_row2", "_row3"))))
colnames(df_A_inv) <- c("Sample_Name", "Array_Size","lm_identical", "row", "col1", "col2", "col3")
for (row in names(slide_IDs)){
  df_A[paste0(row, "_row1"),] <- c(row, array_size[[row]],json_check[[row]]$slide2image$A_identical, "row1", json_check[[row]]$slide2image$A[1,])
  df_A[paste0(row, "_row2"),] <- c(row, array_size[[row]],json_check[[row]]$slide2image$A_identical, "row2",json_check[[row]]$slide2image$A[2,])
  df_A[paste0(row, "_row3"),] <- c(row, array_size[[row]],json_check[[row]]$slide2image$A_identical, "row3",json_check[[row]]$slide2image$A[3,])
  df_A_inv[paste0(row, "_row1"),] <- c(row, array_size[[row]],json_check[[row]]$image2slide$A_inv_identical, "row1",json_check[[row]]$image2slide$A_inv[1,])
  df_A_inv[paste0(row, "_row2"),] <- c(row, array_size[[row]], json_check[[row]]$image2slide$A_inv_identical,"row2",json_check[[row]]$image2slide$A_inv[2,])
  df_A_inv[paste0(row, "_row3"),] <- c(row, array_size[[row]], json_check[[row]]$image2slide$A_inv_identical,"row3",json_check[[row]]$image2slide$A_inv[3,])}

utils::write.csv(df_A, file.path(out_dir, "slide2image_transformation_matrix_A.csv"), row.names = FALSE)
utils::write.csv(df_A_inv, file.path(out_dir, "image2slide_transformation_matrix_A_inv.csv"), row.names = FALSE)

#--------------------------------------------------------------#

# Slide Coordinates
slide_coordinates <- list()
for (sample_id in names(json_check)){
  json_data <- json_check[[sample_id]]$json_file
  slide_coordinates[[sample_id]] <- derive_slide_coordinates_json(json_file = json_data)}

slide_df <- data.frame(matrix(nrow=length(slide_IDs), ncol=11), row.names = names(slide_IDs))
colnames(slide_df) <- c("Sample_Name", "Array_Size","dx", "dy", "xy_ratio",
                        "range_cols",  "range_rows", "range_x", "range_y", "range_imageX", "range_imageY")
for (sample_id in names(slide_IDs)){
  df <- json_check[[sample_id]]$json_file$oligo
  slide_df[sample_id,] <- c(sample_id, array_size[[sample_id]], slide_coordinates[[sample_id]]$oligo[c("dx", "dy", "xy_ratio")],
                      paste(range(df$col), collapse = "-"), paste(range(df$row), collapse = "-"), 
                      paste(range(df$x), collapse = "-"),paste(range(df$y), collapse = "-"),
                      paste(range(df$imageX), collapse = "-"), paste(range(df$imageY), collapse = "-"))}

utils::write.csv(slide_df, file.path(out_dir, "slide_coordinates.csv"), row.names = FALSE)

#--------------------------------------------------------------#
# Fiducial Symbol Positions in tissue image upload
#--------------------------------------------------------------#

# JSON and Tissue position combination
symbol_list <- list()
for (sample_id in names(json_check)){
  json_data <- json_check[[sample_id]]$json_file
  symbol_list[[sample_id]] <- fiducial_symbols_json(json_file=json_data, sample_id=sample_id)}

corner_positions <- c("Upper_left", "Lower_left", "Lower_right", "Upper_right")
symbol_df <- data.frame(matrix(nrow=length(slide_IDs), ncol=12), row.names = names(slide_IDs))
colnames(symbol_df) <- c("Sample_Name", "Array_Size",  corner_positions, "dx", "dy", "xy_ratio",
                         "range_cols",  "range_rows", "range_x", "range_y", "range_imageX", "range_imageY")
for (sample_id in names(json_check)){
  df <- symbol_list[[sample_id]]$symbol_df
  d_oligo <- unique(json_check[[sample_id]]$json_file$oligo$dia)
  d_fid <- unique(json_check[[sample_id]]$json_file$fiducial$dia)
  stopifnot(length(d_oligo)==1L && length(d_fid)==1L)
  symbol_df[sample_id,] <- c(sample_id, array_size[[sample_id]], 
                             unname(unlist(symbol_list[[sample_id]]$symbol_corners)[match(corner_positions, names(unlist(symbol_list[[sample_id]]$symbol_corners)))]),
                             slide_coordinates[[row]]$fiducials[c("dx", "dy", "xy_ratio")],
                             paste(range(df$col), collapse = "-"), paste(range(df$row), collapse = "-"), 
                             paste(range(df$x), collapse = "-"),paste(range(df$y), collapse = "-"),
                             paste(range(df$imageX), collapse = "-"), paste(range(df$imageY), collapse = "-"))}

utils::write.csv(symbol_df, file.path(out_dir, "symbol_df.csv"), row.names = FALSE)

#--------------------------------------------------------------#
# Relation between tissue positions and the manual 10X Loupe Manual Alignment file
#--------------------------------------------------------------#

# JSON and Tissue position combination
tp_list <- list()
for (sample_id in names(json_check)){
  json_data <- json_check[[sample_id]]$json_file
  tp_list[[sample_id]] <- combine_tp_json(json_file=json_data, base_dir=base_dir, sample_id=sample_id)}

n_spots <- data.frame(matrix(nrow=length(slide_IDs), ncol=5), row.names = names(slide_IDs))
colnames(n_spots) <- c("Sample_Name", "Array_Size", "n_spots_tp", "n_spots_json", "n_spots_combined")
for (sample_id in names(tp_list)){
  n_spots[sample_id,] <- c(sample_id, array_size[[sample_id]], unname(unlist(tp_list[[sample_id]]$n_spots)))}
utils::write.csv(n_spots, file.path(out_dir, "n_spots_combined.csv"), row.names = FALSE)

#--------------------------------------------------------------#

# Diameter information
diameter_cols <- c("diameter_oligo", "diameter_fiducials")
diameter_df <- data.frame(matrix(nrow=length(slide_IDs), ncol=10), row.names = names(slide_IDs))
colnames(diameter_df) <- c("Sample_Name", "Array_Size",  diameter_cols, "det_A", "d_oligo/d_fid", "d_fid*det_A", "(d_oligo/d_fid)/(d_fid*det_A)",
                           "sf_tp_oligo", "sf_tp_fiducials")

for (sample_id in names(json_check)){
  d_oligo <- unique(json_check[[sample_id]]$json_file$oligo$dia)
  d_fid <- unique(json_check[[sample_id]]$json_file$fiducial$dia)
  det_A <- det(json_check[[sample_id]]$slide2image$A)
  sf <- tp_list[[sample_id]]$scalefactors
  stopifnot(length(d_oligo)==1L && length(d_fid)==1L && det_A != 0)
  diameter_df[sample_id,] <- c(sample_id, array_size[[sample_id]], d_oligo,  d_fid, det_A, d_oligo/d_fid, d_fid*det_A, (d_oligo/d_fid)/(d_fid*det_A),
                               sf$spot_diameter_fullres, sf$fiducial_diameter_fullres)}

utils::write.csv(diameter_df, file.path(out_dir, "diameter_df.csv"), row.names = FALSE)


#--------------------------------------------------------------#

#  Corresponding Tissue Positions
selected_spots <- list("spatial_Human_Ovarian_Cancer" = "V10L13-020-D1_selected_spots",
                       "spatial_Human_Ovarian_Carcinoma" ="V52Y09-003-B1_selected_spots")

check_selection <- list()
for (sample_id in names(selected_spots)){
  json_dir <- file.path(base_dir, sample_id)
  check_selection[[sample_id]] <- check_transformation_matrix_json(json_dir= json_dir, slide_id = selected_spots[[sample_id]])}

tp_info <- list()
for (sample_id in names(selected_spots)){
  json_data <- check_selection[[sample_id]]$json_file
  tp_info[[sample_id]] <- combine_tp_json(json_file=json_data, base_dir=base_dir, sample_id=sample_id)}

for (sample_id in names(tp_info)){
  df <- as.data.frame(table(tp_info[[sample_id]]$tp_json[,c("in_tissue", "tissue")], useNA="ifany"))
  colnames(df) <- c("tp_in_tissue", "json_in_tissue", "n")
  utils::write.csv(df, file.path(out_dir, paste0("tissue_position_check_",sample_id,".csv")), row.names = FALSE)}


#--------------------------------------------------------------#

# pixel conversion Matrices
tp_json_conversion <- list()
for (sample_id in names(slide_IDs)){
  tp_json <- tp_list[[sample_id]]$tp_json
  tp_json_conversion[[sample_id]] <- pixel_conversion_json_tp(tp_json= tp_json, sample_id= sample_id)}

# transformation matrix
df_B <- data.frame(matrix(nrow=length(slide_IDs)*3, ncol=6), 
                   row.names =c(paste0(rep(names(slide_IDs), each=3),c("_row1","_row2", "_row3"))))
colnames(df_B) <- c("Sample_Name", "Array_Size", "row", "col1", "col2", "col3")
df_B_inv <- data.frame(matrix(nrow=length(slide_IDs)*3, ncol=7), 
                       row.names =c(paste0(rep(names(slide_IDs), each=3),c("_row1","_row2", "_row3"))))
colnames(df_B_inv) <- c("Sample_Name", "Array_Size","lm_identical", "row", "col1", "col2", "col3")
for (row in names(slide_IDs)){
  df_B[paste0(row, "_row1"),] <- c(row, array_size[[row]], "row1", tp_json_conversion[[row]]$JSON2TP$B[1,])
  df_B[paste0(row, "_row2"),] <- c(row, array_size[[row]], "row2",tp_json_conversion[[row]]$JSON2TP$B[2,])
  df_B[paste0(row, "_row3"),] <- c(row, array_size[[row]], "row3",tp_json_conversion[[row]]$JSON2TP$B[3,])
  df_B_inv[paste0(row, "_row1"),] <- c(row, array_size[[row]],tp_json_conversion[[row]]$TP2JSON$B_inv_identical, "row1",tp_json_conversion[[row]]$TP2JSON$B_inv[1,])
  df_B_inv[paste0(row, "_row2"),] <- c(row, array_size[[row]], tp_json_conversion[[row]]$TP2JSON$B_inv_identical,"row2",tp_json_conversion[[row]]$TP2JSON$B_inv[2,])
  df_B_inv[paste0(row, "_row3"),] <- c(row, array_size[[row]], tp_json_conversion[[row]]$TP2JSON$B_inv_identical,"row3",tp_json_conversion[[row]]$TP2JSON$B_inv[3,])}

utils::write.csv(df_B, file.path(out_dir, "JSON2TP_transformation_matrix_B.csv"), row.names = FALSE)
utils::write.csv(df_B_inv, file.path(out_dir, "TP2JSON_transformation_matrix_B_inv.csv"), row.names = FALSE)


#--------------------------------------------------------------#

# conversion error
# B: JSON2TP
# B_inv: TP2JSON
df_B <- data.frame(matrix(nrow=length(slide_IDs), ncol=4), row.names = names(slide_IDs))
df_B_inv <- df_B
colnames(df_B) <- c("Sample_Name", "Array_Size", "max_abs_pxl_row", "max_abs_pxl_col")
colnames(df_B_inv) <- c("Sample_Name", "Array_Size", "max_abs_imageX", "max_abs_imageY")
for (row in names(slide_IDs)){
  df_B[row,] <- c(row, array_size[[row]], tp_json_conversion[[row]]$JSON2TP$error_pixels)
  df_B_inv[row,] <- c(row, array_size[[row]], tp_json_conversion[[row]]$TP2JSON$error_image_inv)}

utils::write.csv(df_B, file.path(out_dir, "JSON2TP_transformation_inaccuracy.csv"), row.names = FALSE)
utils::write.csv(df_B_inv, file.path(out_dir, "TP2JSON_transformation_inaccuracy.csv"), row.names = FALSE)

