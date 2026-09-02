
# Download Vignette Files

# download 6.5mm files
mkdir D:\Visium_6500
cd D:\Visium_6500
$baseUrl = "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer"

$files = @(
    "Visium_FFPE_Human_Ovarian_Cancer_possorted_genome_bam.bam",
    "Visium_FFPE_Human_Ovarian_Cancer_possorted_genome_bam.bam.bai",
    "Visium_FFPE_Human_Ovarian_Cancer_molecule_info.h5",
    "Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5",
    "Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.tar.gz",
    "Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.h5",
    "Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.tar.gz",
    "Visium_FFPE_Human_Ovarian_Cancer_analysis.tar.gz",
    "Visium_FFPE_Human_Ovarian_Cancer_spatial.tar.gz",
    "Visium_FFPE_Human_Ovarian_Cancer_spatial_enrichment.csv",
    "Visium_FFPE_Human_Ovarian_Cancer_metrics_summary.csv",
    "Visium_FFPE_Human_Ovarian_Cancer_web_summary.html",
    "Visium_FFPE_Human_Ovarian_Cancer_cloupe.cloupe"
)

foreach ($file in $files) {
    Write-Host "Downloading $file"
    curl.exe -L -C - -o $file "$baseUrl/$file"

    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Download failed for $file"
    }
}

# inspect tar files
tar -tzf Visium_FFPE_Human_Ovarian_Cancer_spatial.tar.gz
# spatial/
# spatial/detected_tissue_image.jpg
# spatial/aligned_fiducials.jpg
# spatial/tissue_positions_list.csv
# spatial/tissue_lowres_image.png
# spatial/scalefactors_json.json
# spatial/tissue_hires_image.png
tar -tzf Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.tar.gz
# filtered_feature_bc_matrix/
# filtered_feature_bc_matrix/barcodes.tsv.gz
# filtered_feature_bc_matrix/features.tsv.gz
# filtered_feature_bc_matrix/matrix.mtx.gz
tar -tzf Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.tar.gz
# raw_feature_bc_matrix/
# raw_feature_bc_matrix/barcodes.tsv.gz
# raw_feature_bc_matrix/features.tsv.gz
# raw_feature_bc_matrix/matrix.mtx.gz

# extract tar files
tar -xzf Visium_FFPE_Human_Ovarian_Cancer_spatial.tar.gz
tar -xzf Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.tar.gz
tar -xzf Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.tar.gz

# size of files in current folder
ls | measure length -s
# Count    : 13
# Average  :
# Sum      : 13072112981
# Maximum  :
# Minimum  :
# Property : length

# Data for vignette
# -a----         8/12/2026   4:13 PM       26548409 Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5
# -a----         8/12/2026   4:13 PM       82878588 Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.tar.gz
# -a----         8/12/2026   4:13 PM       40881630 Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.h5
# -a----         8/12/2026   4:13 PM      110985967 Visium_FFPE_Human_Ovarian_Cancer_raw_feature_bc_matrix.tar.gz
# ----                 -------------         ------ ----
# -a----          4/7/2022   9:52 PM        1817752 detected_tissue_image.jpg
# -a----          4/7/2022   9:52 PM         186470 tissue_positions_list.csv
# -a----          4/7/2022   9:52 PM        1045184 tissue_lowres_image.png
# -a----          4/7/2022   9:52 PM       11409022 tissue_hires_image.png

# -------------------------------------------------------------------- #

# download 11mm files
mkdir D:\Visium_11000
cd D:\Visium_11000
$baseUrl = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma"

$files = @(
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_possorted_genome_bam.bam",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_possorted_genome_bam.bam.bai",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_molecule_info.h5",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.h5",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.h5",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.tar.gz",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_analysis.tar.gz",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_metrics_summary.csv",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_web_summary.html",
    "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_cloupe.cloupe"
)

foreach ($file in $files) {
    Write-Host "Downloading $file"
    curl.exe -L -C - -o $file "$baseUrl/$file"

    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Download failed for $file"
    }
}


# inspect tar files
tar -tzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz
# spatial/
# spatial/cytassist_image.tiff
# spatial/detected_tissue_image.jpg
# spatial/scalefactors_json.json
# spatial/tissue_hires_image.png
# spatial/tissue_lowres_image.png
# spatial/tissue_positions.csv
# spatial/aligned_tissue_image.jpg
# spatial/spatial_enrichment.csv
# spatial/aligned_fiducials.jpg

tar -tzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz
# filtered_feature_bc_matrix/
# filtered_feature_bc_matrix/features.tsv.gz
# filtered_feature_bc_matrix/matrix.mtx.gz
# filtered_feature_bc_matrix/barcodes.tsv.gz

tar -tzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.tar.gz
# raw_feature_bc_matrix/
# raw_feature_bc_matrix/barcodes.tsv.gz
# raw_feature_bc_matrix/matrix.mtx.gz
# raw_feature_bc_matrix/features.tsv.gz

# extract tar files
tar -xzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz
tar -xzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz
tar -xzf CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.tar.gz

# size of files in current folder
ls | measure length -s
# Count    : 12
# Average  :
# Sum      : 13771460755
# Maximum  :
# Minimum  :
# Property : length

# Data for vignette
# -a----         8/12/2026   5:08 PM       37939439 CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.h5
# -a----         8/12/2026   5:09 PM      122922633 CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz
# -a----         8/12/2026   5:09 PM       55712353 CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.h5
# -a----         8/12/2026   5:09 PM      153825900 CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_raw_feature_bc_matrix.tar.gz
# ----                 -------------         ------ ----
# -a----         7/14/2022   6:34 PM        1749037 detected_tissue_image.jpg
# -a----         7/14/2022   6:34 PM        3831187 tissue_hires_image.png
# -a----         7/14/2022   6:34 PM         318324 tissue_lowres_image.png
# -a----         7/14/2022   6:34 PM         559837 tissue_positions.csv

# -------------------------------------------------------------------- #

# size of folders and subfolders
ls -r | measure length -s

# -------------------------------------------------------------------- #
# Validation of spatial geometry
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# 11mm spatial files
# -------------------------------------------------------------------- #

# download 11mm files
cd mnt/d
mkdir Visium_spatial
cd mnt/d/Visium_spatial
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz
tar -tzf CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz
#spatial/
#spatial/cytassist_image.tiff
#spatial/spatial_enrichment.csv
#spatial/detected_tissue_image.jpg
#spatial/tissue_positions.csv
#spatial/scalefactors_json.json
#spatial/tissue_hires_image.png
#spatial/aligned_tissue_image.jpg
#spatial/aligned_fiducials.jpg
#spatial/tissue_lowres_image.png

# extract tar files
tar -xzf CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz
# rename spatial folder to avoid overwriting
mv spatial spatial_Glioblastoma

# CytAssist_11mm_FFPE_Human_Lung_Cancer
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz
tar -tzf CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz
tar -xzf CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz
mv spatial spatial_Lung_Cancer

# -------------------------------------------------------------------- #
# 6.5mm spatial files
# -------------------------------------------------------------------- #

# download 6.5mm files
# CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz
tar -tzf CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz
tar -xzf CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz
mv spatial spatial_Lung_Squamous_Cell_Carcinoma

wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz
tar -tzf CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz
tar -xzf CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz
mv spatial spatial_Skin_Melanoma