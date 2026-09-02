# define directories
code_dir <- file.path(getwd(), "demo_spatial")
dir.create(file.path(code_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(code_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# load functions
source(file.path(code_dir, "functions_spatial_resolution.R"))

# ============================================================

# Demo Array with small margins
x <- "array_col"
y <- "array_row"
demo_rows <- 0:8
demo_cols <- 0:13
n_rows <- length(demo_rows)
n_cols <- length(demo_cols)

# Obtain small Visium demo grid
even_spots <- expand.grid("array_row"=seq(min(demo_rows),max(demo_rows),by=2),
                          "array_col"=seq(min(demo_cols),max(demo_cols),by=2))
even_spots$type <- "even"
# Obtain odd spot positions
odd_spots <- expand.grid("array_row"=seq(min(demo_rows)+1,max(demo_rows),by=2),
                         "array_col"=seq(min(demo_cols)+1,max(demo_cols),by=2))
odd_spots$type <- "odd"
array_grid <- rbind(even_spots, odd_spots)

# spatial conversion
dx <- 50
dy <- sqrt(3)*50
spot_diameter <- 55
r <- spot_diameter/2

# define array boundary
x_min <- -dx
x_max <- (max(demo_cols)+1)*dx
y_min <- -dy
y_max <- (max(demo_rows)+1)*dy

spatial_grid <- array_grid
spatial_grid$spatial_x <- array_grid$array_col*dx
spatial_grid$spatial_y <- array_grid$array_row*dy

# create spatial sf object
spatial_sf <- sf::st_as_sf(spatial_grid, coords = c("spatial_x", "spatial_y"), remove = FALSE)
spatial_sf <- sf::st_buffer(spatial_sf, dist = r)

# png image
grDevices::png(filename = file.path(code_dir, "figures", "demo_array_positions.png"), width=10, height=13, units = "cm", res=600)

plot(NULL, xlim=c(x_min, x_max), xlab=x, xaxt="n", 
     ylim=c(y_max, y_min), ylab=y,  yaxt="n", main="Demo Array with small margins", asp=1, cex.main=0.9, cex.axis=0.9)
rect(xleft = x_min, ybottom = y_min, xright = x_max, ytop= y_max, col="skyblue")
graphics::axis(1, at = demo_cols*dx, cex.lab=0.7, labels = demo_cols)
graphics::axis(2, at = demo_rows*dy, las=2, cex.lab=0.7, labels = demo_rows)
graphics::abline(v=demo_cols*dx, col=c("black", "blue"))
graphics::abline(h=demo_rows*dy, col=c("black","blue"))

# Add diagonal lines
# positive-slope diagonals
diag_pos <-sort(unique(array_grid$array_row - array_grid$array_col))
intercept_pos <- diag_pos*dy
for(i in intercept_pos){
  graphics::abline(a=i, b=dy/dx, col="darkblue", lty=3, lwd=1)}

# negative-slope diagonals
diag_neg <-sort(unique(array_grid$array_row + array_grid$array_col))
intercept_neg <- diag_neg*dy
for(i in intercept_neg){
  graphics::abline(a=i, b=-dy/dx, col="darkgrey", lty=3, lwd=1)}

# add points
plot(sf::st_geometry(spatial_sf), pch = 19, col =ifelse(spatial_grid$type == "odd", "blue", "black"), border=NA, add = TRUE)

grDevices::dev.off()

# ============================================================

# calculate outer hexagon bounding radius from point distances
R <-  dx*2/sqrt(3)

# Apply hexagon builder
hex_list <- lapply(1:nrow(spatial_grid), function(i) {
  make_hexagon(spatial_grid$spatial_x[i], spatial_grid$spatial_y[i], R, upper_geom = "corner")})

# Convert polygon collection into an sf object
spatial_hex_sf <- sf::st_sf(type = spatial_grid$type, geometry = sf::st_sfc(hex_list))

# create plot with hexagon boundaries
grDevices::png(filename = file.path(code_dir, "figures", "demo_hexagonal_grid.png"), width=10, height=10, units = "cm", res=600)
plot(sf::st_geometry(spatial_hex_sf), 
     border = ifelse(spatial_hex_sf$type == "odd", "blue", "black"), 
     main = "Hexagonal Grid Layout", 
     lwd = 2)

# add spatial resolved points as filled circles
plot(sf::st_geometry(spatial_sf), pch = 19, col =ifelse(spatial_grid$type == "odd", "blue", "black"), border= NA, add = TRUE)

grDevices::dev.off()

# ============================================================

# Candidate spatial reconstructions for Visium arrays
n_cases <- 6

# Array specifications
array_info <- data.frame(
  array_id = c("6.5mm", "11mm"),
  L = c(6500, 11000),
  n_rows = c(78L, 128L),
  n_cols = c(128L, 224L),
  spot_diameter = c(55, 55))

# Evaluate all arrays and cases
all_results <- vector(mode = "list", length = nrow(array_info) * n_cases)

result_index <- 1L

for (array_index in seq_len(nrow(array_info))) {
  spec <- array_info[array_index, ]
  for (case_id in seq_len(n_cases)) {
    current_result <- calculate_spatial_case(
      case_id = case_id,
      n_cases = n_cases,
      L = spec$L,
      n_rows = spec$n_rows,
      n_cols = spec$n_cols,
      spot_diameter = spec$spot_diameter)
    
    current_result$array_id <- spec$array_id
    
    all_results[[result_index]] <- current_result
    result_index <- result_index + 1L}}

df_cases <- do.call(rbind, all_results)

rownames(df_cases) <- NULL

# reorder columns
df_cases <- df_cases[ , c("array_id", "case_id","case_name", setdiff(colnames(df_cases), c("array_id", "case_id", "case_name")))]

# ------------------------------------------------------------

# Generate a 6-panel figure per array
for (array_index in seq_len(nrow(array_info))) {
  
  spec <- array_info[array_index, ]
  
  array_results <- df_cases[df_cases$array_id == spec$array_id, , drop = FALSE]
  
  # png image
  grDevices::png(filename= file.path(code_dir, "figures", paste0("spatially_resolved_",array_info$L[array_index],"_micrometer_array_positions.png")), 
      width=30, height=22, units = "cm", res=800)
  graphics::par(mfrow=c(2,3), oma=c(0,0,2,0))
  
  # Get active row and column counts
  current_grid <- par("mfrow")
  
  # Multiply the active dimensions and compare to z
  if(current_grid[1] * current_grid[2] != n_cases){
    warning("Plot grid doesn't match number of cases.")}
  
  point_cex <- if (spec$L == 6500) {
    0.26
  } else if (spec$L == 11000) {
    0.1}
  
  for (case_id in seq_len(n_cases)) {
    
    current_result <- array_results[array_results$case_id == case_id, , drop = FALSE]
    
    plot_spatial_case(result = current_result, point_cex = point_cex)}
  
  graphics::mtext(paste0("Spatially resolved Spot Positions for the ",spec$array_id," Visium array"),
        outer = TRUE, side = 3, line = -0.5, cex = 1.2)
  
  grDevices::dev.off()}

# ------------------------------------------------------------

# Save df
utils::write.csv(df_cases,file.path(code_dir, "tables","visium_candidate_spatial_reconstructions.csv"), row.names = FALSE)
