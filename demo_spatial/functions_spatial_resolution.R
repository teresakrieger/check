

# Helper function to generate 6 outer vertices for a single center point
make_hexagon <- function(cx, cy, r_size, upper_geom="line") {
  
  upper_geom <- match.arg(upper_geom, choices = c("line", "corner"))
  
  # Generate angles for the first 6 points only
  angles <- seq(0, length.out = 6) * (2 * pi / 6)
  
  if (upper_geom == "corner"){
    x_hex <- cx + r_size*sin(angles)
    y_hex <- cy + r_size*cos(angles)}

  if (upper_geom == "line"){
  x_hex <- cx + r_size*cos(angles)
  y_hex <- cy + r_size*sin(angles)}
  
  # Append first point to the end to force absolute closure
  x_hex <- c(x_hex, x_hex[1])
  y_hex <- c(y_hex, y_hex[1])
  
  return(sf::st_polygon(list(matrix(c(x_hex, y_hex), ncol = 2))))}

# ------------------------------------------------------------

# Calculate geometry for each case
calculate_spatial_case <- function(case_id, L, n_rows, n_cols, n_cases =6, spot_diameter = 55, tolerance = 1e-3) {
  
  # sanity check
  stopifnot(length(case_id) == 1L, case_id %in% seq_len(n_cases),
    is.numeric(L), length(L) == 1L, is.finite(L), L > 0,
    is.numeric(n_rows), length(n_rows) == 1L, n_rows == floor(n_rows), n_rows >= 2L,
    is.numeric(n_cols), length(n_cols) == 1L, n_cols == floor(n_cols), n_cols >= 2L,
    is.numeric(spot_diameter), length(spot_diameter) == 1L, is.finite(spot_diameter), spot_diameter > 0,
    is.numeric(tolerance), length(tolerance) == 1L, is.finite(tolerance), tolerance >= 0)
  
  r <- spot_diameter / 2
  # dx: physical increment per array_col index
  # dy: physical increment per array_row index
  # m_x: margin along the spatial-col/x axis
  # m_y: margin along the spatial-row/y axis
  if (case_id == 1L) {
    # case 1: 100µm vertical & horizontal distance
    case_name <- paste0("100µm column & row distance between array units")
    dx <- 100
    dy <- 100
    m_y <- (L - (n_rows - 1L) * dx) / 2
    m_x <- (L - (n_cols - 1L) * dy) / 2
  } else if (case_id == 2L) {
    # case 1: 50µm vertical & horizontal distance
    case_name <- paste0("50µm column & row distance between array units")
    dx <- 50
    dy <- 50
    m_y <- (L - (n_rows - 1L) * dx) / 2
    m_x <- (L - (n_cols - 1L) * dy) / 2
  } else if (case_id == 3L) {
    # case 2: 100µm diagonal distance & equal margins
    case_name <- paste0("100µm column & diagonal spot distance in hexagonal grid")
    dx <- 50
    dy <- 50*sqrt(3)
    m_y <- (L - (n_rows - 1L) * dy) / 2
    m_x <- (L - (n_cols - 1L) * dx) / 2
  } else if (case_id == 4L) {
    # case 4: 50µm vertical distance & equal margins
    case_name <- paste0("50µm column distance between array units & equal margins")
    dx <- 50 
    m_x <- (L - (n_cols - 1L) * dx) / 2 
    m_y <- m_x 
    dy <- (L - 2 * m_y) / (n_rows - 1L)
  } else if (case_id == 5L){
    # case 5:  Minimal margins
    case_name <- paste0("Minimal margins (m = r)")
    m_x <- r
    m_y <- m_x
    dy <- ( L - 2*r)/(n_rows-1)
    dx <- (L - 2*r)/(n_cols-1)
  }else if (case_id == 6L){
    # case 8: m_axis = d_axis
    case_name <- paste0("Axis-distance-margin relation of 1:1")
    dy <- L / (n_rows + 1L)
    dx <- L / (n_cols + 1L)
    m_x <- dx
    m_y <- dy}
  
  # calculate diagonal distance
  d_diag <- sqrt(dx^2 + dy^2)

  # upper centre-coordinate limits
  y_max_center <- m_y + (n_rows - 1L) * dy
  x_max_center <- m_x + (n_cols - 1L) * dx
  
  # Complete occupied limits, including spot radius
  x_min <- m_x -r
  x_max <- x_max_center + r
  
  y_min <- m_y - r
  y_max <- y_max_center + r
  
  df <- data.frame(case_id = case_id, 
             case_name = case_name,
             L = L,
             n_rows = n_rows,
             n_columns = n_cols,
             n_spots = n_rows * n_cols / 2,
             spot_diameter = spot_diameter,
            
             d_row  = dy,
             d_column = dx,
             d_diagonal = d_diag,
            
             min_center_row = m_y,
             max_center_row = y_max_center,
             min_center_column= m_x,
             max_center_column = x_max_center,
                        
             min_spot_row = y_min,
             max_spot_row = y_max,
             min_spot_column = x_min,
             max_spot_column = x_max,
             
             identical_margins = isTRUE(all.equal(m_x, m_y, tolerance = tolerance)),
            
             # Candidate interpretations of the reported 100 um distance
             matches_tested_100 = any(
               abs(d_diag - 100) <= tolerance,
               abs(dx - 100) <= tolerance,
               abs(dy - 100) <= tolerance,
               abs(2 * dx - 100) <= tolerance,
               abs(2 * dy - 100) <= tolerance),
            
             # physical feasibility
             # no spot overlaps
             spots_rows_nonoverlap = dy >= r - tolerance,
             spots_columns_nonoverlap= dx >= r - tolerance,
             spots_diagonal_nonoverlap = d_diag >= spot_diameter - tolerance,
             spots_nonoverlapping = dy >= r - tolerance && dx >= r - tolerance && d_diag >= spot_diameter - tolerance,
             # within array
             spots_rows_inside =  y_min >= - tolerance && y_max <= L + tolerance,
             spots_columns_inside = x_min >= - tolerance && x_max <= L + tolerance,
             case_possible = (y_min >= - tolerance && y_max <= L + tolerance) &&
                             (x_min >= - tolerance && x_max <= L + tolerance) &&
               (dx >= r - tolerance && dy >= r - tolerance && d_diag >= spot_diameter - tolerance))
  
  return(df)}


# ------------------------------------------------------------

# Generate parity-compatible spot coordinates
generate_spot_coordinates <- function(n_rows, n_cols, d_row, d_col, m_row, m_col) {
  
  spots <- expand.grid(array_row = 0:(n_rows - 1L), array_col = 0:(n_cols - 1L))
  
  # Orange-crate packing:
  # valid spots have matching row/column parity
  spots <- spots[spots$array_row %% 2L == spots$array_col %% 2L, , drop = FALSE]
  
  spots$parity <- ifelse(spots$array_row %% 2L == 0L, "even", "odd")
  
  spots$spatial_row <- (m_row + spots$array_row * d_row)
  
  spots$spatial_col <- (m_col + spots$array_col * d_col)
  
  rownames(spots) <- NULL
  return(spots)}

# ------------------------------------------------------------

# Plot geometry for each case
plot_spatial_case <- function(result, point_cex = 0.35, axis_step = NULL) {
  
  stopifnot(nrow(result) == 1L)
  
  L <- result$L
  
  if (is.null(axis_step)) {
    axis_step <- if (L <= 6500) 500 else 1000}
  
  spots <- generate_spot_coordinates(
    n_rows = result$n_rows,
    n_cols = result$n_columns,
    d_row = result$d_row,
    d_col = result$d_column,
    m_row = result$min_center_row,
    m_col = result$min_center_column)
  
  stopifnot(nrow(spots) == result$n_spots)
  
  # add padding in µm
  plot_padding <- max(100,
    -result$min_spot_row, result$max_spot_row -L,
    -result$min_spot_column, result$max_spot_column -L,
    na.rm = TRUE) + 10
  
  plot(NA, 
       xlim = c(-plot_padding, L + plot_padding), 
       ylim = c(L + plot_padding, -plot_padding),
       xlab = "Spatial column (µm)", 
       ylab = "Spatial row (µm)",
       xaxt = "n", yaxt = "n", asp = 1,
       main = paste0(result$case_id,") ", result$case_name),
       cex.main=0.95)
  
  axis(side = 1, at = seq(0, L, by = axis_step), las = 2, cex.axis = 0.85)
  axis(side = 2, at = seq(0,L, by = axis_step), las = 2, cex.axis = 0.85)
  
  # Add capture area
  rect(xleft = 0, ybottom = 0, xright = L, ytop = L, col = "skyblue", border = "black")
  
  # Add spots
  even <- spots$parity == "even"
  points(x = spots$spatial_col[even], y = spots$spatial_row[even], pch = 19, cex = point_cex, col = "black")
  odd <- spots$parity == "odd"
  points(x = spots$spatial_col[odd], y = spots$spatial_row[odd], pch = 19, cex = point_cex, col = "blue")
  
  invisible(spots)}

