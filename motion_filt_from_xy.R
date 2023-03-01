# custom motion filter from x-y position vector
motion_filt_from_xy <- function(x_pos, y_pos,
                                center_around_ith_pos = NaN, # which position should we center the matrix around?
                                w_pos = NULL, # weights for each position
                                spatial_pad_size = 10, 
                                spatial_resolution = 1, 
                                pointwidth = 1, # how large (in the defined grid) is a point?
                                diagnostic_plots = FALSE
) {
  require(pracma)
  require(assertthat)
  n_samples <- length(x_pos) # how many samples have we got?
  are_equal(n_samples, length(y_pos))
  assert_that(is.null(w_pos) || length(w_pos)==length(y_pos) )
  
  # global parameters
  critical_dist <- sqrt(pointwidth^2 + pointwidth^2)
  if (is.na(center_around_ith_pos)) {
    center_around_ith_pos <- n_samples
  }
  
  # create the pixel matrix with dimensions according to the positions provided
  # this code also makes sure that the matrix has odd dimensions
  x_span <- (ceiling(max(x_pos))+spatial_pad_size) - (floor(min(x_pos))-spatial_pad_size)
  mat_size_x = c(head(seq(x_pos[center_around_ith_pos]-x_span, x_pos[center_around_ith_pos], 
                          by = spatial_resolution), -1), 
                 x_pos[center_around_ith_pos],
                 seq(x_pos[center_around_ith_pos], x_pos[center_around_ith_pos]+x_span, 
                     by = spatial_resolution)[-1] )# columns, x coordinate
  y_span <- (ceiling(max(y_pos))+spatial_pad_size) - (floor(min(y_pos))-spatial_pad_size)
  mat_size_y = c(head(seq(y_pos[center_around_ith_pos]-y_span, y_pos[center_around_ith_pos], 
                          by = spatial_resolution), -1), 
                 y_pos[center_around_ith_pos],
                 seq(y_pos[center_around_ith_pos], y_pos[center_around_ith_pos]+y_span, 
                     by = spatial_resolution)[-1] ) # row, y coordinate
  mat_basic <- matrix(data = 0, nrow = length(mat_size_y), ncol = length(mat_size_x),
                      dimnames = list(mat_size_y, mat_size_x) )
  # create the meshgrid that we'll use to locate points within the matrix
  mesh_basic <- meshgrid(x = mat_size_x, y = mat_size_y)
  # fill the matrix with the positions
  mat_final <- mat_basic # this is where the sum of all position weights will be stored
  for (sample_i in (1:n_samples)) {
    x_now <- x_pos[sample_i] # our current positions
    y_now <- y_pos[sample_i] 
    if (is.null(w_pos)) { w_now <- 1 } else { w_now <- w_pos[sample_i] }
    dist_to_x <- abs(mesh_basic$X - x_now)
    dist_to_y <- abs(mesh_basic$Y - y_now)
    dist_to_xy <- sqrt(dist_to_x^2 + dist_to_y^2)
    # find the indeces in the grid that are closest to the current position
    index_pixels_closest <- which(dist_to_x < pointwidth & dist_to_y < pointwidth)
    dist_pixels_closest <- dist_to_xy[index_pixels_closest]
    # the weights of the pixels should be inversely proportional to their distance to the position
    weight_pixels_closest <- critical_dist - dist_pixels_closest
    norm_weight_pixels_closest <- w_now * weight_pixels_closest / sum(weight_pixels_closest)
    # this is the weight matrix
    weight_mat <- mat_basic
    weight_mat[index_pixels_closest] <- norm_weight_pixels_closest
    # if (diagnostic_plots) {
    #   dimnames(weight_mat) <- list(mat_size_y, mat_size_x)
    #   p <- plot_heatmap(weight_mat, do_print = FALSE) 
    #   p <- p + geom_point(aes(x = x_now, y = y_now), color = "red")
    #   print(p)
    # }
    # add to the existing weights
    mat_final <- mat_final + weight_mat
  }
  # normalize for convolution
  mat_final <- mat_final / sum(mat_final)
  # final plot
  if (diagnostic_plots) {
    if (!is.null(w_pos)) {
      df_pos <- data.frame(x = x_pos, y = y_pos, z = w_pos)
    } else {
      df_pos <- data.frame(x = x_pos, y = y_pos, z = 1)
    }
    p <- plot_heatmap(mat_final, do_print = FALSE) 
    # p <- p + geom_point(data = df_pos, aes(x = x, y = y, color = z), fill = NA, shape = 0) +
    #   scale_color_viridis_c(option = "inferno")
    p <- p + labs(x = "x position", y = "y position", color = "position weights", 
                  fill = "resulting filter value")
    print(p)
  }
  return(mat_final)
}