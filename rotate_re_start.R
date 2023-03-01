# function to rotate so that x component coincides with the x axis
rotate_re_start <- function(sac_x, sac_y, do_plot=FALSE) {
  require(assertthat)
  require(pracma)
  are_equal(length(sac_x), length(sac_y))
  # rotate the data so that the origin is the first sample and the last sample lies on the same plane
  # see: https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
  # 1. find the first and last samples
  sac_x_start <- sac_x[1]
  sac_y_start <- sac_y[1]
  sac_x_end <- sac_x[length(sac_x)]
  sac_y_end <- sac_y[length(sac_y)]
  # 2. compute the angle
  sac_angle <- rad2deg(atan2(sac_y_end-sac_y_start, sac_x_end-sac_x_start))
  alpha <- deg2rad(sac_angle * (-1))
  # 3. rotate
  sac_x_rotated <- (sac_x-sac_x_start)*cos(alpha) - (sac_y-sac_y_start)*sin(alpha) + sac_x_start
  sac_y_rotated <- (sac_x-sac_x_start)*sin(alpha) + (sac_y-sac_y_start)*cos(alpha) + sac_y_start
  # make sure this worked
  are_equal(round(sac_y_rotated[1], 3), round(sac_y_rotated[length(sac_y_rotated)], 3))
  # 4. plot?
  if (do_plot) {
    plot(sac_x, sac_y, 
         xlim = c(min(c(sac_x, sac_x_rotated)), max(c(sac_x, sac_x_rotated))), 
         ylim = c(min(c(sac_y, sac_y_rotated)), max(c(sac_y, sac_y_rotated))) )
    points(sac_x_rotated, sac_y_rotated, col = "red")
  }
  return(list(sac_x_rotated, sac_y_rotated, sac_angle))
}