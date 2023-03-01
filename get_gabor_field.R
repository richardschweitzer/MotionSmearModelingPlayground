get_gabor_field <- function(rf_freq_dva, rf_width_dva=NaN, rf_ori, rf_amp=0.5, 
                            scr.ppd, ppd_scaler = 1, # the ppd scaler can be applied to produce more highres gabor filters
                            create_aperture=TRUE, gaussian_aperture=TRUE) {
  # CREATES A GABOR FILTER OF A CERTAIN SF AND SIZE
  # in case a gaussian aperture is used (with SD = rf_width_dva/6), then rf_width_dva should be 3 times RF width
  # rf_amp: the color range is from 0 to 1, so the amp should cover half of the entire range
  require(pracma)
  require(assertthat)
  # if needed, compute RF size based on SF
  rf_freq <- rf_freq_dva / scr.ppd
  if (is.na(rf_width_dva) | is.null(rf_width_dva)) {
    rf_width_dva <- exp(-0.5 * log(rf_freq_dva) + log(0.8))
  } 
  rf_width <- rf_width_dva * scr.ppd
  rf_mesh_width <- 3 * rf_width
  # create the grating
  mesh_resolution <- 1 / ppd_scaler
  rf_seq <- seq(-(rf_mesh_width/2), (rf_mesh_width/2), by = mesh_resolution)
  meshlist <- meshgrid(x = rf_seq, y = rf_seq) 
  rf_grating_zero <- rf_amp*cos(meshlist$X*(2*pi*rf_freq*cos(rf_ori)) + meshlist$Y*(2*pi*rf_freq*sin(rf_ori)) + 0*pi)
  rf_grating_half <- rf_amp*cos(meshlist$X*(2*pi*rf_freq*cos(rf_ori)) + meshlist$Y*(2*pi*rf_freq*sin(rf_ori)) + 1/2*pi)
  # create the aperture
  rf_field_zero <- rf_grating_zero
  rf_field_half <- rf_grating_half
  if (create_aperture) {
    if (gaussian_aperture) {
      gaussian_aperture_sd <- rf_width / 2 # Aperture width is twice the space constant of the vignetting Gaussian (Anderson & Burr, 1987)
      rf_gaussian <- exp(-(meshlist$X^2+meshlist$Y^2)/(2*gaussian_aperture_sd^2)) # compute gaussian aperture of stimulus
      rf_field_zero <- rf_field_zero * rf_gaussian
      rf_field_half <- rf_field_half * rf_gaussian
    } else {
      rf_circle <- meshlist$X^2 + meshlist$Y^2 <= round(rf_width/2)^2
      rf_field_zero[rf_circle==FALSE] <- 0
      rf_field_half[rf_circle==FALSE] <- 0
      # zap matrix back to size of circle to reduce computation time
      rows_nonzero <- apply(rf_field_zero, 1, function(x) { !all(x==0) })
      cols_nonzero <- apply(rf_field_zero, 2, function(x) { !all(x==0) })
      rf_field_zero <- rf_field_zero[rows_nonzero, cols_nonzero]
      rf_field_half <- rf_field_half[rows_nonzero, cols_nonzero]
    }
  }
  # return data
  assert_that(all(dim(rf_field_half)==dim(rf_field_zero)) & dim(rf_field_zero)[1]==dim(rf_field_zero)[2] )
  rf_info <- c(rf_freq_dva = rf_freq_dva, rf_width_dva = rf_width_dva, rf_width = rf_width, 
               rf_mesh_width = rf_mesh_width, 
               rf_ori = rf_ori, 
               scr.ppd = scr.ppd, ppd_scaler = ppd_scaler, mesh_resolution = mesh_resolution )
  return(list(rf_field_zero, rf_field_half, rf_info))
}