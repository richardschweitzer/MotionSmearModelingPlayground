get_fft2 <- function(stim_image_gray, scr.ppd, var.name="", 
                     use_torch = FALSE, 
                     downsample_factor = 2) {
  require(data.table)
  require(gsignal)
  require(reshape2)
  require(pracma)
  if (use_torch) { require(torch) }
  dimnames(stim_image_gray) <- NULL
  if (use_torch) {
    do_fft <- gsignal::fftshift(as_array(torch_fft_fft(torch_fft_fft(self = stim_image_gray, 
                                                                     dim = 1), 
                                                       dim = 2)$cpu() ), 
                                MARGIN = c(1,2))
  } else {
    do_fft <- gsignal::fftshift(fft(stim_image_gray), MARGIN = c(1,2))
  }
  max_fft <- do_fft[which.max(abs(do_fft))]
  # plot_heatmap(log(abs(do_fft)))
  # what's the size of the matrix
  fft_size_x = dim(do_fft)[1] # number of rows, this is quite unintuitive, I agree
  fft_size_y = dim(do_fft)[2] # number of columns
  # get the descriptives for plotting etc, according to: 
  # https://de.mathworks.com/matlabcentral/answers/13896-fftshift-of-an-image#answer_19151
  Fs = scr.ppd # pixels per dva
  d = 1/Fs     # dva per pixel
  x = d*seq(0, fft_size_x-1) # size in dva
  y = d*seq(0, fft_size_y-1)
  dFx = Fs/fft_size_x              # cycles per dva
  dFy = Fs/fft_size_y
  Fx = seq(-Fs/2+dFx/2, Fs/2-dFx/2, by = dFx) # Fx = (-Fs/2+dFx/2:dFx:Fs/2-dFx/2)'  
  Fy = seq(-Fs/2+dFy/2, Fs/2-dFy/2, by = dFy) # Fy = (-Fs/2+dFy/2:dFy:Fs/2-dFy/2)'
  if (downsample_factor > 1) {
    require(fields)
    image_obj <- list(x = Fx, y = Fy, z = abs(do_fft)) 
    # which portion of the image should we filter?
    select_x <- seq(-Fs/2+dFx/2, Fs/2-dFx/2, by = dFx * downsample_factor ) # rows
    select_y <- seq(-Fs/2+dFy/2, Fs/2-dFy/2, by = dFy * downsample_factor ) # cols
    query_points <- expand.grid(rows = select_x, cols = select_y)
    # select current subimage via linear interpolation (requires package fields)
    resampled_fft <- interp.surface( obj = image_obj, loc = query_points )
    resampled_fft <- matrix(data = resampled_fft, nrow = length(select_x), ncol = length(select_y), 
                            byrow = FALSE, dimnames = NULL )
    # meshgrid with downsampled dimensions
    meshi = meshgrid(select_y, select_x)
    return_fft <- resampled_fft
  } else {
    meshi = meshgrid(Fy, Fx)
    return_fft <- abs(do_fft)
  }
  # radial spatial frequency here
  SF_im = sqrt(meshi$X^2+meshi$Y^2)
  Ori_im = atan2(meshi$Y, meshi$X)*180/pi
  # now let the model predict
  SF_Ori_df <- merge.data.table(x = reshape2::melt(data = SF_im, value.name = "rf_freq_dva"), 
                                y = reshape2::melt(data = Ori_im, value.name = "rf_ori"), 
                                by = c("Var1", "Var2"))
  SF_Ori_df <- merge.data.table(x = SF_Ori_df, 
                                y = reshape2::melt(data = return_fft, value.name = paste0("power", var.name) ), 
                                by = c("Var1", "Var2") )
  # add the dimension names to the returned matrices
  if (downsample_factor > 1) {
    dimnames(Ori_im) = list(select_x, select_y)
    dimnames(SF_im) = list(select_x, select_y)
    dimnames(return_fft) = list(select_x, select_y) 
  } else {
    dimnames(Ori_im) = list(Fx, Fy)
    dimnames(SF_im) = list(Fx, Fy)
    dimnames(return_fft) = list(Fx, Fy) 
  }
  # also return the imaginary and real part
  return_fft_imag <- Im(do_fft)
  dimnames(return_fft_imag) = list(Fx, Fy) 
  return_fft_real <- Re(do_fft)
  dimnames(return_fft_real) = list(Fx, Fy) 
  return(list(SF_Ori_df, SF_im, Ori_im, return_fft, return_fft_real, return_fft_imag))
}