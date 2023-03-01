get_gabor_filter_bank <- function(all_SF, all_Ori, scr_ppd, ppd_scaler = 1, 
                                  do_on_GPU=FALSE, # transfer Gabor filters to GPU memory?
                                  rf_width_override = NaN,
                                  make_gaussian_aperture = FALSE, 
                                  do_plot=FALSE) {
  if (do_on_GPU) {
    require(torch) # we can use torch to save these gabors on the GPU memory
    if (cuda_is_available()) {
      devi <- 'cuda'
    } else {
      devi <- 'cpu'
    }
  }
  # CREATE A BANK OF GABOR FILTERS
  resulting_gabors <- vector(mode = "list", length = length(all_SF)*length(all_Ori))
  result_i <- 0
  for (now_SF in all_SF) {
    for (now_Ori in all_Ori) {
      result_i <- result_i + 1
      this_gabor <- get_gabor_field(rf_freq_dva = now_SF, rf_ori = now_Ori, 
                                    create_aperture = TRUE, 
                                    gaussian_aperture = make_gaussian_aperture,
                                    rf_width_dva = rf_width_override, # TO DO: different eccentricities? -> magnification function: RFarea = A * ecc^(-B)
                                    scr.ppd = scr_ppd, ppd_scaler = ppd_scaler)
      # do we want to perform this operation on the GPU? If so, transfer to GPU using torch
      if (do_on_GPU) {
        this_gabor[[1]] <- torch_tensor(this_gabor[[1]], device = devi)
        this_gabor[[2]] <- torch_tensor(this_gabor[[2]], device = devi)
      }
      resulting_gabors[[result_i]] <- this_gabor
      if (do_plot & !do_on_GPU) {
        p_gabor <- plot_heatmap(resulting_gabors[[result_i]][[1]])
        print(p_gabor)
      }
    }
  }
  return(resulting_gabors)
}