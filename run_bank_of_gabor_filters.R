run_bank_of_gabor_filters <- function(eyepos_x, eyepos_y, 
                                      rf_ecc_x=0, rf_ecc_y=0, 
                                      stim_image, stim_image_midpoint=0.5, 
                                      gabor_list, 
                                      do_on_GPU=FALSE, # perform dot product on GPU?
                                      compute_lum_rms=TRUE, # compute local luminance and rms contrast
                                      normalize_brightness=TRUE, # subtract the gabor filter response to a equiluminant field
                                      normalize_size=TRUE, # divide by number of pixels to normalize across different RF sizes
                                      do_plot=FALSE) {
  # packages required
  require(data.table)
  require(assertthat)
  require(fields)
  if (do_on_GPU) {
    require(torch) # we can use torch to save these gabors on the GPU memory
    if (cuda_is_available()) {
      devi <- 'cuda'
    } else {
      devi <- 'cpu'
    }
  }
  # gather some info on the stim_image
  colnames_stim_image <- as.numeric(colnames(stim_image))
  rownames_stim_image <- as.numeric(rownames(stim_image))
  assert_that(!is.null(colnames_stim_image) & !is.null(rownames_stim_image))
  image_obj <- list(x = rownames_stim_image, y = colnames_stim_image, 
                    z = stim_image - stim_image_midpoint) # here we subtract the midpoint
  # RUNS THE BANK OF GABOR FILTERS ON THE IMAGE (GIVEN A CERTAIN EYE POSITION)
  gabor_responses <- vector(mode = "list", length = length(gabor_list))
  if (do_plot) {
    target_images <- vector(mode = "list", length = length(gabor_list))
  }
  target_image <- NULL
  for (gabor_i in 1:length(gabor_list)) {
    # select the right Gabor filter
    gabor_now <- gabor_list[[gabor_i]]
    # check whether we have to extract a new target image according to the size of the RF
    if (is.null(target_image) || !all(dim(target_image)==dim(gabor_now[[1]]))) {
      # # find valid indeces of current path
      # current_row_index <- (round(eyepos_y-dim(gabor_now[[1]])[1]/2)):(round(eyepos_y-dim(gabor_now[[1]])[1]/2)+dim(gabor_now[[1]])[1]-1)
      # current_row_index_valid <- current_row_index>=1 & current_row_index<=nrow(stim_image)
      # current_col_index <- (round(eyepos_x-dim(gabor_now[[1]])[2]/2)):(round(eyepos_x-dim(gabor_now[[1]])[2]/2)+dim(gabor_now[[1]])[2]-1)
      # current_col_index_valid <- current_col_index>=1 & current_col_index<=ncol(stim_image)
      # target_image <- stim_image[current_row_index[current_row_index_valid], 
      #                            current_col_index[current_col_index_valid]]
      # target_image <- target_image - stim_image_midpoint
      # plot_heatmap(target_image)
      
      # which portion of the image should we filter?
      select_col <- seq(from = eyepos_x + rf_ecc_x - as.numeric(gabor_now[[3]]['rf_mesh_width'])/2, 
                        to = eyepos_x + rf_ecc_x + as.numeric(gabor_now[[3]]['rf_mesh_width'])/2, 
                        by = as.numeric(gabor_now[[3]]['mesh_resolution']) )
      assert_that(length(select_col)==ncol(gabor_now[[1]]) & length(select_col)==ncol(gabor_now[[2]]) )
      select_row <- seq(from = eyepos_y + rf_ecc_y - as.numeric(gabor_now[[3]]['rf_mesh_width'])/2, 
                        to = eyepos_y + rf_ecc_y + as.numeric(gabor_now[[3]]['rf_mesh_width'])/2, 
                        by = as.numeric(gabor_now[[3]]['mesh_resolution']) )
      assert_that(length(select_row)==nrow(gabor_now[[1]]) & length(select_row)==nrow(gabor_now[[2]]) )
      assert_that(length(select_row)==length(select_col))
      query_points <- expand.grid(rows = select_row, cols = select_col)
      # select current subimage via linear interpolation (requires package fields)
      target_image <- interp.surface( obj = image_obj, loc = query_points )
      target_image <- matrix(data = target_image, nrow = length(select_row), ncol = length(select_col), 
                             byrow = FALSE, dimnames = list(select_row, select_col))
      # substitute NAs with zeros (it's dark)
      target_image[is.na(target_image)] <- 0
      # plot_heatmap(target_image)
      if (do_plot) {
        target_images[[gabor_i]] <- target_image
      }
      # check whether all of this worked
      assert_that(all(dim(target_image)==dim(gabor_now[[1]])))
      assert_that(all(dim(gabor_now[[2]])==dim(gabor_now[[1]])))
      # compute the mean of the image, in case we want the brightness regulation?
      if (do_on_GPU) { # crucially, should this target image be on the GPU?
        target_image <- torch_tensor(target_image, device = devi)
        if (normalize_brightness) {
          target_image_mean <- torch_full(size = dim(target_image), 
                                          fill_value = torch_mean(target_image))
        } else {
          target_image_mean <- NULL
        }
      } else {
        if (normalize_brightness) {
          target_image_mean <- target_image
          target_image_mean[ , ] <- mean(target_image_mean)
        } else {
          target_image_mean <- NULL
        }
      }
    } # end of checking whether we have to extract a new target image
    # extract the filter energy
    gabor_response <- get_gabor_filter_resp(target_image = target_image, gabor_pair = gabor_now, 
                                            target_image_mean = target_image_mean, 
                                            do_normalize = compute_lum_rms, 
                                            divide_by_size = normalize_size, 
                                            use_GPU = do_on_GPU )
    gabor_info <- data.table(t(gabor_now[[3]]))
    gabor_info[ , ':=' (gabor_i = gabor_i, rf_size_pix = dim(target_image)[1], 
                        rf_ecc_x_pix = rf_ecc_x, rf_ecc_y_pix = rf_ecc_y, 
                        resp_zero = gabor_response[[1]], resp_half = gabor_response[[2]], 
                        L = gabor_response[[3]], C = gabor_response[[4]]) ]
    # gabor_info[ , gabor_i := gabor_i ]
    # gabor_info[ , rf_size_pix := dim(target_image)[1]]
    # gabor_info[ , rf_ecc_x_pix := rf_ecc_x]
    # gabor_info[ , rf_ecc_y_pix := rf_ecc_y]
    # gabor_info[ , resp_zero := gabor_response[[1]]]
    # gabor_info[ , resp_half := gabor_response[[2]]]
    # gabor_info[ , L := gabor_response[[3]]]
    # gabor_info[ , C := gabor_response[[4]]]
    gabor_responses[[gabor_i]] <- gabor_info
  }
  gabor_responses <- rbindlist(gabor_responses)
  if (do_plot) {
    return(list(gabor_responses, target_images))
  } else {
    return(gabor_responses)
  }
}