# this function computes the gabor filter response
get_gabor_filter_resp <- function(target_image, gabor_pair, 
                                  target_image_mean, 
                                  divide_by_size=TRUE, 
                                  do_normalize=TRUE, 
                                  use_GPU=FALSE) {
  # EXTRACTS FILTER ENERGY FROM TARGET_IMAGE
  if (use_GPU) {
    require(torch)
  }
  # perform dot product here
  if (!is.null(target_image_mean)) { # there is a mean image
    if (use_GPU) {
      dot_zero <- torch_tensordot(target_image, gabor_pair[[1]]) - torch_tensordot(target_image_mean, gabor_pair[[1]])
      dot_half <- torch_tensordot(target_image, gabor_pair[[2]]) - torch_tensordot(target_image_mean, gabor_pair[[2]])
    } else {
      dot_zero <- sum(target_image*gabor_pair[[1]], na.rm = TRUE) - sum(target_image_mean*gabor_pair[[1]], na.rm = TRUE)
      dot_half <- sum(target_image*gabor_pair[[2]], na.rm = TRUE) - sum(target_image_mean*gabor_pair[[2]], na.rm = TRUE)
    }
  } else { # mean image provided
    if (use_GPU) {
      dot_zero <- torch_tensordot(target_image, gabor_pair[[1]])
      dot_half <- torch_tensordot(target_image, gabor_pair[[2]])
    } else {
      dot_zero <- sum(target_image*gabor_pair[[1]], na.rm = TRUE)
      dot_half <- sum(target_image*gabor_pair[[2]], na.rm = TRUE)
    }
  }
  # size normalization?
  if (divide_by_size==TRUE) {
    if (use_GPU) {
      # this is very inefficient and thus commented out
      #nonnan_len <- target_image$numel() - torch_sum(target_image$isnan())$cuda()
      #resp_zero <- dot_zero / nonnan_len 
      #resp_half <- dot_half / nonnan_len
      # we need to copy those back to the CPU
      resp_zero <- as.numeric(dot_zero$cpu()) / length(target_image) # if one cell is NA, resp will be NA anyway
      resp_half <- as.numeric(dot_half$cpu()) / length(target_image)
    } else {
      resp_zero <- dot_zero / length(target_image[!is.na(target_image)])
      resp_half <- dot_half / length(target_image[!is.na(target_image)])
    }
  } else {
    if (use_GPU) {
      resp_zero <- as.numeric(dot_zero$cpu())
      resp_half <- as.numeric(dot_half$cpu())
    } else {
      resp_zero <- dot_zero
      resp_half <- dot_half
    }
  }
  
  # do normalization, as described by: https://www.sciencedirect.com/science/article/pii/S0042698905005559
  if (do_normalize==TRUE) {
    w <- matrix(NaN, nrow = dim(target_image)[1], ncol = dim(target_image)[2])
    siz <- dim(target_image)[1] / 3 
    for (row_i in 1:dim(w)[1]) {
      for (col_i in 1:dim(w)[2]) {
        w[row_i, col_i] = 0.5 * (cos(pi/(siz/2) * sqrt((row_i - siz/2)^2 + (col_i - siz/2)^2)) + 1)
      }
    }
    w <- w / sum(w)
    L <- sum(w*target_image) / sum(w) # local luminance
    C <- sqrt(sum(w*(((target_image-L)^2)/L^2))) # local RMS contrast
    # TO DO: adjust resp amp according to L and/or C ??? BUT HOW ???
    
  } else {
    L <- NaN
    C <- NaN
  }
  return(list(resp_zero, resp_half, L, C))
}