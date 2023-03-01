# function to extract a subportion of the an image according to its coordinates
# performs 2D bilinear interpolation as implemented in the fields package
# by Richard Schweitzer

extract_target_image <- function(eyepos_x, eyepos_y, 
                                 rf_ecc_x=0, rf_ecc_y=0, 
                                 stim_image, 
                                 target_image_size_x, target_image_size_y, target_image_res) {
  # packages required
  require(data.table)
  require(assertthat)
  require(fields)
  # gather some info on the stim_image
  colnames_stim_image <- as.numeric(colnames(stim_image))
  rownames_stim_image <- as.numeric(rownames(stim_image))
  assert_that(!is.null(colnames_stim_image) & !is.null(rownames_stim_image))
  image_obj <- list(x = rownames_stim_image, y = colnames_stim_image, 
                    z = stim_image) # here we subtract the midpoint
  # which portion of the image should we filter?
  select_col <- seq(from = eyepos_x + rf_ecc_x - target_image_size_x/2, 
                    to = eyepos_x + rf_ecc_x + target_image_size_x/2, 
                    by = target_image_res )
  select_row <- seq(from = eyepos_y + rf_ecc_y - target_image_size_y/2, 
                    to = eyepos_y + rf_ecc_y + target_image_size_y/2, 
                    by = target_image_res )
  query_points <- expand.grid(rows = select_row, cols = select_col)
  # select current subimage via linear interpolation (requires package fields)
  target_image <- interp.surface( obj = image_obj, loc = query_points )
  target_image <- matrix(data = target_image, nrow = length(select_row), ncol = length(select_col), 
                         byrow = FALSE, dimnames = list(select_row, select_col))
  # substitute NAs with zeros (it's dark)
  target_image[is.na(target_image)] <- 0
  return(target_image)
}