# this is the combined function for orientation, spatial frequency in a global motion display
# by Richard Schweitzer

gabor_ori_vel <- function(gabor_ori, global_dir, global_vel) {
  # gabor_ori is the individual grating orientation
  # global_dir is the overall direction of the movement at global_vel
  require(pracma)
  # according to the equation (Scarfe & Johnston, 2011, Eq. 1): Sc = Sg * cos(Oc - Og)
  # here, however, we need sine because we want velocity to be zero when differences of angles are zero
  gabor_individual_drift_vel = global_vel * sin(deg2rad(gabor_ori) - deg2rad(global_dir) )
  return(gabor_individual_drift_vel)
}

kelly_sens_ori_vel_sf <- function(gabor_ori, gabor_sf, global_dir, global_vel, 
                                  global_orth_vel = NaN, return_final_vel = FALSE) {
  # if global_dir is set to NaN, then all orientations get the same velocity, which is global_vel
  if (any(is.na(global_dir))) { # set velocity constant for all orientations
    final_gabor_vel <- global_vel
  } else { # compute velocities for varying orientations
    # given a given global motion, how fast does a Gabor of a certain orientation move?
    # think of seeing a grating of a certain orientation move in a certain direction at a certain velocity
    gabor_vel <- gabor_ori_vel(gabor_ori = gabor_ori, global_dir = global_dir, global_vel = global_vel)
    # assuming that this is not an ideal straight motion, we can think of some (very small) orthogonal component,
    # such as saccade curvature. We should compute that to make sure that velocity is never zero
    if (any(is.na(global_orth_vel))) {
      global_orth_vel <- global_vel / 40
    }
    gabor_orth_vel <- gabor_ori_vel(gabor_ori = gabor_ori, global_dir = global_dir-90, global_vel = global_orth_vel)
    # combine these velocities
    final_gabor_vel <- sqrt(gabor_vel^2 + gabor_orth_vel^2)
  }
  # now we compute the contrast sensitivity to that certain gabor grating given its velocity and SF
  # we add this baseline_vel to make sure that velocity is never zero (which may happen when parallel)
  if (return_final_vel) {
    gabor_sens <- final_gabor_vel
  } else {
    gabor_sens <- kelly_vel(sf = gabor_sf, v = final_gabor_vel)
  }
  return(gabor_sens)
}