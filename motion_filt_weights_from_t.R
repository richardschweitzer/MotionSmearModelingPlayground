# Functions to compute time-dependent weights
# by Richard Schweitzer

## the Gamma temporal response function, for fitting, and with amplitude
gamma_special <- function(x, latency, shape, scale, norm_amp = NaN) {
  if (any(shape <= 0) | any(scale <= 0)) {
    gamma_res <- dgamma(x = x - latency, 
                        shape = 0.1, scale = 0.1) 
  } else {
    gamma_res <- dgamma(x = x - latency, 
                        shape = shape, scale = scale) 
  }
  if (all(!is.na(norm_amp))) { # normalize amplitude to one?
    if (max(gamma_res) != 0) {
      gamma_res <- gamma_res / max(gamma_res)
    }
    gamma_res <- gamma_res * norm_amp
  }
  return(gamma_res)
}

## function to assign weights:
motion_filt_weights_from_t <- function(t_ms, gamma_shap=1.72, gamma_sca=20.94) {
  w <- vector(mode = "numeric", length = length(t_ms))
  t_ms <- t_ms - t_ms[1]
  max_t_ms <- NaN
  for (it in length(t_ms):1) {
    # find out when the maximum at the last time point occurs
    if (is.na(max_t_ms) && it==length(t_ms)) {
      last_times <- seq(t_ms[it], t_ms[it]+500, by = 0.1)
      last_gamma <- gamma_special(x = last_times, 
                                  latency = t_ms[it], norm_amp = NaN, 
                                  shape = gamma_shap, scale = gamma_sca)
      max_t_ms <- last_times[which.max(last_gamma)]
    }
    # now find out what the expected contribution for stimuli occurring at different
    # latencies should be at that time
    w[it] <- gamma_special(x = max_t_ms, latency = t_ms[it], norm_amp = NaN, 
                           shape = gamma_shap, scale = gamma_sca)
  }
  # normalize so that last value (largest value) is 1
  w <- w / max(w)
  return(w)
}