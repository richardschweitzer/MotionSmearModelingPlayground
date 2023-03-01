
# we could fit thresholds over time re saccade onset with this: 
compress_fun <- function(t, baseline, sd0, t0, a0) {
  res <- baseline / (1 + ((a0-1) * (exp( -((t-t0)^2) / (2*sd0^2))) ))
  return(res)
}

# use this two-piece function for a bit more flexibility
double_compress_fun <- function(t, baseline, sd0_pre, sd0_post, t0, a0) {
  res_1 <- compress_fun(t, baseline, sd0_pre, t0, a0)
  res_2 <- compress_fun(t, baseline, sd0_post, t0, a0)
  res <- c(res_1[t<t0], res_2[t>=t0])
  return(res)
}

# now we can combine the SF- and time-dependence into one model of saccadic suppression
suppression_fun <- function(baseline_sens=NaN, 
                            time, time_onset=-4, time_sd_1=20, time_sd_2=50,
                            sf, sf_log_slope=-0.43, sf_log_intercept=0.857) {
  # compute the suppression ratio (the size of suppression, as measured in fixation vs saccade)
  ratios <- exp(sf_log_intercept + sf_log_slope * log(sf))
  # use this suppression ratio as the amplitude for the time course
  if (any(is.na(baseline_sens))) {
    suppression <- double_compress_fun(t = time, t0 = time_onset, baseline = 1, 
                                       a0 = ratios, 
                                       sd0_pre = time_sd_1, sd0_post = time_sd_2)
  } else {
    suppression <- double_compress_fun(t = time, t0 = time_onset, baseline = baseline_sens, 
                                       a0 = ratios, 
                                       sd0_pre = time_sd_1, sd0_post = time_sd_2)
  }
  # now we have a suppression time course for a certain SF
  return(suppression)
}

