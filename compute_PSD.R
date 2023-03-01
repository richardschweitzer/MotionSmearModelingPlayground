# compute power spectrum for a time series
compute_PSD <- function(time_ms, signal, 
                        pad_n_zeros=5000,
                        return_these_freq=c(),
                        use_method="multitaper", # can be "spectrum", "psd", "spec.ar", "multitaper", "welch"
                        do_plot=FALSE) {
  require(assertthat)
  if (use_method=="psd") {
    require(psd)
  } else if (use_method=="multitaper") {
    require(multitaper)
  } else if (use_method=="welch") {
    require(bspec)
  }
  # from: https://math.mcmaster.ca/~bolker/eeid/2010/Ecology/Spectral.pdf
  # see also: https://vru.vibrationresearch.com/lesson/what-is-the-psd/
  del <- unique(round(diff(time_ms/1000), 10)) # sampling interval
  assert_that(is.sorted(time_ms), msg = "time_ms is not sorted! Check the signal!")
  if (length(del)!=1) {
    warning('Sampling does not seem to be uniform!')
    del <- median(diff(time_ms/1000))
  }
  if (pad_n_zeros>0) {
    signal_padded <- c(signal, rep(0, pad_n_zeros))
  } else {
    signal_padded <- signal
  }
  if (use_method=="spectrum" | use_method=="spec.ar") { # use standard periodogram included in R
    if (use_method=="spectrum") {
      x.spec <- spectrum(signal_padded, log="no", span=10, plot=FALSE)
    } else {
      x.spec <- spec.ar(signal_padded, plot=FALSE) # this produces much cleaner results, but is parametric
    }
    spx <- x.spec$freq/del
    spy <- 2*x.spec$spec
  } else if (use_method=="psd") { # use psd package: https://cran.r-project.org/web/packages/psd/vignettes/psd_overview.pdf
    x.spec1 <- pspectrum(x = signal_padded, verbose = FALSE, x.frqsamp = 1/del)
    spx <- x.spec1$freq
    spy <- x.spec1$spec
  } else if (use_method=="multitaper") { # use multitaper package: https://rdrr.io/cran/multitaper/man/spec.mtm.html
    timeser <- ts(data = signal, frequency = 1/del, start = min(time_ms/1000))
    x.spec <- spec.mtm(timeSeries = timeser, nFFT = pad_n_zeros + length(timeser), plot = FALSE )
    spx <- x.spec$freq
    spy <- x.spec$spec
  } else if (use_method=="welch") { # use WelchPSD: https://www.rdocumentation.org/packages/bspec/versions/1.6/topics/welchPSD
    timeser <- ts(data = signal_padded, frequency = 1/del, start = min(time_ms/1000))
    x.spec <- welchPSD(timeser, seglength=0.004)
    spx <- x.spec$frequency
    spy <- x.spec$power
  }
  # interpolate to extract very specific frequencies
  if (is.null(return_these_freq)) {
    return_spx <- spx
    return_spy <- spy
  } else { # estimate in the specified range
    interpo_fun <- approxfun(x = spx, y = spy, method = "linear", yleft = NaN, yright = NaN)
    return_spx <- return_these_freq
    return_spy <- interpo_fun(return_spx)
  }
  if (do_plot) {
    plot(return_spx, return_spy, xlab="frequency", ylab="spectral density", type="l")
  }
  return(list(tfreq = return_spx, PSD = return_spy))
}

