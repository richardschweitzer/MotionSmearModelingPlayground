# D.H. Kelly's velocity function
# from: Kelly, D. H. (1979). Motion and vision. II. Stabilized spatio-temporal threshold surface. Josa, 69(10), 1340-1349.
kelly_vel <- function(sf, v) {
  alpha <- sf  * 2*pi
  # G <- (6.1 + 7.3 * (abs(log10(v/3))^3)) * v * sf^2 * exp(-2*sf*(v+2)/45.9)
  k <- 6.1 + 7.3 * abs(log10(v/3))^3
  alpha_max <- 45.9 / (v+2)
  G <- k*v*(alpha^2) * exp(-2*alpha/alpha_max)
  return(G)
}