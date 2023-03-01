# Linear Motion filter, translated from Matlab's fspecial('motion')
# by Richard Schweitzer
motion_filt <- function(filter_len, filter_ang) {
  require(pracma) # contains sooooooo many Matlab functions
  len = max(c(1,filter_len))
  half = (len-1)/2  # rotate half length around center
  phi = mod(filter_ang,180)/180*pi
  cosphi = cos(phi)
  sinphi = sin(phi)
  xsign = sign(cosphi)
  linewdt = 1
  # define mesh for the half matrix, eps takes care of the right size for 0 & 90 rotation
  sx = Fix(half*cosphi + linewdt*xsign - len*eps())
  sy = Fix(half*sinphi + linewdt - len*eps())
  meshi = meshgrid(seq(0, sx, by = xsign), 0:sy)
  x = meshi$X
  y = meshi$Y
  # define shortest distance from a pixel to the rotated line
  dist2line = (y*cosphi-x*sinphi) # distance perpendicular to the line
  rad = sqrt(x^2 + y^2)
  # find points beyond the line's end-point but within the line width
  lastpix = which((rad >= half) & (abs(dist2line)<=linewdt))
  # distance to the line's end-point parallel to the line
  x2lastpix = half - abs((x[lastpix] + dist2line[lastpix]*sinphi)/cosphi)
  dist2line[lastpix] = sqrt(dist2line[lastpix]^2 + x2lastpix^2)
  dist2line = linewdt + eps() - abs(dist2line)
  dist2line[dist2line<0] <- 0  # zero out anything beyond line width
  # unfold half-matrix to the full size
  h = rot90(dist2line,2)
  # this Matlab command doesn't work in R: h(end+(1:end)-1,end+(1:end)-1) = dist2line;
  # therefore, we do this:
  h2 <- matrix(data = NaN, nrow = size(h, 1)+size(h, 1)-1, ncol = size(h, 2)+size(h, 2)-1)
  h2[1:size(h, 1), 1:size(h, 2)] <- h
  h2[size(h, 1)+(1:size(h, 1))-1, size(h, 2)+(1:size(h, 2))-1] = dist2line 
  h2[is.na(h2)] <- 0
  h2 = h2/(sum(h2) + eps()*len*len)
  # get the phi right
  if (cosphi>0 && nrow(h2)>1) { h2 = flipud(h2) }
  return(h2)
}