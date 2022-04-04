######################################
# To create an outline of ellipsoid
######################################

ellipse.outline <- function(n, a, b, c=0, xc=0, yc=0, zc=0) {a <- 20
   x <- 2*a/n*(0:n)-a
   y <- 2*b/n*(0:n)-b
   pts <- {}

   for (X in 1:length(x)) {
      Z <- 1 - (x[X]/a)^2 - (y/b)^2
      Z[abs(Z)< 10^(-12)] <- 0
      z <- c*sqrt(Z[Z>=0])
      pts <- rbind(pts, cbind(x[X],y[Z>=0],z))
      pts <- rbind(pts, cbind(x[X],y[Z>=0],-z))
   }
   dimnames(pts) <- list(NULL, c("x","y","z"))
   pts <- t(t(pts) + c(xc, yc, zc))
   pts
} 

