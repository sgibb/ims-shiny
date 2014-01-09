plotSlice <- function(im, roirange, pixel, ...) {
  colRamp <- colorRamp(c("black", "blue", "green", "yellow", "red", "purple"))

  i <- MALDIquant:::.which.closest(roirange, as.double(colnames(im)))

  mI <- max(im, na.rm=TRUE)
  rM <- rowMeans(im[, i], na.rm=TRUE)
  rM[is.nan(rM)] <- 0
  col <- rgb(colRamp(rM/mI), maxColorValue=255)

  p <- matrix(NA, nrow=nrow(pixel), ncol=ncol(pixel), byrow=TRUE)

  p[pixel] <- col

  plot(NA, type="n", xlim=c(1, ncol(p)), ylim=c(1, nrow(p)),
       xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE, asp=1, ...)
  rasterImage(as.raster(p), xleft=1, xright=ncol(p), ybottom=1, ytop=nrow(p),
              interpolate=FALSE)
}
