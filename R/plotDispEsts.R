plotDispEsts <- function( cds, ymin, pch=".", pointcol = "black", linecol="#ff000080", ... )
{
  px = rowMeans( counts( cds, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]

  py = fitInfo(cds)$perGeneDispEsts[sel]
  ##stopifnot(all(py>=0, na.rm=TRUE))
  if(!all(py>=0, na.rm=TRUE))
      warning(sprintf("How peculiar, some of the per-gene dispersion estimates are negative (minimum: %5.3g)", min(py, na.rm=TRUE)))
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
  py = pmax(py, ymin)

  plot(px, py, xlab="mean normalised counts", ylab="dispersion",
    log="xy", pch=pch, col=pointcol, ... )
  xg <- 10^seq( -.5, 5, length.out=100 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col=linecol, lwd=4)
}
