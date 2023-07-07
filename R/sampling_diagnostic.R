#' sampling_diagnostic
#'
#' Diagnostic plots to evaluate the potential impact of sampling effort on SDM predictions
#'
#' @param x a lazySDM object
#' @param s a sampling effort SpatRaster object (see example)
#' @param w spatial window to evaluate local effects (default=5)
#' @param pal a color palette (see example)
#'
#' @return A four panel plot, showing: a) global sampling effort, b)background modeled areas, c) habitability, and c) local correlation between the sampling effort and habitability.
#' A SpatRaster object of local correlations is also provided.
#' @export
#'
#' @examples
#'
#' data(whiteshark)
#' cmip1=terra::rast(system.file("extdata/SSP1_2.6_1985_2014.tif", package="simpleSDM"))
#' cmip2=terra::rast(system.file("extdata/SSP1_2.6_2070_2099.tif", package="simpleSDM"))
#' out=lazySDM(whiteshark,buff=500,cmip1,cmip2)
#'
#' # Load a global raster of the sampling effort, i.e., number of total ocurrences in the OBIS database at one-degree resolution
#' # Data re-analyzed from Chaudhary et al. 2021 (https://doi.org/10.1073/pnas.2015094118)
#' obis=terra::rast(system.file("extdata/OBIS_sampling_effort.tif", package="simpleSDM"))
#' pal1=colorRampPalette(c("blue", "yellow", "red")) # a simple color palette

#' effortbias=sampling_diagnostic(x=out,s=log10(obis+1),w=5,pal=pal1)
#'
#' # An histogram showing the distribution of local correlation values
#' hist(effortbias)
#' global(effortbias,"mean",na.rm=T)
#'
sampling_diagnostic=function(x,s,w=5,pal=pal1)
{
  r_stack=c(x$prob.time1,s)
  fc=focalCor(r_stack, w = w,cor)
  land=terra::rast(system.file("extdata/landmass.tif", package="simpleSDM"))[[1]]
  par(mfrow=c(2,2),mar=c(1,1,1,5))
  plot(land,col="tan",legend=F,axes=F)
  plot(s,add=T,col=pal(100))
  mtext("sampling effort",side=3,adj=0)

  plot(land,col="tan",legend=F)
  plot(x$target.area,col="gray90",legend=F,axes=F,add=T)
  mtext("background area",side=3,adj=0)

  plot(land,col="tan",legend=F,axes=F)
  plot(x$prob.time1,col=pal(100),range=c(0,1),add=T)
  mtext("Habitability",side=3,adj=0)

  plot(land,col="tan",legend=F,axes=F)
  plot(x$target.area,col="gray90",legend=F,axes=F,add=T)
  plot(fc,col=pal(100),range=c(-1,1),add=T)
  mtext("local correlation-sampling effort vs habitability",side=3,adj=0)

  return(fc)
}
