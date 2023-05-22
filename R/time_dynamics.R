#' Title
#'
#'Temporal shifts across species' ranges. Compares two rasters of habitability in time 1 and time 2.
#'
#' @param x a lazySDM object
#'
#' @return A plot showing the species' range dynamics. A SpatRaster object with four numeric classes
#' 0=unsuitable areas
#' 1=colonization areas
#' 2=extinction areas
#' 3=stable areas
#' @export
#'
#' @examples
#' data(whiteshark) # load Great White shark occurrences
#' # Now load oceanographic rasters downloaded from the NOAA (https://psl.noaa.gov/ipcc/cmip6/)
#' cmip1=terra::rast(system.file("extdata/SSP1_2.6_1985_2014.tif", package="simpleSDM")) # Historical climatology
#' cmip2=terra::rast(system.file("extdata/SSP1_2.6_2070_2099.tif", package="simpleSDM")) # soft-core future scenario
#' cmip3=terra::rast(system.file("extdata/SSP5_8.5_2070_2099.tif", package="simpleSDM")) # Mad-Max future scenario
#'
#' out1=lazySDM(whiteshark,buff=500,cmip1,cmip2)
#' out2=lazySDM(whiteshark,buff=500,cmip1,cmip3)
#'
#' par(mfrow=c(2,1))
#' time_dynamics(out1)
#' time_dynamics(out2)
#'
time_dynamics=function(x)
{
  a=x$occupancy.time2+2*(x$occupancy.time1)
  land=terra::rast(system.file("extdata/landmass.tif", package="simpleSDM"))[[1]]
  terra::plot(land,col="tan",legend=F,axes=F)
  terra::plot(a,type="classes",col=c("gray90","blue","red","green"),
              levels=c("unsuitable","colonization","extinction","stable"),
              plg=list(x="bottom",horiz=T),add=T)
  a=list(a)
  names(a)=c("dynamic")
  return(a)
}
