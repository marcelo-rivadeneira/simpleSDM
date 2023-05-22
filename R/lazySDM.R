#' LazySDM
#'
#' Black-box correlational SDM, with minimum specification steps. Calibrates SDM in time1 and projects onto time2.
#' @param dato data.frame containing geographic coordinates of a species presence
#' @param buff geographic buffer (in km) where to generate pseudo-absence. It is also the area where the model is tested.
#' @param stck1 Raster stack (terra class) with environmental variables for time 1
#' @param stck2 Raster stack (terra class) with environmental variables for time 2
#'
#' @return Four objects
#'      summary: Main diagnostic statistics
#'      target.area: raster of the buffer zone
#'      occupancy.time1: raster of estimated habitability for time 1
#'      occupancy.time2: raster of estimated habitability for time 2
#' @import terra
#' @import sp
#' @import SDMtune
#' @import ENMeval
#' @import usdm
#' @import PresenceAbsence
#' @import modEvA
#' @export
#'
#' @examples
#' data(whiteshark)
#' cmip1=terra::rast(system.file("extdata/SSP1_2.6_1985_2014.tif", package="simpleSDM"))
#' cmip2=terra::rast(system.file("extdata/SSP1_2.6_2070_2099.tif", package="simpleSDM"))
#'
#' out=lazySDM(whiteshark,buff=500,cmip1,cmip2)
#' out$summary
#
lazySDM=function(dato,buff=500,stck1,stck2)
{
  options(warn=-1)
  x=dato$x
  y=dato$y
  coords=data.frame(x,y)
  coordinates(coords)= ~x + y
  crs(coords)="+proj=longlat +datum=WGS84"
  buf=vect(buffer(coords,(buff*1e3))) # default, 500 km buffer
  mask=stck1[[1]]  # Create a mask with target resolution and extent from climate layers
  values(mask)[!is.na(values(mask))]=1
  maskedbuf=mask(mask,buf) # Set all raster cells outside the buffer to NA.
  ext.coods=intersect(ext(coords), ext(maskedbuf))
  maskedbuf2=crop(maskedbuf, ext.coods)
  bg_dat=spatSample(maskedbuf2,size=length(x),na.rm=T,as.points=T)
  bg_dat=crds(bg_dat)
  pr_dat=data.frame(x,y)
  puntos=rbind(pr_dat,bg_dat)

  ## Prepare the tunning
  totuneo=prepareSWD(species="bla",p = pr_dat,a = bg_dat, env = stck1)

  ## spatial folding (4 checkerboard) for cross-validation
  check_folds=get.checkerboard2(occs = pr_dat,env=stck1,bg = bg_dat,aggregation.factor = 4)
  folds.p=check_folds$occs.grp
  block_folds_formatted1=matrix(ncol=length(unique(folds.p)),nrow=length(folds.p),F)
  block_folds_formatted2=matrix(ncol=length(unique(folds.p)),nrow=length(folds.p),T)
  for(i in 1:nrow(block_folds_formatted1))
  {
    block_folds_formatted1[i,folds.p[i]]=T
    block_folds_formatted2[i,folds.p[i]]=F
  }
  block_folds_formatted=list(block_folds_formatted1,block_folds_formatted2)
  names(block_folds_formatted)=c("train","test")

  # Cross-validated trained RF model with hyperparameter tunning
  randfo.model=train("RF", data = totuneo, folds = block_folds_formatted,verbose=F)
  randfo.auc.test=SDMtune::auc(randfo.model,test=T)


  # Thresholding
  rf.predicted=predict(randfo.model, data = totuneo)
  spatcoord=totuneo@coords
  coordinates(spatcoord)=~X+Y
  umbral=cbind(1:length(totuneo@pa),sobs=totuneo@pa,
               avg.prob=rf.predicted)
  av.thr=optimal.thresholds(umbral)[3,2]

  # temporal collinearity shift between predictors
  messy.stck1=data.frame(totuneo@data)
  messy.stck2=terra::extract(x=stck2,y=vect(spatcoord,crs="+proj=longlat +datum=WGS84"),ID=F)
  vif.stck1=stats::cor(messy.stck1)
  vif.stck2=cor(messy.stck2)
  cor.present.stck1=cor(as.dist(vif.stck1),as.dist(vif.stck2))

  # VIF of predictors
  vif1=vif(messy.stck1)
  vif2=vif(messy.stck2)
  VIFvar.time1=vif1$Variables[which(vif1$VIF==max(vif1$VIF))]
  VIFvalue.time1=vif1$VIF[which(vif1$VIF==max(vif1$VIF))]
  VIFvar.time2=vif2$Variables[which(vif2$VIF==max(vif2$VIF))]
  VIFvalue.time2=vif2$VIF[which(vif2$VIF==max(vif2$VIF))]

  #MESS
  mess.p1=modEvA::MESS(V = totuneo@data[totuneo@pa==1,], P = totuneo@data[totuneo@pa==0,],verbosity=0)
  mess1=sum(ifelse(mess.p1$TOTAL<0,1,0))/nrow(mess.p1)

  mess.p2=modEvA::MESS(V = totuneo@data[totuneo@pa==1,], P = messy.stck2,verbosity=0)
  mess2=sum(ifelse(mess.p2$TOTAL<0,1,0))/nrow(mess.p2)

  #target area
  target1=crop(stck1, maskedbuf,mask=T)
  target2=crop(stck2, maskedbuf,mask=T)

  ## Predicted distribution in the target area for present and future scenarios
  hab.stck1=predict(randfo.model, data = target1)
  hab.prob1=hab.stck1
  hab.stck1[hab.stck1>av.thr]=1;hab.stck1[hab.stck1<1]=0

  hab.stck2=predict(randfo.model, data = target2)
  hab.prob2=hab.stck2
  hab.stck2[hab.stck2>av.thr]=1;hab.stck2[hab.stck2<1]=0

  ## Predicted area of occupancy (AOO) for each scenario
  aoo.stck1=round(global(cellSize(hab.stck1,unit="km")*hab.stck1, "sum",na.rm=T)$sum,1)
  aoo.stck2=round(global(cellSize(hab.stck2,unit="km")*hab.stck2, "sum",na.rm=T)$sum,1)

  ## Output summary
  outta=data.frame(pres=nrow(pr_dat), abs=nrow(bg_dat),buff,auc.test=randfo.auc.test,
                   VIFvar.time1,VIFvalue.time2,VIFvar.time2,VIFvalue.time2, cor.present.stck1,
                   mess1, mess2, aoo.t1=aoo.stck1,aoo.t2=aoo.stck2)
  outta=data.frame(t(outta))

  salida=list(outta,maskedbuf,hab.stck1,hab.stck2,hab.prob1,hab.prob2)
  names(salida)=c("summary","target.area","occupancy.time1","occupancy.time2",
                  "prob.time1","prob.time2")
  return(salida)
  options(warn=0)
}

