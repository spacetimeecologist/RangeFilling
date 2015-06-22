climateEnvelope2 <- function(test, train, p.abs, pres, envi){
  workers <- makeCluster(4)
  registerDoParallel(workers)
  clusterEvalQ(workers, c(library(dismo), library(rgdal), library(sp)))
  clusterExport(workers, varlist=c('envi', 'train', 'pres', 'test', 'p.abs')
                , envir=environment())
  
  bc <- llply(.data=c(1:length(pres)), 
              .fun=function(X){X=names(pres[X])
                               bioclim(envi, train[[X]])}, .parallel=T)
  names(bc) <- names(pres)
  print('Model created')
  clusterExport(workers, varlist=c('envi', 'pres', 'bc', 'test', 'p.abs'), envir=environment())
  bc.e <- llply(.data=c(1:length(pres)), 
                .fun=function(X){X=names(pres[X])
                                 evaluate(p=test[[X]], a=p.abs[[X]], model=bc[[X]], x=envi)}, 
                .parallel=T)
  names(bc.e) <- names(pres)
  print('Model evaluated')
  bc.tr <- llply(.data=c(1:length(pres)), 
                 .fun=function(X){X=names(pres[X])
                                  threshold(bc.e[[X]], stat='spec_sens')}, .progress='text')
  names(bc.tr) <- names(pres)
  print('Threshold determined')
  
  bc.pred <- llply(.data=c(1:length(pres)),
                   .fun=function(X){X=names(pres[X])
                                    predict(bc[[X]], envi)},.parallel=T)
  names(bc.pred) <- names(pres)
  print('Predicted habitat')
  
  bc.prab <- llply(.data=c(1:length(pres)),
                   .fun=function(X){X=names(pres[X])
                                    bc.pred[[X]] > bc.tr[[X]]}, .progress='text')
  names(bc.prab) <- names(bc.pred)
  print('Threshold applied')
  return(bc.prab)
  stopCluster(workers)
}
