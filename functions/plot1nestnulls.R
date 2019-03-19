plot1.nestnulls=function(x){
  ## parameters ##
  equi.null=x$parameters$equi.null
  prop.null=x$parameters$prop.null
  index=x$parameters$index
  binindex=x$parameters$binindex
  wprob=x$parameters$wprob
  AS=paste("samp",x$nsampled,sep="")
  ## plot 1: density plot ####
    if(equi.null){deneq=density(x$equinull[[AS]],na.rm=T)}else{deneq=list(y=NA)}
    if(prop.null){denprop=density(x$propnull[[AS]],na.rm=T)}else{denprop=list(y=NA)}
    if(is.element(index,c("MD","SR","binSR"))){
      XLIM=c(min(c(equinull[[AS]],propnull[[AS]],x$observed),na.rm=T),max(c(equinull[[AS]],propnull[[AS]],x$observed),na.rm=T))
    }else{XLIM=c(0,1)}
    YLIM=c(0,max(c(deneq$y,denprop$y),na.rm=T))
    plot(NA,ylim=YLIM,xlim=XLIM, lwd=2,main=NA, xlab=index, ylab="density")
    if(equi.null){lines(deneq$x,deneq$y,lwd=2,col="#C65C75FF")}
    if(prop.null){
      if(binindex&!wprob){
        lines(denprop$x,denprop$y,lwd=2,col="#d7bd64FF")}
      else{lines(denprop$x,denprop$y,lwd=2,col="#565196FF")}
    }
    abline(v = x$observed, lwd=2,col="#205501ff")
}