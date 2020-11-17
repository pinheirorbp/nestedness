plot2.nestmodels=function(x){
  ## parameters ##
  equi.model=x$parameters$equi.model
  prop.model=x$parameters$prop.model
  index=x$parameters$index
  binindex=x$parameters$binindex
  wprob=x$parameters$wprob
  AS=paste("samp",x$nsampled,sep="")
  calc.at=x$parameters$calc.at
  if(length(calc.at)>1){
    if(is.element(index,c("MD","SR","binSR"))){
      YLIM=c(min(c(unlist(x$summary.equi),unlist(x$summary.prop),x$observed),na.rm=T),max(c(unlist(x$summary.equi),unlist(x$summary.prop),x$observed),na.rm=T))
    }else{YLIM=c(0,1)}
    plot(NA,ylim=YLIM,xlim=c(min(calc.at),max(calc.at)),xlab="sampling",ylab=index)
    # observed
    points(x$nsampled,x$observed,pch=16,cex=1.5,col="#205501FF")
    # equiprobable
    if(equi.model){
      mean.equi=numeric()
      plus95.equi=numeric()
      minus95.equi=numeric()
      for (i in 1:length(calc.at)){
        mean.equi[i]=x$summary.equi[[i]]$mean.equi
        plus95.equi[i]=x$summary.equi[[i]]$plus95.equi
        minus95.equi[i]=x$summary.equi[[i]]$minus95.equi
      }
      polygon(x=c(calc.at,rev(calc.at)),y=c(plus95.equi,rev(minus95.equi)), border=F,col="#C65C7540")
      lines(calc.at,mean.equi,lwd=2,col="#C65C75FF")
    }
    #proportional
    if(prop.model){
      mean.prop=numeric()
      plus95.prop=numeric()
      minus95.prop=numeric()
      for (i in 1:length(calc.at)){
        mean.prop[i]=x$summary.prop[[i]]$mean.prop
        plus95.prop[i]=x$summary.prop[[i]]$plus95.prop
        minus95.prop[i]=x$summary.prop[[i]]$minus95.prop
      }
      polygon(x=c(calc.at,rev(calc.at)),y=c(plus95.prop,rev(minus95.prop)), border=F,col="#56519640")
      if(binindex&!wprob){lines(calc.at,mean.prop,lwd=2,col="#d7bd64FF")
      } else{lines(calc.at,mean.prop,lwd=2,col="#565196FF")}
    }
  }
}