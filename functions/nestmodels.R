nestmodels=function(M, index="wnoda",equi.model=T, prop.model=T, n.model=1000, calc.at=NULL, print.at.each=NULL, wprob=F,wsamp=F, density.plot=T,sampling.plot=F){
  # functions ####
  #function for one tailed one sample z test
  z1samp= function (x, mean, sd){
    z= (x-mean)/sd
    p= pnorm(-abs(z))
    out= list(ZScore=z,Pvalue=p)
    return(out)
  }
  if(index=="SR"|index=="MD"|index=="binSR"){
    # The code used to calculate Spectral Radius and Manhattan Distance was adapted from the FALCON package developed by Stephen Beckett. 
    # FALCON (Framework of Adaptive ensembLes for the Comparison Of Nestedness) package: https://github.com/sjbeckett/FALCON
    
    SPECTRAL_RADIUS <- function(MATRIX){
      #The spectral radius method is described in Staniczenko et al., 2013.
      
      #PPA Staniczenko, JC Kopp, S Allesina. 2013.
      #The ghost of nestedness in ecological networks
      #Nature Communications 4(1391). (http://dx.doi.org/10.1038/ncomms2422)
      
      
      r <- dim(MATRIX)[1]
      c <- dim(MATRIX)[2]
      
      tot <- r+c
      
      #Build adjacency matrix
      ADJ <- matrix(0,tot,tot)
      ADJ[1:r, (r + 1):tot] <- MATRIX
      ADJ <- ADJ+t(ADJ)
      
      #return the real eigenvalue with largest magnitude
      MEASURE <- Re(eigen(ADJ)$value[1])
      
      return(MEASURE)
    }
    # Part of the FALCON (Framework of Adaptive ensembLes for the Comparison Of
    # Nestedness) package: https://github.com/sjbeckett/FALCON
    # Last updated: 10th March 2014
    
    MANHATTAN_DISTANCE <- function(MATRIX) {
      #The Tau-temperature method based on the Manhattan Distance is described in
      #Corso & Britton, 2012. This function returns the Manhattan distance of a
      #matrix. The Tau-temperature is calculated as the Manhattan distance of the
      #test matrix divided by the expected Manhattan distance as calculated from
      #the mean of an ensemble. In FALCON the Tau-temperature is shown as the
      #NormalisedTemperature in the outputs.
      
      #G Corso, NF Britton. 2012.
      #Nestedness and ??-temperature in ecological networks
      #Ecological Complexity 11: 137 - 143. (http://dx.doi.org/10.1016/j.ecocom.2012.05.003)
      
      
      
      rows <- dim(MATRIX)[1]
      cols <- dim(MATRIX)[2]
      
      #make sure binary.
      MATRIX <- 1*(MATRIX>0)
      
      #Find distance matrix for row and column indices
      Rnumbers <- rep(1:rows,cols)
      distancematrix_R <- matrix(Rnumbers,nrow=rows)
      
      #Find the sum of the distances for each index from the input matrix.
      Cnumbers <- rep(1:cols, each=rows)
      distancematrix_C <- matrix(Cnumbers,nrow=rows)
      
      
      dx <- sum(sum(distancematrix_R*MATRIX)) ;
      dy <- sum(sum(distancematrix_C*MATRIX)) ;
      
      DISTANCE <- dx + dy;
      
      
      return(DISTANCE)
    }
    ####
    
    }
  if(index=="nodf"|index=="wnodf"){library(vegan)}
  if(index=="wnoda"){ library (bipartite)}
  if(index=="temperature"){library(bipartite)}
  if(index=="betaNES"|index=="betaNES2"){library(betapart)}
  # binary index (F / T) ####
  if(is.element(index,c("nodf","temperature","MD","binSR", "betaNES", "betaNES2"))){
    binindex=T}else{binindex=F}
  # interaction matrix ####
  #removing empty rows / columns
  M=M[rowSums(M, na.rm=T)>0,colSums(M, na.rm=T)>0]
  net=M
  if(binindex){
    net[M>0]=1
    if(!wsamp){calc.at=sum(net, na.rm=T)}
  }
  # probability matrix ####
  if(prop.model&!wprob&binindex){
    PM=matrix(NA,nrow(net),ncol(net))
    RS=rowSums(net, na.rm=T)/ncol(net)
    CS=colSums(net,na.rm=T)/nrow(net)
    for (r in 1:nrow(net)){
      for (c in 1:ncol(net)){
        PM[r,c]=(RS[r]+CS[c])/2
      }}
  }else if(prop.model){
    PM=matrix(NA,nrow(M),ncol(M))
    RS=rowSums(M, na.rm=T)/sum(M, na.rm=T)
    CS=colSums(M, na.rm=T)/sum(M, na.rm=T)
    for (r in 1:nrow(M)){
      for (c in 1:ncol(M)){
        PM[r,c]=RS[r]*CS[c]
      }}}
  # others ####
  if (is.null(calc.at)){calc.at=sum(M, na.rm=T)}
  if ((!binindex|wsamp)&!is.element(sum(M, na.rm=T),calc.at)){
    calc.at=c(calc.at,sum(M, na.rm=T))
    calc.at=sort(calc.at)
  }
  # observed values ####
  if(index=="wnodf"){
    observed=nestednodf(net,weighted = T)[[3]][[3]]/100
  }else if(index=="wnoda"){
    observed=nest.smdm(net,weighted = T,decreasing = "abund")$WNODAmatrix/100
  }else if(index=="SR"|index=="binSR"){
    observed=SPECTRAL_RADIUS(net)
  }else if(index=="MD"){
    observed=MANHATTAN_DISTANCE(net)
  }else if(index=="nodf"){
    observed=nestednodf(net)[[3]][[3]]/100
  }else if(index=="temperature"){
    observed=bipartite::networklevel(net, index="nestedness")/100
  }else if(index=="betaNES"){
    observed=beta.multi(net)$beta.SNE
  }else if(index=="betaNES2"){
    bm=beta.multi(net)
    observed=bm$beta.SNE/bm$beta.SOR
  }
  # loop for randomized models ####
  equimodel=list()
  propmodel=list()
  for (N in 1:n.model){
    # equiprobable ####
    if(equi.model){
      if(binindex&!wsamp){INT=sample(1:length(net),size = calc.at,replace = F)
      }else{INT=sample(1:length(M),size = max(calc.at),replace = T)}
      for(CA in calc.at){
        TAB.INT= table(INT[1:CA])
        NM=matrix(0,nrow(net),ncol(net))
        NM[as.numeric(rownames(TAB.INT))]=as.numeric(TAB.INT)
        NM=NM[rowSums(NM)>0,colSums(NM)>0]
        if(binindex&wsamp){NM[NM>0]=1}
        samp=paste("samp",CA,sep="")
        # indices
        if(index=="wnodf"){
          equimodel[[samp]][N]=nestednodf(NM,weighted = T)[[3]][[3]]/100
        }else if(index=="wnoda"){
          equimodel[[samp]][N]=nest.smdm(NM,weighted = T,decreasing = "abund")$WNODAmatrix/100
        }else if(index=="SR"|index=="binSR"){
          equimodel[[samp]][N]=SPECTRAL_RADIUS(NM)
        }else if(index=="MD"){
          equimodel[[samp]][N]=MANHATTAN_DISTANCE(NM)
        }else if(index=="nodf"){
          equimodel[[samp]][N]=nestednodf(NM)[[3]][[3]]/100
        }else if(index=="temperature"){
          equimodel[[samp]][N]=bipartite::networklevel(NM, index="nestedness")/100
        }else if(index=="betaNES"){
          equimodel[[samp]][N]=beta.multi(NM)$beta.SNE
        }else if(index=="betaNES2"){
          bm=beta.multi(NM)
          equimodel[[samp]][N]=bm$beta.SNE/bm$beta.SOR
        }
      }
    }else{equimodel=NULL}
    # proportional ####
    if(prop.model){
      if(binindex&!wsamp){INT=sample(1:length(net),size = calc.at,replace = F,prob = as.numeric(PM))
      }else{INT=sample(1:length(M),size = max(calc.at),replace = T,prob = as.numeric(PM))}
      for(CA in calc.at){
        TAB.INT= table(INT[1:CA])
        NM=matrix(0,nrow(net),ncol(net))
        NM[as.numeric(rownames(TAB.INT))]=as.numeric(TAB.INT)
        NM=NM[rowSums(NM)>0,colSums(NM)>0]
        if(binindex&wsamp){NM[NM>0]=1}
        samp=paste("samp",CA,sep="")
        # indices
        if(index=="wnodf"){
          propmodel[[samp]][N]=nestednodf(NM,weighted = T)[[3]][[3]]/100
        }else if(index=="wnoda"){
          propmodel[[samp]][N]=nest.smdm(NM,weighted = T,decreasing = "abund")$WNODAmatrix/100
        }else if(index=="SR"|index=="binSR"){
          propmodel[[samp]][N]=SPECTRAL_RADIUS(NM)
        }else if(index=="MD"){
          propmodel[[samp]][N]=MANHATTAN_DISTANCE(NM)
        }else if(index=="nodf"){
          propmodel[[samp]][N]=nestednodf(NM)[[3]][[3]]/100
        }else if(index=="temperature"){
          propmodel[[samp]][N]=bipartite::networklevel(NM, index="nestedness")/100
        }else if(index=="betaNES"){
          propmodel[[samp]][N]=beta.multi(NM)$beta.SNE
        }else if(index=="betaNES2"){
          bm=beta.multi(NM)
          propmodel[[samp]][N]=bm$beta.SNE/bm$beta.SOR
        }
      }
    }else{propmodel=NULL}
    if(!is.null(print.at.each)){
      if(N%%as.integer(print.at.each)==0){print(N)}}
  }
  ## summaries ####
  # for the observed matrix
  if(binindex&!wsamp){
    as0=sum(net, na.rm=T)
  }else{as0=sum(M, na.rm=T)}
  AS=paste("samp",as0,sep="")
  significance=list()
  if(equi.model){
    PZ=z1samp(observed,mean(equimodel[[AS]],na.rm=T),sd(equimodel[[AS]],na.rm=T))
    significance$pvalue.equi=PZ$Pvalue
    significance$Zscore.equi=PZ$ZScore
    }
  if(prop.model){
    PZ=z1samp(observed,mean(propmodel[[AS]],na.rm=T),sd(propmodel[[AS]],na.rm=T))
    significance$pvalue.prop=PZ$Pvalue
    significance$Zscore.prop=PZ$ZScore
  }
  # for each sampling level
  summary.equi=list()
  if(equi.model){
    for (i in 1:length(calc.at)){
      SAM=names(equimodel)[i]
      summary.equi[[SAM]]$nnull=sum(!is.na(equimodel[[i]]))
      summary.equi[[SAM]]$mean.equi=mean(equimodel[[i]], na.rm = T)
      summary.equi[[SAM]]$sd.equi=sd(equimodel[[i]], na.rm = T)
      summary.equi[[SAM]]$plus95.equi=quantile(equimodel[[i]], na.rm = T, probs = 0.975)[[1]]
      summary.equi[[SAM]]$minus95.equi=quantile(equimodel[[i]], na.rm = T, probs = 0.025)[[1]]
    }
  }
  summary.prop=list()
  if(prop.model){
    for (i in 1:length(calc.at)){
      SAM=names(propmodel)[i]
      summary.prop[[SAM]]$mean.prop=mean(propmodel[[i]], na.rm = T)
      summary.prop[[SAM]]$sd.prop=sd(propmodel[[i]], na.rm = T)
      summary.prop[[SAM]]$plus95.prop=quantile(propmodel[[i]], na.rm = T, probs = 0.975)[[1]]
      summary.prop[[SAM]]$minus95.prop=quantile(propmodel[[i]], na.rm = T, probs = 0.025)[[1]]
    }
  }
  ## plot 1: density plot ####
  if(density.plot){
    if(equi.model){deneq=density(equimodel[[AS]],na.rm=T)}else{deneq=list(y=NA)}
    if(prop.model){denprop=density(propmodel[[AS]],na.rm=T)}else{denprop=list(y=NA)}
    if(is.element(index,c("MD","SR","binSR"))){
      XLIM=c(min(c(equimodel[[AS]],propmodel[[AS]],observed),na.rm=T),max(c(equimodel[[AS]],propmodel[[AS]],observed),na.rm=T))
    }else{XLIM=c(0,1)}
    YLIM=c(0,max(c(deneq$y,denprop$y),na.rm=T))
    plot(NA,ylim=YLIM,xlim=XLIM, lwd=2,main=NA, xlab=index, ylab="density")
    if(equi.model){lines(deneq$x,deneq$y,lwd=2,col="#C65C75FF")}
    if(prop.model){
      if(binindex&!wprob){
        lines(denprop$x,denprop$y,lwd=2,col="#d7bd64FF")}
      else{lines(denprop$x,denprop$y,lwd=2,col="#565196FF")}
      }
    abline(v = observed, lwd=2,col="#205501ff")
  }
  ## plot 2: sampling plot ####
  if(sampling.plot&length(calc.at)>1){
    if(is.element(index,c("MD","SR","binSR"))){
      YLIM=c(min(c(unlist(summary.equi),unlist(summary.prop),observed),na.rm=T),max(c(unlist(summary.equi),unlist(summary.prop),observed),na.rm=T))
    }else{YLIM=c(0,1)}
    plot(NA,ylim=YLIM,xlim=c(min(calc.at),max(calc.at)),xlab="sampling",ylab=index)
    # observed
    points(sum(M),observed,pch=16,cex=1.5,col="#205501FF")
    # equiprobable
    if(equi.model){
      mean.equi=numeric()
      plus95.equi=numeric()
      minus95.equi=numeric()
      for (i in 1:length(calc.at)){
        mean.equi[i]=summary.equi[[i]]$mean.equi
        plus95.equi[i]=summary.equi[[i]]$plus95.equi
        minus95.equi[i]=summary.equi[[i]]$minus95.equi
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
        mean.prop[i]=summary.prop[[i]]$mean.prop
        plus95.prop[i]=summary.prop[[i]]$plus95.prop
        minus95.prop[i]=summary.prop[[i]]$minus95.prop
      }
        polygon(x=c(calc.at,rev(calc.at)),y=c(plus95.prop,rev(minus95.prop)), border=F,col="#56519640")
        if(binindex&!wprob){lines(calc.at,mean.prop,lwd=2,col="#d7bd64FF")
        } else{lines(calc.at,mean.prop,lwd=2,col="#565196FF")}
    }
  }
  ## output ####
  parameters=list(index=index,equi.model=equi.model, prop.model=prop.model, n.model=n.model, calc.at=calc.at, print.at.each=print.at.each, wprob=wprob, wsamp=wsamp, binindex=binindex)
  return(list(observed=observed,significance=significance,summary.equi=summary.equi,summary.prop=summary.prop,equimodel=equimodel, propmodel=propmodel, nsampled=as0, parameters=parameters))
}
