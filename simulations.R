### Authors ####

#Rafael Barros Pereira Pinheiro 
#Carsten F. Dormann 
#Gabriel Moreira F?lix Ferreira
#Marco Aurelio Ribeiro Mello 

# Contact: rafael-bpp@hotmail.com

### Information ####

# The code used to calculate Spectral Radius and Manhattan Distance was adapted from the FALCON package developed by Stephen Beckett. 
# FALCON (Framework of Adaptive ensembLes for the Comparison Of Nestedness) package: https://github.com/sjbeckett/FALCON

## We used parallel processing to produce the simulations.
## The easiest way to run this code without paralelizing is by registring only 1 core [registerDoParallel(cores=1)], but be prepared because it may take a LONG time.
## On the other hand, if you are interested in applying these methods to a given system, the function nestnull (file= nestnull.R) is probable the best choice, not this script.

## This script run the simulations. 
## To produce the plots, use the R Markdown files: Appendix_S2.Rmd and Appendix_S3.Rmd

#### Packages / functions ####
library(sads)
library(bipartite)
library(betapart)
library(doParallel)
source("functions/nest.smdm.R")
source("functions/SR_MD_Becket.R")

### Probability Matrices ####
## 20x20
M=matrix(NA,10000,20)
for (i in 1:10000){
  M[i,]=sort(rpoilog(20,mu = 1,sig = 0.6),decreasing = T)
}
ProbLognorm=colMeans(M)
ProbLognorm=ProbLognorm/sum(ProbLognorm)
ProbLinDec=20:1 
ProbLinDec=ProbLinDec/sum(ProbLinDec)
rm(i,M)
PMLogn= matrix(NA,20,20)
for (i in 1:20){
  for (j in 1:20){
    PMLogn[i,j]=ProbLognorm[i]*ProbLognorm[j]
  }
}
PMLin= matrix(NA,20,20)
for (i in 1:20){
  for (j in 1:20){
    PMLin[i,j]=ProbLinDec[i]*ProbLinDec[j]
  }
}
PMEqui= matrix(0.0025,20,20)
rm(i,j)
PM20=list(Logn=PMLogn,Lin=PMLin,Equi=PMEqui)
## 10x10
M=matrix(NA,10000,10)
for (i in 1:10000){
  M[i,]=sort(rpoilog(10,mu = 1,sig = 0.6),decreasing = T)
}
ProbLognorm=colMeans(M)
ProbLognorm=ProbLognorm/sum(ProbLognorm)
ProbLinDec=10:1 
ProbLinDec=ProbLinDec/sum(ProbLinDec)
rm(i,M)
PMLogn= matrix(NA,10,10)
for (i in 1:10){
  for (j in 1:10){
    PMLogn[i,j]=ProbLognorm[i]*ProbLognorm[j]
  }
}
PMLin= matrix(NA,10,10)
for (i in 1:10){
  for (j in 1:10){
    PMLin[i,j]=ProbLinDec[i]*ProbLinDec[j]
  }
}
PMEqui= matrix(0.01,10,10)
rm(i,j)
PM10=list(Logn=PMLogn,Lin=PMLin,Equi=PMEqui)
## 5x5
M=matrix(NA,10000,5)
for (i in 1:10000){
  M[i,]=sort(rpoilog(5,mu = 1,sig = 0.6),decreasing = T)
}
ProbLognorm=colMeans(M)
ProbLognorm=ProbLognorm/sum(ProbLognorm)
ProbLinDec=5:1 
ProbLinDec=ProbLinDec/sum(ProbLinDec)
rm(i,M)
PMLogn= matrix(NA,5,5)
for (i in 1:5){
  for (j in 1:5){
    PMLogn[i,j]=ProbLognorm[i]*ProbLognorm[j]
  }
}
PMLin= matrix(NA,5,5)
for (i in 1:5){
  for (j in 1:5){
    PMLin[i,j]=ProbLinDec[i]*ProbLinDec[j]
  }
}
PMEqui= matrix(0.04,5,5)
rm(i,j)
PM5=list(Logn=PMLogn,Lin=PMLin,Equi=PMEqui)
#
PM=list(mod1=list(PM5=PM5,PM10=PM10,PM20=PM20))
save(PM, file="files/PM.RData")
### Simulations 1 ####
#Defined by the total weights
load("files/PM.RData")
## general information
matprob= rep(c("Logn","Lin","Equi"),each=80,times=3)
sampling= rep(c(50,100,200,400,800,1600,3200,6400),each=10,times=9)
size= rep(c(20,10,5),each=240)
#
LOG=character()
save(LOG, file="log.RData")
registerDoParallel(cores=10)
simulations=foreach(N=1:720,.packages =c("bipartite","betapart","vegan"))%dopar%{
  ### list for results ###
  results=list()
  ## Probability matrix ##
  PROB=as.numeric(PM[[paste("PM",size[[N]],sep="")]][[matprob[[N]]]])
  for (I in 1:1000){
  ## Simulated matrix ##
  INT=sample(1:(size[[N]]^2),size =sampling[[N]], prob = PROB,replace = T)
  TAB.INT= table(INT)
  MAT=matrix(0,size[[N]],size[[N]])
  MAT[as.numeric(rownames(TAB.INT))]=as.numeric(TAB.INT)
  net=MAT
  net=net[rowSums(net)>0,colSums(net)>0]
  results$net[[I]]=net #0
  netBin=net
  netBin[net>=1]=1
  ## indices ##
  results$nrow[I]=nrow(net)#1
  results$ncol[I]=ncol(net)#2
  results$connectance[I]=bipartite::networklevel(netBin,index = "connectance")#3
  # binary
  results$betasim[I]=beta.multi(netBin)$beta.SIM#4
  results$betasne[I]=beta.multi(netBin)$beta.SNE#5
  results$betasor[I]=beta.multi(netBin)$beta.SOR#6
  results$temp[I]=1-bipartite::networklevel(netBin, index="nestedness")/100#7
  results$nodf[I]=vegan::nestednodf(netBin)$statistic[[3]]/100#8
  # weighted
  results$wnodf[I]=nest.smdm(x = net,weights = T,decreasing = "fill",sort = T)$WNODFmatrix[1]/100#12
  results$wnoda[I]=nest.smdm(x = net,weights = T,decreasing = "abund",sort = T)$WNODAmatrix[1]/100#13
  results$SR[I]=SPECTRAL_RADIUS(net)#16
  }
  ## Saving results ##
  fileID=paste("simulations/part1/simsID",N,".RData",sep = "")
  save(results, file = fileID)
  ##### LOG ####
  load("log.RData")
  LOG=c(LOG, paste(N,Sys.time()))
  save(LOG, file="log.RData")
  rm(LOG)
  ###
  output="done"
  output
}
rm(sampling,matprob, PM, size, simulations)
### Simulations 2 ####
# Defined by connectance (only binary indices) ###
load("files/PM.RData")
# # general information
matprob= rep(c("Logn","Lin","Equi"),each=80,times=3)
connectance= rep(c(.3,.5,.7,.9,.92,.94,.96,.98),each=10,times=9)
size= rep(c(20,10,5),each=240)
#
registerDoParallel(cores=10)
simulations=foreach(N=1:720,.packages =c("bipartite","betapart","vegan"))%dopar%{
  ### list for results ###
  results=list()
  PROB=as.numeric(PM[[paste("PM",size[[N]],sep="")]][[matprob[[N]]]])
  for (I in 1:1000){
    ## Simulated matrix ##
    ndiscarded=0
    SIZE_CONSERVED=F
    while(!SIZE_CONSERVED){
      INT=sample(1:(size[[N]]^2),size =connectance[[N]]*(size[[N]]^2), prob = PROB,replace = F)
      TAB.INT= table(INT)
      MAT=matrix(0,size[[N]],size[[N]])
      MAT[as.numeric(rownames(TAB.INT))]=as.numeric(TAB.INT)
      net=MAT
      if(sum(rowSums(net)>0)==size[[N]]&sum(colSums(net)>0)==size[[N]]){SIZE_CONSERVED=T}
      ndiscarded=ndiscarded+1
    }
  # network
  results$net[[I]]=net #0
  results$ndiscarded[I]=ndiscarded#0.1
  ## indices ##
  results$nrow[I]=nrow(net)#1
  results$ncol[I]=ncol(net)#2
  results$connectance[I]=bipartite::networklevel(net,index = "connectance")#3
  # binary
  results$betasim[I]=beta.multi(net)$beta.SIM#4
  results$betasne[I]=beta.multi(net)$beta.SNE#5
  results$betasor[I]=beta.multi(net)$beta.SOR#6
  results$temp[I]=1-bipartite::networklevel(net, index="nestedness")/100#7
  results$nodf[I]=vegan::nestednodf(net)$statistic[[3]]/100#8
  results$MD[I]=MANHATTAN_DISTANCE(net)#10
  results$binSR[I]=SPECTRAL_RADIUS(net)#11
  }
  ## Saving results ##
  fileID=paste("simulations/part2/simsID",N,".RData",sep = "")
  save(results, file = fileID)
  ###
  output="done"
  output
}
## removing objects ##
rm(connectance,matprob, PM, size, simulations)
### Simulations 3 ####
# Proportional null model defined by node degrees in simulations 2#
# general information
connectance= rep(c(.3,.5,.7,.9,.92,.94,.96,.98),each=10,times=9)
size= rep(c(20,10,5),each=240)
#
registerDoParallel(cores=10)
simulations=foreach(N=1:720,.packages =c("bipartite","betapart","vegan"))%dopar%{
  # loading binary matrices (simulations 2) #
  loadfileID=paste("simulations/part2/simsID",N,".RData",sep = "")
  load(loadfileID)
  NETS=results$net
  rm(results)
  ### list for results ###
  results=list()
  for (I in 1:1000){
    # Probability matrix based on node degrees #
    RS=rowSums(NETS[[I]])/size[N]
    CS=colSums(NETS[[I]])/size[N]
    PM= matrix(NA,size[N],size[N])
    for (i in 1:size[N]){
      for (j in 1:size[N]){
        PM[i,j]=(RS[i]+CS[j])/2
      }
    }
    PROB=as.numeric(PM)
    ## Simulated matrix ##
    ndiscarded=0
    SIZE_CONSERVED=F
    while(!SIZE_CONSERVED){
      INT= rbinom((size[[N]]^2),1,prob = PROB)
      MAT=matrix(INT,size[[N]],size[[N]])
      net=MAT
      if(sum(rowSums(net)>0)==size[[N]]&sum(colSums(net)>0)==size[[N]]){SIZE_CONSERVED=T}
      ndiscarded=ndiscarded+1
    }
    # network
    results$net[[I]]=net #0
    results$ndiscarded[I]=ndiscarded#0.1
    ## indices ##
    results$nrow[I]=nrow(net)#1
    results$ncol[I]=ncol(net)#2
    results$connectance[I]=bipartite::networklevel(net,index = "connectance")#3
    # binary
    results$betasim[I]=beta.multi(net)$beta.SIM#4
    results$betasne[I]=beta.multi(net)$beta.SNE#5
    results$betasor[I]=beta.multi(net)$beta.SOR#6
    results$temp[I]=1-bipartite::networklevel(net, index="nestedness")/100#7
    results$nodf[I]=vegan::nestednodf(net)$statistic[[3]]/100#8
    results$MD[I]=MANHATTAN_DISTANCE(net)#10
    results$binSR[I]=SPECTRAL_RADIUS(net)#11
  }
  ## Saving results ##
  fileID=paste("simulations/part3/simsID",N,".RData",sep = "")
  save(results, file = fileID)
  ###
  output="done"
  output
}
## removing objects ##
rm(connectance,matprob, PM, size, simulations)

### Summary - Simulations 1 ####
## general information
SIM1= list()
SIM1$PM=rep(c("Logn","Lin","Equi"),each=8,times=3)
SIM1$sampling=rep(c(50,100,200,400,800,1600,3200,6400),times=9)
SIM1$size= rep(c(20,10,5),each=24)
 for (i in 1:72){
  # Indices for the 10000 matrices with each setup
  for (j in 1:10){
    SIMID=paste("simulations/part1/simsID",((i-1)*10)+j,".RData",sep="")
    load(SIMID)
    if(j==1){INDICES=results[-1]}else{
      for (k in 1:13){
        INDICES[[k]][(((j-1)*1000)+1):(j*1000)]=results[[k+1]]
      }}}
  rm(results,j,SIMID,k)
  # Summaries for the values for each setup
  for(k in 1:14){
    if(k==14){
      # proportional betanes
      options(warn=-1)
      INDICE=as.numeric(INDICES$betasne/INDICES$betasor)
      options(warn=0)
      INDNAME="betasne_prop"
    }else{
      options(warn=-1)
      INDICE=as.numeric(INDICES[[k]])
      options(warn=0)
      INDNAME=names(INDICES)[k]}
    SIM1[[paste("na.",INDNAME,sep = "")]][i]=sum(is.na(INDICE))
    SIM1[[paste("mean.",INDNAME,sep = "")]][i]=mean(INDICE,na.rm=T)
    SIM1[[paste("sd.",INDNAME,sep = "")]][i]=sd(INDICE,na.rm=T)
    SIM1[[paste("median.",INDNAME,sep = "")]][i]=median(INDICE,na.rm=T)
    SIM1[[paste("q025.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.025,na.rm = T)[[1]]
    SIM1[[paste("q05.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.05,na.rm = T)[[1]]
    SIM1[[paste("q1.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.1,na.rm = T)[[1]]
    SIM1[[paste("q9.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.9,na.rm = T)[[1]]
    SIM1[[paste("q95.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.95,na.rm = T)[[1]]
    SIM1[[paste("q975.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.975,na.rm = T)[[1]]
    SIM1[[paste("max.",INDNAME,sep = "")]][i]=max(INDICE,na.rm = T)
    SIM1[[paste("min.",INDNAME,sep = "")]][i]=min(INDICE,na.rm = T)
  }
print(i)
}
save(SIM1, file="results/simulations1.RData")
rm(i, INDICES, INDNAME,k, INDICE,SIM1)
### Summary - Simulations 2 ####
## general information
SIM2=list()
SIM2$PM=rep(c("Logn","Lin","Equi"),each=8,times=3)
SIM2$connectance= rep(c(.3,.5,.7,.9,.92,.94,.96,.98),times=9)
SIM2$size= rep(c(20,10,5),each=24)
for (i in 1:72){
  # Indices for the 10000 matrices with each setup
  for (j in 1:10){
    SIMID=paste("simulations/part2/simsID",((i-1)*10)+j,".RData",sep="")
    load(SIMID)
    if(j==1){INDICES=results[-1]}else{
      for (k in 1:11){
        INDICES[[k]][(((j-1)*1000)+1):(j*1000)]=results[[k+1]]
      }}}
  rm(results,j,SIMID,k)
  # Summaries for the values for each setup
  for(k in 1:12){
    if(k==12){
      # proportional betanes
      options(warn=-1)
      INDICE=as.numeric(INDICES$betasne/INDICES$betasor)
      options(warn=0)
      INDNAME="betasne_prop"
    }else{
      options(warn=-1)
      INDICE=as.numeric(INDICES[[k]])
      options(warn=0)
      INDNAME=names(INDICES)[k]}
    SIM2[[paste("na.",INDNAME,sep = "")]][i]=sum(is.na(INDICE))
    SIM2[[paste("mean.",INDNAME,sep = "")]][i]=mean(INDICE,na.rm=T)
    SIM2[[paste("sd.",INDNAME,sep = "")]][i]=sd(INDICE,na.rm=T)
    SIM2[[paste("median.",INDNAME,sep = "")]][i]=median(INDICE,na.rm=T)
    SIM2[[paste("q025.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.025,na.rm = T)[[1]]
    SIM2[[paste("q05.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.05,na.rm = T)[[1]]
    SIM2[[paste("q1.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.1,na.rm = T)[[1]]
    SIM2[[paste("q9.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.9,na.rm = T)[[1]]
    SIM2[[paste("q95.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.95,na.rm = T)[[1]]
    SIM2[[paste("q975.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.975,na.rm = T)[[1]]
    SIM2[[paste("max.",INDNAME,sep = "")]][i]=max(INDICE,na.rm = T)
    SIM2[[paste("min.",INDNAME,sep = "")]][i]=min(INDICE,na.rm = T)
  }
  print(i)
}
save(SIM2, file="results/simulations2.RData")
rm(i, INDICES, INDNAME,k, INDICE,SIM2)
### Summary - Simulations 3 ####
## general information
SIM3=list()
SIM3$PM=rep(c("Logn","Lin","Equi"),each=8,times=3)
SIM3$connectance= rep(c(.3,.5,.7,.9,.92,.94,.96,.98),times=9)
SIM3$size= rep(c(20,10,5),each=24)
for (i in 1:72){
  # Indices for the 10000 matrices with each setup
  for (j in 1:10){
    SIMID=paste("simulations/part3/simsID",((i-1)*10)+j,".RData",sep="")
    load(SIMID)
    if(j==1){INDICES=results[-1]}else{
      for (k in 1:11){
        INDICES[[k]][(((j-1)*1000)+1):(j*1000)]=results[[k+1]]
      }}}
  rm(results,j,SIMID,k)
  # Summaries for the values for each setup
  for(k in 1:12){
    if(k==12){
      # proportional betanes
      options(warn=-1)
      INDICE=as.numeric(INDICES$betasne/INDICES$betasor)
      options(warn=0)
      INDNAME="betasne_prop"
    }else{
      options(warn=-1)
      INDICE=as.numeric(INDICES[[k]])
      options(warn=0)
      INDNAME=names(INDICES)[k]}
    SIM3[[paste("na.",INDNAME,sep = "")]][i]=sum(is.na(INDICE))
    SIM3[[paste("mean.",INDNAME,sep = "")]][i]=mean(INDICE,na.rm=T)
    SIM3[[paste("sd.",INDNAME,sep = "")]][i]=sd(INDICE,na.rm=T)
    SIM3[[paste("median.",INDNAME,sep = "")]][i]=median(INDICE,na.rm=T)
    SIM3[[paste("q025.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.025,na.rm = T)[[1]]
    SIM3[[paste("q05.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.05,na.rm = T)[[1]]
    SIM3[[paste("q1.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.1,na.rm = T)[[1]]
    SIM3[[paste("q9.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.9,na.rm = T)[[1]]
    SIM3[[paste("q95.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.95,na.rm = T)[[1]]
    SIM3[[paste("q975.",INDNAME,sep = "")]][i]=quantile(INDICE,probs =.975,na.rm = T)[[1]]
    SIM3[[paste("max.",INDNAME,sep = "")]][i]=max(INDICE,na.rm = T)
    SIM3[[paste("min.",INDNAME,sep = "")]][i]=min(INDICE,na.rm = T)
  }
  print(i)
}
save(SIM3, file="results/simulations3.RData")
rm(i, INDICES, INDNAME,k, INDICE,SIM3)