### Colors ####
col1="#7DBC58FF"
col2="#D7BD64FF"
col3="#C65C75FF"
col4="#565196FF"
col6="#614b01FF"
col5="#205501FF"
col7="#590115FF"
col8="#0e0a44FF"

### FIGURE 4 ####
load("results/simulations1.RData")
#define the index to plot
#general info: nrow,ncol,connectance
#binary: betasim,betasne,betasor,temp,nodf,ORBoverlap,ORBdeviance,ORBnulldeviance,MD,binSR
#quantitative: wnodf,wnoda,ORLoverlap,ORLdeviance,ORLnulldeviance,ORPoverlap,ORPdeviance,ORPnulldeviance,SR
indice="betasne" 
#define the summaries information to plot
#mean,median,sd,min,max,q025,q05,q1,q9,q95,q975
summ=c("mean","q025","q975")
# define the size of the matrices
S=10
################# #
N=paste(summ,".",indice,sep = "")
SL=SIM1$size==S
plot(0,xlim = c(0.5,8),ylim=c(0,1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',pch=NA, xaxs="i", yaxs="i")
axis(side = 1,lwd=1,at = 1:8,labels = c("50","100","200","400","800","1600","3200","6400"))
axis(side=2,lwd=1)
for (i in 1:length(N)){
  LTY=c(1,2,2)[i]
  lines(SIM1[[N[i]]][SIM1$PM=="Logn"&SL],col=col1,lty=LTY)
  lines(SIM1[[N[i]]][SIM1$PM=="Lin"&SL],col=col4,lty=LTY)
  lines(SIM1[[N[i]]][SIM1$PM=="Equi"&SL],col=col3,lty=LTY)
}
############ SR
indice="SR" 
######## #
#define the summaries information to plot
#mean,median,sd,min,max,q025,q05,q1,q9,q95,q975
summ=c("mean","q025","q975")
# define the size of the matrices
S=10
######## #
N=paste(summ,".",indice,sep = "")
SL=SIM1$size==S
plot(0,xlim = c(0.5,8),ylim=c(0.075,0.25), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',pch=NA, xaxs="i", yaxs="i")
axis(side = 1,lwd=1,at = 1:8,labels = c("50","100","200","400","800","1600","3200","6400"))
axis(side=2,lwd=1)
for (i in 1:length(N)){
  LTY=c(1,2,2)[i]
  lines(SIM1[[N[i]]][SIM1$PM=="Logn"&SL]/c(50,100,200,400,800,1600,3200,6400),col=col1,lty=LTY)
  lines(SIM1[[N[i]]][SIM1$PM=="Lin"&SL]/c(50,100,200,400,800,1600,3200,6400),col=col4,lty=LTY)
  lines(SIM1[[N[i]]][SIM1$PM=="Equi"&SL]/c(50,100,200,400,800,1600,3200,6400),col=col3,lty=LTY)
}


### FIGURE 5 ####
load("results/simulations2.RData")
load("results/simulations3.RData")
#define the index to plot
#general info: nrow,ncol
#binary: betasim,betasne,betasor,betasne_prop,temp,nodf,ORBoverlap,ORBdeviance,ORBnulldeviance,MD,binSR
indice="betasne"
#define the summaries information to plot
#mean,median,sd,min,max,q025,q05,q1,q9,q95,q975
summ=c("mean","q025","q975")
# define the size of the matrices
S=10
################# #
N=paste(summ,".",indice,sep = "")
SL=SIM2$size==S
plot(0,xlim = c(0.2,1),ylim=c(0,1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',pch=NA, xaxs="i", yaxs="i")
axis(side = 1,lwd=1,at = c(0.3,0.5,0.7,0.9,0.92,0.94,0.96,0.98))
axis(side=2,lwd=1)
for (i in 1:length(N)){
  LTY=c(1,2,2)[i]
  lines(c(0.3,0.5,0.7,0.9,0.92,0.94,0.96,0.98),SIM2[[N[i]]][SIM2$PM=="Logn"&SL],col=col1,lty=LTY)
  lines(c(0.3,0.5,0.7,0.9,0.92,0.94,0.96,0.98),SIM3[[N[i]]][SIM2$PM=="Logn"&SL],col=col2,lty=LTY)
  lines(c(0.3,0.5,0.7,0.9,0.92,0.94,0.96,0.98),SIM2[[N[i]]][SIM2$PM=="Equi"&SL],col=col3,lty=LTY)
}

