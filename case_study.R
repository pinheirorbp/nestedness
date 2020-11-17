### Authors ####

#Rafael Barros Pereira Pinheiro 
#Carsten F. Dormann 
#Gabriel Moreira Felix Ferreira
#Marco Aurelio Ribeiro Mello 

# Contact: rafael-bpp@hotmail.com

### Case study ###
# The empirical network analyzed was downloaded from the Web of Life database (www.web-of-life.es, network: M_PL_060_17).
# Originally described by Kaiser-Bunbury 2010.
#Kaiser-Bunbury CN, Muff S, Memmott J, M?ller CB, Caflisch A. 2010. The robustness of pollination networks to the loss of species and interactions: a quantitative approach incorporating pollinator behaviour. Ecol. Lett. 13, 442-452. (doi:10.1111/j.1461-0248.2009.01437.x).
net1=as.matrix(read.table("data/net.txt"))
### Requirements ####
library(bipartite)
source("functions/nestmodels.R")
### Plot of the matrix ####
plotmatrix(net1, binary=F, within_color = "#0e0a44FF",base_color = "#d7d6e9ff")
### WNODF ####
WNODF=nestmodels(M=net1,index="wnodf",prop.model = T,n.model = 1000,calc.at = c(50,100,200,300,400,444), print.at.each = 100,density.plot = T,sampling.plot = T)
WNODF$observed
WNODF$significance
### WNODA ####
WNODA=nestmodels(M=net1,index="wnoda",prop.model = T,n.model = 1000,print.at.each = 10,density.plot = T,sampling.plot = T,calc.at = c(50,100,200,300,400,444))
WNODA$observed
WNODA$significance
### NODF ####
# Accounting for weighted marginal sums
NODF1=nestmodels(M=net1,index="nodf",prop.model = T,n.model = 1000, print.at.each = 10,density.plot = T,equi.model = T,wprob = T,wsamp=F)
NODF1$observed
NODF1$significance
# Only binary information
NODF2=nestmodels(M=net1,index="nodf",prop.model = T,n.model = 1000, print.at.each = 10,density.plot = T,equi.model = T,wprob = F,wsamp=F)
NODF2$observed
NODF2$significance

