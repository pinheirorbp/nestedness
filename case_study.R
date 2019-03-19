### Authors ####

#Rafael Barros Pereira Pinheiro 
#Carsten F. Dormann 
#Gabriel Moreira F?lix Ferreira
#Marco Aurelio Ribeiro Mello 

# Contact: rafael-bpp@hotmail.com

### Case study ###
# The empirical network analyzed was downloaded from the Web of Life database (www.web-of-life.es, network: M_PL_060_17).
# Originally described by Kaiser-Bunbury 2010.
#Kaiser-Bunbury CN, Muff S, Memmott J, Müller CB, Caflisch A. 2010. The robustness of pollination networks to the loss of species and interactions: a quantitative approach incorporating pollinator behaviour. Ecol. Lett. 13, 442-452. (doi:10.1111/j.1461-0248.2009.01437.x).
### Requirements ####
net1=as.matrix(read.table("data/net.txt"))
library(bipartite)
source("functions/nestnulls.R")
source("functions/plotmatrix.R")
source("functions/nest.smdm.R")
### Plot of the matrix ####
plotmatrix(net1, binary=F, within_color = "#0e0a44FF",base_color = "#d7d6e9ff")
### WNODF ####
WNODF=nestnulls(M=net1,index="wnodf",prop.null = T,n.null = 10000,calc.at = c(50,100,200,300,400,444), print.at.each = 100,density.plot = T)
### WNODA ####
WNODA=nestnulls(M=net1,index="wnoda",prop.null = T,n.null = 10000,print.at.each = 100,density.plot = T,sampling.plot = T)
### NODF ####
NODF=nestnulls(M=net1,index="nodf",prop.null = T,n.null = 10000, print.at.each = 1000,density.plot = T,equi.null = T,wprob = T,wsamp=F)
NODF2=nestnulls(M=net1,index="nodf",prop.null = T,n.null = 10000, print.at.each = 100,density.plot = T,equi.null = T,wprob = F, wsamp=F)
