# SPECTRAL_RADIUS.R
# Part of the FALCON (Framework of Adaptive ensembLes for the Comparison Of
# Nestedness) package: https://github.com/sjbeckett/FALCON
# Last updated: 12th March 2014

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


# MANHATTAN_DISTANCE.R
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
