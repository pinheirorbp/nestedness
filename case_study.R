################################################################################

# Supplementary code and data from the study: Pinheiro, R. B. P., Dormann, 
# C. F., Felix, G. M. F., & Mello M. A. R. Linking null models to hypotheses
# to improve nestedness analysis. Submitted.

################################################################################


# Authors:

# Rafael R. B. Pinheiro 
# Carsten F. Dormann 
# Gabriel M. Felix
# Marco A. R. Mello 

# Contact: rafael-bpp@hotmail.com

# Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com.

# See README for the complete user guide to our new protocol:
# https://github.com/pinheirorbp/nestedness/blob/master/readme.pdf

# Disclaimer: You may use this script freely for commercial or non-commercial 
# purposes at your own risk. We assume no responsibility or liability for the 
# use of this software, convey no license or title under any patent, copyright,
# or mask work right to the product. We reserve the right to make changes in 
# the software without notification. We also make no representation or warranty
# that such application will be suitable for the specified use without further
# testing or modification. If this script helps you produce any academic work
# (paper, book, chapter, dissertation, monograph, thesis, report, talk, keynote,
# lecture, and similar), please acknowledge the authors and cite the source.


##### CASE STUDY #####


# Use the following code to reproduce the case study presented in our paper.
# The file net.txt represents the network analyzed by Kaiser-Bunbury et al. 
# (2010). If you want to apply our new protocol to your own data, just replace
# net.txt with another network of your choice, following the same format.

# The empirical network analyzed in our paper was downloaded from the Web of 
# Life database (http://www.web-of-life.es, network: M_PL_060_17). It was
# originally analyzed by: 

# Kaiser-Bunbury CN, Muff S, Memmott J, M?ller CB, Caflisch A. 2010. The 
# robustness of pollination networks to the loss of species and interactions: 
# a quantitative approach incorporating pollinator behaviour. Ecol. Lett. 13,
# 442-452. (http://dx.doi.org/10.1111/j.1461-0248.2009.01437.x).

# Create a new object with the network to be analyzed.
net1 <- as.matrix(read.table("data/net.txt"))

# Load the required packages
library(bipartite)

# Load our custom-made function to calculate nestedness and run the null models.
# Check out the comments in the respective R file to obtain detailed information
# about the function, including its arguments and options.
source("functions/nestmodels.R")

# Plot the matrix to get a feeling about its structure.
plotmatrix(net1, 
           binary = F, 
           within_color = "#0e0a44FF",
           base_color = "#d7d6e9ff")

# The simulations included in our protocol are quite resource-consuming and may
# require several minutes, depending on your computer's processing power. 
# Therefore, choose the number of randomized matrices (iterations) carefully.
# In our paper, we ran 10,000 iterations, but we suggest experimenting with 
# smaller numbers first.

matrices <- 10000 #This number will be the same for all functions

# Generate randomized matrices and calculate their degree of nestedness using
# different metrics.

# WNODF
WNODF=nestmodels(M=net1,
                 index="wnodf",
                 n.model = matrices,
                 prop.model = T,
                 print.at.each = 100,
                 density.plot = T,
                 sampling.plot = T,
                 calc.at = c(50,100,200,300,400,444))
WNODF$observed
WNODF$significance

# WNODA
WNODA=nestmodels(M=net1,
                 index="wnoda",
                 n.model = matrices,
                 prop.model = T,
                 print.at.each = 100,
                 density.plot = T,
                 sampling.plot = T,
                 calc.at = c(50,100,200,300,400,444))
WNODA$observed
WNODA$significance

# NODF
# Accounting for weighted marginal sums
NODF1=nestmodels(M=net1,
                 index="nodf",
                 n.model = matrices, 
                 prop.model = T,
                 equi.model = T,
                 print.at.each = 100,
                 density.plot = T,
                 wprob = T,
                 wsamp=F)
NODF1$observed
NODF1$significance

# Using only binary information
NODF2=nestmodels(M=net1,
                 index="nodf",
                 n.model = matrices, 
                 prop.model = T,
                 equi.model = T,
                 print.at.each = 100,
                 density.plot = T,
                 wprob = F,
                 wsamp=F)
NODF2$observed
NODF2$significance

