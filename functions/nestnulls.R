nestnulls=function(M, index="wnoda",equi.null=T, prop.null=T, n.null=1000, calc.at=NULL, print.at.each=NULL, wprob=F,wsamp=F, density.plot=T,sampling.plot=T){
  # functions ####
  z1samp= function (x, mean, sd){
    z= (x-mean)/sd
    p= 2*pnorm(-abs(z))
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
  if(index=="wnoda"){nest.smdm =function(x,constrains=NULL, weights=F, decreasing="fill", sort=T){
      ### Checking inputs ####
      if (!is.null(constrains)&length(unique(constrains))==1){
        warning("Only one module. Nestedness calculated only for the entire matrix")
        constrains=NULL}
      if(is.element(NA,constrains)|is.element(NaN,constrains)){
        warning("NA or NaN in constrains. Nestedness calculated only for the entire matrix")
        constrains=NULL
      }
      if (!is.null(constrains)&length(constrains)!=nrow(x)+ncol(x)){
        stop("constrains vector is not of the same length that network vertices")
      }
      if (weights==F&any(x!=0&x!=1)){
        x[x>0]=1
        warning ("binary metric applied")
      }
      if (decreasing!="fill"&decreasing!="abund"){
        stop("decreasing should be fill or abund")
      }
      if (!is.null(constrains)){constrains=as.character(constrains)}
      if(is.null(rownames(x))){
        xrnames=paste("R",1:nrow(x),"")
        rownames(x)<-xrnames
      }
      if(is.null(colnames(x))){
        xcnames=paste("C",1:ncol(x),"")
        colnames(x)<-xcnames
      }
      ### Unweighted NODF Function ####
      unweightednodf=function (x,constrains){
        # Sorting matrix order by row and collumn sums
        if (sort==T){tab0=x[sort(rowSums(x), index=T, decreasing=TRUE)$ix,
                            sort(colSums(x), index=T, decreasing=TRUE)$ix]}
        else {tab0=x}
        
        # N for rows
        MTrow= rowSums(tab0)
        Nrow= matrix(rep(NA,times=nrow(tab0)^2),nrow(tab0),nrow(tab0))
        dimnames(Nrow)=list(rownames(tab0),rownames(tab0))
        
        for (jrow in 2:nrow(tab0)){
          for (irow in 1:(jrow-1)){
            if (MTrow[jrow]>=MTrow[irow]){Nrow[jrow,irow]=0} 
            else {
              S=0
              for(i in 1:ncol(tab0)){
                if (tab0[jrow,i]==1&tab0[jrow,i]==tab0[irow,i]) {
                  S=S+1
                }
              }
              Nrow[jrow,irow]=S*100/MTrow[jrow]
            }
          }
          
        }      
        Nrow=Nrow[rownames(x), rownames(x)]
        
        # NODF for rows
        NODFrow= mean(Nrow,na.rm = T)
        
        # N for collumns
        
        MTcol= colSums(tab0)
        Ncol= matrix(rep(NA,times=ncol(tab0)^2),ncol(tab0),ncol(tab0))
        dimnames(Ncol)=list(colnames(tab0),colnames(tab0))
        
        for (jcol in 2:ncol(tab0)){
          for (icol in 1:(jcol-1)){
            if (MTcol[jcol]>=MTcol[icol]){Ncol[jcol,icol]=0} 
            else {
              S=0
              for(i in 1:nrow(tab0)){
                if (tab0[i,jcol]==1&tab0[i,jcol]==tab0[i,icol]) {
                  S=S+1
                }
              }
              Ncol[jcol,icol]=S*100/MTcol[jcol]
            }
          }
          
        }      
        Ncol=Ncol[colnames(x),colnames(x)]
        
        # NODF for rows
        NODFcol= mean(Ncol,na.rm = T)
        
        # NODF for the entire matrix
        NODFmatrix= mean(c(Ncol,Nrow),na.rm=T)
        
        #### NODF SM/DM ###
        if (!is.null(constrains)){
          # Constrains for rows
          
          rowcons=cbind (rownames(x),constrains[1:nrow(x)])
          tabrcons=table(rowcons[,1],rowcons[,2])
          distrcons= dist(tabrcons,method = "binary")
          distrcons= as.matrix (distrcons)
          distrcons=distrcons[rownames(x),rownames(x)]
          rm(rowcons,tabrcons)
          
          # NODF SM/DM for rows
          SM_Nrow=0
          SM_nrow=0
          DM_Nrow=0
          DM_nrow=0
          for (i in 1:nrow(x)){
            for (j in 1:nrow(x)){
              if (!is.na(Nrow[i,j])){
                if(distrcons[i,j]==0){
                  SM_Nrow=SM_Nrow+Nrow[i,j]
                  SM_nrow=SM_nrow+1
                }
                else{
                  DM_Nrow=DM_Nrow+Nrow[i,j]
                  DM_nrow=DM_nrow+1
                }
              }
            }
          }
          NODF_SM_row= SM_Nrow/SM_nrow
          NODF_DM_row= DM_Nrow/DM_nrow
          
          # Constrains for collumns
          
          colcons=cbind (colnames(x),constrains[(nrow(x)+1):length(constrains)])
          tabccons=table(colcons[,1],colcons[,2])
          distccons= dist(tabccons,method = "binary")
          distccons= as.matrix (distccons)
          distccons=distccons[colnames(x),colnames(x)]
          rm(colcons,tabccons)
          
          # NODF SM/DM for collumns
          SM_Ncol=0
          SM_ncol=0
          DM_Ncol=0
          DM_ncol=0
          for (i in 1:ncol(x)){
            for (j in 1:ncol(x)){
              if (!is.na(Ncol[i,j])){
                if(distccons[i,j]==0){
                  SM_Ncol=SM_Ncol+Ncol[i,j]
                  SM_ncol=SM_ncol+1
                }
                else{
                  DM_Ncol=DM_Ncol+Ncol[i,j]
                  DM_ncol=DM_ncol+1
                }
              }
            }
          }
          NODF_SM_col= SM_Ncol/SM_ncol
          NODF_DM_col= DM_Ncol/DM_ncol
          
          # NODF SM/DM for matrix
          
          NODF_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
          NODF_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
          # return
          return(list(NODFrow=NODFrow,NODFcol=NODFcol, NODFmatrix=NODFmatrix,
                      NODF_SM_row= NODF_SM_row, NODF_DM_row=NODF_DM_row, 
                      NODF_SM_col= NODF_SM_col, NODF_DM_col=NODF_DM_col,
                      NODF_SM_matrix= NODF_SM_matrix, NODF_DM_matrix=NODF_DM_matrix))
          
        }
        else {
          return(list(NODFrow=NODFrow,NODFcol=NODFcol, NODFmatrix=NODFmatrix))}
      }
      ### Weighted NODF function ####
      weightednodf=function (x,constrains){
        # Sorting matrix order by row and collumn sums
        if(sort==T){tab0=x[sort(rowSums(x!=0), index=T, decreasing=TRUE)$ix,
                           sort(colSums(x!=0), index=T, decreasing=TRUE)$ix]}
        else{tab0=x}
        
        # N for rows
        MTrow= rowSums(tab0)
        Frow= rowSums(tab0!=0)
        Nrow= matrix(rep(NA,times=nrow(tab0)^2),nrow(tab0),nrow(tab0))
        dimnames(Nrow)=list(rownames(tab0),rownames(tab0))
        
        for (jrow in 2:nrow(tab0)){
          for (irow in 1:(jrow-1)){
            if (Frow[jrow]>=Frow[irow]){Nrow[jrow,irow]=0} 
            else {
              S=0
              for(i in 1:ncol(tab0)){
                if (tab0[jrow,i]!=0&tab0[jrow,i]<tab0[irow,i]) {
                  S=S+1
                }
              }
              Nrow[jrow,irow]=S*100/Frow[jrow]
            }
          }
          
        }      
        Nrow=Nrow[rownames(x), rownames(x)]
        
        # WNODF for rows
        NODFrow= mean(Nrow,na.rm = T)
        
        # N for collumns
        
        MTcol= colSums(tab0)
        Fcol= colSums(tab0!=0)
        Ncol= matrix(rep(NA,times=ncol(tab0)^2),ncol(tab0),ncol(tab0))
        dimnames(Ncol)=list(colnames(tab0),colnames(tab0))
        
        for (jcol in 2:ncol(tab0)){
          for (icol in 1:(jcol-1)){
            if (Fcol[jcol]>=Fcol[icol]){Ncol[jcol,icol]=0}
            else {
              S=0
              for(i in 1:nrow(tab0)){
                if (tab0[i,jcol]!=0&tab0[i,jcol]<tab0[i,icol]) {
                  S=S+1
                }
              }
              Ncol[jcol,icol]=S*100/Fcol[jcol]
            }
          }
          
        }      
        Ncol=Ncol[colnames(x),colnames(x)]
        
        # WNODF for rows
        NODFcol= mean(Ncol,na.rm = T)
        
        # WNODF for the entire matrix
        NODFmatrix= mean(c(Ncol,Nrow),na.rm=T)
        
        #### WNODF SM/DM ###
        if (!is.null(constrains)){
          # Constrains for rows
          
          rowcons=cbind (rownames(x),constrains[1:nrow(x)])
          tabrcons=table(rowcons[,1],rowcons[,2])
          distrcons= dist(tabrcons,method = "binary")
          distrcons= as.matrix (distrcons)
          distrcons=distrcons[rownames(x),rownames(x)]
          rm(rowcons,tabrcons)
          
          # WNODF SM/DM for rows
          SM_Nrow=0
          SM_nrow=0
          DM_Nrow=0
          DM_nrow=0
          for (i in 1:nrow(x)){
            for (j in 1:nrow(x)){
              if (!is.na(Nrow[i,j])){
                if(distrcons[i,j]==0){
                  SM_Nrow=SM_Nrow+Nrow[i,j]
                  SM_nrow=SM_nrow+1
                }
                else{
                  DM_Nrow=DM_Nrow+Nrow[i,j]
                  DM_nrow=DM_nrow+1
                }
              }
            }
          }
          NODF_SM_row= SM_Nrow/SM_nrow
          NODF_DM_row= DM_Nrow/DM_nrow
          
          # Constrains for collumns
          
          colcons=cbind (colnames(x),constrains[(nrow(x)+1):length(constrains)])
          tabccons=table(colcons[,1],colcons[,2])
          distccons= dist(tabccons,method = "binary")
          distccons= as.matrix (distccons)
          distccons=distccons[colnames(x),colnames(x)]
          rm(colcons,tabccons)
          
          # WNODF SM/DM for collumns
          SM_Ncol=0
          SM_ncol=0
          DM_Ncol=0
          DM_ncol=0
          for (i in 1:ncol(x)){
            for (j in 1:ncol(x)){
              if (!is.na(Ncol[i,j])){
                if(distccons[i,j]==0){
                  SM_Ncol=SM_Ncol+Ncol[i,j]
                  SM_ncol=SM_ncol+1
                }
                else{
                  DM_Ncol=DM_Ncol+Ncol[i,j]
                  DM_ncol=DM_ncol+1
                }
              }
            }
          }
          NODF_SM_col= SM_Ncol/SM_ncol
          NODF_DM_col= DM_Ncol/DM_ncol
          
          # WNODF SM/DM for matrix
          NODF_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
          NODF_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
          # return
          return(list(WNODFrow=NODFrow,WNODFcol=NODFcol, WNODFmatrix=NODFmatrix,WNODF_SM_row= NODF_SM_row, WNODF_DM_row=NODF_DM_row,WNODF_SM_col= NODF_SM_col, WNODF_DM_col=NODF_DM_col,WNODF_SM_matrix= NODF_SM_matrix, WNODF_DM_matrix=NODF_DM_matrix))
          
        }
        else {
          return(list(WNODFrow=NODFrow,WNODFcol=NODFcol, WNODFmatrix=NODFmatrix))}
      }
      ### Weighted NODA funcion ####
      weightednoda=function (x,constrains){
        # Sorting matrix order by row and collumn sums
        if(sort==T){tab0=x[sort(rowSums(x), index=T, decreasing=TRUE)$ix,
                           sort(colSums(x), index=T, decreasing=TRUE)$ix]}
        else{tab0=x}
        
        # N for rows
        MTrow= rowSums(tab0)
        Frow= rowSums(tab0!=0)
        Nrow= matrix(rep(NA,times=nrow(tab0)^2),nrow(tab0),nrow(tab0))
        dimnames(Nrow)=list(rownames(tab0),rownames(tab0))
        
        for (jrow in 2:nrow(tab0)){
          for (irow in 1:(jrow-1)){
            if (MTrow[jrow]>=MTrow[irow]){Nrow[jrow,irow]=0}
            else {
              S=0
              for(i in 1:ncol(tab0)){
                if (tab0[jrow,i]!=0&tab0[jrow,i]<tab0[irow,i]) {
                  S=S+1
                }
              }
              Nrow[jrow,irow]=S*100/Frow[jrow]
            }
          }
          
        }      
        Nrow=Nrow[rownames(x), rownames(x)]
        
        # WNODA for rows
        NODArow= mean(Nrow,na.rm = T)
        
        # N for collumns
        
        MTcol= colSums(tab0)
        Fcol= colSums(tab0!=0)
        Ncol= matrix(rep(NA,times=ncol(tab0)^2),ncol(tab0),ncol(tab0))
        dimnames(Ncol)=list(colnames(tab0),colnames(tab0))
        
        for (jcol in 2:ncol(tab0)){
          for (icol in 1:(jcol-1)){
            if (MTcol[jcol]>=MTcol[icol]){Ncol[jcol,icol]=0}
            else {
              S=0
              for(i in 1:nrow(tab0)){
                if (tab0[i,jcol]!=0&tab0[i,jcol]<tab0[i,icol]) {
                  S=S+1
                }
              }
              Ncol[jcol,icol]=S*100/Fcol[jcol]
            }
          }
          
        }      
        Ncol=Ncol[colnames(x),colnames(x)]
        
        # NODA for rows
        NODAcol= mean(Ncol,na.rm = T)
        
        # NODA for the entire matrix
        NODAmatrix= mean(c(Ncol,Nrow),na.rm=T)
        
        #### WNODA SM/DM ###
        if (!is.null(constrains)){
          
          # Constrains for rows
          rowcons=cbind (rownames(x),constrains[1:nrow(x)])
          tabrcons=table(rowcons[,1],rowcons[,2])
          distrcons= dist(tabrcons,method = "binary")
          distrcons= as.matrix (distrcons)
          distrcons=distrcons[rownames(x),rownames(x)]
          rm(rowcons,tabrcons)
          
          # WNODA SM/DM for rows
          SM_Nrow=0
          SM_nrow=0
          DM_Nrow=0
          DM_nrow=0
          for (i in 1:nrow(x)){
            for (j in 1:nrow(x)){
              if (!is.na(Nrow[i,j])){
                if(distrcons[i,j]==0){
                  SM_Nrow=SM_Nrow+Nrow[i,j]
                  SM_nrow=SM_nrow+1
                }
                else{
                  DM_Nrow=DM_Nrow+Nrow[i,j]
                  DM_nrow=DM_nrow+1
                }
              }
            }
          }
          NODA_SM_row= SM_Nrow/SM_nrow
          NODA_DM_row= DM_Nrow/DM_nrow
          
          # Constrains for collumns
          
          colcons=cbind (colnames(x),constrains[(nrow(x)+1):length(constrains)])
          tabccons=table(colcons[,1],colcons[,2])
          distccons= dist(tabccons,method = "binary")
          distccons= as.matrix (distccons)
          distccons=distccons[colnames(x),colnames(x)]
          rm(colcons,tabccons)
          
          # WNODA SM/DM for collumns
          SM_Ncol=0
          SM_ncol=0
          DM_Ncol=0
          DM_ncol=0
          for (i in 1:ncol(x)){
            for (j in 1:ncol(x)){
              if (!is.na(Ncol[i,j])){
                if(distccons[i,j]==0){
                  SM_Ncol=SM_Ncol+Ncol[i,j]
                  SM_ncol=SM_ncol+1
                }
                else{
                  DM_Ncol=DM_Ncol+Ncol[i,j]
                  DM_ncol=DM_ncol+1
                }
              }
            }
          }
          NODA_SM_col= SM_Ncol/SM_ncol
          NODA_DM_col= DM_Ncol/DM_ncol
          
          # WNODA SM/DM for matrix
          
          NODA_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
          NODA_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
          # return
          return(list(WNODArow=NODArow,WNODAcol=NODAcol, WNODAmatrix=NODAmatrix,
                      WNODA_SM_row= NODA_SM_row, WNODA_DM_row=NODA_DM_row, 
                      WNODA_SM_col= NODA_SM_col, WNODA_DM_col=NODA_DM_col,
                      WNODA_SM_matrix= NODA_SM_matrix, WNODA_DM_matrix=NODA_DM_matrix))
          
        }
        else {
          return(list(WNODArow=NODArow,WNODAcol=NODAcol, WNODAmatrix=NODAmatrix))}
      }
      
      ### Using functions ####
      if(decreasing=="abund"){
        return(weightednoda(x,constrains))
      }
      if (decreasing=="fill"){
        if (weights==F){
          return(unweightednodf(x,constrains))
        }
        if (weights==T){
          return(weightednodf(x,constrains))
        }
      }
    }}
  if(index=="temperature"){library(bipartite)}
  if(index=="betaNES"|index=="betaNES2"){library(betapart)}
  # binary index (F / T) ####
  if(is.element(index,c("nodf","temperature","MD","binSR", "betaNES", "betaNES2"))){
    binindex=T}else{binindex=F}
  # network ####
  M=M[rowSums(M, na.rm=T)>0,colSums(M, na.rm=T)>0]
  net=M
  if(binindex){
    net[M>0]=1
    if(!wsamp){calc.at=sum(net, na.rm=T)}
  }
  # probability matrix ####
  if(prop.null&!wprob&binindex){
    PM=matrix(NA,nrow(net),ncol(net))
    RS=rowSums(net, na.rm=T)/ncol(net)
    CS=colSums(net,na.rm=T)/nrow(net)
    for (r in 1:nrow(net)){
      for (c in 1:ncol(net)){
        PM[r,c]=(RS[r]+CS[c])/2
      }}
  }else if(prop.null){
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
    observed=nest.smdm(net,weights = T,decreasing = "abund")$WNODAmatrix/100
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
  # loop for null models ####
  equinull=list()
  propnull=list()
  for (N in 1:n.null){
    # equiprobable ####
    if(equi.null){
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
          equinull[[samp]][N]=nestednodf(NM,weighted = T)[[3]][[3]]/100
        }else if(index=="wnoda"){
          equinull[[samp]][N]=nest.smdm(NM,weights = T,decreasing = "abund")$WNODAmatrix/100
        }else if(index=="SR"|index=="binSR"){
          equinull[[samp]][N]=SPECTRAL_RADIUS(NM)
        }else if(index=="MD"){
          equinull[[samp]][N]=MANHATTAN_DISTANCE(NM)
        }else if(index=="nodf"){
          equinull[[samp]][N]=nestednodf(NM)[[3]][[3]]/100
        }else if(index=="temperature"){
          equinull[[samp]][N]=bipartite::networklevel(NM, index="nestedness")/100
        }else if(index=="betaNES"){
          equinull[[samp]][N]=beta.multi(NM)$beta.SNE
        }else if(index=="betaNES2"){
          bm=beta.multi(NM)
          equinull[[samp]][N]=bm$beta.SNE/bm$beta.SOR
        }
      }
    }else{equinull=NULL}
    # proportional ####
    if(prop.null){
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
          propnull[[samp]][N]=nestednodf(NM,weighted = T)[[3]][[3]]/100
        }else if(index=="wnoda"){
          propnull[[samp]][N]=nest.smdm(NM,weights = T,decreasing = "abund")$WNODAmatrix/100
        }else if(index=="SR"|index=="binSR"){
          propnull[[samp]][N]=SPECTRAL_RADIUS(NM)
        }else if(index=="MD"){
          propnull[[samp]][N]=MANHATTAN_DISTANCE(NM)
        }else if(index=="nodf"){
          propnull[[samp]][N]=nestednodf(NM)[[3]][[3]]/100
        }else if(index=="temperature"){
          propnull[[samp]][N]=bipartite::networklevel(NM, index="nestedness")/100
        }else if(index=="betaNES"){
          propnull[[samp]][N]=beta.multi(NM)$beta.SNE
        }else if(index=="betaNES2"){
          bm=beta.multi(NM)
          propnull[[samp]][N]=bm$beta.SNE/bm$beta.SOR
        }
      }
    }else{propnull=NULL}
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
  if(equi.null){
    PZ=z1samp(observed,mean(equinull[[AS]],na.rm=T),sd(equinull[[AS]],na.rm=T))
    significance$pvalue.equi=PZ$Pvalue
    significance$Zscore.equi=PZ$ZScore
    }
  if(prop.null){
    PZ=z1samp(observed,mean(propnull[[AS]],na.rm=T),sd(propnull[[AS]],na.rm=T))
    significance$pvalue.prop=PZ$Pvalue
    significance$Zscore.prop=PZ$ZScore
  }
  # for each sampling level
  summary.equi=list()
  if(equi.null){
    for (i in 1:length(calc.at)){
      SAM=names(equinull)[i]
      summary.equi[[SAM]]$nnull=sum(!is.na(equinull[[i]]))
      summary.equi[[SAM]]$mean.equi=mean(equinull[[i]], na.rm = T)
      summary.equi[[SAM]]$sd.equi=sd(equinull[[i]], na.rm = T)
      summary.equi[[SAM]]$plus95.equi=quantile(equinull[[i]], na.rm = T, probs = 0.975)[[1]]
      summary.equi[[SAM]]$minus95.equi=quantile(equinull[[i]], na.rm = T, probs = 0.025)[[1]]
    }
  }
  summary.prop=list()
  if(prop.null){
    for (i in 1:length(calc.at)){
      SAM=names(propnull)[i]
      summary.prop[[SAM]]$mean.prop=mean(propnull[[i]], na.rm = T)
      summary.prop[[SAM]]$sd.prop=sd(propnull[[i]], na.rm = T)
      summary.prop[[SAM]]$plus95.prop=quantile(propnull[[i]], na.rm = T, probs = 0.975)[[1]]
      summary.prop[[SAM]]$minus95.prop=quantile(propnull[[i]], na.rm = T, probs = 0.025)[[1]]
    }
  }
  ## plot 1: density plot ####
  if(density.plot){
    if(equi.null){deneq=density(equinull[[AS]],na.rm=T)}else{deneq=list(y=NA)}
    if(prop.null){denprop=density(propnull[[AS]],na.rm=T)}else{denprop=list(y=NA)}
    if(is.element(index,c("MD","SR","binSR"))){
      XLIM=c(min(c(equinull[[AS]],propnull[[AS]],observed),na.rm=T),max(c(equinull[[AS]],propnull[[AS]],observed),na.rm=T))
    }else{XLIM=c(0,1)}
    YLIM=c(0,max(c(deneq$y,denprop$y),na.rm=T))
    plot(NA,ylim=YLIM,xlim=XLIM, lwd=2,main=NA, xlab=index, ylab="density")
    if(equi.null){lines(deneq$x,deneq$y,lwd=2,col="#C65C75FF")}
    if(prop.null){
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
    if(equi.null){
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
    if(prop.null){
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
  parameters=list(index=index,equi.null=equi.null, prop.null=prop.null, n.null=n.null, calc.at=calc.at, print.at.each=print.at.each, wprob=wprob, wsamp=wsamp, binindex=binindex)
  return(list(observed=observed,significance=significance,summary.equi=summary.equi,summary.prop=summary.prop,equinull=equinull, propnull=propnull, nsampled=as0, parameters=parameters))
}