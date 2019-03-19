plotmatrix= function(matrix, background_color="white", within_color = "black", between_color="gray", modules_colors=NULL, border = F,row_partitions=NULL,col_partitions=NULL, border_color="black", base_color="white", binary=T, plot_labels=F,xlab=NA,ylab=NA, offset = 0.4,...){
  
  matrix <- as.matrix(matrix)
  
  if(!is.null(row_partitions)){
    index=unique(row_partitions)
    row_partitions=match(row_partitions,index)
    col_partitions=match(col_partitions,index)
  }
  if(!is.null(row_partitions)){
    row_partitions=as.numeric(row_partitions)
    col_partitions=as.numeric(col_partitions)
  }
  ## Checking inputs
  if (border==T&is.null(row_partitions)){
    border=F
    warning("borders cannot be ploted if partitions are not defined")
  }
  if(is.null(row_partitions)){modules_colors=NULL}
  if (!is.null(modules_colors)&length(unique(row_partitions))!=length(modules_colors)){
    stop("wrong number of modules colors")
  }
  #if(!is.null(row_partitions)|!is.null(col_partitions)){
    #if(!identical(sort(unique(row_partitions)),sort(unique(col_partitions)))){
      #stop("different numbers of col and row partitions")
   #}
    #if(length(row_partitions)!=nrow(matrix)|length(col_partitions)!=ncol(matrix)){
      #stop ("partitions with inapropriate length")
    #}
  #}
  ####
  #Defining colors
  if(!is.null(modules_colors)){within_color=NULL}
  colors=c(between_color,within_color,modules_colors)
  #matrix2 = interaction matrices in which numbers defines the color of the rectangle
  matrix2=matrix
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
      if(!is.null(row_partitions)){
        if(is.null(modules_colors)){
          if(matrix[i,j]==0){matrix2[i,j]=0} else{
            matrix2[i,j]=ifelse(row_partitions[i]==col_partitions[j],2,1)
          }}
        if(!is.null(modules_colors)){
          if(matrix[i,j]==0){matrix2[i,j]=0} else{
            matrix2[i,j]=ifelse(row_partitions[i]!=col_partitions[j],1,row_partitions[i]+1)}
        }
      }
    }
  }
image(y=1:nrow(matrix2),x=1:ncol(matrix2),z=0*t(matrix2), col = background_color,xaxt = "n", yaxt = "n",xlab=xlab,ylab=ylab,...)
  if (plot_labels){axis(1,at = 1:ncol(matrix2), labels = colnames(matrix), ...)}
  if (plot_labels){axis(2,at = 1:nrow(matrix2), labels = rownames(matrix)[length(rownames(matrix)):1], ...)}
  box()
  
  TRmatrix2 <- t(matrix2[ nrow(matrix2):1, ] )
  
  if(binary==F){
    TRmatrix <- t(matrix[ nrow(matrix):1, ] )
    TRmatrix=TRmatrix/max(TRmatrix)
  }
  
  offset=offset
  
  #plotting filled positions as rectangles
  for (a in 1:nrow(TRmatrix2)) {
    for (b in 1:ncol(TRmatrix2)) { 
      if (TRmatrix2[a,b] != 0){
        if (is.null(row_partitions)){
          if(binary==T){rect(a-offset,b-offset,a+offset,b+offset,col=within_color,border='NA')}
          if(binary==F){
            MAX=max(matrix)
            funcol=colorRamp(c(base_color,within_color))
            rect_color=funcol(TRmatrix[a,b])
            rect(a-offset,b-offset,a+offset,b+offset,col=rgb(rect_color[1],rect_color[2],rect_color[3],255,maxColorValue = 255),border='NA')}
          
        } else { 
          if(binary==T){rect(a-offset,b-offset,a+offset,b+offset,
               col=colors[TRmatrix2[a,b]],border='NA')}
          if(binary==F){
            MAX=max(matrix)
            funcol=colorRamp(c(base_color,colors[TRmatrix2[a,b]]))
            rect_color=funcol(TRmatrix[a,b])
            rect(a-offset,b-offset,a+offset,b+offset,col=rgb(rect_color[1],rect_color[2],rect_color[3],255,maxColorValue = 255),border='NA')
          }
        }
      }
    }
  }
  # plotting borders
  if (border==T) {
    c=0
    for (i in 1:length(matrix)){
      if(i%%nrow(matrix)==1){c=c+1}
      r=i%%nrow(matrix)
      if(r==0){r=nrow(matrix)}
      mod=row_partitions[r]==col_partitions[c]
      if(r==1){modtop=F
      }else{modtop=row_partitions[r-1]==col_partitions[c]}
      if(c==1){modleft=F
      }else{modleft=row_partitions[r]==col_partitions[c-1]}
      
      if(mod&!modtop&!modleft){
        xleft=c-0.5
        ytop=nrow(matrix)-r+1.5
        endmod=F
        j=c
        while(!endmod){
          j=j+1
          if((j-1)==ncol(matrix)){endmod=T
          }else{
            endmod=col_partitions[j]!=col_partitions[j-1]}
        }
        xright=j-0.5
        endmod=F
        j=r
        while(!endmod){
          j=j+1
          if((j-1)==nrow(matrix)){endmod=T
          }else{
            endmod=row_partitions[j]!=row_partitions[j-1]}
        }
        ybottom=nrow(matrix)-j+1.5
        rect(xleft = xleft,xright = xright,ybottom = ybottom,ytop = ytop,col="NA",border=border_color)
      }
    }
    
  }
}
