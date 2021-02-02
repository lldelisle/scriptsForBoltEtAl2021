# From v16: bug solved when normalizing with extension
library(RColorBrewer)
colp<-colorRampPalette(c("white","#ffe1aa","#ffbc41","#ff7f00","#ff3000","#bc0000","#5a0000"))

isValidColor <- function(colorname){
  if(is.numeric(colorname)){
    return(TRUE)
  }
  if(colorname %in% colors()){
    return(TRUE)
  }
  if(is.character(colorname)){
    if(nchar(colorname)==7 || nchar(colorname)==9){
      if(substr(colorname,1,1)=="#"){
        #I should do other checks
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

library(plotrix)
bigFunctionDrawTADs<-function(start,end,bin,matrixDF,colForNA="grey",domainsDF=NULL,giveM=F,drawTAD=F,
                   title="",threshold=NULL,annotation=NULL,shift=0,extension="",
                   TADline=NULL,TADlinecolor="blue",start2=NULL,end2=NULL,matrixDF2=NULL, extension2="",
                   operation="",maxim=NULL,forceRange=NULL,maximl2R=NULL,forceRangel2R=NULL,wholeScaleBar=F,normalize=F){
  if(!operation%in%c("","dif","l2r")){
    stop("This operation is not allowed.")
  }
  list1<-getSubMatAndRanges(matrix = matrixDF,start = start,end = end,bin = bin,start2 = start2,end2 = end2,extension = extension)
  m1<-list1[[1]]
  range1<-list1[[2]]
  range2<-list1[[3]]
  if(!is.null(matrixDF2)){
    list2<-getSubMatAndRanges(matrix = matrixDF2,start = start,end = end,bin = bin,start2 = start2,end2 = end2,extension = extension2)
    m2<-list2[[1]]
    range1_from2<-list2[[2]]
    range2_from2<-list2[[3]]
    #In case the ranges are not the same:
    range1<-intersect(range1,range1_from2)
    range2<-intersect(range2,range2_from2)
    m1<-getSubMatAndRanges(matrix = matrixDF,start = startFromRange(range1,bin),end = endFromRange(range1,bin),bin = bin,start2 = startFromRange(range2,bin),end2 = endFromRange(range2,bin),extension = extension)[[1]]
    m2<-getSubMatAndRanges(matrix = matrixDF2,start = startFromRange(range1,bin),end = endFromRange(range1,bin),bin = bin,start2 = startFromRange(range2,bin),end2 = endFromRange(range2,bin),extension = extension)[[1]]
    if(normalize){
      # I use to normalize a matrix 3 times wider than the one which is drawn
      bigRange<-min(range1,range2):max(range1,range2)
      new.start<-startFromRange(bigRange,bin)-length(bigRange)*bin
      new.end<-endFromRange(bigRange,bin)+length(bigRange)*bin
      list1norm<-getSubMatAndRanges(matrix = matrixDF,start = new.start,end = new.end,bin = bin, extension = extension)
      list2norm<-getSubMatAndRanges(matrix = matrixDF2,start = new.start,end = new.end,bin = bin, extension = extension2)
      new.start<-startFromRange(intersect(list1norm[[2]],list2norm[[2]]), bin)
      new.end<-endFromRange(intersect(list1norm[[2]],list2norm[[2]]), bin)
      norm1Mat<-getSubMatAndRanges(matrix = matrixDF,start = new.start,end = new.end,bin = bin, extension = extension)[[1]]
      norm2Mat<-getSubMatAndRanges(matrix = matrixDF2,start = new.start,end = new.end,bin = bin, extension = extension2)[[1]]
      m1<-m1/mean(as.matrix(norm1Mat),na.rm=T)
      m2<-m2/mean(as.matrix(norm2Mat),na.rm=T)
    }
    if(operation=="dif"){
      m<-m1-m2
    } else if(operation=="l2r"){
      m<-log2(m1/m2)
    }
  } else{
    if(normalize){
      if(min(m1,na.rm=T)>=0){
        # I use to normalize a matrix 3 times wider than the one which is drawn
        bigRange<-min(range1,range2):max(range1,range2)
        new.start<-startFromRange(bigRange,bin)-length(bigRange)*bin
        new.end<-endFromRange(bigRange,bin)+length(bigRange)*bin
        norm1Mat<-getSubMatAndRanges(matrix = matrixDF,start = new.start,end = new.end,bin = bin)[[1]]
        m1<-m1/mean(as.matrix(norm1Mat),na.rm=T)
      } else{
        cat("There is no sense to normalize a matrix with negative values\nNo normalization applied.\n")
      }
    }
    m<-m1
  }
  copym<-m
  #Saturation
  if(operation==""){
    if(is.null(threshold)){
      threshold<-quantile(as.vector(as.matrix(m)),probs = 0.95,na.rm=T)
      print(threshold)
    }
    m<-m1/threshold
    m[m>1]<-1
    #Testing
    minm<-min(as.vector(m),na.rm = T)
    maxm<-max(as.vector(m),na.rm=T)
    minc<-1+ceiling(minm*98)
    maxc<-1+floor(maxm*98)
    colorForMat<-colp(99)[minc:maxc]
    #Before:
    #colorForMat<-colp(200)
  } else if(operation%in%c("dif","l2r")){
    if(operation=="l2r"){
      maxim<-maximl2R
      forceRange<-forceRangel2R
    }
    if(!is.null(maxim)){
      m[m>maxim]<-maxim
      m[m<(-maxim)]<- -maxim
    }
    if(!is.null(forceRange)){
      minr<-forceRange[1]
      maxr<-forceRange[2]
      m[m>maxr]<-maxr
      m[m<minr]<- minr
      minm<-min(as.vector(m),na.rm = T)
      maxm<-max(as.vector(m),na.rm=T)
      minc<-1+ceiling((minm-minr)/(maxr-minr)*98)
      maxc<-1+floor((maxm-minr)/(maxr-minr)*98)
    } else{
      minm<-min(as.vector(m),na.rm = T)
      maxm<-max(as.vector(m),na.rm = T)
      if(minm>=0){
        minc<-1
        maxc<-99
      } else{
        if ((-minm)<maxm){
          maxc<-99
          minc<-ceiling(50+(minm/maxm)*49)
        } else{
          minc<-1
          maxc<-floor(50+(maxm/(-minm))*49)
        }
      }
    }
    colorForMat<-colorRampPalette(c("blue","white","red"))(99)[minc:maxc]
  }
  layout(matrix(c(1,2), nrow=2, ncol=1), widths=c(lcm(ncol(m)/8+4)), heights=c(lcm(nrow(m)/8+4),lcm(2.5)))
  par(mai=rep(2/2.54,4))
  image(1:ncol(m), 1:nrow(m), as.matrix(t(m[nrow(m):1,])),
        col = colorForMat,axes=F,xlab = "",ylab="")
  polygon(x=c(0.5,0.5,ncol(m)+0.5,ncol(m)+0.5),y=c(0.5,nrow(m)+0.5,nrow(m)+0.5,0.5),col=colForNA,border =NA)
  image(1:ncol(m), 1:nrow(m), as.matrix(t(m[nrow(m):1,])),add=T,
        col = colorForMat,axes=F,xlab = "",ylab="")  
  if(start==start2 && end==end2){
    title(main=paste(title,paste0((floor(start/bin))*bin/1e6,"Mb -",ceiling(end/bin)*bin/1e6,"Mb"),sep="\n"))
    yaxis=pretty(c(ceiling(start/bin)*bin,floor(end/bin)*bin))
    ylabels=paste(yaxis/1000000,"Mb")
    axis(2, at = yaxis/bin-(floor(start/bin))+0.5,labels=rev(ylabels))
  } else {
    title(main=paste(title,paste0((floor(start/bin))*bin/1e6,"Mb -",ceiling(end/bin)*bin/1e6,"Mb x",
                                  (floor(start2/bin))*bin/1e6,"Mb -",ceiling(end2/bin)*bin/1e6,"Mb"),sep="\n"))
    yaxis=pretty(c(ceiling(start/bin)*bin,floor(end/bin)*bin))
    ylabels=paste(rev(yaxis)/1000000,"Mb")
    axis(2, at = yaxis/bin-(floor(start/bin))+0.5,labels=ylabels)   
  }
  xaxis=pretty(c(floor(start2/bin)*bin,ceiling(end2/bin)*bin))
  xlabels=paste(xaxis/1000000,"Mb")
  axis(1, at = xaxis/bin-(floor(start2/bin))+0.5,labels=xlabels)
  bottom<-tail(range1,1)
  left<-range2[1]
  if(drawTAD){
    vec<-subset(domainsDF,from.id<max(tail(range1,1),tail(range2,1))&to.id>min(range1[1],range2[1]))
    if(is.null(TADline)){
      apply(vec,1,function(x){
        rect(as.numeric(x["from.id"])-left+.5,bottom-as.numeric(x["to.id"])+.5,
             as.numeric(x["to.id"])-left+1.5,bottom-as.numeric(x["from.id"])+1.5,
             border=TADlinecolor,lwd=100/max(length(range1),length(range2)))
      })
    } else {
      apply(vec,1,function(x){
        rect(as.numeric(x["from.id"])-left+.5,bottom-as.numeric(x["to.id"])+.5,
             as.numeric(x["to.id"])-left+1.5,bottom-as.numeric(x["from.id"])+1.5,
             border=TADlinecolor,lwd=TADline)
      })
    }
    
  }
  if(start==start2 && end==end2){
    polygon(x=c(0.5,0.5,nrow(m)+0.5),y=c(0.5,nrow(m)+0.5,0.5),col="white",border = "white")
  }
  if(!is.null(annotation)){
    newStart<-sapply(annotation[,2],function(x){0.5+(x-floor(start2/bin)*bin)/bin})
    newEnd<-sapply(annotation[,3],function(x){0.5+(x-floor(start2/bin)*bin)/bin})
    if(start==start2 && end==end2){
      if(shift==0){
        for (i in 1:nrow(annotation)){
          segments(x0=newStart[i], y0=nrow(m)+1-newStart[i], x1 = newEnd[i],
                   y1 = nrow(m)+1-newEnd[i],col = "black",lwd=3,lend=1)
        }
      } else{
        for (i in 1:nrow(annotation)){
          if(i==1 | all(newStart[i]>newEnd[1:(i-1)])){
            polygon(x=c(newStart[i]-2*shift,newStart[i]-shift,newEnd[i]-shift,newEnd[i]-2*shift),
                    y=c(nrow(m)+1-newStart[i]-2*shift,nrow(m)+1-newStart[i]-shift,nrow(m)+1-newEnd[i]-shift,nrow(m)+1-newEnd[i]-2*shift),
                    border=NULL,col="black",lty=0)
            
          } else{
            polygon(x=c(newStart[i]-3*shift,newStart[i]-2*shift,newEnd[i]-2*shift,newEnd[i]-3*shift),
                    y=c(nrow(m)+1-newStart[i]-3*shift,nrow(m)+1-newStart[i]-2*shift,nrow(m)+1-newEnd[i]-2*shift,nrow(m)+1-newEnd[i]-3*shift),
                    border=NULL,col="black",lty=0)
            
          }
        }
      }
    } else {
      par(xpd=TRUE)
      newStartY<-sapply(annotation[,2],function(x){0.5+(x-floor(start/bin)*bin)/bin})
      newEndY<-sapply(annotation[,3],function(x){0.5+(x-floor(start/bin)*bin)/bin})
      for (i in 1:nrow(annotation)){
        if(newStart[i]<ncol(m) & newEnd[i]>0){
          segments(x0=max(0,newStart[i]), y0=nrow(m)+1, x1 = min(newEnd[i],ncol(m)+1),col = "black",lwd=3,lend=1)
        }
        if(newStartY[i]<nrow(m) & newEndY[i]>0){
          segments(x0=ncol(m)+1, y0=min(nrow(m)+1,nrow(m)+1-newStartY[i]), y1 = max(0,nrow(m)+1-newEndY[i]),col = "black",lwd=3,lend=1)
        }
      }
      par(xpd=FALSE)
    }
  }
  # TODO: check if it is correct, if we could add a margin
  par(mai=c(0,0,0,0))
  if(operation==""){
    image.scale(matrix(c(0,threshold)),col=colp(99))
  } else if(wholeScaleBar & exists("minr")){
    image.scale(matrix(c(minr,maxr)), col=colorRampPalette(c("blue","white","red"))(99))
  } else {
    image.scale(m, col=colorForMat)
  }
  box()
  par(mar=c(2,2,2,2)+0.1)
  par(mfrow=c(1,1))
  if(giveM){
    return(copym)
  }
}

drawTADs<-bigFunctionDrawTADs

drawdifTADs<-function(start,end,bin,matrixDF1,matrixDF2,domainsDF=NULL,giveM=F,drawTAD=F,title="",
                      normalize=F,maxim=NULL,forceRange=NULL,annotation=NULL,shift=0,extension1="",extension2="",TADline=NULL,wholeScaleBar=F){
  bigFunctionDrawTADs(start=start,end=end,bin=bin,matrixDF=matrixDF1,domainsDF=domainsDF,giveM=giveM,drawTAD=drawTAD,
                      title=title,threshold=threshold,annotation=annotation,shift=shift,extension=extension1,
                      TADline=TADline,start2=NULL,end2=NULL,matrixDF2=matrixDF2, extension2=extension2,
                      operation="dif",maxim=maxim,forceRange=forceRange,wholeScaleBar=wholeScaleBar,normalize=normalize)
}

drawl2rTADs<-function(start,end,bin,matrixDF1,matrixDF2,domainsDF=NULL,giveM=F,drawTAD=F,title="",
                      normalize=F,maxim=NULL,forceRange=NULL,annotation=NULL,shift=0,extension1="",extension2="",TADline=NULL,wholeScaleBar=F){
  bigFunctionDrawTADs(start=start,end=end,bin=bin,matrixDF=matrixDF1,domainsDF=domainsDF,giveM=giveM,drawTAD=drawTAD,
                      title=title,threshold=threshold,annotation=annotation,shift=shift,extension=extension1,
                      TADline=TADline,start2=NULL,end2=NULL,matrixDF2=matrixDF2, extension2=extension2,
                      operation="l2r",maxim=maxim,forceRange=forceRange,wholeScaleBar=wholeScaleBar,normalize=normalize)  
}

#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "axis.pos" argument
#defines the side of the axis. The "add.axis" argument defines
#whether the axis is added (default: TRUE)or not (FALSE).
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}

#This function take a vector with rownames or colnames are return different parameters in a list
#It assumes the names are : id|assembly|chr:start-end
getPara<-function(vecOfNames){
  splitNames<-strsplit(vecOfNames,"\\|")
  if(length(splitNames[[1]])==1){
    splitNames<-splitNames[2:length(splitNames)]
  }
  if(!all(lapply(splitNames,length)==3)){
    return(NULL)
  }
  assembly<-splitNames[[1]][2]
  startPos<-splitNames[[1]][3]
  chrName<-strsplit(startPos,":")[[1]][1]
  startPos<-as.numeric(strsplit(strsplit(startPos,":")[[1]][2],"-")[[1]][1])
  nextStartPos<-as.numeric(strsplit(strsplit(splitNames[[2]][3],":")[[1]][2],"-")[[1]][1])
  endPos<-as.numeric(strsplit(strsplit(splitNames[[length(splitNames)]][3],":")[[1]][2],"-")[[1]][2])
  return(list(chr=chrName,bin=nextStartPos-startPos,start=startPos,end=endPos,assembly=assembly))
}

library(tools)
#This function take a vector with file names and return both matrices and rownames/colnames:
getMatAndMeta<-function(filesValues){
  #First check if all files are here
  for (i in 1:length(filesValues)){
    if(!file.exists(filesValues[i])){
      stop(paste("This file does not exist:",filesValues[i]))
    }
  }
  #Then download the matrix and arrange the size
  m<-list()
  rownameOfMat<-list()
  colnameOfMat<-list()
  for (i in 1:length(filesValues)){
    cat("File",filesValues[i],"\n")
    if (file_ext(filesValues[i])=="gz"){
      m[[i]]<-read.delim(gzfile(filesValues[i]),h=F,comment.char = "#",na.strings = c("NA","nan","Na","na","Nan"))
    }else{
      m[[i]]<-read.delim(filesValues[i],h=F,comment.char = "#",na.strings = c("NA","nan","Na","na","Nan"))
    }
    cat("Downloaded.\n")
    m[[i]]<-tryCatch(apply(m[[i]],2,as.numeric),warning = function(w) m[[i]])
    while(!(all(is.numeric(m[[i]])) & nrow(m[[i]])==ncol(m[[i]]))){
      if(ncol(m[[i]])>nrow(m[[i]])){
        rownameOfMat[[i]]<-m[[i]][,1]
        m[[i]]<-data.frame(m[[i]][,-c(1)])
      } else if(ncol(m[[i]])<nrow(m[[i]])){
        colnameOfMat[[i]]<-as.vector(t(m[[i]])[,1])
        m[[i]]<-data.frame(m[[i]][-c(1),])
      } else{
        rownameOfMat[[i]]<-m[[i]][,1]
        colnameOfMat[[i]]<-as.vector(t(m[[i]])[,1])
        m[[i]]<-data.frame(m[[i]][-c(1),-c(1)])
      }
      m[[i]]<-tryCatch(apply(m[[i]],2,as.numeric),warning = function(w) m[[i]])
    }
    cat(paste("Adjusted to",nrow(m[[i]]),"rows and columns.\n"))
  }
  return(list(matrices=m,rownameOfMat=rownameOfMat,colnameOfMat=colnameOfMat))
}

getSubMatAndRanges<-function(matrix,start,end,bin,start2=NULL,end2=NULL,extension=""){
  range1<-(floor(start/bin)+1):ceiling(end/bin)
  if(is.null(start2)){
    start2=start
  }
  if(is.null(end2)){
    end2=end
  }  
  range2<-(floor(start2/bin)+1):ceiling(end2/bin)
  if(extension==""){
    return(list(matrix[range1,range2],range1,range2))
  } else{
    borders<-strsplit(extension,"__")[[1]][2]
    borders<-as.numeric(strsplit(borders,"-")[[1]])*1000000
    #To be able to do the intersection
    bordersIndex<-c(round(borders[1]/bin),ceiling(borders[2]/bin))
    range1<-intersect(range1,bordersIndex[1]:bordersIndex[2])
    range2<-intersect(range2,bordersIndex[1]:bordersIndex[2])
    return(list(matrix[range1-bordersIndex[1],range2-bordersIndex[1]],range1,range2))
  }
}

startFromRange<-function(range,bin){
  return((range[1]-1)*bin)
}
endFromRange<-function(range,bin){
  return(tail(range,1)*bin)
}