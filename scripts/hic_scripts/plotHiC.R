options(scipen=999)
options(stringsAsFactors=F)
rm(list=ls())
if(! "devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

library(tools)
#Install the libraries needed if they are not already present
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("plotrix")
if(length(commandArgs(TRUE))>0){
  f<-commandArgs(TRUE)[1]
} else{
  #Ask for the config file
  f <- file.choose()
}

#Check the necessary options
if(!file.exists(f)){
  stop("This file does not exist.")
}
source(f)
if(!exists("coordinatesToPlot")){
  stop("The config file do not have start definition.")
} else{
  coordinatesToPlotDF<-as.data.frame(matrix(unlist(lapply(coordinatesToPlot,function(s){v<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(s),":")),"-")),"__"))
  if(length(v)==3){return(c(v,v))}else if(length(v)==6){return(v)}else{stop("Invalid coordinate, expected chrN:start-end or chrN1:start1-end1__chrN2:start2-end2")}
  })),ncol=6,byrow = T))
  colnames(coordinatesToPlotDF)<-c("chr1","start1","end1","chr2","start2","end2")
  coordinatesToPlotDF[,c(2,3,5,6)]<-apply(coordinatesToPlotDF[,c(2,3,5,6)],2,gsub,pattern=",",replacement="")
  coordinatesToPlotDF[,c(2,3,5,6)]<-apply(coordinatesToPlotDF[,c(2,3,5,6)],2,as.numeric)
  if(!all(c(coordinatesToPlotDF$chr1,coordinatesToPlotDF$chr2)==coordinatesToPlotDF$chr1[1])){
    stop("For the moment it is not possible to plot 2 regions with different chromosomes.")
  }
}
if(exists("useColNames")){
  if(!is.logical(useColNames)){
    cat("useColNames is not boolean and is not taken into account.\n The col names will not be used.\n")
    useColNames<-F
  }
} else{
  useColNames<-F
}

if(!exists("bin") & !useColNames){
  stop("The config file do not have bin size definition.")
}
if(!exists("pathForFunctionsHiC")){
  stop("No path for functions.")
} else if(!file.exists(pathForFunctionsHiC)){
  stop(paste("This file does not exists:",pathForFunctionsHiC))
}
source(pathForFunctionsHiC)
if(!exists("bigFunctionDrawTADs")){
  stop("This version of plotHiC require functionsHiCv11 or older.")
}
files<-setdiff(ls()[grepl("file",tolower(ls()))],c("files","filesValues"))
if(length(files)<1){
  stop("There is no matrix put in the configFile the name of the matrix may contain file")
}

#Download the matrix
cat(paste("There are",length(files),"to download.\n"))
filesValues<-sapply(files,function(n){eval(parse(text=n))})
bigList<-getMatAndMeta(filesValues)
m<-bigList[["matrices"]]
rownameOfMat<-bigList[["rownameOfMat"]]
colnameOfMat<-bigList[["colnameOfMat"]]
nameFiles<- sapply(files,function(f){strsplit(f,"file")[[1]][2]})

if(useColNames){
  chrInNames<-list()
  binsInNames<-list()
  extensionInNames<-list()
  for(i in 1:length(filesValues)){
    rowname<-tryCatch(rownameOfMat[[i]],error= function(e) NULL)
    colname<-tryCatch(colnameOfMat[[i]],error= function(e) NULL)
    if(is.null(rowname) & is.null(colname)){
      useColNames<-F
    } else if(! is.null(rowname)){
      para1<-getPara(rowname)
      if(is.null(colname)){
        para<-para1
      } else{
        para2<-getPara(colname)
        if(all(unlist(lapply(1:length(para1),function(j){para1[[j]]==para2[[j]]})))){
          para<-para1
        } else{
          cat("For the moment it is not possible to plot matrices with rowbins and colbins different as in ")
          cat(files[i])
          cat(". The information from the header is not taken into account.\n")
          para<-NULL
          useColNames<-F
        }
      }
    } else{
      para<-getPara(colname)
    }
    if(useColNames){
      chrInNames[[i]]<-para[["chr"]]
      binsInNames[[i]]<-para[["bin"]]
      extensionInNames[[i]]<-paste0("__",para[["start"]]/1000000,"-",para[["end"]]/1000000)
    }
  }
  if(useColNames){
    if(all(chrInNames==chrInNames[[1]])){
      chr<-chrInNames[[1]]
    }
    if(all(binsInNames==binsInNames[[1]])){
      bin<-binsInNames[[1]]
    }
  }
}
if(!useColNames & !exists("bin")){
  stop("The config file do not have bin size definition and the matrix "+files+" does have neither rownames nor colnames.")
}
if(!useColNames){
  binsInNames<-as.list(rep(bin,length(m)))
  extensionInNames<-as.list(rep("",length(m)))
}
####Fit the optional parameters
#common to both plots
if(exists("annotBed")){
  annot<-read.delim(annotBed,h=F,nrow=1,comment.char = "#")
  if(ncol(annot)==1){
    annot<-read.delim(annotBed,h=F,skip=1,comment.char = "#")
  } else{
    annot<-read.delim(annotBed,h=F,comment.char = "#")
  }
  ###This may be changed
  if(length(unique(annot[,1]))>1 & !exists("chr")){
    cat("There is more than one chromosome in the annotation bed file and there is no chromosome defined.\nThe annotations will not be plotted.\n")
    annot<-NULL
  } else{
    if(exists("chr")){
      annot<-subset(annot,V1%in%c(chr))
      if(nrow(annot)<1){
        annot<-NULL
      }
    }
  }
  ####
  if(!is.null(annot)){
    if(ncol(annot)<3 | !all(is.numeric(c(annot$V2,annot$V3))) ){
      cat("This is not a correct annotation bed file. The annotations will not be plotted.\n")
      annot<-NULL
    }
  }
} else{
  annot<-NULL
}

if(!is.null(annot)){
  if(exists("shift")){
    if(!is.numeric(shift)){
      cat("The shift put in the config file is not numeric.\n")
      shift<-1
    }
  } else{
    shift<-1
  }
} else{
  shift<-1
}

if(exists("domainsbed")){
  domains<-read.delim(domainsbed,h=F,nrow=1,comment.char = "#")
  if(ncol(domains)==1){
    domains<-read.delim(domainsbed,h=F,skip=1,comment.char = "#")
  } else{
    domains<-read.delim(domainsbed,h=F,comment.char = "#")
  }
  if(length(unique(domains[,1]))>1 & !exists("chr")){
    cat("There is more than one chromosome in the domains bed file and there is no chromosome defined.\nThe domains will not be plotted.\n")
    domains<-NULL
  } else{
    if(exists("chr")){
      domains<-subset(domains,V1%in%c(chr))
    }
  }
  if(!is.null(domains)){
    if(ncol(domains)<3 | !all(is.numeric(c(domains$V2,domains$V3))) ){
      cat("This is not a correct annotation bed file. The annotations will not be plotted.\n")
      domains<-NULL
    } else{
      if(ncol(domains)>=4 ){
        if(sum(domains$V4%in%c("domain"))>0){
          domains<-subset(domains,V4=="domain")
        }
      }
      domains$from.id<-sapply(domains$V2,function(start){1+floor(start/bin)})
      domains$to.id<-sapply(domains$V3,function(end){ceiling(end/bin)})
    }
  }
} else{
  domains<-NULL
}

if(!is.null(domains)){
  if(exists("TADline")){
    if(!is.numeric(TADline)){
      cat("The TADline put in the config file is not numeric.\n")
      TADline<-NULL
    }
  } else {
    TADline<-NULL
  }
  if(exists("TADlinecolor")){
    if(!isValidColor(TADlinecolor)){
      cat("The color specified as TADlinecolor cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThe TAD lines will be blue.\n")
      TADlinecolor<-"blue"
    }
  } else{
    TADlinecolor<-"blue"
  }
  if(exists("plotWithAndWithoutTADs")){
    if(!is.logical(plotWithAndWithoutTADs)){
      cat("plotWithAndWithoutTADs is not boolean and is not taken into account.\n The plots will not be drawn without TADs.\n")
      plotWithAndWithoutTADs<-F
    }
  } else{
    plotWithAndWithoutTADs<-F
  }
} else{
  TADline<-NULL
  TADlinecolor<-"blue"
  plotWithAndWithoutTADs<-F
}

if(exists("usePng")){
  if(!is.logical(usePng)){
    cat("usePng is not boolean and is not taken into account.\n The output will be pdf.\n")
    usePng<-F
  }
} else{
  usePng<-F
}

if(usePng){
  if(! exists("pngRes")){
    cat("You need to specify a resolution for png files.(pngRes<-96 or whatever).\n By default it will be 96")
    pngRes<-96
  } else{
    if(! is.numeric(pngRes)){
      cat("The resolution for png files need to be a number, for example pngRes<-96 or whatever.\n By default it will be 96")
      pngRes<-96
    }
  }
}


if(!exists("outputPath")){
  outputPath<-paste0(dirname(filesValues[1]),"/",gsub(" ","_",Sys.time()))
  cat("The plots will be in :")
  cat(outputPath)
  cat(".pdf\n")
} else{
  if(!dir.exists(dirname(outputPath))){
    dir.create(dirname(outputPath))
  }
}

if(exists("colForNA")){
  if(!isValidColor(colForNA)){
    cat("The color specified as colForNA cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\nThe NA bins will be grey.\n")
    colForNA<-"grey"
  }
} else{
  colForNA<-"grey"
}

#For the classical plots
if(exists("threshold")){
  if(!is.numeric(threshold)){
    cat("The threshold put in the config file is not numeric.\n")
    threshold<-NULL
  }
} else{
  threshold<-NULL
}

#For the operation plots
allOp<-ls()[grepl("operation",tolower(ls()))]
summaryAllValidOp<-NULL
if(length(allOp)>0){
  for(k in 1:length(allOp)){
    currentValueOfOp<-eval(parse(text=allOp[k]))
    #Check if it is a difference:
    mNames<-strsplit(currentValueOfOp,"-")[[1]]
    if(length(mNames)==2){
      currentOp<-"dif"
      currentName<-allOp[k]
      currentName1<-mNames[1]
      currentName2<-mNames[2]
      currentNameForPngFile<-gsub("file","",currentValueOfOp)
      currentTitle<-gsub("file","",paste0("Substraction ",mNames[1],"-",mNames[2],"\nBlue means higher in ",mNames[2]))
    } else{
      #Check if it is a ratio:
      mNames<-strsplit(currentValueOfOp,"/")[[1]]
      if(length(mNames)==2){
        currentOp<-"l2r"
        currentName<-allOp[k]
        currentName1<-mNames[1]
        currentName2<-mNames[2]
        currentNameForPngFile<-gsub("file","",paste0("log2",mNames[1],"over",mNames[2]))
        currentTitle<-gsub("file","",paste0("Log 2 of ",mNames[1],"/",mNames[2],"\nBlue means higher in ",mNames[2]))
      } else{
        cat(paste0("Cannot plot ",currentName,"=",currentValueOfOp," because it is neither of type file1-file2 nor file1/file2.\n"))
        next
      }
    }
    currentFile1<-which(files==currentName1)
    if(length(currentFile1)==0){
      cat(paste0("Cannot plot ",currentName,"=",currentValueOfOp," because ",currentName1," is not defined.\n"))
      next       
    }
    currentFile2<-which(files==currentName2)
    if(length(currentFile2)==0){
      cat(paste0("Cannot plot ",currentName,"=",currentValueOfOp," because ",currentName2," is not defined.\n"))
      next       
    }
    if(!exists("bin")){
      #I need to check that the bins are compatible
      if(binsInNames[[currentFile1]]!=binsInNames[[currentFile2]]){
        cat(paste0("Cannot plot ",currentName,"=",currentValueOfOp," because the files are not with the same bin size.\n"))
        next       
      }
    }
    if(!exists("chr")){
      #I need to check that the chrs are compatible
      if(chrInNames[[currentFile1]]!=chrInNames[[currentFile2]]){
        cat(paste0("Cannot plot ",currentName,"=",currentValueOfOp," because the files are not with the same chr.\n"))
        next       
      }
    }
    #Everything have been checked:
    summaryAllValidOp<-rbind(summaryAllValidOp,c(currentOp,currentFile1,currentFile2,currentNameForPngFile,currentTitle))
  }
  if("dif"%in%summaryAllValidOp[,1]){
    if(exists("maxim")){
      if(!is.numeric(maxim)){
        cat("The maxim put in the config file is not numeric.\n")
        maxim<-NULL
      }
    } else {
      maxim<-NULL
    }
    if(exists("forceRange")){
      if(!is.numeric(forceRange)){
        cat("The forceRange put in the config file is not numeric.\n")
        forceRange<-NULL
      }else if(length(forceRange)!=2){
        cat("The forceRange do not have 2 values. It is necessary to have one minimum and one maximum.\n")
        forceRange<-NULL
      }else if (forceRange[1]>forceRange[2]){
        cat("The first value of forceRange is greater than the second. The values are inverted.\n")
        forceRange<-forceRange[2:1]
      } 
    } else {
      forceRange<-NULL
    }
    if(!is.null(maxim) & !is.null(forceRange)){
      if(maxim<forceRange[2]){
        cat("This is not good to use a maxim smaller than the range of color.\n")
        maxim<-NULL
      }
    }
    if(exists("normDifMat")){
      cat("normDifMat does not exists anymore. It has been merged with normL2RMat to normBeforeDiffOrl2R.\n")
    }
  }
  if("l2r"%in%summaryAllValidOp[,1]){
    if(exists("maximl2R")){
      if(!is.numeric(maximl2R)){
        cat("The maximl2R put in the config file is not numeric.\n")
        maximl2R<-NULL
      }
    } else {
      maximl2R<-NULL
    }
    if(exists("forceRangel2R")){
      if(!is.numeric(forceRangel2R)){
        cat("The forceRangel2R put in the config file is not numeric.\n")
        forceRangel2R<-NULL
      }else if(length(forceRangel2R)!=2){
        cat("The forceRangel2R do not have 2 values. It is necessary to have one minimum and one maximum.\n")
        forceRangel2R<-NULL
      }else if (forceRangel2R[1]>forceRangel2R[2]){
        cat("The first value of forceRangel2R is greater than the second. The values are inverted.\n")
        forceRangel2R<-forceRangel2R[2:1]
      } 
    } else {
      forceRangel2R<-NULL
    }
    if(!is.null(maximl2R) & !is.null(forceRangel2R)){
      if(maximl2R<forceRangel2R[2]){
        cat("This is not good to use a maximl2R smaller than the range of color.\n")
        maximl2R<-NULL
      }
    }
    if(exists("normL2RMat")){
      cat("normL2RMat does not exists anymore. It has been merged with normDifMat to normBeforeDiffOrl2R.\n")
    }
  }
  if(!is.null(summaryAllValidOp)){
    summaryAllValidOp<-as.data.frame(summaryAllValidOp)
    colnames(summaryAllValidOp)<-c("operation","index1","index2","nameForPngFile","title")
    summaryAllValidOp[,grep("index",colnames(summaryAllValidOp))]<-apply(summaryAllValidOp[,grep("index",colnames(summaryAllValidOp))],2,as.numeric)
    if(exists("plotWholeScaleBar")){
      if(!is.logical(plotWholeScaleBar)){
        cat("plotWholeScaleBar is not boolean and is not taken into account.\n The scale bar will not be plotted to the extreme values.\n")
        plotWholeScaleBar<-F
      }
    }else{
      plotWholeScaleBar<-F
    }
    if(exists("normBeforeDiffOrl2R")){
      if(!is.logical(normBeforeDiffOrl2R)){
        cat("normBeforeDiffOrl2R is not boolean and is not taken into account.\n The matrix will be normalized.\n")
        normBeforeDiffOrl2R<-T
      }
    }else{
      normBeforeDiffOrl2R<-T
    }
  }
}

#Do the plot
if(TRUE){
cat("\n\nCommand lines that will be launched\n\n")
cat("if (!usePng){\n")
cat("  dimOfMat<-max(sapply(binsInNames,function(bin){max(apply(coordinatesToPlotDF,1,function(x){max(as.numeric(x[\"end1\"])-as.numeric(x[\"start1\"]),
                                                                                        as.numeric(x[\"end2\"])-as.numeric(x[\"start2\"]))/bin}))}))\n\n")
cat("  pdf(paste0(outputPath,\".pdf\"),width = max(7,0.075*dimOfMat),height = max(7,0.075*dimOfMat+2),title=basename(f))\n")
cat("}\n")
cat("for (i in 1:length(files)){\n")
cat("  for (j in 1:nrow(coordinatesToPlotDF)){\n")
cat("    if (usePng){\n")
cat("      dimOfCurrentMat<-max(coordinatesToPlotDF$end1[j]-coordinatesToPlotDF$start1[j],
                       coordinatesToPlotDF$end2[j]-coordinatesToPlotDF$start2[j])/binsInNames[[i]]\n")
cat("      png(paste0(outputPath,\"_\",nameFiles[files[i]],\"_\",coordinatesToPlotDF$start1[j]/1e6,\"-\",coordinatesToPlotDF$end1[j]/1e6,\"_\",coordinatesToPlotDF$start2[j]/1e6,\"-\",coordinatesToPlotDF$end2[j]/1e6,\".png\"),
        width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units=\"in\",res=pngRes)\n")
cat("    }\n")
cat("\n    bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i]],m[[i]],colForNA=colForNA,domainsDF=domains,
        giveM=F,drawTAD=!is.null(domains),title=nameFiles[files[i]],threshold=threshold,annotation=annot,shift=shift,extension=extensionInNames[[i]],
        TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j])\n")
cat("    if (plotWithAndWithoutTADs){\n")
cat("      if (usePng){\n")
cat("        dev.off()\n")
cat("        png(paste0(outputPath,\"_\",nameFiles[files[i]],\"_\",coordinatesToPlotDF$start1[j]/1e6,\"-\",coordinatesToPlotDF$end1[j]/1e6,\"_\",coordinatesToPlotDF$start2[j]/1e6,\"-\",coordinatesToPlotDF$end2[j]/1e6,\"_withoutTAD.png\"),
          width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units=\"in\",res=pngRes)\n")
cat("      }\n")
cat("\n      bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i]],m[[i]],colForNA=colForNA,domainsDF=domains,
          giveM=F,drawTAD=F,title=nameFiles[files[i]],threshold=threshold,annotation=annot,shift=shift,extension=extensionInNames[[i]],
          TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j])\n")
cat("    }\n")
cat("    if (usePng){\n")
cat("      dev.off()\n")
cat("    }\n")
cat("  }\n}\n")
if(!is.null(summaryAllValidOp)){
  cat("for(k in 1:nrow(summaryAllValidOp)){\n")
  cat("  i1<-summaryAllValidOp$index1[k]\n")
  cat("  i2<-summaryAllValidOp$index2[k]\n")
  cat("  for(j in 1:nrow(coordinatesToPlotDF)){\n")
  cat("    if (usePng){\n")
  cat("      dimOfCurrentMat<-max(coordinatesToPlotDF$end1[j]-coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end2[j]-coordinatesToPlotDF$start2[j])/binsInNames[[i1]]\n")
  cat("      png(paste0(outputPath,\"_\",summaryAllValidOp$nameForPngFile[k],\"_\",coordinatesToPlotDF$start1[j]/1e6,\"-\",coordinatesToPlotDF$end1[j]/1e6,\"_\",coordinatesToPlotDF$start2[j]/1e6,\"-\",coordinatesToPlotDF$end2[j]/1e6,\".png\"),
      width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units=\"in\",res=pngRes)\n")
  cat("    }\n")
  cat("\n    bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i1]],m[[i1]],colForNA=colForNA,matrixDF2=m[[i2]],
      domainsDF=domains,giveM=F,drawTAD=!is.null(domains),title=summaryAllValidOp$title[k],operation=summaryAllValidOp$operation[k],
      normalize=normBeforeDiffOrl2R,maxim=maxim,forceRange=forceRange,maximl2R=maximl2R,forceRangel2R=forceRangel2R,
      annotation=annot,shift=shift,extension=extensionInNames[[i1]],extension2=extensionInNames[[i2]],
      TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j],wholeScaleBar=plotWholeScaleBar)\n")
  cat("    if (plotWithAndWithoutTADs){\n")
  cat("      if (usePng){\n")
  cat("        dev.off()\n")
  cat("        png(paste0(outputPath,\"_\",summaryAllValidOp$nameForPngFile[k],\"_\",coordinatesToPlotDF$start1[j]/1e6,\"-\",coordinatesToPlotDF$end1[j]/1e6,\"_\",coordinatesToPlotDF$start2[j]/1e6,\"-\",coordinatesToPlotDF$end2[j]/1e6,\"_withoutTAD.png\"),
        width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units=\"in\",res=pngRes)\n")
  cat("      }\n")
  cat("\n      bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i1]],m[[i1]],colForNA=colForNA,matrixDF2=m[[i2]],
        domainsDF=domains,giveM=F,drawTAD=F,title=summaryAllValidOp$title[k],operation=summaryAllValidOp$operation[k],
        normalize=normBeforeDiffOrl2R,maxim=maxim,forceRange=forceRange,maximl2R=maximl2R,forceRangel2R=forceRangel2R,
        annotation=annot,shift=shift,extension=extensionInNames[[i1]],extension2=extensionInNames[[i2]],
        TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j],wholeScaleBar=plotWholeScaleBar)\n")
  cat("    }\n")
  cat("    if (usePng){\n")
  cat("      dev.off()\n")
  cat("    }\n")
  cat("  }\n}\n")
}
cat("if (!usePng){\n")
cat("  dev.off()\n")
cat("}\n")
}
if (!usePng){
  dimOfMat<-max(sapply(binsInNames,function(bin){max(apply(coordinatesToPlotDF,1,function(x){max(as.numeric(x["end1"])-as.numeric(x["start1"]),
                                                                                                 as.numeric(x["end2"])-as.numeric(x["start2"]))/bin}))}))
  
  pdf(paste0(outputPath,".pdf"),width = max(7,0.075*dimOfMat),height = max(7,0.075*dimOfMat+2),title=basename(f))
}
for (i in 1:length(files)){
  for (j in 1:nrow(coordinatesToPlotDF)){
    if (usePng){
      dimOfCurrentMat<-max(coordinatesToPlotDF$end1[j]-coordinatesToPlotDF$start1[j],
                           coordinatesToPlotDF$end2[j]-coordinatesToPlotDF$start2[j])/binsInNames[[i]]
      png(paste0(outputPath,"_",nameFiles[files[i]],"_",coordinatesToPlotDF$start1[j]/1e6,"-",coordinatesToPlotDF$end1[j]/1e6,"_",coordinatesToPlotDF$start2[j]/1e6,"-",coordinatesToPlotDF$end2[j]/1e6,".png"),
          width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units="in",res=pngRes)
    }
    
    bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i]],m[[i]],colForNA=colForNA,domainsDF=domains,
                        giveM=F,drawTAD=!is.null(domains),title=nameFiles[files[i]],threshold=threshold,annotation=annot,shift=shift,extension=extensionInNames[[i]],
                        TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j])
    if (plotWithAndWithoutTADs){
      if (usePng){
        dev.off()
        png(paste0(outputPath,"_",nameFiles[files[i]],"_",coordinatesToPlotDF$start1[j]/1e6,"-",coordinatesToPlotDF$end1[j]/1e6,"_",coordinatesToPlotDF$start2[j]/1e6,"-",coordinatesToPlotDF$end2[j]/1e6,"_withoutTAD.png"),
            width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units="in",res=pngRes)
      }
      
      bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i]],m[[i]],colForNA=colForNA,domainsDF=domains,
                          giveM=F,drawTAD=F,title=nameFiles[files[i]],threshold=threshold,annotation=annot,shift=shift,extension=extensionInNames[[i]],
                          TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j])
    }
    if (usePng){
      dev.off()
    }
  }
}
if(!is.null(summaryAllValidOp)){
  for(k in 1:nrow(summaryAllValidOp)){
    i1<-summaryAllValidOp$index1[k]
    i2<-summaryAllValidOp$index2[k]
    for(j in 1:nrow(coordinatesToPlotDF)){
      if (usePng){
        dimOfCurrentMat<-max(coordinatesToPlotDF$end1[j]-coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end2[j]-coordinatesToPlotDF$start2[j])/binsInNames[[i1]]
        png(paste0(outputPath,"_",summaryAllValidOp$nameForPngFile[k],"_",coordinatesToPlotDF$start1[j]/1e6,"-",coordinatesToPlotDF$end1[j]/1e6,"_",coordinatesToPlotDF$start2[j]/1e6,"-",coordinatesToPlotDF$end2[j]/1e6,".png"),
            width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units="in",res=pngRes)
      }
      
      bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i1]],m[[i1]],colForNA=colForNA,matrixDF2=m[[i2]],
                          domainsDF=domains,giveM=F,drawTAD=!is.null(domains),title=summaryAllValidOp$title[k],operation=summaryAllValidOp$operation[k],
                          normalize=normBeforeDiffOrl2R,maxim=maxim,forceRange=forceRange,maximl2R=maximl2R,forceRangel2R=forceRangel2R,
                          annotation=annot,shift=shift,extension=extensionInNames[[i1]],extension2=extensionInNames[[i2]],
                          TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j],wholeScaleBar=plotWholeScaleBar)
      if (plotWithAndWithoutTADs){
        if (usePng){
          dev.off()
          png(paste0(outputPath,"_",summaryAllValidOp$nameForPngFile[k],"_",coordinatesToPlotDF$start1[j]/1e6,"-",coordinatesToPlotDF$end1[j]/1e6,"_",coordinatesToPlotDF$start2[j]/1e6,"-",coordinatesToPlotDF$end2[j]/1e6,"_withoutTAD.png"),
              width = max(7,0.075*dimOfCurrentMat),height = max(7,0.075*dimOfCurrentMat+2),units="in",res=pngRes)
        }
        
        bigFunctionDrawTADs(coordinatesToPlotDF$start1[j],coordinatesToPlotDF$end1[j],binsInNames[[i1]],m[[i1]],colForNA=colForNA,matrixDF2=m[[i2]],
                            domainsDF=domains,giveM=F,drawTAD=F,title=summaryAllValidOp$title[k],operation=summaryAllValidOp$operation[k],
                            normalize=normBeforeDiffOrl2R,maxim=maxim,forceRange=forceRange,maximl2R=maximl2R,forceRangel2R=forceRangel2R,
                            annotation=annot,shift=shift,extension=extensionInNames[[i1]],extension2=extensionInNames[[i2]],
                            TADline=TADline,TADlinecolor=TADlinecolor,start2=coordinatesToPlotDF$start2[j],end2=coordinatesToPlotDF$end2[j],wholeScaleBar=plotWholeScaleBar)
      }
      if (usePng){
        dev.off()
      }
    }
  }
}
if (!usePng){
  dev.off()
}

