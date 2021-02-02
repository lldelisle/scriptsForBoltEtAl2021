###Put the name of the matrices you want to ues in your plots
###The name MUST begin by file
###Every line which have a # in front of it will be ignored
filewtDFL="E12.5_wt_DFL_CHIC.txt"
fileinv2DFL="E12.5_inv2_DFL_CHIC.txt"

###Required parameters
coordinatesToPlot=c("chr2:73900700-74703000__chr2:74638000-74783000")#If only one region is provided, this means you want the classical triangle. If you want to see specific interactions between 2 regions, you need to separate them by __
pathForFunctionsHiC="/home/ldelisle/softwares/scriptsForBoltEtAl2021/scripts/hic_scripts/functionsHiC.R"

###For the bin size and the chromosome
###Either you use the rownames or colnames if they exists:
useColNames=T
###Or you define them explicitely this means that you have the whole chromosome
#chr="chr2"
#bin=40000

###If you want to compute differences/log2ratio:
###The name of the operation MUST begin by operation
operation1="filewtDFL-fileinv2DFL" #To compute a difference put the name of the file as defined above separated by -
#operation2="fileNPC-fileES" #To compute a difference put the name of the file as defined above separated by -
#operation3="fileCN/fileES" #To compute a l2r put the name of the file as defined above separated by /


### Optional parameters (Uncomment the lines you want)
##Common for all plots
#Output location
outputPath="/scratch/ldelisle/BoltPlot/figure2Sb_heatmap" #The output file (without extension)
#Output format
usePng=F #If you want to use png replace F by T
pngRes=300 #This is the resolution of the png file.
#NA color
colForNA="grey" #If you want to change the color for NA values, put another color

#Features
annotBed="/home/ldelisle/softwares/scriptsForBoltEtAl2021/scripts/annotations/HoxD_Elements_mm10_Minimal.bed" # If you want to have black box under the plot with the annotations
shift=5 #the distance between annotations and the matrix 
#domainsbed="/home/ldelisle/Dropbox/scripts/toShare/drawMatrix/demo/es_10kb_ws24.bed" # If you want to have the TADs (for example) ploted as blue triangles.
#TADline=5 #If you want to change manually the width of the line drawing the TADs
#TADlinecolor="black" #Default is blue
#plotWithAndWithoutTADs=T # If you want to have both heatmaps with TADs and without.


#For the plot of files (classical heatmaps)
threshold=0.01 #if you want to saturate the plot to have more contrast decrease this value.

#Common for diff and l2r
plotWholeScaleBar=T #If you want to see in the scale bar both values of forceRange, even if they are not reached.
normBeforeDiffOrl2R=F #If you do not want to normalize the matrix before doing the difference or the l2R (for example when you are using matrices normalized the same way) remove the # at the beginning of the line.

#For the plot of diff
forceRange=c(-0.005,0.005) #if you want to saturate the plot or have the same scale for all plots. This forces the values for Blue and Red.

#For the plot of l2r
#forceRangel2R=c(-3,3) #if you want to saturate the plot or have the same scale for all plots. This forces the values for Blue and Red.

