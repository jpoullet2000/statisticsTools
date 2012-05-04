## @input groups-samples file, feature-length file, relative abundance file, comparisons
## @param abundFile string abundance file
## @param featureLengthFile string feature length file 
## @param groupSampFile string groups-samples file
## @param compar string comparisons to be made
## @output 

## ---------------------------------------------------------------
## ABUNDANCE TAB
## ---------------------------------------------------------------
## FeatureID S1   S2   S3
## CDS1     2.3  3.9  3.8
## ...
## CDSN     3.8  2.1  3.0


## groupSampFile is a matrix (rows = groups, columns = samples)
## ---------------------------------------------------------------
## GROUP TAB
## ---------------------------------------------------------------
## GroupID  S1   S2   S3   S4
## Group1   1    0    0    0 
## Group2   0    1    1    0
## Group3   0    0    0    1

## Comparisons should be of the form:
## G1-G3;G2-G6-G3 (comparisons of samples of G1 vs samples of G3 ; and samples of G2 vs G6 vs G3) => ";-," cannot be used in group names  
## RMK: Multiple groups comparisons is not yet implemented. If more than 2 groups are given, all combinations of 2 groups are taken as comparisons  


## @usage R CMD BATCH '--args groupSampFile abundFile featureLengthFile compar methodName' differentialExpressionAnalysis.R


## ------ ---------------- FUNCTIONS ------------------------------------------------
## @brief pairwiseComparisonsDifferentialExpression
## dataObject ===  list data object with the following structure:
## dataObject -- "groupLabels" -- string list of group labels
##            -- "groupSamples" -- group1 -- sample list (ex:c(sample1,sample5))
##                              -- ...
##                              -- groupN -- sample list (ex:c(sample2,sample8))
##            -- "data"         -- data.frame (rows = features, column = samples)
pairCompDiffExp <- function(
                            dataObject,
                            methodName, # "DESeq","KW" 
                            withplot = 0 # generate plots
                            ){
  sampleNames <- colnames(dataObject$"data")
  conds <- c()
  resAll<- c()
  ## Get allcomb samples from groups
  for (nam in sampleNames){
    print(nam)
    for (group in dataObject$"groupLabels"){
      if (nam %in% dataObject$"groupSamples"[[group]]){
        conds <- c(conds,group)
        break
      }
    }
  }
  if (length(conds)!=length(sampleNames)){
    stop("Some of the sample names in the group table are not found in the abundance table")
  }
  ## going through all pairwise combinations
  mystring<-dataObject$"groupLabels"
  allpaircomb<-unlist(sapply(mystring[-length(mystring)],function(x) paste(x,mystring[(grep(x,mystring)+1):length(mystring)],sep=",")))
  ## check if biological replicates exist
  for (group in dataObject$"groupLabels"){
    if (length(dataObject$"groupSamples"[[group]]) == 1){
      flagBlind <- 1
      break
    }else{
      flagBlind <- 0
    }
  }
  for (pair in allpaircomb){
    if (methodName=="DESeq"){
      g <- strsplit(pair,",")[[1]]
      cds <- newCountDataSet(dataObject$"data",conds)
      cds <- estimateSizeFactors(cds)
      ## cds <- estimateVarianceFunctions(cds)
      if (flagBlind == 1){
        #cds <- estimateVarianceFunctions(cds,method="blind")#,locfit_extra_args=list(maxk=200))
	cds <- estimateDispersions(cds,method="blind",sharingMode="fit-only",fitType="local")
		
#cds <- estimateDispersions(cds,method="pooled")
      }else{
        #cds <- estimateVarianceFunctions(cds)
	cds <- estimateDispersions(cds)
      }
      res <- nbinomTest(cds,g[1],g[2])
      resAll<- c(resAll,res)
      ## Diagnostic plot to check the fit of the variance functions 
      if (withplot){
        pdf(paste("checkVarFitPlot",".pdf",sep="",collapse=""))
        for (group in dataObject$"groupLabels"){
          diagForGroup<-varianceFitDiagnostics(cds,group)
          smoothScatter(log10(diagForGroup$baseMean),log10(diagForGroup$baseVar))
          lines(log10(fittedBaseVar)~log10(baseMean),diagForGroup[order(diagForGroup$baseMean),],col="red")
          title(paste("Fit of variance - group ",group,sep="",collapse=""))
          dev.off()
        }
      }
    }else if (methodName == "KW"){
      print(conds)
      print(colnames(dataObject$"data"))
      print(dataObject$"data" %in% conds)
      myd <- dataObject$"data"#[,which(colnames(dataObject$"data") %in% conds)]
      ## colnames(myd)[which(colnames(myd)=="Deinococcus-Thermus_uncertain")]<-"Deinococcus_Thermus_uncertain"
      ## colnames(myd)[which(colnames(myd)=="Deinococcus-Thermus_uncertain")]<-"Deinococcus_Thermus_uncertain"
      myd <- rbind(myd,conds)
      rownames(myd)[nrow(myd)]<-"group"
      myd<-t(myd)
      ## Check whether the first character is an integer and change the column name if this is the case
      nbrList <- as.character(0:9)
      if (substring(colnames(myd)[1],1,1) %in% nbrList){
        colnames(myd)[1:(ncol(myd)-1)]<-unlist(lapply(colnames(myd)[1:(ncol(myd)-1)],function(x) paste("C",x,sep="",collapse="")))
      }
      resAll<-c()
      resAll$names<-c()
      resAll$pvalue<-c()
      names<-c()
      pvalue<-c()
      for (i in 1:ncol(myd)){
        print(i)
        res<-kruskal.test(eval(parse(text=colnames(myd)[i])) ~ group,data=myd)
        #res<-kruskal.test(colnames(myd)[i] ~ group,data=myd)
        res$names<-colnames(myd)[i]
        res$pvalue<-res$p.value
        #resAll <- c(resAll,res)
        names <- c(names,res$names)
        pvalue <- c(pvalue,res$pvalue)
        ## resAll$names <- c(resAll$names,res$names)
        ## resAll$pvalue <- c(resAll$pvalue,res$pvalue)
      }
      print(names)
      print(pvalue)
      
      resAll$names <- names
      resAll$pvalue <- pvalue
    }
  }
  ret <- list(resAll=resAll, allpaircomb=allpaircomb)
}
## ----------------------------------------------------------------------------------

# @brief Generate table and figure for the DESeq method
DESeqGenFigAndTab <- function(
                              com,#comparison list
                              resAll,#result object 
                              padj, #p-value
                              allpaircomb, # all pair combinations
                              dirout # output directory 
                              ){
  bitmap(file.path(dirout,paste("meanVSvarPlot-",com[i],".pdf",sep="",collapse="")),type = "png256")
  print(length(resAll))
  print(class(resAll))
  for (ir in 1:length(resAll)){
    plot(
         resAll[[ir]]$baseMean,
         resAll[[ir]]$log2FoldChange,
         log="x",pch=20,cex=0.1,
         col = ifelse( resAll[[ir]]$padj < padj, "red", "black" ) ## 
         )
    title(paste("Pair ",allpaircomb[ir],sep="",collapse=""))
    dev.off()
  }
  
  
  ## writing results
  for (ir in 1:length(resAll)){
    #write.table(resAll[[ir]],paste("diffExp-",com[i],"-",allpaircomb[ir],".csv",sep="",collapse=""))
    write.table(resAll[[ir]],file.path(dirout,paste("diffExp-",com[i],".csv",sep="",collapse="")),row.names=F)
    dat <- resAll[[ir]]
    newTab <- dat[which(dat[,"padj"]<padj),c("id","padj")]
    write.table(newTab,file.path(dirout,paste("summ",com[i],".csv",sep="",collapse="")),row.names=F)
  }
}

## ----------------------------------------------------------------------------------
# @brief Generate table and figure for the DESeq method
KWGenFigAndTab <- function(
                   com,#comparison list
                   resAll,#result object 
                   padj, #p-value
                   allpaircomb, # all pair combinations
                   outdir # output directory 
                   ){
  for (ir in 1:length(resAll)){
    ## toW<-c()
    ## for (i in nrow(resAll[[ir]])){
    ##   toW<-rbind(toW,cbind(resAll[[ir]]$names,resAll[[ir]]$pvalue))
    ##   print(toW)
    ## }
    write.table(resAll[[ir]],file.path(outdir,paste("diffExp-",com[ir],".csv",sep="",collapse="")),row.names=F)
  }
}



## ----------------------------- MAIN -----------------------------------------------
## loading arguments
args <- commandArgs(TRUE)
abundFile <- args[1]
featureLengthFile <- args[2]
groupSampFile <- args[3]
compar <- args[4]
outdir <- args[5] # output directory 
methodName <- args[6]


## Loading libraries
library(DESeq)

## Method names
metNames = c("DESeq","KW")

## making group comparisons
com <- strsplit(compar,";") # comparisons (ex: c("G1-G3","G2-G6-G3))
nbcom <-length(com) # nb of comparisons

## loading groups-sample table + loading relative abundance table
## nbrTot <- system(paste("head -n 2 ",featureLengthFile,sep="",collapse=""),intern=TRUE) 
## nbrTot <- strsplit(nbrTot[2],"\t")[[1]] #list of total numbers of reads 

abtab <- read.csv(abundFile,header=T,sep="\t",check.names=F)
rownames(abtab) <- abtab[,1] # 1st column contains the rownames
abtab <- abtab[,-1]

grtab <- read.csv(groupSampFile,header=T,sep="\t",check.names=F)
rownames(grtab) <- grtab[,1] # 1st column contains the rownames 
grtab <- grtab[,-1] 


## Statistics
allcomb<-list()

for (i in 1:length(com)){ # going through comparisons
  allcomb[[i]] <- list()
  grps <-strsplit(com[[i]],"-")[[1]] # sample groups (ex: c("G1","G3") )
  allcomb[[i]]$"groupLabels" <- grps
  allcomb[[i]]$"groupSamples" <- list()
  sampList<-c()
  for (grp in grps){ # going through groups
    grpline <- grtab[which(rownames(grtab)==grp),] #
    samps <- colnames(grtab)[which(grpline==1)] #sample list
    allcomb[[i]]$"groupSamples"[[grp]] <- samps
    sampList<-c(sampList,samps)
  }
  allcomb[[i]]$"data" <- abtab[,which(colnames(abtab) %in% sampList)]
  res <- pairCompDiffExp(allcomb[[i]],methodName,withplot=0)
  resAll <- res$resAll
  if (length(allcomb)==1){
    tmp<-resAll
    rm(resAll)
    resAll<-list()
    resAll[[1]] <- tmp
  }
  allpaircomb <- res$allpaircomb
  ## Plotting files
  padj = 0.1 ## p-value after correcting for multiple testing using Benjamini-Hochberg procedure which controls false discovery rate (FDR)
  if (methodName == "DESeq"){
    DESeqGenFigAndTab(
                      com,#comparison list
                      resAll,#result object 
                      padj, #p-value
                      allpaircomb, # all pair combinations
                      outdir # output directory 
                      )
  }else if (methodName == "KW"){
    KWGenFigAndTab(
                   com,#comparison list
                   resAll,#result object 
                   padj, #p-value
                   allpaircomb, # all pair combinations
                   outdir # output directory 
                   )
  }else{
    stop(paste("Method name must be in: ",metNames))  
  }
}



