############################################################
##Author: Elmer A. Fernández
## Institution : CIDIE - CONICET - UCC
## Version : 0.1
## Date: 13/04/2020
## Last Changes:
##it is on GitHub
##CHANGES: 
##########################################################


.debug <- TRUE
.version <- 1.0

is.debug <- function(.stop = FALSE){
  if(.debug){
    cat("ESTAMOS EN MODO DEBUG")
    if(.stop) {
      stop()
    }
  }
}

#dependencies
library(e1071)
library(parallel)
#library(preprocessCore)
library(stringr)
#extras
library(nnls)
library(MASS)
#for table manipulation
library(plyr)
library(abind)#for multi-arrays manipulation
#graphics
library(ggplot2)
#excel files handler
library(openxlsx)


####----


#### Main CORE function: This function should not be called by end user
#' .MIXER
#' this is the main core function, it should not be called by users
#' @param expressionMatrix a GxS gene expression matrix. with row names genes ID as in the signature Matrix. Genes in rows, samples in columns
#' @param signatureMatrix a NxL gene signature matrix for L cell types
#' @param functionMixture   : the deconvolution function: cibersort, nu.svm.robust.RFE, ls.rfe.abbas (tiped without quotation marks)
#' @param useCores          : (integer) the number of cpus to use.

.MIXER <- function(expressionMatrix, signatureMatrix, functionMixture, useCores = 1L, verbose =FALSE){
  
  ## OS independent multi-core platform (single machine)
  bp.param <- BiocParallel::bpparam()
  
  .ncores = ifelse(useCores > parallel::detectCores(), parallel::detectCores()-1, useCores)
  
  BiocParallel::bpworkers(bp.param) <- .ncores
  # print(bp.param)
  X <- data.matrix(signatureMatrix)
  Yorig <- data.matrix(expressionMatrix)
  ##No Idea why
  Ymedian <- max(median(Yorig),1)
  
  nx.genes <- nrow(X)
  id.xy <- rownames(X) %in% rownames(Yorig)
  id.yx <- rownames(Yorig) %in% rownames(X)
  
  ##Normalizing the signature Matrix
  ##this follows the normalization by Newman, overall mean and overall sd
  X <- (X - mean(X))/ sd(as.vector(X))
  
  #Select common genes
  X <- X[which(id.xy),]
  Y <- Yorig[which(id.yx),,drop=FALSE]
  X <- data.matrix(X[order(rownames(X)),,drop=FALSE])
  Y <- data.matrix(Y[order(rownames(Y)),,drop=FALSE])
  
  if(verbose){
    cat("--------------\n")
    cat(paste("We found ", nrow(Y), 
              "(", round(100*nrow(Y)/nx.genes,0), "%) genes from the ", nx.genes, " signature genes "))
    cat(paste("\nAnalyzing ", ncol(Y)," subjects \n"))
    cat("---------------")
    
  }
  
  
  ##Normalization of data
  ##truly single subject
  y.median.by.subj <- apply(Y, 2, median)
  Yn <- scale(Y, center = TRUE, scale = TRUE)
  
  ##perform mixture analysis
  
  out <- BiocParallel::bplapply(1:ncol(Yn), function(yi) {
  #  print(paste("id ",yi))
    functionMixture(X,y = as.numeric(Yn[,yi]))
  }  , BPPARAM = bp.param)
  
  ##building output
  mix.mat.abs <- do.call(rbind, lapply(out, function(x) x$Wa))
  rownames(mix.mat.abs) <- colnames(Yorig)
  colnames(mix.mat.abs) <- colnames(X)
  mix.mat.prop <- do.call(rbind, lapply(out, function(x) x$Wp))
  rownames(mix.mat.prop) <- colnames(Yorig)
  colnames(mix.mat.prop) <- colnames(X)
  mat.res <- do.call(rbind,lapply(out, function(x) c( x$RMSEa, x$RMSEp, x$Ra, x$Rp, x$BestParams, x$Iter)))
  colnames(mat.res) <- c("RMSEa", "RMSEp", "Ra", "Rp", "BestParams","Iter")
  rownames(mat.res) <- rownames(mix.mat.prop)
  
  
  return(list(MIXabs = mix.mat.abs, MIXprop = mix.mat.prop, ACCmetric = mat.res))
  
}

###---------
##PermutationNUll---------
##Main function
#' MIXTURE
#' 
#' the MIXTURE algorithm solves Y = X*B estimating the B coefficients or cell type proportions (B>0)
#' @param expressionMatrix a GxS gene expression matrix. with row names genes ID as in the signature Matrix. Genes in rows, samples in columns
#' @param signatureMatrix a NxL gene signature matrix for L cell types
#' @param iter integer, number of iterations to estimate the p values
#' @param functionMixture   : (default nu.svm.robust.RFE) the deconvolution function: cibersort, nu.svm.robust.RFE, ls.rfe.abbas (typed without quotation marks)
#' @param useCores          : (integer) the number of cpus to use.
#' @param verbose           : (logical) if TRUE, msgs on the screen
#' @param nullDist          : (character) select the null distribution mode for p.valu estimation of the regression estimate
#'                            one of the followings : "none" , "PopulationBased" use the whole expressionMatrix to draw "iter" random samples at once. 
#'                                         Then all the models are compared against this null distribution
#'@param filesave         : (character) the name of the output EXCEL(r) file in xlsx format
#'
#'@export
#'
#'@return   an S3 object of class MIXTURE  (list) with the following slots: 
#'Subjects        : a list with the following slots:
#'        MIXabs    : a data.frame with S rows, and L columns with the absoluite regression coefficients
#'        MIXprop   : a data.frame with S rows, and L columns with the proportions
#'        ACCmetric : a data.frame with S rows, and RMSEa (Root mean squared error from absolute coeficientes),
#'        RMSEp (RMSE with proportions coefficientes)
#'        Ra (correlation absolute coeffs)
#'        Rp (correlation proportion coeffs)
#'        BestParams (Best model parameters)Iter (number of itereation of the RFE to reach the best model)
#'        
#'  PermutedMatrix  : the simulated wxpression matrix for null distribution estimation (if nullDist == "none", PermutedMatrix = NULL)
#'  p.values        : The p.values for each subject and for each of the ACC metrix columns (if nullDist == "none", p.value = NULL)
#'  method          : the null distribution method (nullDist)
#'  usedGenes       : the intersection gene list between expressionMatrix and signatureMatrix
#'
#'
MIXTURE <- function(expressionMatrix , signatureMatrix, iter = 100, functionMixture = nu.svm.robust.RFE  , useCores =1L,
                             verbose = FALSE, nullDist = c("none","SingleSubjectTest","PopulationBased"),
                    fileSave){
#This function perform the decovolution of the signatureMatrix over the gene expression subject matrix.

## Returns:
  bp.param <- BiocParallel::bpparam()
  cat(paste0("\nMIXTURE version ",.version))
  .ncores = ifelse(useCores > parallel::detectCores(), parallel::detectCores()-1, useCores)
  
  BiocParallel::bpworkers(bp.param) <- .ncores
  ###bla bla bla
  nullDist <- match.arg(nullDist)
  if(verbose) {
    cat("\nRunning...\n")
  }
  if(missing(signatureMatrix)){
    cat(paste0("\nLoading LM22 signature matrix"))
    data("LM22")
    signatureMatrix <- LM22
  }else{
    cat(paste0("\nProvided signature matrix"))
  }
  cat(paste0(" ",nrow(signatureMatrix)," genes and ", ncol(signatureMatrix)," cell types"))
  .list.of.genes <- rownames(expressionMatrix)[which((rownames(expressionMatrix) %in% rownames(signatureMatrix)))]
  
  if(length(.list.of.genes)<20){
    stop(paste0("ERRORless than 20 genes in common with signatureMatrix"))
  }
  cat(paste0("\nOverlapping genes : ", length(.list.of.genes)," (",round(100*length(.list.of.genes)/nrow(signatureMatrix),2),"%)"))
  
  ##compute the deconvolution onto the original samples
  Orig <- .MIXER(expressionMatrix , signatureMatrix, functionMixture , useCores = .ncores)
  
  ##compute the total "immune content"
  total.median.by.subj <- apply(data.matrix(expressionMatrix), 2, median, na.rm=T)
  total.median.by.subj.lm22genes <- apply(data.matrix(expressionMatrix[.list.of.genes,]), 2, median, na.rm=T)
  max.median.full.cohort <- max(median(data.matrix(expressionMatrix)),1)
  
  Orig$ACCmetric <- cbind(Orig$ACCmetric, IscBySbj = total.median.by.subj.lm22genes/total.median.by.subj, 
                          IscPob = total.median.by.subj.lm22genes/max.median.full.cohort)
  
  if(verbose){
    cat("Original Samples run\n")
  }
  
  if (nullDist == "none") {
    out.list <- list(Subjects = Orig,PermutedMetrix = NULL, p.values = NULL, method = "none", usedGenes = .list.of.genes)#, pvalue = pvalues))   
  }
  #id.iter <- lapply(1:iter, function(x) sample(nrow(expressionMatrix)))
  if (nullDist == "SingleSubjectTest") {
    M.O <- expressionMatrix[rownames(signatureMatrix),]
    out.l <- lapply(1:iter, function(i) {
      if (verbose) {
        cat(paste("%",round(100*i/iter), " / "))  
      }
      nsamp <- sample(nrow(M.O)) 
      M.aux <- M.O[nsamp,]
      rownames(M.aux) <- rownames(M.O)
      return(.MIXER(M.aux , signatureMatrix , functionMixture , useCores = .ncores )$ACCmetric)
    } )
    
    if(verbose) cat("\nfinish\n")
    metrix <- do.call(function(...) abind(along = 3,... ), lapply(out.l, function(x) x))
    
    out.list <- list(Subjects = Orig,PermutedMetrix = metrix, method = "SingleSubjectTest", usedGenes = .list.of.genes)
    #return(list(Subjects = Orig,PermutedMetrix = metrix, method = "SingleSubjectTest"))#, pvalue = pvalues))  
  }
  
  if (nullDist == "PopulationBased") {
    expressionMatrix <- data.matrix(expressionMatrix)
    if(verbose) cat("\nPopulation based null distribution\n")
    if(verbose) cat("\nBuilding random population\n")
        M.aux <- do.call(cbind, BiocParallel::bplapply(1:iter, function(i) {
        as.vector(expressionMatrix)[sample(nrow(expressionMatrix)*ncol(expressionMatrix), size = nrow(expressionMatrix))]
      }, BPPARAM = bp.param  ))
      
  
   rownames(M.aux) <- rownames(expressionMatrix)
   if(verbose) cat("\nBuilding null distribution\n")  
    out.mix <- .MIXER(M.aux , signatureMatrix , functionMixture , useCores = .ncores )$ACCmetric
    if(verbose)  cat("\nfinish\n")
    
    
    out.list <- list(Subjects = Orig,PermutedMetrix = out.mix, method = "PopulationBased", usedGenes = .list.of.genes)
    out.list$p.values <- GetPvalues(out.list)
    #return(list(Subjects = Orig,PermutedMetrix = out.mix, method = "PopulationBased"))#, pvalue = pvalues))  
  }
  ## save to excel file
  if( !missing(fileSave) ){
    SaveExcel(out.list, file = fileSave)
    cat(paste("\n",fileSave,"....OK"))
  }
  class(out.list) <- c("MIXTURE",class(out.list))
  return(out.list)
}

##Access functions----
#' GetMixture
#' returns the regresssion coefficient matrix B as absolute values or as proportion (normalized) values
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @param type a character string "proportion" or "absolute". short word is allowed
#' @export
#' @return a SxK data matrix with K cell types (according to the signature matrix) and S subjects
#' 
GetMixture <- function(obj, type = c("proportion","absolute") ){
##This function returns the mixture coefficientes estimated by MIXTURE
##Args:
  ## obj: the list object returned by MIXTURE
  ## type : (character) the type of the coefficients: "proportion" or "absolute" (default: "proportion")
##Returns:
  ## a SxL coefficient matrix (N: Number of subjects, L: number of cell types in the signature matrix)
  type <- match.arg(type)
  switch(type,
         absolute = matrix.type <- "MIXabs",
         proportion = matrix.type <- "MIXprop")
  
  M.ret <- t(apply(obj$Subjects[[matrix.type]] ,1, function(x) as.numeric(x)))
  colnames(M.ret) <- colnames(obj$Subjects[[matrix.type]])
  rownames(M.ret) <- rownames(obj$Subjects[[matrix.type]])
  return(M.ret)
}
#' GetRMSE
#' returns the Root Mean Square Error between y - XB, with unnormalized B>0 (absolute) or normalized \code{sum(B)==1} (proportions)
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @param type a character string "proportion" or "absolute". short word is allowed (default: "proportion")
#' @export
#' @return a SxK data matrix with K cell types (according to the signature matrix) and S subjects
#'
GetRMSE <- function(obj, type = c("proportion","absolute") ){
##This function returns the RMSE estimated by MIXTURE
  ##Args:
  ## obj: the list object returned by MIXTURE
  ## type : (character) the type of the coefficients: "proportion" or "absolute"
  ##Returns:
  ## a vector of length N (number of subjects) with the corresponding RMSE
  
    type <- match.arg(type)
  switch(type,
         absolute = rmsetype <- "RMSEa",
         proportion = rmsetype <- "RMSEp")
  
  unlist(obj$Subjects$ACCmetric[,rmsetype])
}
#' GetCorrelation
#' returns the Correlation \code{cor(Y,XB)} unnormalized B>0 (absolute) or normalized \code{sum(B)==1} (proportions)
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @param type a character string "proportion" or "absolute". short word is allowed (default: "proportion")
#' @export
#' @return a SxK data matrix with K cell types (according to the signature matrix) and S subjects
#'

GetCorrelation <- function(obj, type = c("proportion","absolute") ){
##This function returns the Correlation estimated by MIXTURE
  ##Args:
  ## obj: the list object returned by MIXTURE
  ## type : (character) the type of the coefficients: "proportion" or "absolute"
  ##Returns:
  ## a vector of length N (number of subjects) with the corresponding Correlation
  
  type <- match.arg(type)
  switch(type,
         absolute = cortype <- "Ra",
         proportion = cortype <- "Rp")
  
  unlist(obj$Subjects$ACCmetric[,cortype])
}

#' GetMetrix
#' returns the Root Mean Square Error between y - XB, with unnormalized B>0 (absolute) or normalized \code{sum(B)==1} (proportions)
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @param metric : (character) "all" all the metrics (type ignored), "RMSE": (see GetRMSE), "R": (see GetCorrelation)
#' @param type a character string "proportion" or "absolute". short word is allowed (default: "proportions")
#' @export
#' @return if "all" : ACCmetric matrix, otherwise see GetRMSE or GetCorrelation
GetMetrics <- function(obj, metric = c("all","RMSE", "R"), type = c("proportion","absolute")){
##This function returns the metrix matrix
  ##Args:
  ## obj: the list object returned by MIXTURE
  ## metric : (character) "all" all the metrics (type ignored), "RMSE": (see GetRMSE), "R": (see GetCorrelation)
  ## type : (character) the type of the coefficients: "proportion" or "absolute"
  ##Returns:
  ## if "all" : ACCmetric matrix, otherwise see GetRMSE or GetCorrelation
  
  metric <- match.arg(metric)
  
  switch(metric,
         all = return(obj$Subject$ACCmetric),
         RMSE = return(GetRMSE(obj, type)),
         R = return(GetCorrelation(obj, type)))
}

#' GetCellTypes
#' returns the Root Mean Square Error between y - XB, with unnormalized B>0 (absolute) or normalized \code{sum(B)==1} (proportions)
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @param delta (float) default 0, the threshold value to define a coefficient as cero
#' @export
#' @return a logical SxL matrix with TRUE if the coefficients > delta , FALSE otherwise.
GetCellTypes <- function(obj, delta = 0){
## This function returns those coefficiente > delta 
##Args:
  ## obj: the list object returned by MIXTURE
##Returns:
  ## a boolean SxL matrix with TRUE if the coefficients > 0 , FALSE otherwise.
  return(GetMixture(obj)>delta)
}

# GetSignificantsSubjects <- function(obj, pval){
# ##This function evaluates the null distribution statistics
# ## Args:
#   ## obj: the list object returned by MIXURE (only if run with nullDist = "PopulationBased")
#   ##pval : p value trheshold for significance
# ##Returns:
#   ## a Sx4 matrix of 
#   if (obj$method == "none"){
#     stop("error: run MIXTURE with nullDist = \"PopulationBased\" ")
#   }
#   
#   if (pval >0 & pval <=1) {
#     if (obj$method == "PopulationBased"){
#       nd <- apply(obj$PermutedMetrix,2,quantile, c(pval,1-pval) )
#       
#       pvalued <- cbind(nd[1,"RMSEa"] > obj$Subjects$ACCmetric[, "RMSEa"],
#                        nd[1,"RMSEp"] > obj$Subjects$ACCmetric[, "RMSEp"],
#                        nd[2,"Ra"] < obj$Subjects$ACCmetric[, "Ra"],
#                        nd[2,"Rp"] > obj$Subjects$ACCmetric[, "Rp"])
#       colnames(pvalued) <- colnames(obj$Subjects$ACCmetric)[1:4]  
#       
#       return(pvalued)  
#     }else{
#       cat("error, SingleSubjectTest not implemented yet")
#     }
#     
#   }else{
#     cat("error, pval not in range 0 < pval <= 1")
#   }
# }

#' GetPvalues
#' returns the p values matrix.
#' @description if not calculated, stop
#' @param obj an object of class MIXTURE (see MIXTURE)
#' @export
#' @return the p value matrix
GetPvalues <- function(obj){
##This function returns the p values of the model evaluatin each ACC Metric
  # ## Args:
  #   ## obj: the list object returned by MIXURE (only if run with nullDist = "PopulationBased")
    # ##Returns:
  #   ## a Sx4 matrix of p values (columns RMSEa, RMSEp, Ra, Rp)
  # 
  if (obj$method == "none"){
    stop("error: run MIXTURE with nullDist = SingleSubjectTest or PopulationBased")
    
  }
  if (obj$method == "PopulationBased"){
    pvalores <- apply(obj$Subjects$ACCmetric, MARGIN = 1, FUN = function(x, nulld) {
      
      c(RMSEa = sum(nulld[,"RMSEa"] < x[1], na.rm=TRUE),
        RMSp = sum(nulld[,"RMSEp"] < x[2], na.rm=TRUE), 
        Ra = sum(nulld[,"Ra"] > x[3], na.rm=TRUE),
        Rp = sum(nulld[,"Rp"] > x[4], na.rm=TRUE))
    }, nulld = obj$PermutedMetrix) / nrow(obj$PermutedMetrix)
    return(t(pvalores))  
  }else{
    cat("SingleSubjectTest not implemented yet")
  }
  
}
#' Predict
#' This function will predict the profile of the molecular signature detected cell types
#' @description 
#' It will provide \code{Y=X*B} with the absolute or proportion values. The absolute mode (default) is suggested
#' @param obj the object provided by MIXTURE
#' @param signatureMatrix the signature matrix used to fit the models.
#' @param type any of this "absolute" or "proportion", default ("absolute")
#' @export
#' @return the NxS predicted gene expression matrix
Predict <- function(obj, signatureMatrix, type = c("absolute","proportion")){
##This function provides the prediction of the fitted model
##Args:
  ##obj: the object provided by MIXTURE
  ## signatureMatrix: the signature matrix used to fit the models.
##Returns:
  
  type <- match.arg(type, choices = c("absolute","proportion"))
  mat <- GetMixture(obj, type)
  apply(mat, 1, function(beta, sm) {
    u <- sweep(sm,MARGIN=2,beta,'*')
    return(apply(u, 1, sum))
  }, sm = signatureMatrix)
}

#' GetUsedGenes
#' return the amount of signature genes in common with the input data
#' @param obj the object provided by MIXTURE
#' @export
#' @return a vector of Gene identifiers (i.e GeneSymbol) 
GetUsedGenes <- function(obj){
##This function returns the gene list intersected between the expressionMatrix and the signatureMatrix
##Args: 
  ##obj: the list object returned by MIXTURE
##Returns:
  ## a vector of GeneSymbols 
  return(obj$usedGenes)
}

# GetMergedCellTypes <- function(obj, ctClassType  = c("Original","merge4","merge5","merge10","merge11")){
GetMergedCellTypes <- function(obj, ctClassType){  
##This function merges and summarize the cell types proportion from the 22 cell-types to the 11 ones defined by Newman et al.  
  # ctClassType should be a vector of 1x22 , whith names.
  # for instance
  # Type	No.ct	B cells naive	B cells memory	Plasma cells	T cells CD8	T cells CD4 naive	T cells CD4 memory resting	T cells CD4 memory activated	T cells follicular helper	T cells regulatory (Tregs)	T cells gamma delta	NK cells resting	NK cells activated	Monocytes	Macrophages M0	Macrophages M1	Macrophages M2	Dendritic cells resting	Dendritic cells activated	Mast cells resting	Mast cells activated	Eosinophils	Neutrophils
  # LM22-original	22	Bn	Bm	Pcs	CD8T	CD4T	CD4Tr	CD4Ta	Tfh	Tregs	gdt	Nkr	Nka	MM	M0	M1	M2	Dcr	Dca	Mcr	Mca	E	N
  # LM22-merge4	4	B	B	B	CD8T	CD4T	CD4T	CD4T	CD4T	CD4T	CD8T	R	R	R	R	R	R	R	R	R	R	R	R
  # then ctClassType = c("B","B","B","CD8T","CD4T","CD4T","CD4T","CD4T","CD4T","CD8T","R","R","R","R","R","R","R","R","R","R","R","R")
  #this represent a vector of 4 innmuno cell types group, as in https://www.nature.com/articles/s41587-019-0114-2
  
  
    colnames(obj) <- ctClassType
    
    ##This function merges and summarize the cell types proportion from the 22 cell-types to the 11 ones defined by Newman et al.  
    # dat <- data.matrix(t(aggregate(data.frame(t(mixt)), by= list(composite), FUN = sum)[,-1]))
    ret<- do.call(cbind, lapply(unique(ctClassType), function(x) rowSums(obj[,ctClassType==x,drop=FALSE],na.rm = T)))
    colnames(ret) <- unique(ctClassType)
    return(ret)
    
}

##Saving - Open Functions----

#' SaveExcel
#' It saves MIXTURE results as an Excel xlsx file
#' @description 
#' It will build an xlsx file with named sheets as:
#' Absolute (with the absolute coefficients)
#' Proportions
#' Metrics
#' Pvalues
#' UsedGenes
#' @param obj the object provided by MIXTURE
#' @param fileSave (character) the xlsx file name
#' @export
#' @seealso \code{\link{LoadMixtureResultsFromExcel}}
#' 
SaveExcel <- function(obj, fileSave){
##This function save the MIXTURE results into an EXCEL(r) file. It is internally called by MIXTURE if specified
##Args:
  ##obj: The list objects provided by MIXTURE
  ## fileSave. (character) the filename with the corresponing path directory
  wb <- createWorkbook()
  addWorksheet(wb, "Absolute")
  x <- GetMixture(obj, "abs")
  x <- cbind(Subjects = rownames(x),x)
  
  writeData(wb, "Absolute", x)
  addWorksheet(wb, "Proportions")
  x <- GetMixture(obj, "prop")
  x <- cbind(Subjects = rownames(x),x)
  
  writeData(wb, "Proportions", x)
  
  addWorksheet(wb, "Metrics")
  x <- GetMetrics(obj, "all")
  x <- cbind(CellTypes = rownames(x),x)
  writeData(wb, "Metrics", x)
  
  addWorksheet(wb, "Pvalues")
  x <- GetPvalues(obj)
  x <- cbind(Subjects = rownames(x),x)
  writeData(wb, "Pvalues", x)
  
  addWorksheet(wb, "UsedGenes")
  x <- c("GeneSymbol",GetUsedGenes(obj))
  writeData(wb, "UsedGenes", x)
  saveWorkbook(wb, fileSave, overwrite = TRUE)          
}
#' LoadMixtureResultsFromExcel
#' This function load a MIXTURE object from an excel file 
#' @param path the file path to the xlsx file
#' @export
#' @seealso \code{\link{SaveExcel}}
LoadMixtureResultsFromExcel <- function(path){
  if ( str_detect(path,"xlsx") == FALSE) return( NULL)
  abs <- read.xlsx(path, sheet = "Absolute")
  rownames(abs) <- abs[,1]
  abs <- abs[,-1]
  prop <- read.xlsx(path, sheet = "Proportions")
  rownames(prop) <- prop[,1]
  prop <- prop[,-1]
  met <- read.xlsx(path, sheet = "Metrics")
  rownames(met) <- met[,1]
  met <- met[,-1]
  Orig <- list(MIXabs = abs, 
               MIXprop = prop, ACCmetric = met)

  pval <- read.xlsx(path, sheet = "Pvalues")
  rownames(pval) <- pval[,1]
  pval <- pval[,-1]
  ugen <- read.xlsx(path, sheet = "UsedGenes")
  list(Subjects = Orig,PermutedMetrix = NULL, p.values = pval, method = "from file", usedGenes = ugen)#, pvalue = pvalues))   
}  
##Predictions

#Graphics functions----------
## Under development!!
#' ProportionPlot
#' a bar plot for each subject displaying the cell type proportions 
#' 
#' @param obj a MIXTURE object, see \code{\link{LoadMixtureResultsFromExcel}}
#' @param sortBy character indicating the cell type used to sort the profiles along the proportion plot (mutually exclusive with subjOrder). It should match one of the cell type names.
#' @param subjOrder an indexing vector to order the subject along the proportion plot (mutually exclusive with sortBy)
#' @seealso \code{\link{LoadMixtureResultsFromExcel}}
#' @export
#' @return it return a ggplot object 
#' see vignette for examples
ProportionPlot <- function(obj , sortBy, subjOrder){
  if(all(c(!missing(sortBy),!missing(subjOrder) ) ) ){
    stop("Error: only sortBy or subjOrder is allowed")
  }
  
  m.mix <- GetMixture(obj)
  
  
  col.cel.types <- c("chocolate1", "brown4", "black",
                     "tan1","green4", "green2", "lawngreen", "olivedrab3", "olivedrab", "chartreuse4",
                     "goldenrod","gold4","yellow","violetred","orangered1","red",
                     "plum4","plum","navy","mediumblue","cyan",
                     "grey28")
  col.cel.types <- col.cel.types[1:ncol(m.mix)]
  colores <- data.frame(CT=as.character(unique(colnames(m.mix))), Colores=col.cel.types)
  # print(ncol(signature$Mat))
  
  if(missing(subjOrder)==FALSE){
    if(length(subjOrder)!=nrow(m.mix)){
      stop("Error: colOrder should have the same length as nrow(GetMixture(obj))")
    }
    m.mix <- m.mix[subjOrder,]
  }
  
  if(missing(sortBy)){
    
    df.test <- data.frame(b = as.vector(t(m.mix)), 
                          ct = factor(rep(colnames(m.mix),nrow(m.mix))),
                          sbj = factor(rep(rownames(m.mix),each=ncol(m.mix)), 
                                       levels = rownames(m.mix)))   
  }else{
    cn <- colnames(m.mix)
    sortBy <- match.arg(sortBy, cn)
    ord <- order(m.mix[,sortBy])
    
    id.sb <- which(stringr::str_detect(cn, sortBy))
    cn <- c(sortBy,cn[-id.sb])
    m.mix <- m.mix[ord,cn]
    df.test <- data.frame(b = as.vector(t(m.mix)), 
                          ct = factor(rep(cn,nrow(m.mix)), levels = unique(cn)),
                          sbj = factor(rep(rownames(m.mix),each=ncol(m.mix)), 
                                       levels = rownames(m.mix)))     
  }
  
  
  
  # colores <- data.frame(CT=as.character(unique(df.test$ct)), Colores=col.cel.types)
  rownames(colores) <- colores[,1]
  
  ggplot(df.test, aes(sbj, b, fill = ct)) +   geom_col() + 
    scale_fill_manual("Cell Type",breaks = colores[levels(df.test$ct),1], values  = colores[levels(df.test$ct),2]) + 
    theme(axis.text.x = element_text(angle = 90))+ 
    xlab("Subjects") + ylab("Proportions") 
}





PlotNullDensity <- function(obj){
  if (obj$method == "none"){
    par(mfrow = c(2,2))
    plot(density(obj$Subjects$ACCmetric[,"RMSEa"]), main = "RMSE absolute", col="red")
    plot(density(obj$Subjects$ACCmetric[,"RMSEp"]), main = "RMSE proportion", col="red")
    plot(density(obj$Subjects$ACCmetric[,"Ra"]), main = "Corr absolute", col="red")
    plot(density(obj$Subjects$ACCmetric[,"Rp"]), main = "Corr proportion", col="red")
  }
  if (obj$method == "PopulationBased" ){
    par(mfrow = c(2,2))
    plot(density(obj$PermutedMetrix[,"RMSEa"]), main = "RMSE absolute")
    lines(density(obj$Subjects$ACCmetric[,"RMSEa"]), col="red")
    plot(density(obj$PermutedMetrix[,"RMSEp"]), main = "RMSE proportion")
    lines(density(obj$Subjects$ACCmetric[,"RMSEp"]),  col="red")
    plot(density(obj$PermutedMetrix[,"Ra"]), main = "Corr absolute")
    lines(density(obj$Subjects$ACCmetric[,"Ra"]), col="red")
    plot(density(obj$PermutedMetrix[,"Rp"]), main = "Corr proportion")
    lines(density(obj$Subjects$ACCmetric[,"Rp"]), col="red")
  }else{
    cat("not available")
  }
} 



##estimation functions-----
## Modelo original CIBERSORT ----
cibersort <- function(X, y, absolute, abs_method){
##This function was downaloded in October 2018 from the CIBERSORT site.
##the name was changed and the return values to meet the MIXTURE output requirements.
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  ##Absolute
  wo<-q
  u <- sweep(X,MARGIN=2,w,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  ##proportions
  w<-q/sum(q)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  
  
  return(list(Wa = wo, Wp = w,  RMSEa = nusvm, RMSEp = nusvmw, Ra = corrv, Rp = corrvw,  BestParams= out[[mn]]$nu , Iter = mn))
  
}



nu.svm.robust.RFE.old <- function(X,y, nu = c(0.25,0.5,0.75), minProp = 1e-3, maxiter = 6){
#this function is not supossed to be directly called 
#Args:
# X : Nxc gene expression data for the "c" molecular signatures with N genes.
# y : normalized (centred) gene expression data from the subject Mixture (Nx1)
# nu: the nu value (values in a verctor)  for the SVR
# minProp : noise upper bound (0.1% of the FULL proportion range)
# maxiter : maximal number of allowed iterations
# Return:
# a list of:
#     Wa = absolute coefficients (1xc)
#     Wp = proportional coefficients (1xc)
#     RMSEa = RMSE for absolute coefficients (numeric)
#     RMSEp= RMSE for proportion coefficients (numeric) 
#     Ra = correlation for abolute coeffs
#     Rp = correlation for proportional coeffs
#     BestParams = optimal paramters for SVR
#     Iter= number or reached iterations
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  colnames(wsel) <- make.names(colnames(X))
  ok<- TRUE
  
  iter <- 0
  while(ok){
    iter <- iter + 1
    model.tune <- tune.svm(X[,which(wsel>0)],y, type = "nu-regression", kernel = "linear", nu = nu, scale = FALSE  )
    model <- model.tune$best.model
    w <- t(model$coefs) %*% model$SV
    w[w<0] <- 0
    w <- w/sum(w)
    if(any(w < minProp)){
      wsel[which(colnames(wsel) %in% colnames(w)[-which(w >= minProp)]) ] <- 0      
    } else{
      ok <- FALSE
      #print(iter)
    }
    if(iter > maxiter) {
      ok=FALSE
    }
    
  }
  
  wo<-w
  wsel[ which(colnames(wsel) %in% colnames(wo))] <- wo
  
  ##wsel es el absolute, sin normalizar por la sum(w)
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  
  ##calculo de fracciones
  rm(w)
  w<-wsel/sum(wsel)#proportions or fractions (abs_method = sig.scores)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  
  
  return(list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw,  BestParams = unlist(model.tune$best.parameters), Iter=iter))
  
}
## Modelo via regresion eps ----
##eps svm regression (under research)
eps.svm.optim.RFE <- function(X,y, maxiter = 5){
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  
  colnames(wsel) <- make.names(colnames(X))
  
  ok<- TRUE
  sdv <- sd(lm.fit(X,y)$residuals)
  iter <- 0
  while(ok){
    iter <- iter + 1
    model.tune <- tune.svm(X[,which(wsel>0)],y, type = "eps-regression", kernel = "linear", epsilon = c(0.1, 0.25, 0.5, 1), scale = FALSE  )
    model <- model.tune$best.model
    w <- t(model$coefs) %*% model$SV
    
    if(any(w < 0)){
      wsel[which(colnames(wsel) %in% colnames(w)[-which(w>0)]) ] <- 0      
    } else{
      ok <- FALSE
    }
    if(iter > maxiter) {
      ok=FALSE
    }
    
  }
  
  wo<-w
  wsel[ which(colnames(wsel) %in% colnames(wo))] <- wo
  #aqui lo hace como cibersort original, via proporciones del total
  ##wsel es el absolute, sin normalizar por la sum(w)
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  
  ##calculo de fracciones
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  
  return(list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams= unlist(model.tune$best.parameters), Iter = iter))
  
}

########
#This is the model implemented in BMC Medical Genomics volume 12, Article number: 169 (2019)
eps.svm <- function(X,y, eps = 0.0){
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  
  colnames(wsel) <- make.names(colnames(X))
  
  
    # model.tune <- tune.svm(X[,which(wsel>0)],y, type = "eps-regression", kernel = "linear", epsilon = c(0.1, 0.25, 0.5, 1), scale = FALSE  )
    model <- svm(X,y, type = "eps-regression", kernel = "linear", epsilon = eps, scale = FALSE  )
      # model.tune$best.model
    w <- t(model$coefs) %*% model$SV
    
    w[ w < 0] <- 0.0
    
  
  wo<-w
  wsel <- wo
  #aqui lo hace como cibersort original, via proporciones del total
  ##wsel es el absolute, sin normalizar por la sum(w)
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  
  ##calculo de fracciones
  w <- w/sum(w,na.rm = T)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  
  return(list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams= model$tot.nSV, Iter = 1))
  
}

## Modelo via minimos cuadrados no negativos
##non negative least squares (under research)
nnls.optim.RFE <- function(X,y, maxiter = 5){
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  
  colnames(wsel) <- make.names(colnames(X))
  colnames(X) <- make.names(colnames(X))
  ok<- TRUE
  eps <- 0.00000001
  iter <- 0
  while(ok){
    iter <- iter + 1
    mod <- nnls(X[,which(wsel>0)],y)
    w <- coef(mod)
    names(w) <- colnames(X[,which(wsel>0)])
    if(any(w < 0)){
      wsel[which(colnames(wsel) %in% names(w)[-which(w>0)]) ] <- 0      
    } else{
      ok <- FALSE
     # print(iter)
    }
    if(iter > maxiter) {
      ok=FALSE
    }
  }  
  
  wo<-w
  if(all(wo< eps)){
    return(list(Wa=wo, Wp = wo, RMSEa = -1, RMSEp= -1 , Ra=-1, Rp=-1, BestParams = 0, Iter = 0 ))
  }
  wsel <- matrix(0, ncol=ncol(X), nrow=1)
  colnames(wsel) <- make.names(colnames(X))
  wsel[ which(colnames(wsel) %in% names(wo))] <- wo
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  ##fractions or proportions
  w<-wsel/sum(wsel)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  
  return(list(Wa=wo, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams = iter, Iter = iter ))
  
}
ls.rfe.abbas <- function(X,y, maxiter = 5){
  ##The ABBAS model, taken from TIMER source code. Modification only to fit MIXTURE output requirements
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  
  colnames(wsel) <- make.names(colnames(X))
  colnames(X) <- make.names(colnames(X))
  ok<- TRUE
  eps <- 0.00000001
  iter <- 0
  while(ok){
    iter <- iter + 1
    mod <- lsfit(X[,which(wsel>0)],y, intercept = FALSE)##abbas model
    w <- coef(mod)
    names(w) <- colnames(X[,which(wsel>0)])
    
    if(any(w < 0)){
      if(length(w)==1) {
        w<-0
        break
      }
      if(all(w<=0)){
        w <- 0*w
        break
      }
      wsel[which(colnames(wsel) %in% names(w)[-which(w>=0)]) ] <- 0      
      
    } else{
      ok <- FALSE
      #  print(iter)
    }
    if(iter > maxiter) {
      ok=FALSE
    }
    
    
  }  
  
  wo<-w
  
  wsel <- matrix(0, ncol=ncol(X), nrow=1)
  colnames(wsel) <- make.names(colnames(X))
  wsel[ which(colnames(wsel) %in% names(wo))] <- wo
  
  
    
  if(all(wsel< 0)){
    return(list(Wa=wsel, Wp = wsel, RMSEa = -1, RMSEp= -1 , Ra=-1, Rp=-1, BestParams = 0, Iter = 0 ))
  }
  
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  ##fractions or proportions
  w<-wsel/sum(wsel)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  ret <- list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams = iter, Iter = iter ) 

  return(ret)
  
}


nc.rfe.rlm <- function(X,y, minProp = 7e-3, maxiter = 5){
  ##The ABBAS model, taken from TIMER source code. Modification only to fit MIXTURE output requirements
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  
  colnames(wsel) <- make.names(colnames(X))
  colnames(X) <- make.names(colnames(X))
  ok<- TRUE
  eps <- 0.00000001
  iter <- 0
  while(ok){
    iter <- iter + 1
    mod <- rlm(X[,which(wsel>=minProp)],y)##abbas model
    w <- coef(mod)
    names(w) <- colnames(X[,which(wsel >= minProp)])
    
    w[w<0] <- 0
    w <- w/sum(w)
    if(all(is.nan(w))){
      return(list(Wa=rep(NA,ncol(X)), Wp = rep(NA,ncol(X)), RMSEa = NA, RMSEp= NA , Ra=NA, Rp=NA,  BestParams = unlist(model$nu), Iter=iter))
    }
    
    if(any(w < minProp) ){
      if(sum(w > 0) == 1) break
      wsel[which(colnames(wsel) %in% colnames(w)[-which(w >= minProp)]) ] <- 0      
    } else{
      ok <- FALSE
      #  print(iter)
    }
    if(iter > maxiter) {
      ok=FALSE
    }
    
    
  }  
  
  wo<-w
  
  wsel <- matrix(0, ncol=ncol(X), nrow=1)
  colnames(wsel) <- make.names(colnames(X))
  wsel[ which(colnames(wsel) %in% names(wo))] <- wo
  
  if(all(wsel< minProp)){
    return(list(Wa=wsel, Wp = wsel, RMSEa = -1, RMSEp= -1 , Ra=-1, Rp=-1, BestParams = 0, Iter = 0 ))
  }
  
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  ##fractions or proportions
  w<-wsel/sum(wsel)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  ret <- list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams = iter, Iter = iter ) 
  return(ret)
  
}




rlm.abis <- function(X,y, maxiter = 5){
  ##The ABIS model, taken from TIMER source code. Modification only to fit MIXTURE output requirements
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  iter <- NA
  colnames(wsel) <- make.names(colnames(X))
  colnames(X) <- make.names(colnames(X))
  
  wo<-coef(rlm(X,y))
  
  wsel <- matrix(0, ncol=ncol(X), nrow=1)
  colnames(wsel) <- make.names(colnames(X))
  wsel[ which(colnames(wsel) %in% names(wo))] <- wo
  
  
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum)
  nusvm <- sqrt((mean((k - y)^2)))
  corrv <- cor(k, y)
  ##fractions or proportions
  w<-wsel/sum(wsel)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum)
  nusvmw <- sqrt((mean((kw - y)^2)))
  corrvw <- cor(kw, y)
  ret <- list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw, BestParams = iter, Iter = iter ) 
  return(ret)
  
}


## My tune
TuneSvmForDeconv <- function(XX,Y, nuseq, delta){
##Internal function, to solve step 3 of the RFE with noise restriction algorithm
  nu.svr<-lapply(nuseq, function(.nu){
     modelo <-svm(XX,Y,type = "nu-regression", kernel = "linear", nu = .nu, scale = FALSE  )
     beta <- t(modelo$coefs) %*% modelo$SV #do we expect a zero intercept?
     beta[beta<0] <-0
     beta<-beta/sum(beta)
     beta[beta < delta] <- 0
     new.W <- sweep(XX,MARGIN=2,beta,'*')
     pred <- apply(new.W, 1, sum, na.rm=T)
     rmse.pred <- sqrt((mean((Y-pred)^2,na.rm=T)))
     return(list(modelo=modelo, RMSEpred = rmse.pred ))
  })
  selector <- which.min(unlist(lapply(nu.svr, function(x) x$RMSEpred)))[1]
  ret.model <- nu.svr[[selector]]$modelo
  colnames(ret.model$SV) <- colnames(XX)
  return(ret.model)
}

#' nu.svm.robust.RFE
#' noise constrained & recursive feature extraction based support vector nu-regression
#' it will solve the y = X*B problem, providing B>0
#' @param X the signature matrix, a NxL matrix of N genes and L cell lines
#' @param y a vector of bulk gene expression sample 
#' @param nu (float) the nu value, default nu = c(0.25,0.5,0.75)
#' @param minProp float. noise threshold level
#' @param maxiter default 6
#' 
#' @export
#' 
#' @author Elmer A. Fernández
nu.svm.robust.RFE <- function(X,y, nu = c(0.25,0.5,0.75), minProp = 7e-3, maxiter = 6){
  #this function is not supossed to be directly called 
  #Args:
  # X : Nxc gene expression data for the "c" molecular signatures with N genes.
  # y : normalized (centred) gene expression data from the subject Mixture (Nx1)
  # nu: the nu value (values in a verctor)  for the SVR
  # minProp : noise upper bound (0.1% of the FULL proportion range)
  # maxiter : maximal number of allowed iterations
  # Return:
  # a list of:
  #     Wa = absolute coefficients (1xc)
  #     Wp = proportional coefficients (1xc)
  #     RMSEa = RMSE for absolute coefficients (numeric)
  #     RMSEp= RMSE for proportion coefficients (numeric) 
  #     Ra = correlation for abolute coeffs
  #     Rp = correlation for proportional coeffs
  #     BestParams = optimal paramters for SVR
  #     Iter= number or reached iterations
  if(all(is.nan(y))) {
    return(list(Wa=rep(NA,ncol(X)), Wp = rep(NA,ncol(X)), RMSEa = NA, RMSEp= NA , Ra=NA, Rp=NA,  BestParams = NA, Iter=NA))
  }
  if(all(is.na(y))) {
    return(list(Wa=rep(NA,ncol(X)), Wp = rep(NA,ncol(X)), RMSEa = NA, RMSEp= NA , Ra=NA, Rp=NA,  BestParams = NA, Iter=NA))
  }
  wsel <- matrix(1, ncol=ncol(X), nrow=1)
  colnames(wsel) <- colnames(X)
  ok<- TRUE
  
  iter <- 0
  while(ok){
    iter <- iter + 1
    model <- TuneSvmForDeconv(XX=X[,which(wsel>0),drop=FALSE],Y = y, nuseq = nu, delta = minProp)
    
    w <- t(model$coefs) %*% model$SV
    w[w<0] <- 0
    w.abs <- w
    w <- w/sum(w,na.rm=T)
   # w[ w < minProp] <- 0
   if(all(is.nan(w))){
     return(list(Wa=rep(NA,ncol(X)), Wp = rep(NA,ncol(X)), RMSEa = NA, RMSEp= NA , Ra=NA, Rp=NA,  BestParams = unlist(model$nu), Iter=iter))
   }
    if(any(w < minProp) ){#normlized test
      
      wsel[which(colnames(wsel) %in% colnames(w)[-which(w >= minProp)]) ] <- 0      
      if(sum(w > 0) == 1) break
    } else{
      ok <- FALSE
     # print(iter)
    }
    if(iter > maxiter) {
      ok=FALSE
    }
    
    
  }
  
  wo <- w.abs
  wsel <- matrix(0, ncol=ncol(X), nrow=1)
  colnames(wsel) <- colnames(X)
  wsel[ which(colnames(wsel) %in% colnames(wo))] <- wo
  
  ##wsel es el absolute, sin normalizar por la sum(w)
  u <- sweep(X,MARGIN=2,wsel,'*')
  k <- apply(u, 1, sum,na.rm=T)
  nusvm <- sqrt((mean((k - y)^2,na.rm=T)))
  corrv <- cor(k, y)
  
  ##calculo de fracciones
  rm(w)
  w<-wsel/sum(wsel)#proportions or fractions (abs_method = sig.scores)
  uw <- sweep(X,MARGIN=2,w,'*')
  kw <- apply(uw, 1, sum,na.rm=T)
  nusvmw <- sqrt((mean((kw - y)^2,,na.rm=T)))
  corrvw <- cor(kw, y)
  
  
  return(list(Wa=wsel, Wp = w, RMSEa = nusvm, RMSEp= nusvmw , Ra=corrv, Rp=corrvw,  BestParams = unlist(model$nu), Iter=iter))
  
}

### PARA REVISAR----
GetMIXTUREfromExcelFile <- function(path){
  #this function recovers the MIXTURE object from the EXCEL(r) results file
  #NOTE: we need to save the null distribution mode!
  #NOTE2: "Metrix" should be changed for "Metrics"
  sn <- getSheetNames(path)
  if( !all(sn %in% c("Absolute","Proportions","Metrix","Pvalues","UsedGenes"))){
    stop("not a MIXTURE EXCEL FILE")
  }
  Subject <- vector("list",3)
  names(Subject) <- c("MIXabs","MIXprop","ACCmetric")
  Subject$MIXabs <- read.xlsx(xlsxFile = path, sheet = "Absolute")
  Subject$MIXprop <- read.xlsx(xlsxFile = path, sheet = "Proportions")
  Subject$ACCmetric <- read.xlsx(xlsxFile = path, sheet = "Metrix")
  mix.list <- vector("list",4) 
  names(mix.list) <- 
    return(
      list(Subjects = Subject,
           PermutedMatrix = NA,
           p.values = read.xlsx(xlsxFile = path, sheet = "Pvalues"),
           usedGenes = read.xlsx(xlsxFile = path, sheet = "UsedGenes"))
    )
}

require(plyr)
ContrastCellTypes <- function(obj, pheno, test = c("wilcoxon","ttest"),
                              type = c("proportion","absolute"),...){
  ##falta controles de tipo
  test <- match.arg(test, choices = c("wilcoxon","ttest"))
  type <- match.arg(type, choices = c("proportion","absolute"))
  cts <- GetMixture(obj, type = type)
  if(test == "wilcoxon"){
    res <- adply(.data = cts, .margins = 2, .fun = function(x, ph,...){
      prueba <- wilcox.test(a~b, data.frame(a=x,b=ph),...)
      return(c(prueba$statistic, prueba$p.value))
    },ph = pheno, ...)
  }else{##assume t.test
    res <- adply(.data = cts, .margins = 2, .fun = function(x, ph,...){
      prueba <- t.test(a~b, data.frame(a=x,b=ph),...)
      return(c(prueba$statistic, prueba$p.value))
    },ph = pheno, ...)    
  }
  colnames(res) <- c("statistic","p.value")
  return(res)
}

###Miscelanias
ReadCibersortWebResults <- function(file, type = "csv", nct = 22){
  type <- match.arg(type, choices = c("csv", "xlsx"))
  if(type == "csv"){
    tabla <- read.csv(file = file, h=T)
    mix.est <- data.matrix(tabla[,2:(nct+1)])
    rownames(mix.est) <- tabla[,1]
    colnames(mix.est) <- unlist(str_replace_all(colnames(mix.est), "\\.", " "))
    return(mix.est)
  }else{
    cat("not implemented yet")
  }
}


FindSurvCutoff <- function(df.fs, n.intervals,plot=FALSE,title=""){
  cutoff <- seq(min(df.fs$CT, na.rm=T),max(df.fs$CT, na.rm=T),length.out = n.intervals)[-c(1,n.intervals)]
  pvals <- unlist(lapply(cutoff, function(i){
    df.fs$Q <- "H"
    df.fs$Q[df.fs$CT <= i] <- "L"
    fit.m <- try(survdiff(Surv(Time, Event) ~ Q, data = df.fs))
    if(class(fit.m) == "try-error") return(NA)
    pchisq(fit.m$chisq, length(fit.m$n)-1, lower.tail = FALSE)
  }))
  
  ic <- cutoff[which.min(pvals)]
  pval <- pvals[which.min(pvals)]
  df.fs$Q <- "H"
  df.fs$Q[df.fs$CT <= ic] <- "L"
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.fs)
  gs <- ggsurvplot(fit.m, df.fs, pval = TRUE, risk.table = T) + ggtitle(paste(round(ic,3),title,sep=" - "))
  
  if(plot) print(gs)
  
  invisible(return(list(gplot=gs, cutoff = ic, pval = pval, 
                        surv.summary = surv_summary(fit.m, df.fs))))
  # p = 0.00226
  # e=0.1
  # z <- qnorm(1-(0.5*p))
  # (4*dnorm(z)/z)+dnorm(z)*(z-(1/z))*log((1-e)*(1-e)/(e*e))
  
}

CoxSurvival <- function(mixtureObj, subjectsDF){
  
}
