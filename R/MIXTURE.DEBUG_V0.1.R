############################################################
##Author: Elmer A. Fernández
## Institution : CIDIE - CONICET - UCC
## Version : 0.1
## Date: 07/01/2019
## Last Changes:
##it is on GitHub
##CHANGES: 
##########################################################


.debug <- TRUE

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
MIXER <- function(expressionMatrix, signatureMatrix, functionMixture, useCores = 1L, verbose =FALSE){
  X <- data.matrix(signatureMatrix)
  Yorig <- data.matrix(expressionMatrix)
  ##No Idea why
  Ymedian <- max(median(Yorig),1)
  
  nx.genes <- nrow(X)
  id.xy <- rownames(X) %in% rownames(Yorig)
  id.yx <- rownames(Yorig) %in% rownames(X)
  
  ##Normalizing the signature Matrix
  
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
  ##trully simgle subject
  y.median.by.subj <- apply(Y, 2, median)
  Yn <- scale(Y, center = TRUE, scale = TRUE)
  
  ##perform mixture analysis
  
  out <- mclapply(1:ncol(Yn), function(yi) {
  #  print(paste("id ",yi))
    functionMixture(X,y = as.numeric(Yn[,yi]))
  }  , mc.cores = useCores)
  
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
  
  
  return(list(MIXabs = mix.mat.abs, MIXprop = mix.mat.prop, ACCmetrix = mat.res))
  
}

###---------
##PermutationNUll---------
##Main function
#' MIXTURE
#' the MIXTURE algorithm. It solves Y = X*B estimating the B coefficients or cell type proportions (B>0)
#' @param expressionMatrix a GxS gene expression matrix. with row names genes ID as in the signature Matrix. Genes in rows, samples in columns
#' @param signatureMatrix a NxL gene signature matrix for L cell types
#' @param iter integer, number of iterations to estimate the p values
#' @param functionMixture   : the deconvolution function: cibersort, nu.svm.robust.RFE, ls.rfe.abbas (tiped without quotation marks)
#' @param useCores          : (integer) the number of cpus to use.
#' @param verbose           : boolean. if TRUE, msgs on the screen
#' @param nullDist          : (character) select the null distribution mode for p.valu estimation of the regression estimate
#'                            one of the followings : "none" , "PopulationBased" use the whole expressionMatrix to draw "iter" random samples at once. 
#'                                         Then all the models are compared against this null distribution
#'@param filesave           : (character) the name of the output EXCEL(r) file in xlsx format
#'
#'@export
#'
#'@return   A list with the following slots: 
#'Subjects        : a list with the following slots:
#'        MIXabs    : a data.frame with S rows, and L columns with the absoluite regression coefficients
#'        MIXprop   : a data.frame with S rows, and L columns with the proportions
#'        ACCmetrix : a data.frame with S rows, and RMSEa (Root mean squared error from absolute coeficientes),
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
MIXTURE <- function(expressionMatrix , signatureMatrix, iter = 100, functionMixture , useCores ,
                             verbose = FALSE, nullDist = c("none","SingleSubjectTest","PopulationBased"),
                    fileSave){
#This function perform the decovolution of the signatureMatrix over the gene expression subject matrix.

## Returns:
  bpParam <- bpparam()
  
  .ncores = ifelse(useCores > parallel::detectCores(), parallel::detectCores()-1, useCores)
  
  
  ###bla bla bla
  nullDist <- match.arg(nullDist)
  if(verbose) {
    cat("\nRunning...\n")
  }
  
  .list.of.genes <- rownames(expressionMatrix)[which((rownames(expressionMatrix) %in% rownames(signatureMatrix)))]
  
  ##compute the deconvolution onto the original samples
  Orig <- MIXER(expressionMatrix , signatureMatrix, functionMixture , useCores = .ncores)
  
  ##compute the total "immune content"
  total.median.by.subj <- apply(data.matrix(expressionMatrix), 2, median, na.rm=T)
  total.median.by.subj.lm22genes <- apply(data.matrix(expressionMatrix[.list.of.genes,]), 2, median, na.rm=T)
  max.median.full.cohort <- max(median(data.matrix(expressionMatrix)),1)
  
  Orig$ACCmetrix <- cbind(Orig$ACCmetrix, IscBySbj = total.median.by.subj.lm22genes/total.median.by.subj, 
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
      return(MIXER(M.aux , signatureMatrix , functionMixture , useCores = .ncores )$ACCmetrix)
    } )
    
    cat("\nfinish\n")
    metrix <- do.call(function(...) abind(along = 3,... ), lapply(out.l, function(x) x))
    
    out.list <- list(Subjects = Orig,PermutedMetrix = metrix, method = "SingleSubjectTest", , usedGenes = .list.of.genes)
    #return(list(Subjects = Orig,PermutedMetrix = metrix, method = "SingleSubjectTest"))#, pvalue = pvalues))  
  }
  
  if (nullDist == "PopulationBased") {
    expressionMatrix <- data.matrix(expressionMatrix)
    cat("\nPopulation based null distribution\n")
    cat("\nBuilding random population\n")
        M.aux <- do.call(cbind, mclapply(1:iter, function(i) {
        as.vector(expressionMatrix)[sample(nrow(expressionMatrix)*ncol(expressionMatrix), size = nrow(expressionMatrix))]
      }, mc.cores = .ncores  ))
      
  
   rownames(M.aux) <- rownames(expressionMatrix)
   cat("\nBuilding null distribution\n")  
    out.mix <- MIXER(M.aux , signatureMatrix , functionMixture , useCores = .ncores )$ACCmetrix
    cat("\nfinish\n")
    
    
    out.list <- list(Subjects = Orig,PermutedMetrix = out.mix, method = "PopulationBased", usedGenes = .list.of.genes)
    out.list$p.values <- GetPvalues(out.list)
    #return(list(Subjects = Orig,PermutedMetrix = out.mix, method = "PopulationBased"))#, pvalue = pvalues))  
  }
  ## save to excel file
  if( !missing(fileSave) ){
    SaveExcel(out.list, file = fileSave)
    cat(paste("\n",fileSave,"....OK"))
  }
  
  return(out.list)
}

##Access functions----
GetMixture <- function(obj, type = c("proportion","absolute") ){
##This function returns the mixture coefficientes estimated by MIXTURE
##Args:
  ## obj: the list object returned by MIXTURE
  ## type : (character) the type of the coefficients: "proportion" or "absolute"
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
  
  unlist(obj$Subjects$ACCmetrix[,rmsetype])
}

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
  
  unlist(obj$Subjects$ACCmetrix[,cortype])
}


GetMetrics <- function(obj, metric = c("all","RMSE", "R"), type = c("proportion","absolute")){
##This function returns the metrix matrix
  ##Args:
  ## obj: the list object returned by MIXTURE
  ## metric : (character) "all" all the metrics (type ignored), "RMSE": (see GetRMSE), "R": (see GetCorrelation)
  ## type : (character) the type of the coefficients: "proportion" or "absolute"
  ##Returns:
  ## if "all" : ACCmetrix matrix, otherwise see GetRMSE or GetCorrelation
  
  metrix <- match.arg(metric)
  
  switch(metric,
         all = return(obj$Subject$ACCmetrix),
         RMSE = return(GetRMSE(obj, type)),
         R = return(GetCorrelation(obj, type)))
}

GetCellTypes <- function(obj, delta = 0){
## This function returns those coefficientes >0 
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
#       pvalued <- cbind(nd[1,"RMSEa"] > obj$Subjects$ACCmetrix[, "RMSEa"],
#                        nd[1,"RMSEp"] > obj$Subjects$ACCmetrix[, "RMSEp"],
#                        nd[2,"Ra"] < obj$Subjects$ACCmetrix[, "Ra"],
#                        nd[2,"Rp"] > obj$Subjects$ACCmetrix[, "Rp"])
#       colnames(pvalued) <- colnames(obj$Subjects$ACCmetrix)[1:4]  
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

GetPvalues <- function(obj){
##This function returns the p values of the model evaluatin each ACC Metric
  # ## Args:
  #   ## obj: the list object returned by MIXURE (only if run with nullDist = "PopulationBased")
    # ##Returns:
  #   ## a Sx4 matrix of p values (columns RMSEa, RMSEp, Ra, Rp)
  # 
  if (obj$method == "none"){
    cat("error: run MIXTURE with nullDist = SingleSubjectTest or PopulationBased")
    return()
  }
  if (obj$method == "PopulationBased"){
    pvalores <- apply(obj$Subjects$ACCmetrix, MARGIN = 1, FUN = function(x, nulld) {
      
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

GetPredictedMixture <- function(obj, signatureMatrix){
##This function provides the prediction of the fitted model
##Args:
  ##obj: the object provided by MIXTURE
  ## signatureMatrix: the signature matrix used to fit the models.
##Returns:
  ##
  apply(obj$Subjects$MIXprop, 1, function(beta, sm) {
    u <- sweep(sm,MARGIN=2,beta,'*')
    return(apply(u, 1, sum))
  }, sm = signatureMatrix)
}

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
               MIXprop = prop, ACCmetrix = met)

  pval <- read.xlsx(path, sheet = "Pvalues")
  rownames(pval) <- pval[,1]
  pval <- pval[,-1]
  ugen <- read.xlsx(path, sheet = "UsedGenes")
  list(Subjects = Orig,PermutedMetrix = NULL, p.values = pval, method = "from file", usedGenes = ugen)#, pvalue = pvalues))   
}  
##Predictions

#Graphics functions----------
## Under development!!
ProportionsBarPlot <- function(obj, type = c("proportion", "absolute") ){
  # type <- match.arg(type)
  # matrix <- GetMixture(obj, type)
  # 
  # if(is.null(colnames(matrix)) | is.null(rownames(matrix)) ) stop("Error col/row NULL")
  # df <- data.frame(p = as.vector(matrix), 
  #                  ctype = rep(colnames(matrix), nrow(matrix)), 
  #                  subj = rep(rownames(matrix), each = ncol(matrix)))
  # ggplot(data=df, aes(x=subj, y=p, fill=ctype)) +
    # geom_bar(stat="identity")
  m.mix <- GetMixture(obj, type)
  # print(ncol(signature$Mat))
  df.test <- data.frame(b = as.vector(t(m.mix)), 
                        ct = rep(colnames(m.mix),nrow(m.mix)),
                        sbj = factor(rep(rownames(m.mix),each=ncol(m.mix)), 
                                     levels = rownames(m.mix))) 
  col.cel.types <- c("chocolate1", "brown4", "black",
                     "tan1","green4", "green2", "lawngreen", "olivedrab3", "olivedrab", "chartreuse4",
                     "goldenrod","gold4","yellow","violetred","orangered1","red",
                     "plum4","plum","navy","mediumblue","cyan",
                     "grey28")
  col.cel.types <- col.cel.types[1:ncol(m.mix)]
  
  colores <- data.frame(CT=as.character(unique(df.test$ct)), Colores=col.cel.types)
  rownames(colores) <- colores[,1]
  
  ggplot(df.test, aes(sbj, b)) +   geom_col(aes(fill=ct)) + 
    scale_fill_manual(values  = as.character(colores[levels(df.test$ct),2])) + 
    theme(axis.text.x = element_text(angle = 90))+ 
    xlab("Subjects") + ylab("Proportions")
  }





PlotNullDensity <- function(obj){
  if (obj$method == "none"){
    par(mfrow = c(2,2))
    plot(density(obj$Subjects$ACCmetrix[,"RMSEa"]), main = "RMSE absolute", col="red")
    plot(density(obj$Subjects$ACCmetrix[,"RMSEp"]), main = "RMSE proportion", col="red")
    plot(density(obj$Subjects$ACCmetrix[,"Ra"]), main = "Corr absolute", col="red")
    plot(density(obj$Subjects$ACCmetrix[,"Rp"]), main = "Corr proportion", col="red")
  }
  if (obj$method == "PopulationBased" ){
    par(mfrow = c(2,2))
    plot(density(obj$PermutedMetrix[,"RMSEa"]), main = "RMSE absolute")
    lines(density(obj$Subjects$ACCmetrix[,"RMSEa"]), col="red")
    plot(density(obj$PermutedMetrix[,"RMSEp"]), main = "RMSE proportion")
    lines(density(obj$Subjects$ACCmetrix[,"RMSEp"]),  col="red")
    plot(density(obj$PermutedMetrix[,"Ra"]), main = "Corr absolute")
    lines(density(obj$Subjects$ACCmetrix[,"Ra"]), col="red")
    plot(density(obj$PermutedMetrix[,"Rp"]), main = "Corr proportion")
    lines(density(obj$Subjects$ACCmetrix[,"Rp"]), col="red")
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
  names(Subject) <- c("MIXabs","MIXprop","ACCMetrix")
  Subject$MIXabs <- read.xlsx(xlsxFile = path, sheet = "Absolute")
  Subject$MIXprop <- read.xlsx(xlsxFile = path, sheet = "Proportions")
  Subject$ACCMetrix <- read.xlsx(xlsxFile = path, sheet = "Metrix")
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
