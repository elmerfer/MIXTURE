##funciones para evaluacion de las firmas

require(ComplexHeatmap)
#' SeflTest 
#' This function performe the Self Test, which means that the model should estimate
#' as pure cell types each cell type (column) present in the signature matrix
#' @param signatureMatrix a NxK molecular signature
#' @export
#' @return an MIXTURE object 
#' @example 
SelfTest <- function(signatureMatrix){
  res <- MIXTURE(expressionMatrix = signatureMatrix,
                 signatureMatrix = signatureMatrix, useCores = 2)
  M <- GetMixture(res)
  M[M==0] <- NA
  hp <- Heatmap(cor(signatureMatrix), cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE,col = c("blue","red"),
          column_title = "Correlation", show_heatmap_legend = FALSE) +
    Heatmap(M, cluster_rows = FALSE, cluster_columns = FALSE,
            column_title = "Pred. Coefs",col = c("blue", "red"),show_heatmap_legend = FALSE)
  print(hp)
  return(invisible(res))
}

#' SimulatedMixtures
#' 
#' It returns the simulated MIxtures from the given molecular signature matrix
#' @description the Beta coefficnets are drown from a uniform distribution between 0.2 to 1
#' @param signatureMatrix the NxK molecular signature matrix
#' @param maxCoefs (integer) the maximum number of simulated coefficients > 0
#' @param maxSamples (integer) the maximum number of simulated samples
#' @param noisy (logical) if TRUE, noisy samples will be provided
#' @export
#' @return a list with a simulated bulk tissues gene expression matrix nrow(signatureMatrix) by maxCoefs columns and their corresponding Beta coefficient matrix.
SimulatedMixtures <- function(signatureMatrix, maxCoefs, maxSamples, noisy = TRUE){
  betas.list <- lapply(1:maxSamples, function(x, Mat, nrep) {

      ns <-  sample(2:nrep,1)
      r <- runif(ns, 0.2,1)
      id <- sample(ncol(Mat),ns)
      betas <- rep(0,ncol(Mat))
      betas[id] <- r/sum(r)
      A <- (Mat %*% betas)
      if(noisy){
        A <- A + + matrix(as.vector(data.matrix(Mat))[sample(1:(nrow(Mat)*ncol(Mat)),nrow(Mat))],ncol=1) 
      }
      list(beta = betas, id = id,A = A)
    }, Mat = data.matrix(signatureMatrix), nrep = maxCoefs)
  M <- do.call(cbind, lapply(betas.list, function(x) x$A)) 
  B <- do.call(rbind, lapply(betas.list, function(x) x$beta))
  colnames(B) <- colnames(signatureMatrix)
  colnames(M) <- paste("Sim",1:ncol(M),sep="")
  return(invisible(list(M=M,B=B)))
}

#' SimulationTest
#' 
#' @param signatureMatrix the NxK molecular signature matrix
#' @param maxCoefs (integer) the maximum number of simulated coefficients > 0
#' @param maxSamples (integer) the maximum number of simulated samples
#' @param noisy (logical) if TRUE, noisy samples will be provided
#' @export
#' @return a list with two slots
#'     MIXTURE (an object from \code{\link{MIXTURE()}}) 
#'     
#'     SimulatedData from \code{\link{SimulatedMixtures()}}

SimulationTest <- function(signatureMatrix, maxCoefs, maxSamples, noisy = TRUE, useCores){
  useCores <- ifelse(missing(useCores), 2L, useCores)
  SM <- SimulatedMixtures(signatureMatrix, maxCoefs, maxSamples, noisy = TRUE)
  res <- MIXTURE(expressionMatrix = SM$M, signatureMatrix = signatureMatrix, useCores = useCores)
  ba <- data.frame(betahat = c( as.numeric(GetMixture(res))),
             betasim = as.numeric(SM$B)
  )
  nbeta.sim <- apply(SM$B ,1, function(x) sum(x>0))
  nbeta.hat <- rowSums(GetMixture(res,"prop")>0)
  df.betas <- data.frame(est = c(nbeta.hat),
                               beta = as.factor(c(nbeta.sim))
                         )
  
  ba$difs <- ba$betahat - ba$betasim
  prom <- mean(ba$difs, na.rm=T)
  sdev <- sd(ba$difs, na.rm=T)
  gp <- ggplot(ba, aes(x=betasim, y = difs)) + geom_point(alpha = 0.5) + 
    geom_smooth(data=subset(ba, betahat > 0), aes(x=betasim, y = difs), colour = "green")+
    geom_hline(yintercept =prom, colour ="blue") +
    geom_hline(yintercept =prom + 1.96*sdev, colour ="red") +
    geom_hline(yintercept =prom - 1.96*sdev, colour ="red") + labs( x = "simulated coefficients", y = "Error")
  
  
  bx1 <- ggplot(df.betas, aes(x=beta, y=est)) +
    geom_boxplot(position = position_dodge(width = 0.7) ) +
    labs( x = "True number of coefficients",
          y = "Estimated number of coefficients") + theme_bw() +
    scale_y_continuous(name ="Estimated number of coefficients", 
                       breaks = c(0:15),labels=as.character(c(0:15)))#+ facet_wrap(~Signature) 
  
  gridExtra::grid.arrange(bx1,ggMarginal(gp,  margins = "y", type = "density", color = "black", size =4), ncol = 1)
  return(invisible(list(MIXTURE=res, SimulatedData = SM)))
}