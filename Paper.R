rm(list=ls())
library(data.table)
library(ComplexHeatmap)
library(ade4)#distancia de jaccard
library(ggplot2)
library(circlize)


## Change this path directory to the one where you download the code
##i.e .../.../my_directory/MIXTURE.R
source('~/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/MIXTURE.R')

##Load LM22 Signature
## Change this path directory to the one where you download the code
##i.e .../.../my_directory/LM22.RData
load("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/LM22.RData")

##Utils
TestBetas <- function(mat, b.list){
  #This function estimate the root mean squared error between the estimated beta and the true simulated one
  #Args :
  # mat : the estimated betas (NxC) where N: number of samples, C: number of cell.types in the gene signature
  # b.list  : the list of simulated betas
  # 
  # Returns:
  # the RMSE for each mixture
  err <- vector("numeric", length(b.list))
  for ( i in 1: length(b.list)){
    err[i] <- sqrt(sum((mat[i,b.list[[i]]$id]-b.list[[i]]$beta[b.list[[i]]$id])^2))        
  }
  return(err)
}


TestExtraBetas <- function(mat, b.list){
  #This function estimate the extra proportion given by th extra coefficients
  #Args :
  # mat : the estimated betas (NxC) where N: number of samples, C: number of cell.types in the gene signature
  # b.list  : the list of simulated betas
  # 
  # Returns:
  # total extra proportion
  err <- vector("numeric", length(b.list))
  for ( i in 1: length(b.list)){
    err[i] <- sum(mat[ i, -b.list[[i]]$id] )        
  }
  return(err)
}


##Scenario a)
M <- LM22
M.c <- M
out.cib <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = 3L)

out.mixture <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

out.abbas <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)


##Amount of estimated cell types
ciber<- apply(GetCellTypes(out.cib),1, sum)
mixture <- apply(GetCellTypes(out.mixture),1, sum)
abbas <- apply(GetCellTypes(out.abbas), 1, sum)

df <- data.frame(id = c(1:length(ciber), 1:length(mixture), 1:length(abbas)),N = c(ciber, mixture, abbas), 
                 Method = c(rep("CIBERSORT",length(ciber)),rep("MIXTURE",length(mixture)), rep("ABBAS",length(abbas))) )
summary(cbind(ciber,mixture, abbas))
# ciber          mixture      abbas      
# Min.   :11.00   Min.   :1   Min.   :6.000  
# 1st Qu.:12.00   1st Qu.:1   1st Qu.:6.000  
# Median :13.00   Median :1   Median :7.000  
# Mean   :13.27   Mean   :1   Mean   :7.091  
# 3rd Qu.:14.00   3rd Qu.:1   3rd Qu.:8.000  
# Max.   :17.00   Max.   :1   Max.   :9.000   

round(100*cumsum(table(ciber))/22) 
# # 11  12  13  14  15  16  17 
# % 18  32  55  82  91  95 100 
round(100*cumsum(table(mixture))/22)
# 1 
# 100 
round(100*cumsum(table(abbas))/22) 
# 6   7   8   9 
# 36  64  91 100

##ploting paired data
b <- runif(nrow(df), -0.1, 0.1)

##Figure of Panel A left
ggplot(df) +
  geom_violin(aes(x = as.numeric(Method), y = N, group = Method, fill= Method), trim =FALSE,
              draw_quantiles = c(0.5))+
  geom_point(aes(x = as.numeric(Method) + b, y = N)) +
  geom_line(aes(x  = as.numeric(Method) + b, y = N, group = id)) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("CIBERSORT", "MIXTURE","ABBAS" ))+
  xlab("Method")

##Figure of Panel A left OPTION 2
ggplot(df) +
  geom_violin(aes(x = as.numeric(Method), y = N, group = Method, fill= Method), trim =FALSE,
               draw_quantiles = c(0.5))+ 
   geom_point(aes(x = as.numeric(Method) + b, y = N)) +
  # geom_line(aes(x  = as.numeric(Method) + b, y = N, group = id)) +
 # scale_x_continuous(breaks = c(1,2,3), labels = c("CIBERSORT", "MIXTURE","ABBAS" ))+
  xlab("Method")



wilcox.test(ciber, mixture, paired = TRUE)

# Wilcoxon signed rank test with continuity correction
# 
# data:  ciber and rfe
# V = 252, p-value = 4.712e-05
# alternative hypothesis: true location shift is not equal to 0
# 
# Warning message:
#   In wilcox.test.default(ciber, rfe, paired = TRUE) :
#   cannot compute exact p-value with ties

##Estimated proportion for each estimated cell type i.e s_k, k!=i
M.cib <- GetMixture(out.cib, "proportion")
M.nu.robust <- GetMixture(out.mixture, "proportion")
M.abbas <- GetMixture(out.abbas, "proportion") 


#autocorrelation of cell-types from the LM22 signature
CorrLM22 <- cor(M.c)


m <- max(cbind(M.cib, M.nu.robust,M.nu.rfe,M.abbas),na.rm=TRUE)
col.map <- colorRamp2(c(0,m), c("blue","red"))
col.map.robust <- colorRamp2(c(0,m), c("grey","red"))
#renaming the cell types for graphical 
cell.types.names <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","FH","Tr","TGD","NKr","NKa","M","M0","M1","M2","Dr","Da","Mr","Ma","E","N")

colnames(M.cib) <- rownames(M.cib) <- cell.types.names
dimnames(M.nu.mixture) <- dimnames(M.cib)
dimnames(M.abbas) <- dimnames(M.cib)
dimnames(CorrLM22) <- dimnames(M.cib)
M.cib[M.cib == 0] <- NA

M.abbas[M.abbas == 0] <- NA

#M.nu.robust[M.nu.robust <= 0] <- NA



##Figure of Panel B
setEPS()
postscript("HeatMap_v3.eps")##we can manage better
Heatmap(M.cib, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = col.map) + 
  Heatmap(M.nu.robust , cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE , column_title = "MIXTURE", name = "MIXTURE",col = col.map.robust)  +
Heatmap(M.abbas , cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE , column_title = "ABBAS", name = "ABBAS",col = col.map)  +
Heatmap(CorrLM22 , cluster_rows = FALSE, cluster_columns = FALSE , column_title = "Cell-types cor", name = "Correlation") 
dev.off()




#Heatmap(M.nnls, cluster_rows = FALSE, cluster_columns = FALSE,  column_title = "NNLS", name = "NNLS")


df.s <- data.frame(id = c(1:length(ciber), 1:length(mixture),1:length(mixture)),
                   N = c(diag(M.cib), diag(M.nu.robust),diag(M.abbas) ), 
                   Method = c(rep("CIBERSORT",22),rep("MIXTURE",22),rep("ABBAS",22)))

df.s <- data.frame(N = c(diag(M.cib), diag(M.nu.robust),diag(M.abbas) ), 
                   Method = rep(c("CIBERSORT","MIXTURE","ABBAS"),each=22))

b <- runif(nrow(dfs), -0.1, 0.1)


summary(cbind(diag(M.cib), diag(M.nu.robust),diag(M.abbas) ))
##Figure Panel A right
ggplot(df.s) +
  geom_violin(aes(x = as.numeric(Method), y = N, group = Method, fill= Method), trim =FALSE,
              draw_quantiles = c(0.5)) +
  xlab("Method")

pairwise.wilcox.test(df.s$N,g=df.s$Method, paired = TRUE)

# Pairwise comparisons using Wilcoxon signed rank test 
# 
# data:  df.s$N and df.s$Method 
# 
#            ABBAS   CIBERSORT
# CIBERSORT 1.4e-06 -        
# MIXTURE   1.4e-06 1.4e-06  
# 
# P value adjustment method: holm 
summary(cbind(CIBERSORT=diag(M.cib), MIXTURE=diag(M.nu.robust),ABBAS=diag(M.abbas) ))
# CIBERSORT         MIXTURE      ABBAS       
# Min.   :0.9993   Min.   :1   Min.   :0.9334  
# 1st Qu.:0.9995   1st Qu.:1   1st Qu.:0.9714  
# Median :0.9997   Median :1   Median :0.9850  
# Mean   :0.9996   Mean   :1   Mean   :0.9794  
# 3rd Qu.:0.9998   3rd Qu.:1   3rd Qu.:0.9938  
# Max.   :0.9999   Max.   :1   Max.   :0.9996

## evaluating the estimation of betas
## we randomly chose from 2 to 8 cell-types from LM22
## then, random proportions are choosen for each drawn from a uniform distritribution [0,1] 
## normalized by the sum of betas 
## then, a new mixture is build by LMM * betas

##Scenario b)
set.seed(123)
betas.list <- lapply(1:1000, function(x, Mat, nrep) {
  ns <-  sample(2:nrep,1)
  r <- runif(ns, 0.2,1)
  id <- sample(22,ns)
  betas <- rep(0,22)
  betas[id] <- r/sum(r)
  A <- Mat %*% betas
  list(beta = betas, id = id,A = A)
}, data.matrix(LM22), nrep = 8)

## assigning to M.c matrix the simulated cell-types mixtures
M.pure <- do.call(cbind, lapply(betas.list, function(x) x$A))
##estimating betas by CIBERSORT
out.betas.cib <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = 3L)
##estimating betas by MIXTURE

out.betas.robust <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

out.betas.abbas <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)


##Getting the estimated proportions (beta hat)
M.r.cib <- GetMixture(out.betas.cib, "prop")
M.r.robust <- GetMixture(out.betas.robust, "prop")
M.r.abbas <- GetMixture(out.betas.abbas, "prop")


Beta <- do.call(rbind, lapply(betas.list, function(x) x$beta))

dim(Beta)
##Root mean squared error between estimated and true simulated Betas
Difs <- sqrt(cbind( CIBERSORT = rowSums((M.r.cib-Beta)^2),
        #MIXTURE = rowSums((M.r.rfe-Beta)^2),
        ROBUST = rowSums((M.r.robust-Beta)^2),
        ABBAS = rowSums((M.r.abbas-Beta)^2)))

par(mfrow=c(1,2))
boxplot(Difs)

summary(Difs)
t.test(Difs[,1], Difs[,2], paired = TRUE)

pairwise.wilcox.test(Difs, g = rep(c("CIBERSORT","MIXTURE","ABBAS"),each=1000), paired = TRUE)
# Pairwise comparisons using Wilcoxon signed rank test 
# 
# data:  Difs and rep(c("CIBERSORT", "MIXTURE", "ABBAS"), each = 1000) 
# 
# ABBAS  CIBERSORT
# CIBERSORT <2e-16 -        
#   MIXTURE   <2e-16 <2e-16   
# 
# P value adjustment method: holm

cib.bet.length <- apply(GetCellTypes(out.betas.cib),1, sum)
rfe.bet.length <- apply(GetCellTypes(out.betas.rfe),1, sum)
robust.bet.length <- apply(GetCellTypes(out.betas.robust),1, sum)
abbas.bet.length <- apply(GetCellTypes(out.betas.abbas),1, sum)

beta.length <- apply(Beta ,1, function(x) sum(x>0))

df.betas <- data.frame(est = c(cib.bet.length,rfe.bet.length, robust.bet.length, abbas.bet.length),
                       beta = as.factor(c(beta.length,beta.length,beta.length,beta.length)),
                       method = c(rep("CIBERSORT", length(cib.bet.length)), 
                                rep("MIXTURE",length(rfe.bet.length)),
                                rep("ROBUST",length(robust.bet.length)),
                                rep("ABBAS",length(rfe.bet.length))))
df.betas <- subset(df.betas, method != "MIXTURE")


#Figure 1.C
ggplot(df.betas, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) + ylim(1,22)

## error between expected betas and simulated ones
err.mixture <- TestBetas(GetMixture(out.betas.rfe, "prop"), betas.list)
err.robust <- TestBetas(GetMixture(out.betas.robust, "prop"), betas.list)
err.cibersort <- TestBetas(GetMixture(out.betas.cib, "prop"), betas.list)
err.abbas <- TestBetas(GetMixture(out.betas.cib, "prop"), betas.list)

df.betas.RMSE <- data.frame(est = c(err.cibersort,err.mixture, err.robust, err.abbas),
                       beta = as.factor(c(beta.length,beta.length,beta.length,beta.length)),
                       method = c(rep("CIBERSORT", length(cib.bet.length)), 
                                  rep("MIXTURE",length(rfe.bet.length)),
                                  rep("ROBUST",length(robust.bet.length)),
                                  rep("ABBAS",length(rfe.bet.length))))
df.betas.RMSE <- subset(df.betas.RMSE, method != "MIXTURE")


#Figure 1.C
ggplot(df.betas.RMSE, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) 



summary(cbind(CIBERSORT= TestExtraBetas(GetMixture(out.betas.cib, "prop"), betas.list),
   RFE= TestExtraBetas(GetMixture(out.betas.rfe, "prop"), betas.list),
   ROBUST=TestExtraBetas(GetMixture(out.betas.robust, "prop"), betas.list),
    ABBAS=TestExtraBetas(GetMixture(out.betas.abbas, "prop"), betas.list) ))
#         CIBERSORT             RFE              ROBUST        ABBAS          
# Min.   :4.219e-05   Min.   :-6.983e-05   Min.   :0   Min.   :7.330e-06  
# 1st Qu.:1.181e-04   1st Qu.: 2.041e-06   1st Qu.:0   1st Qu.:3.625e-03  
# Median :1.496e-04   Median : 2.369e-05   Median :0   Median :9.198e-03  
# Mean   :1.610e-04   Mean   : 3.566e-05   Mean   :0   Mean   :1.353e-02  
# 3rd Qu.:1.914e-04   3rd Qu.: 5.790e-05   3rd Qu.:0   3rd Qu.:2.017e-02  
# Max.   :5.223e-04   Max.   : 2.520e-04   Max.   :0   Max.   :7.705e-02 





## evaluating the estimation of betas + noise
## ## evaluating the estimation of betas
## we randomly chose from 2 to 8 cell-types from LM22,
## then, random proportions are choosen for each drawn from a uniform distritribution [0,1] 
## normalized by the sum of betas 
## random selection of 547 gene expression for the whole signature LM22 matrix
## then, a new mixture is build by LMM * betas + the random gene expression of the 547 genes

##Scenario c)
set.seed(124)
betas.noise.list <- lapply(1:1000, function(x, Mat, nrep) {
  ns <-  sample(2:nrep,1)
  r <- runif(ns, 0.2,1)
  id <- sample(22,ns)
  betas <- rep(0,22)
  betas[id] <- r/sum(r)
  A <- (Mat %*% betas) + matrix(as.vector(data.matrix(Mat))[sample(1:(nrow(Mat)*ncol(Mat)),nrow(Mat))],ncol=1) 
  list(beta = betas, id = id,A = A)
}, data.matrix(LM22), nrep = 8)



M.pure.noise <- do.call(cbind, lapply(betas.noise.list, function(x) x$A))
dim(M.pure.noise)
#M <- riaz[rownames(riaz) %in% rownames(LM22),]

out.betas.noise.cib <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = 10L)


out.robust.noise <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 11L)

out.betas.noise.abbas <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)

##Comparison between simulated and estimated betas

Beta.noise<- do.call(rbind, lapply(betas.noise.list, function(x) x$beta))



beta.length.noise <- apply(Beta.noise ,1, function(x) sum(x>0))


cib.bet.length.noise <- rowSums(GetMixture(out.betas.noise.cib,"prop")>0)

robust.bet.length.noise2 <- rowSums(GetMixture(out.robust2.noise,"prop")>0)

abbas.bet.length.noise <- rowSums(GetMixture(out.betas.noise.abbas,"prop")>0)


df.betas.noise <- data.frame(est = c(cib.bet.length.noise,
                                     robust.bet.length.noise2, 
                                     robust.bet.length.noise,
                                     abbas.bet.length.noise),
                             beta = as.factor(c(beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise)),
                             method = c(rep("CIBERSORT", length(cib.bet.length.noise)), 
                                      rep("MIXTURE",length(rfe.bet.length.noise)),
                                      rep("ROBUST",length(robust.bet.length.noise)),
                                      rep("ABBAS",length(abbas.bet.length.noise))))
df.betas.noise <- subset(df.betas.noise, method != "MIXTURE")
#Figure 1.C
ggplot(df.betas.noise, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) + ylim(1,22)

do.call(rbind,lapply(levels(df.betas.noise$beta), function(x,dd) {
  pairwise.wilcox.test(subset(dd,beta == x)$est, subset(dd,beta == x)$method, paired =TRUE)$p.value
  }, dd = df.betas.noise)) < 0.01

#            ABBAS CIBERSORT
# 2
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 3
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 4
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 5
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 6
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 7
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# 8
# CIBERSORT  TRUE        NA
# ROBUST    FALSE      TRUE

err.mixture.noise <- TestBetas(GetMixture(out.betas.noise.rfe,"prop"), betas.noise.list)
err.robust.noise <- TestBetas(GetMixture(out.betas.noise.robust,"prop"), betas.noise.list)
err.cibersort.noise <- TestBetas(GetMixture(out.betas.noise.cib,"prop"), betas.noise.list)
err.abbas.noise <- TestBetas(GetMixture(out.betas.noise.abbas,"prop"), betas.noise.list)


df.betas.noise.RMSE <- data.frame(est = c(err.cibersort.noise,err.mixture.noise, err.robust.noise, err.abbas.noise),
                            beta = as.factor(c(beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise)),
                            method = c(rep("CIBERSORT", length(cib.bet.length.noise)), 
                                       rep("MIXTURE",length(rfe.bet.length.noise)),
                                       rep("ROBUST",length(robust.bet.length.noise)),
                                       rep("ABBAS",length(abbas.bet.length.noise))))
df.betas.noise.RMSE <- subset(df.betas.noise.RMSE, method != "MIXTURE")


#Figure 1.C
ggplot(df.betas.noise.RMSE, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) 



summary(cbind( CIBERSORT=TestExtraBetas(GetMixture(out.betas.noise.cib, "prop"), betas.noise.list),
               MIXTURE=TestExtraBetas(GetMixture(out.betas.noise.rfe, "prop"), betas.noise.list),
           ROBUST=TestExtraBetas(GetMixture(out.betas.noise.robust, "prop"), betas.noise.list),
               ABBAS=TestExtraBetas(GetMixture(out.betas.noise.abbas, "prop"), betas.noise.list)
))


# CIBERSORT           MIXTURE             ROBUST             ABBAS       
# Min.   :0.003904   Min.   :0.000000   Min.   :0.000000   Min.   :0.0000  
# 1st Qu.:0.034637   1st Qu.:0.006012   1st Qu.:0.004714   1st Qu.:0.1081  
# Median :0.050252   Median :0.010999   Median :0.009874   Median :0.1810  
# Mean   :0.056418   Mean   :0.015469   Mean   :0.014215   Mean   :0.1973  
# 3rd Qu.:0.070378   3rd Qu.:0.019237   3rd Qu.:0.017898   3rd Qu.:0.2669  
# Max.   :0.303740   Max.   :0.170585   Max.   :0.170585   Max.   :0.6902  


df.extra.betas.noise <- data.frame(est = c(TestExtraBetas(GetMixture(out.betas.noise.cib, "prop"), betas.noise.list),
#                                     TestExtraBetas(GetMixture(out.betas.noise.rfe, "prop"), betas.noise.list),
                                     TestExtraBetas(GetMixture(out.betas.noise.robust, "prop"), betas.noise.list),
                                     TestExtraBetas(GetMixture(out.betas.noise.abbas, "prop"), betas.noise.list)),
                             beta = as.factor(c(beta.length.noise,beta.length.noise,beta.length.noise)),
                       method = rep(c("CIBERSORT", "MIXTURE","ABBAS"), each = 1000))
df.extra.betas.noise <- subset(df.extra.betas.noise, method != "MIXTURE")

#Figure 1.C

ggplot(df.extra.betas.noise, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) 


quantile(TestExtraBetas(GetMixture(out.betas.noise.robust, "prop"), betas.noise.list), c(0.9,0.99))
summary(TestExtraBetas(GetMixture(out.betas.noise.robust, "prop"), betas.noise.list))

summary(TestExtraBetas(GetMixture(out.betas.noise.abbas, "prop"), betas.noise.list))

df.extra.betas.noise$Method <- df.extra.betas.noise$method

ggplot(df.extra.betas.noise, aes(x=Method, y=est, fill=Method)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + ylab("Total extra proportion estimates")

do.call(rbind,lapply(levels(df.extra.betas.noise$beta), function(x,dd) {
  pairwise.wilcox.test(subset(dd,beta == x)$est, subset(dd,beta == x)$method, paired =TRUE)$p.value
}, dd = df.extra.betas.noise)) < 0.01
# ABBAS CIBERSORT
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE
# CIBERSORT  TRUE        NA
# ROBUST     TRUE      TRUE




exess<- cbind(CIBERSORT=summary((cib.bet.length.noise - beta.length.noise)), 
      MIXTURE=summary((rfe.bet.length.noise - beta.length.noise)),
      ROBUST = summary((robust.bet.length.noise - beta.length.noise)),
      ABBAS = summary((abbas.bet.length.noise - beta.length.noise)))
summary(exess)
# CIBERSORT         MIXTURE           ROBUST          ABBAS       
# Min.   : 3.000   Min.   : 0.000   Min.   :0.000   Min.   :-2.000  
# 1st Qu.: 7.250   1st Qu.: 4.239   1st Qu.:2.196   1st Qu.: 3.250  
# Median : 8.204   Median : 4.979   Median :2.893   Median : 4.090  
# Mean   : 8.402   Mean   : 5.160   Mean   :3.298   Mean   : 4.030  
# 3rd Qu.: 9.602   3rd Qu.: 5.750   3rd Qu.:3.750   3rd Qu.: 4.795  
# Max.   :14.000   Max.   :11.000   Max.   :8.000   Max.   :10.000 

list(CIBERSORT=(table((cib.bet.length.noise - beta.length.noise)))/10, 
     MIXTURE=(table((rfe.bet.length.noise - beta.length.noise)))/10,
     ROBUST=(table((robust.bet.length.noise - beta.length.noise)))/10,
     ABBAS = (table((abbas.bet.length.noise - beta.length.noise)))/10)
# $CIBERSORT
# 
# 3    4    5    6    7    8    9   10   11   12   13   14 
# 0.3  1.3  3.3  8.6 17.9 20.6 22.1 13.5  7.5  3.6  1.2  0.1 
# 
# $MIXTURE
# 
# 0    1    2    3    4    5    6    7    8    9   10   11 
# 0.3  0.9  6.7 14.7 19.2 21.6 16.9 10.6  6.1  2.0  0.9  0.1 
# 
# $ROBUST
# 
# 0    1    2    3    4    5    6    7    8 
# 4.4 17.3 25.0 22.7 16.5  8.7  4.3  0.9  0.2 
# 
# $ABBAS
# 
# -2   -1    0    1    2    3    4    5    6    7    8    9   10 
# 0.1  0.2  1.1  4.3 11.2 18.3 24.0 19.7 10.9  5.9  3.3  0.7  0.3 

list(CIBERSORT=cumsum(table((cib.bet.length.noise - beta.length.noise)))/10, 
                    MIXTURE=cumsum(table((rfe.bet.length.noise - beta.length.noise)))/10,
     ROBUST=cumsum(table((robust.bet.length.noise - beta.length.noise)))/10,
                    ABBAS = cumsum(table((abbas.bet.length.noise - beta.length.noise)))/10)
# $CIBERSORT
# 3     4     5     6     7     8     9    10    11    12    13    14 
# 0.3   1.6   4.9  13.5  31.4  52.0  74.1  87.6  95.1  98.7  99.9 100.0 
# 
# $MIXTURE
# 0     1     2     3     4     5     6     7     8     9    10    11 
# 0.3   1.2   7.9  22.6  41.8  63.4  80.3  90.9  97.0  99.0  99.9 100.0 
# 
# $ROBUST
# 0     1     2     3     4     5     6     7     8 
# 4.4  21.7  46.7  69.4  85.9  94.6  98.9  99.8 100.0 
# 
# $ABBAS
# -2    -1     0     1     2     3     4     5     6     7     8     9    10 
# 0.1   0.3   1.4   5.7  16.9  35.2  59.2  78.9  89.8  95.7  99.0  99.7 100.0 

#### Analisis de permutaciones ver mas adelante
##Analisis de permutacion-------
## evaluo un modelo aleatorio, pero permutando filas
## en el caso de cibersort, construye un nuevo individuo mezclando de todos los
## pacientes.
## Guardo solo los genes que tienen LM22




