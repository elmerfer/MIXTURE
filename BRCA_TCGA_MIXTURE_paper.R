#We test here the TCGA with MIXER and CIBERSORT.
#

library(edgeR)
library(gplots)
##change the directory to your own directory!!!
#the BRCA RNAseq data can be downloaded from https://www.dropbox.com/s/zki1gkx5mq1quah/BRCA_rna.rds?dl=0
brca <- readRDS("/home/elmer/Dropbox/Doctorandos/DarioRocha/BRCA/processed_data/BRCA_rna.rds")
TNBC <- apply(brca$targets[,c("er","pgr","her2")],1, FUN = function(x) all(x == "negative"))
brca$targets$TNBC <- TNBC ##defining Triple Negative BRCA


##normalize brca counts
 dge <- DGEList(counts = brca$E)
 dge <- calcNormFactors(dge)
 brca.norm <- brca
 brca.norm$E <- cpm(dge$counts)
# 
 M.brca.n <- brca.norm$E
# 
 brca.cib.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = ncores2use, verbose  =  TRUE,
                       iter = 1000, nullDist = "PopulationBased")
# 
 brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
                       iter = 1000, nullDist = "PopulationBased")
# 
# 
 brca.abbas.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture = ls.rfe.abbas, useCores = ncores2use, verbose  =  TRUE,
                         iter = 1000, nullDist = "PopulationBased")
# 


#no norm                      
df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.cib.n),1, sum), 
                              apply(GetCellTypes(brca.mix.n),1, sum),
                              apply(GetCellTypes(brca.abbas.n),1, sum)),
                       method = rep(c("CIBERSORT","MIXTURE","ABBAS"),each = ncol(M.brca.n)))


##Ploting the amount of estimated cell - types per subject for each method
ggplot(df.brca, aes(x=method, y=Nb, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) + ggtitle("BRCA")

summary(cbind(CIBERSORT=apply(GetCellTypes(brca.cib.n),1, sum), 
          MIXTURE=apply(GetCellTypes(brca.mix.n),1, sum),
          ABBAS=apply(GetCellTypes(brca.abbas.n),1, sum)))                       

# CIBERSORT        MIXTURE          ABBAS       
# Min.   : 7.00   Min.   : 2.00   Min.   : 0.000  
# 1st Qu.:11.00   1st Qu.: 6.00   1st Qu.: 4.000  
# Median :12.00   Median : 7.00   Median : 5.000  
# Mean   :12.07   Mean   : 7.46   Mean   : 5.266  
# 3rd Qu.:13.00   3rd Qu.: 9.00   3rd Qu.: 7.000  
# Max.   :17.00   Max.   :13.00   Max.   :11.000 


CT.brca.cib <- GetMixture(brca.cib.n,"proportion")
CT.brca.mix <- GetMixture(brca.mix.n,"proportion")
CT.brca.abbas <- GetMixture(brca.abbas.n,"proportion")

## amount of subjects for each cell-type (beta_cell-type > 0)
cbind(CIBERSORT=colSums(CT.brca.cib>0),MIXTURE=colSums(CT.brca.mix>0),ABBAS=colSums(CT.brca.abbas>0))

# CIBERSORT MIXTURE ABBAS
# B cells naive                     1128     549   150
# B cells memory                     210      90    40
# Plasma cells                      1004     839   846
# T cells CD8                        876     416    49
# T cells CD4 naive                   22       4    46
# T cells CD4 memory resting        1173     519   512
# T cells CD4 memory activated       204      49     4
# T cells follicular helper         1203    1144   604
# T cells regulatory (Tregs)        1019     703   512
# T cells gamma delta                 79      16     3
# NK cells resting                   937     285    62
# NK cells activated                 542      63     5
# Monocytes                          696     180   121
# Macrophages M0                    1006     910   943
# Macrophages M1                    1111     965   659
# Macrophages M2                    1207    1193   834
# Dendritic cells resting            339     142   330
# Dendritic cells activated          438     266   138
# Mast cells resting                1174     686   348
# Mast cells activated                57      14    40
# Eosinophils                         86       2   144
# Neutrophils                        150      29     8



## Pvalores de BRCA normalizado
p.cib.n <- GetPvalues(brca.cib.n)
p.mix.n <- GetPvalues(brca.mix.n)
p.abbas.n <- GetPvalues(brca.abbas.n)

##amount of subjects with significant correlation between observed and predicted 
cbind(CIBERSORT=table(p.cib.n[,"Rp"] < 0.01), MIXTURE=table(p.mix.n[,"Rp"] < 0.01),ABBAS=table(p.abbas.n[,"Rp"] < 0.01))
#         CIBERSORT MIXTURE ABBAS
# FALSE       945     929  1094
# TRUE        270     286   121

##overlapp between significant subjects
venn(list(CIBERSORT = rownames(p.cib.n)[p.cib.n[,"Rp"] < 0.01], 
          MIXTURE = rownames(p.mix.n)[p.mix.n[,"Rp"] < 0.01],
          ABBAS=  rownames(p.abbas.n)[p.abbas.n[,"Rp"] < 0.01] ))


##Survival analysis.  Short survival (S): Those subjects who lives less than 3 years
##                    Long Survival (L): Those subjects who lives more than 6 years

brca$targets$SurvClass <- NA

brca$targets$SurvClass[ brca$targets$survival.days >= 365*6 ] <- "L"
brca$targets$SurvClass[ brca$targets$survival.days <= 365*3 & brca$targets$vital.status == 1] <- "S"
brca$targets$SurvClass <- as.factor(brca$targets$SurvClass)

##Survival by PAM50 intrinsic sutype (defined by PBCMC - Fresno et al. https://doi.org/10.1093/bioinformatics/btw704)
clasif <- as.factor(brca$targets$pam50.pbcmc):brca$targets$SurvClass

#Short vs Long survival on Macrophages M2
df.m2 <- data.frame(Clasif = rep(brca$targets$SurvClass, 3),
                    CT = c(GetMixture(brca.cib.n,"prop")[,"Macrophages M2"],
                           GetMixture(brca.mix.n,"prop")[,"Macrophages M2"],
                           GetMixture(brca.abbas.n,"prop")[,"Macrophages M2"]),
                 Method = rep(c("CIBERSORT","MIXTURE", "ABBAS"), each = ncol(brca$E))
                            )

ggplot(na.omit(df.m2)) +
  geom_violin(aes(x = Clasif, y = CT, fill= Method), trim =FALSE,
              draw_quantiles = c(0.5))  + xlab("(survival")

wilcox.test(CT~Clasif, subset(df.m2, Method == "CIBERSORT"))
# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 8116, p-value = 0.0245
# alternative hypothesis: true location shift is not equal to 0  
wilcox.test(CT~Clasif, subset(df.m2, Method == "MIXTURE"))

# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 6984, p-value = 0.0001235
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(CT~Clasif, subset(df.m2, Method == "ABBAS"))
# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 9144, p-value = 0.412
# alternative hypothesis: true location shift is not equal to 0




saveRDS(list(brca.cib, brca.mix2,brca.abbas), "MIXTURE_ABBAS_Results.rds")
