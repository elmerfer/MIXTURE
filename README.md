![MIXTURE](https://github.com/elmerfer/MIXTURE.App/blob/master/www/Logo_B_1.pdf.png)

A noise constrained Recursive Feature Extraction algorithm for robust deconvolution of cell-types mixture from molecular signatures

Since the significant impact of immunotherapy in cancer, the estimation of the immune cell-type proportions present in a tumor becomes crucial. Currently, the deconvolution of the cell mixture content of a tumor is carried out by different analytic tools, yet the accuracy of inferred cell type proportions has room for improvement. We improve tumor immune environment characterization developing MIXTURE, an analytical method based on a noise constrained recursive variable selection for a support vector regression
[Get manuscript](https://www.biorxiv.org/content/10.1101/726562v1)

# NOTE: 

The [MIXTURE shiny App](https://github.com/elmerfer/MIXTURE.App) is now available

The RUN_MIXTURE code was tested on Linux, Windows and Mac. 

# New! [MIXTURE in Python](https://github.com/MsMatias/MixturePy)


## Getting Started


## Installation
```
install.packages("devtools")
library(devtools)
install_github("elmerfer/MIXTURE")
```

## Running MIXTURE

This example tends to estimate the same pure cell-types from LM22 signature matrix from [Newman et al.](http://www.nature.com/nmeth/journal/v12/n5/abs/nmeth.3337.html). The LM22 matrix was stored as an RData obj for convinience and provided as data
```
library(MIXTURE)
##Load signature matrix
data(LM22)
##  Run the self test on LM22 signature
mix.test <- MIXTURE(expressionMatrix = LM22,          #N x ncol(signatureMatrix) gene expresion matrix to evaluate 
                                                      ##rownames(M) should be the GeneSymbols
              signatureMatrix = LM22,                 #the gene signature matrix (W) such that M = W*betas' 
                                                      #(i.e the LM22 from Newman et al)
              iter = 1000,                            #iterations for the statistical test (null distribution)
              functionMixture = nu.svm.robust.RFE,    #cibersort, nu.svm.robust.rfe, ls.rfe.abbas, 
              useCores = 10L,                         #cores for parallel processing/ if using windows set to 1
              verbose = TRUE,                         #TRUE or FALSE messages  
              nullDist = "PopulationBased",           #"none" or "PopulationBased" if the statistical test should
                                                      #be performed
              fileSave = "MIXTURE_FILE_LM22.xlsx")    #EXCEL file name to stare the results 

save(mix.test, file = "MIXTURE_FILE_LM22.RData") #save full list as an RData object.

```

### Documentation
[MIXTURE vignette](https://github.com/elmerfer/MIXTURE/blob/master/vignettes/MIXTURE.pdf)

[Wiki](https://github.com/elmerfer/MIXTURE/wiki)

## Authors

* **Elmer Andrés Fernández** - *Initial work* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CIDIE]- [CONICET](http://www.conicet.gov.ar) - [UCC](http://www.ucc.edu.ar)


## Collaborators
* **Dario Rocha** -- - *Testing and application on TCGA Data
* **Yamil Mahamoud** -- IBYME-CONICET. Application on Cancer immunotherapy
* **Joaquin Merlo** -- IBYME-CONICET. Application on Cancer immunotherapy

## How to cite

MIXTURE: an improved algorithm for immune tumor microenvironment estimation based on gene expression data. Elmer A. Fernández,  Yamil D. Mahmoud, Florencia Veigas, Darío Rocha,  Mónica Balzarini,  Hugo D. Lujan,  Gabriel A. Rabinovich,  M. Romina Girotti doi: https://doi.org/10.1101/726562

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/elmerfer/MIXTURE.App/blob/master/LICENSE) file for details

