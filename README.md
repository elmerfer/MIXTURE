# MIXTURE

A recursive Feature extraction algorithm for robust deconvolution of cell-types mixture from molecular signatures

Since the significant impact of immunotherapy in cancer, the estimation of the immune cell-type proportions present in the tumor becomes crucial. The tumor cell mixture content of a patient has been dissected by statistical and machine learning algorithms, yet the accuracy of inferred proportions is limited by the underlying model. We overcome this drawback with a noise constrained recursive feature extraction added to a ν -SVR algorithm.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The current "functional like" version of the software requires the following libraries

* [data.table](https://cran.r-project.org/web/packages/data.table/)
* [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
* [ade4](https://cran.r-project.org/web/packages/ade4/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
* [e1071](https://cran.r-project.org/web/packages/e1071/index.html)
* [preprocessCore](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)
* [nnls](https://cran.r-project.org/web/packages/nnls/index.html)
* [plyr](https://cran.r-project.org/web/packages/plyr/index.html)
* [abind](https://cran.r-project.org/web/packages/abind/index.html)
* [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)

### Installing

Download the file MIXTURE.R and install them in your favorite directory (i.e ../myFavorite/MIXTURE.R)
Download the file LM22.RData and install them in the same directory (i.e ../myFavorite/LM22.RData) [Newman et al.](http://www.nature.com/nmeth/journal/v12/n5/abs/nmeth.3337.html)


## Running the tests

This example tends to estimate the same pure cell-types from LM22 signature matrix from [Newman et al.](http://www.nature.com/nmeth/journal/v12/n5/abs/nmeth.3337.html)
```
source('~/myFavorite/MIXTURE.R')
##Load signature matrix
load("~/myFavorite/LM22.RData")
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


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Elmer Andrés Fernández** - *Initial work* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

