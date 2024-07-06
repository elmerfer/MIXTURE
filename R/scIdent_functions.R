#' scIdent: Deconvolve Pseudobulk Samples from Single-Cell RNA-seq Data with MIXTURE
#'
#'
#' @param SeuObj A Seurat object containing single-cell RNA-seq data.
#' @param clusters_metadata Column name in `SeuObj@meta.data` with cluster assignments for each cell. Default is `NULL`, which uses `seurat_clusters` variable in `meta.data`.
#' @param pseudobulk Number of pseudobulk samples to generate for each cluster. Default is 1.
#' @param pct If 'pseudobulk' is greater than 1, percentage of cells to randomly select for each pseudobulk sample, between 0.1 and 0.9. Default is 1.
#' @param ms_treshold Minimum similarity threshold for marker gene selection. Default is 0.09.
#' @param sgmtx A matrix where rows are genes and columns are cell types, used as the molecular signature for deconvolution. Default is `LM22`.
#' @param cores Number of CPU cores to use for computation. Default is 10. If using Windows, it must be set to 1.
#' @param SeuratAssay The assay to use from the Seurat object. Default is "RNA".
#' 
#' @return A list containing the following dataframes:
#' \describe{
#'   \item{clust_idents}{A data frame with one row per analyzed cluster. Includes a column "MS_cluster" (1 if the cluster is composed of cell types in the molecular signature, 0 otherwise), and columns "ident_1", "ident_2", "ident_3" for the top three absolute coefficients of identified cell types in the cluster.}
#'   \item{clust_abs}{Absolute coefficients estimated for each cluster. If `pseudobulk` > 1, retrieves the median of calculated absolute coefficients.}
#'   \item{clust_props}{Normalized coefficients, interpretable as proportions.}
#'   \item{sigOverlap}{Percentage of genes in the molecular signature with more than 0 counts in each pseudobulk sample.}
#'   \item{statAbs}{If `pseudobulk` > 1, a data frame holding the median, IQR, first and third quartile, minimum and maximum absolute coefficients for cell types with a median higher than `ms_treshold`.}
#' }
#' 
#' @details
#' `scIdent` aggregates raw counts of cells within each cluster to create pseudobulk samples, and deconvolves them with MIXTURE. Users can specify the number of pseudobulk samples and the percentage of cells to include in each sample. Then, deconvolution with MIXTURE is performed, paired to any molecular signature, with the default being `LM22`. Multiple pseudobulk samples can be generated to assess the reliability of deconvolution estimations, allowing users to analyze variations among predictions.
#' 
#' @export

scIdent <- function(SeuObj, clusters_metadata = NULL, pseudobulk = 1, pct = 1, ms_treshold = 0.09, sgmtx = LM22, cores = 10L, SeuratAssay = "RNA") {
  if (!require("tidyverse", character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
  
  if (!require("Matrix", character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
  
  if (!require("rlist", character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
  
  if (!require("reshape2", character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
  
  options(dplyr.summarise.inform = FALSE)
  
  if (pseudobulk > 1 & pct == 1) {
    stop("If 'pseudobulk' > 1, 'pct' must range between 0.1 and 0.9")
  }
  
  if (Sys.info()['sysname'] == "Windows" & cores != 1) {
    stop("If using Windows, 'cores' must be set to 1.")
  }
  
  if (is.null(clusters_metadata)) {
    clusters <- SeuObj$seurat_clusters
  } else {
    clusters <- SeuObj@meta.data[,names(SeuObj@meta.data) == clusters_metadata]
  }
  
  if (sum(names(SeuObj@meta.data) == clusters_metadata) == 0) {
    stop(paste0("The column '",clusters_metadata,"' is not present in the meta.data of the Seurat object."))
  }
  
  data <- SeuObj[[SeuratAssay]]$counts
  
  if (nrow(SeuObj@meta.data) != length(clusters)) {
    stop("'clusters' does not equal the number of cells in the dataset.")
  }
  clusters <- as.character(clusters)
  clusters[is.na(clusters)] <- "NA"
  n <- nrow(data)
  cluster_names <- unique(clusters)
  results <- list()
  cat(paste0("Generating ",pseudobulk," pseudobulk sample(s) * ", length(cluster_names)," clusters..."))
  for (cluster in cluster_names) {
    cluster_indices <- which(clusters == cluster)
    cluster_size <- length(cluster_indices)
    if (cluster_size == 1) {
      stop("A cluster is composed by only 1 cell, pseudobulk can't be performed.")}
    cluster_pseudobulk <- as.data.frame(matrix(0, nrow = n, ncol = pseudobulk))
    rownames(cluster_pseudobulk) <- rownames(data)
    names(cluster_pseudobulk) <- paste0("*sim",1:pseudobulk)
    
    for (i in 1:pseudobulk) {
      selected <- sample(cluster_indices,round(cluster_size * pct))
      bulk_subjects <- Matrix::rowSums(data[,unique(selected)])
      cluster_pseudobulk[,i] <- bulk_subjects
    }
    results[[paste0("Cluster_", cluster)]] <- as.data.frame(cluster_pseudobulk)
  }
  
  if (pseudobulk == 1) {
    exp_pseudobulk <- as.data.frame(rlist::list.cbind(results))
    names(exp_pseudobulk) <- paste0(names(results),".",names(exp_pseudobulk))
    sigOverlap <- as.data.frame(exp_pseudobulk[rownames(exp_pseudobulk) %in% rownames(LM22),])
    names(sigOverlap) <- cluster_names
    sigOverlap[sigOverlap > 0] <- 1
    sigOverlap <- as.data.frame(colSums(sigOverlap)/length(rownames(sgmtx)))
    names(sigOverlap)[1] <- "sigOverlap"
    sigOverlap$sigOverlap <- round(sigOverlap$sigOverlap,3)
    sigOverlap$cluster <- rownames(sigOverlap)
    sigOverlap <- sigOverlap[,c(2,1)]
  } else {
    exp_pseudobulk <- as.data.frame(rlist::list.cbind(results))
    sigOverlap <- as.data.frame(lapply(results, FUN = function(x) {x <- as.data.frame(x[rownames(x) %in% rownames(sgmtx),]);
    x <- x %>% 
      mutate_if(is.numeric, ~1 * (. > 0));
    round(mean(colSums(x))/length(rownames(sgmtx)),3)}))
    names(sigOverlap) <- cluster_names
    sigOverlap <- as.data.frame(t(sigOverlap))
    names(sigOverlap)[1] <- "sigOverlap"
    sigOverlap$cluster <- rownames(sigOverlap)
    sigOverlap <- sigOverlap[,c(2,1)]
  }
  
  deconv_rs <- MIXTURE(expressionMatrix = exp_pseudobulk,
                       signatureMatrix = sgmtx,
                       iter = 0L,
                       functionMixture = nu.svm.robust.RFE,
                       useCores = cores)
  
  if (pseudobulk != 1) {
    mix_abs <- as.data.frame(deconv_rs$Subjects$MIXabs) %>% 
      mutate(cluster = gsub("\\*.*","",rownames(.))) %>% 
      mutate(cluster = gsub("Cluster_","",cluster)) %>% 
      mutate(cluster = substr(cluster, 1, nchar(cluster)-1)) %>% 
      reshape2::melt(id.vars="cluster") %>%
      mutate(value = ifelse(is.na(value),0,value)) %>% 
      group_by(cluster,variable) %>% 
      dplyr::summarize(median = median(value),
                       IQR = IQR(value)) %>% 
      ungroup() 
    
    clust_idents <- mix_abs  %>% 
      group_by(cluster) %>%
      dplyr::mutate(median_max = max(median)) %>%
      mutate(median = ifelse(median_max < ms_treshold,NA,median),
             IQR = ifelse(median_max < ms_treshold,NA,IQR)) %>% 
      dplyr::summarize(MS_cluster = ifelse(is.na(max(median)),0,1)) %>% 
      ungroup() %>% 
      left_join(., mix_abs  %>% 
                  group_by(cluster) %>%
                  dplyr::mutate(median_max = max(median)) %>%
                  filter(median_max >= ms_treshold) %>% 
                  filter(median != 0) %>% 
                  filter(median >= ms_treshold) %>% 
                  dplyr::summarize(ident_1 = variable[order(median, decreasing = TRUE)][1],
                                   ident_2 = variable[order(median, decreasing = TRUE)][2],
                                   ident_3 = variable[order(median, decreasing = TRUE)][3]),
                by = "cluster")
    
    cluster_abs <- as.data.frame(deconv_rs$Subjects$MIXabs) %>% 
      mutate(cluster = gsub("\\*.*","",rownames(.))) %>% 
      mutate(cluster = gsub("Cluster_","",cluster)) %>% 
      mutate(cluster = substr(cluster, 1, nchar(cluster)-1)) %>% 
      reshape2::melt(id.vars="cluster") %>%
      mutate(value = ifelse(is.na(value),0,value)) %>% 
      group_by(cluster,variable) %>% 
      dplyr::summarize(median = median(value)) %>% 
      ungroup() %>%
      pivot_wider(names_from = variable, values_from = median)
    cluster_abs[,-1][cluster_abs[,-1] < ms_treshold] <- 0
  } else {
    mix_abs <- as.data.frame(deconv_rs$Subjects$MIXabs) %>% 
      mutate(cluster = gsub("\\*.*","",rownames(.))) %>% 
      mutate(cluster = gsub("Cluster_","",cluster)) %>% 
      mutate(cluster = substr(cluster, 1, nchar(cluster)-1)) %>% 
      reshape2::melt(id.vars="cluster") %>%
      mutate(value = ifelse(is.na(value),0,value)) 
    
    clust_idents <- mix_abs  %>% 
      group_by(cluster) %>%
      dplyr::mutate(value_max = max(value)) %>%
      mutate(value = ifelse(value_max < ms_treshold,NA,value)) %>% 
      dplyr::summarize(MS_cluster = ifelse(is.na(max(value)),0,1)) %>% 
      ungroup() %>% 
      left_join(., mix_abs  %>% 
                  group_by(cluster) %>%
                  dplyr::mutate(value_max = max(value)) %>%
                  filter(value_max >= ms_treshold) %>% 
                  filter(value != 0) %>% 
                  filter(value >= ms_treshold) %>% 
                  dplyr::summarize(ident_1 = variable[order(value, decreasing = TRUE)][1],
                                   ident_2 = variable[order(value, decreasing = TRUE)][2],
                                   ident_3 = variable[order(value, decreasing = TRUE)][3]),
                by = "cluster")
    
    cluster_abs <- as.data.frame(deconv_rs$Subjects$MIXabs) %>% 
      mutate(cluster = gsub("\\*.*","",rownames(.))) %>% 
      mutate(cluster = gsub("Cluster_","",cluster)) %>% 
      mutate(cluster = substr(cluster, 1, nchar(cluster)-1)) %>% 
      reshape2::melt(id.vars="cluster") %>%
      mutate(value = ifelse(is.na(value),0,value)) %>% 
      pivot_wider(names_from = variable, values_from = value) 
    cluster_abs[,-1][cluster_abs[,-1] < ms_treshold] <- 0
  }
  
  cat(paste0("\nComplete! ",sum(clust_idents$MS_cluster)," clusters contain MS cell types."))
  
  statAbs <- if (pseudobulk != 1) {as.data.frame(deconv_rs$Subjects$MIXabs) %>% 
      mutate(cluster = gsub("\\*.*","",rownames(.))) %>% 
      mutate(cluster = gsub("Cluster_","",cluster)) %>% 
      mutate(cluster = substr(cluster, 1, nchar(cluster)-1)) %>% 
      reshape2::melt(id.vars="cluster",variable.name="celltypes") %>%
      mutate(value = ifelse(is.na(value),0,value)) %>% 
      group_by(cluster,celltypes) %>% 
      dplyr::summarize(median = median(value),
                       IQR = IQR(value),
                       q1 = quantile(value,0.25),
                       q3 = quantile(value,0.75),
                       min = min(value),
                       max = max(value)) %>% 
      filter(median != 0) %>% 
      arrange(cluster,dplyr::desc(median)) %>% 
      filter(median >= ms_treshold)
    
  } else {NULL}
  normalize_nonzero <- function(x) {
    nonzero_indices <- x != 0
    sum_nonzero <- sum(x[nonzero_indices])
    if (sum_nonzero == 0) {
      return(x)
    }
    x[nonzero_indices] <- x[nonzero_indices] / sum_nonzero
    return(x)
  }
  
  apply_normalize_nonzero <- function(df) {
    as.data.frame(t(apply(df, 1, normalize_nonzero)))
  }
  
  cluster_props <- cluster_abs
  cluster_props[,-1] <- apply_normalize_nonzero(cluster_props[,-1])
  
  return(list(clust_idents = clust_idents,
              cluster_abs = cluster_abs,
              cluster_props = cluster_props,
              sigOverlap = sigOverlap,
              statAbs = statAbs,
              pseudobulk = pseudobulk,
              pct = pct,
              ms_treshold = ms_treshold,
              clusters_metadata = clusters_metadata))
}

#' clustHeatmap: Generate Heatmap of Cluster Coefficients
#'
#' This function generates a heatmap of cluster coefficients (absolute or normalized) from the output of `scIdent`.
#'
#' @param scIdentObj The output object from the `scIdent` function.
#' @param coef A string indicating the type of coefficients to plot. Can be "abs" for absolute coefficients or "prop" for normalized coefficients. Default is "abs".
#' 
#' @return A ggplot object representing the heatmap of cluster coefficients.
#' 
#' @details
#' `clustHeatmap` creates a heatmap using ggplot2, showing the coefficients for each cluster. The function can plot either absolute or normalized coefficients. 
#' 
#' @examples
#' \dontrun{
#' clustHeatmap(scIdentObj, coef = "abs")
#' }
#' 
#' @export
clustHeatmap <- function(scIdentObj, coef = "abs") {
  load_package <- function(package) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
  
  load_package("hrbrthemes")
  load_package("ggplot2")
  load_package("viridis")
  load_package("reshape2")
  load_package("tidyverse")
  
  if (!(coef %in% c("abs", "prop"))) {
    stop("Invalid 'coef'. Choose either 'abs' or 'prop'.")
  }
  
  data <- if (coef == "abs") {
    scIdentObj$cluster_abs
  } else {
    scIdentObj$cluster_props
  }
  
  name <- if (coef == "abs") "Absolute\ncoefficients" else "Normalized\ncoefficients"
  
  data <- data %>%
    reshape2::melt(id.vars = "cluster") %>%
    filter(value > 0) 
  
  data <-  rbind(data, data.frame(cluster = scIdentObj$clust_idents$cluster[which(scIdentObj$clust_idents$MS_cluster == 0)],
                                  variable = rep(data$variable[1],length(scIdentObj$clust_idents$cluster[which(scIdentObj$clust_idents$MS_cluster == 0)])),
                                  value = rep(NA,length(scIdentObj$clust_idents$cluster[which(scIdentObj$clust_idents$MS_cluster == 0)]))
  ))
  
  data <- data %>%
    mutate(value = if (coef == "abs") round(value, 2) else round(value, 3)) %>%
    mutate(variable = as.character(variable)) %>%
    complete(cluster, variable)
  
  p <- ggplot(data, aes(variable, cluster, fill = value)) + 
    geom_tile(color = "black") + 
    scale_fill_gradient2(mid = "black", high = "red", na.value = "grey", name = name) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(face = "bold", color = "black"),    
          axis.title.y = element_text(face = "bold", color = "black")) +
    geom_label(aes(label = if (coef == "abs") value else gsub("NA%", NA, paste0(value * 100, "%"))), color = "black", fill = "grey") +
    coord_fixed() + xlab("MS cell types") + ylab("Cluster")
  
  suppressWarnings(print(p))
}

#' PlotDimCoef: Plot Dimensional Reduction with Cluster Coefficients
#'
#' This function plots dimensional reduction on a Seurat object with cluster coefficients (absolute or normalized).
#'
#' @param SeuratObj A Seurat object containing single-cell RNA-seq data.
#' @param scIdentObj The output object from the `scIdent` function, containing deconvolution results.
#' @param ms_celltypes A vector of cell types from the molecular signature to plot. Default is `NULL`, which plots all identified cell types.
#' @param reduction The type of dimensional reduction to use (e.g., "pca", "tsne", "umap").
#' @param coef A string indicating the type of coefficients to plot. Can be "abs" for absolute coefficients or "prop" for normalized coefficients. Default is "abs".
#' @param ncol Number of columns for the facet wrap. Default is `NULL`.
#' @param nrow Number of rows for the facet wrap. Default is `NULL`.
#' 
#' @return Plots of dimensional reduction with cluster coefficients for each celltype detected by MIXTURE.
#' 
#' @export
PlotDimCoef <- function(SeuratObj, scIdentObj, ms_celltypes = NULL, reduction, coef = "abs", ncol = NULL, nrow = NULL) {
  load_package <- function(package) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
  
  load_package("dplyr")
  load_package("ggplot2")
  load_package("patchwork")
  load_package("reshape2")
  
  if (!(coef %in% c("abs", "prop"))) {
    stop("Invalid 'coef'. Choose either 'abs' or 'prop'.")
  }
  
  clust_label = scIdentObj$clusters_metadata
  cellnames <- rownames(SeuratObj@meta.data)
  
  ms_celltypes <- if (is.null(ms_celltypes)) {unique((reshape2::melt(scIdentObj$cluster_abs,id.vars="cluster") %>% 
                                                        filter(value > 0) %>% 
                                                        mutate(variable = as.character(variable)))$variable)
  } else {ms_celltypes}
  
  cluster_data <- if (coef == "abs") {
    scIdentObj$cluster_abs[, names(scIdentObj$cluster_abs) %in% c("cluster", ms_celltypes)]
  } else {
    scIdentObj$cluster_props[, names(scIdentObj$cluster_props) %in% c("cluster", ms_celltypes)]
  }
  
  ms_coef_max <- reshape2::melt(cluster_data, id.vars = "cluster")
  ms_coef_max <- max(ms_coef_max$value)
  
  SeuratObj@meta.data <- left_join(SeuratObj@meta.data,
                                   cluster_data,
                                   by = setNames("cluster", clust_label))
  rownames(SeuratObj@meta.data) <- cellnames
  
  p1 <- FeaturePlot(object = SeuratObj, reduction = reduction, features = ms_celltypes, combine = FALSE)
  fix.sc <- scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, ms_coef_max))
  p2 <- lapply(p1, function(x) suppressMessages(x + fix.sc))
  
  patchwork::wrap_plots(p2, ncol = ncol, nrow = nrow)
}
