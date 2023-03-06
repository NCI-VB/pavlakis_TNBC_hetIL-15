# Recluster Filtered Seurat Object [scRNA-seq][CCBR] (576fe688-c445-48c6-a40a-033da171c149): v6
reclusterd_so <- function(macrophage_filtered, macrophage_filtered_metadata) {
    #image: png

    ## Libraries.
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    suppressMessages(library(cowplot))
    
    ## Read input SO.
    SO = macrophage_filtered$value
    meta <- macrophage_filtered_metadata
    
    ## Save original clusters.
    prepend_text = "old"
    idents_to_save <- c("SCT_snn_res_0_6")

    if(length(idents_to_save) > 0 & !all(idents_to_save %in% colnames(SO[[]]))){
        colnames(SO[[]]) <- gsub("\\.","_",colnames(SO[[]]))
        if(!all(idents_to_save %in% colnames(SO[[]]))){
            stop("Could not find requested metadata columns!")
        }
    }

    for(i in seq_along(idents_to_save)){
        old_column_name <- idents_to_save[i]
        new_colume_name <- paste(prepend_text,old_column_name, sep="_")
        SO@meta.data[[new_colume_name]] <- SO@meta.data[[old_column_name]]
    }

    ## Remove original cluster columns because they are inaccurate
    columns_to_remove <- c("seurat_clusters",grep("^SCT_snn_res",colnames(SO[[]]), value=TRUE))
    SO@meta.data %>% select(-one_of(columns_to_remove)) -> SO@meta.data
    
    npcs = 0

    ## Find new clusters.
    SO <- FindVariableFeatures(object = SO, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
    
    # Method for automatically picking NPCS
    if(npcs == 0){
        SO <- RunPCA(object = SO, npcs = 30, verbose = FALSE,seed.use = 42) # initial run
        sumpcsd = sum(SO@reductions$pca@stdev)
        pcvar = (SO@reductions$pca@stdev/sumpcsd)*100
        cumu <- cumsum(pcvar)
        co1 <- which(cumu > 80 & pcvar < 5)[1]
        co2 <- sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 0.1), decreasing = T)[1] + 1
        npcs = min(co1,co2)
        print(npcs)
    }

    ## Dimensionality reduction
    SO <- RunPCA(object = SO, npcs = npcs, verbose = FALSE,seed.use = 42)
    SO <- RunUMAP(object = SO, reduction = "pca", dims = 1:npcs, seed.use=42)
    SO <- RunTSNE(object = SO, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
    SO <- FindNeighbors(SO, dims = 1:npcs)
    
    ## Find Clusters.
    resolutions <- seq(0.2,1.2,0.2)
    for (r in resolutions) {
        SO <- FindClusters(SO, resolution = r, algorithm = 1)
    }
    print("Clustering successful!")

    ## Fix orig_ident back to orig.ident.
    colnames(SO@meta.data)[colnames(SO@meta.data) == "orig_ident"] <- "orig.ident"

    ncol = ceiling(length(resolutions)^0.5)
    nrow = ceiling(length(resolutions) / ncol)
    imageWidth = 2000*ncol
    imageHeight = 2000*nrow
    dpi = 300

    png(
      filename="reclusterd_so.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    
    reduction = "umap"
    plot.list <- lapply(paste0("SCT_snn_res.",resolutions),function(z) DimPlot(SO, reduction=reduction, group.by=z) +
                                                                                labs(title=z)+
                                                                                theme(plot.title = element_text(hjust = 0.5)))

    g <- plot_grid(plotlist = plot.list, nrow = nrow, ncol = ncol)
    print(g)

    dev.off()

    ## Return SO.
    return(list(value=SO))
}

print("template_function_reclusterd_so.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_macrophage_filtered<-readRDS(paste0(rds_output,"/var_macrophage_filtered.rds"))
Input_is_Seurat_count <- 0
for(item in var_macrophage_filtered){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_macrophage_filtered<-as.data.frame(var_macrophage_filtered)}else{var_macrophage_filtered <- var_macrophage_filtered}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_macrophage_filtered_metadata<-readRDS(paste0(rds_output,"/var_macrophage_filtered_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_macrophage_filtered_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_macrophage_filtered_metadata<-as.data.frame(var_macrophage_filtered_metadata)}else{var_macrophage_filtered_metadata <- var_macrophage_filtered_metadata}
invisible(graphics.off())
var_reclusterd_so<-reclusterd_so(var_macrophage_filtered,var_macrophage_filtered_metadata)
invisible(graphics.off())
saveRDS(var_reclusterd_so, paste0(rds_output,"/var_reclusterd_so.rds"))
