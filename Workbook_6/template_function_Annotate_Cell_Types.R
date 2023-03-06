# Annotate Cell Types (SingleR) [scRNA-seq] [CCBR] (ba316b6a-e424-49ab-a8af-32af8c85529f): v95
Annotate_Cell_Types <- function(Merged_SO, source_files) {
    #image: png
    
    imageType = "png"
    #Need to copy and paste into console to run:
    #BiocManager::install("GenomeInfoDbData", destdir="/scratch")

    #
    if(!"GenomeInfoDbData" %in% loadedNamespaces()){
        fakepackage("GenomeInfoDbData")
    }
    suppressMessages(library(Seurat))
    suppressMessages(library(SingleR))
    suppressMessages(library(gridExtra))
    suppressMessages(library(tools))
    suppressMessages(library(grid))
    suppressMessages(library(gridBase))
    suppressMessages(library(cowplot))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(magrittr))
    

    so <- Merged_SO$value    
    species <- "Mouse"
    doFineTuning <- FALSE
    useSpark <-FALSE
    clusterList <-so@meta.data$seurat_clusters

    annotations <- function(so) {
        counts = GetAssayData(object = so)[, colnames(x = so)]
        if (species == "Human") {
            singler = CreateSinglerObject(counts = counts, project.name = "projectDesc", annot = NULL, min.genes = 0, technology = '10x',
                species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = doFineTuning,
                numCores = 4,reduce.file.size = T, clusters = so@meta.data$seurat_clusters
, do.signatures = T
            )
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="HPCA_main")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="HPCA")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="BP_encode_main")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="BP_encode")
        }
        if (species == "Mouse") {
            singler = CreateSinglerObject(counts = counts, project.name = "projectDesc" , annot = NULL, min.genes = 0, technology = '10x',
                species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = doFineTuning,
                numCores = 4, reduce.file.size = T ,clusters = so@meta.data$seurat_clusters
, do.signatures = T
            )
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="immgen_main")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="immgen")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="mouseRNAseq_main")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="mouseRNAseq")

        }
        return(so)
    }

    annotations_with_spark <- function(so) {
        N <- 400
        counts = GetAssayData(object = so)[, colnames(x = so)]
        rownames(counts) <- sub(".*?-", "", rownames(counts)) 
        n = ncol(counts)
        s = seq(1, n, by=N)
        doSingleRPerPartition <- function(i) {
           #BiocManager::install("GenomeInfoDbData")
            library(SingleR)
            library(tools)
            A = seq(i,min(i+N-1,n))
            singler = CreateSinglerObject(counts = counts[,A], project.name = "projectDesc", annot = NULL, min.genes = 0, technology = '10x',
                species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = FALSE,
                numCores = 1,reduce.file.size = T, clusters = NULL, do.signatures = T
            )
        }
        if (species == "Human") {
            singler.objects <- spark.lapply(s, doSingleRPerPartition)
            singler = SingleR.Combine(singler.objects)
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="HPCA_main")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="HPCA")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="BP_encode_main")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="BP_encode")

        }
        if (species == "Mouse") {
            singler.objects <- spark.lapply(s, doSingleRPerPartition)
            singler = SingleR.Combine(singler.objects)
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="immgen_main")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
            so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="immgen")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="mouseRNAseq_main")
            so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="mouseRNAseq")

        }

        return(so)
    }

    if (useSpark) {
        so <- annotations_with_spark(so)
        print("done")
    } else {
        so <- annotations(so)
        print("done")
    }

    if (species == "Human") {
        numColors = max(length(unique(so@meta.data$BP_encode_main)),length(unique(so@meta.data$HPCA_main)))
    } else {
        numColors = max(length(unique(so@meta.data$mouseRNAseq_main)),length(unique(so@meta.data$immgen_main)))
    }
    colpaired = colorRampPalette(brewer.pal(12,"Paired"))
    cols=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075",colpaired(numColors))
    
    imageWidth = 5000
    imageHeight = 3000
    dpi = 300

    if (imageType == 'png') {
    png(
      filename="Annotate_Cell_Types.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    } else {
        library(svglite)
        svglite::svglite(
        file="Annotate_Cell_Types.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    if (species =="Human") {
         p1 = DimPlot(so, reduction="tsne", group.by="HPCA_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle("HPCA Main Cell Type Annotations")
         p2 = DimPlot(so, reduction="tsne", group.by="BP_encode_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4))+ ggtitle("BP Encode Main Cell Type Annotations")
    } else {
        p1 = DimPlot(so, reduction="tsne", group.by="immgen_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle("Immgen Main Cell Type Annotations")
        p2 = DimPlot(so, reduction="tsne", group.by="mouseRNAseq_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4))+ ggtitle("Mouse RNAseq Main Cell Type Annotations")
    }
    
    print(plot_grid(p1,p2,nrow=1))
    so@meta.data$Barcode <- rownames(so@meta.data)
    so@meta.data$sample_name <-so@meta.data$orig.ident
    so@meta.data$sample_name <- gsub("-","_",so@meta.data$sample_name)
    so@meta.data$annot <- NULL
            so@meta.data$seurat_clusters <- NULL
    rownames(so@meta.data) <-so@meta.data$Barcode

    slot(so,"commands") <- list()
    cat("\nSingleR Object Checksum:\n")
    print(digest::digest(so))

    return(list(value=so))
}

fakepackage <- function(name, exported = NULL, unexported = NULL, attach = TRUE) {
  # fetch and eval call to create `makeNamespace`
  eval(body(loadNamespace)[[c(8, 4, 4)]])
  # create an empty namespace
  ns <- makeNamespace(name)
  # makethis namespace the closure env of our input functions
  exported <- lapply(exported, `environment<-`, ns)
  unexported <- lapply(unexported, `environment<-`, ns)
  # place these in the namespace
  list2env(exported, ns)
  list2env(unexported, ns)
  # export relevant functions
  namespaceExport(ns, names(exported))
  if(attach) {
    # copy exported funs to "package:pkg" envir of the search path
    attach(exported, name = paste0("package:", name))
  invisible()
}

print("template_function_Annotate_Cell_Types.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Merged_SO<-readRDS(paste0(rds_output,"/var_Merged_SO.rds"))
Input_is_Seurat_count <- 0
for(item in var_Merged_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Merged_SO<-as.data.frame(var_Merged_SO)}else{var_Merged_SO <- var_Merged_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_source_files<-readRDS(paste0(rds_output,"/var_source_files.rds"))
Input_is_Seurat_count <- 0
for(item in var_source_files){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_source_files<-as.data.frame(var_source_files)}else{var_source_files <- var_source_files}
invisible(graphics.off())
var_Annotate_Cell_Types<-Annotate_Cell_Types(var_Merged_SO,var_source_files)
invisible(graphics.off())
saveRDS(var_Annotate_Cell_Types, paste0(rds_output,"/var_Annotate_Cell_Types.rds"))
