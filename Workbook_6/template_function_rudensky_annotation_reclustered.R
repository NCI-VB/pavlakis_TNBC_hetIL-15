rudensky_annotation_reclustered <- function(reclusterd_so, Rudensky_se, format_DC_RNA_seq) {
    #image: png
    
    imageType = "png"
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
    suppressMessages(library(SingleCellExperiment))
    suppressMessages(library(cowplot))
    
    so <- reclusterd_so$value
    ref <- Rudensky_se$value
    ref2 <- format_DC_RNA_seq$value
    print(ref)

    annotations <- function(so,ref) {

        sce <- as.SingleCellExperiment(so)
        cat("Doing main...\n")
        pred.main <- SingleR(method = "single", assay(sce,"logcounts"), assay(ref,"logcounts"), ref$rudensky.main, fine.tune = TRUE)
        so <- AddMetaData(so, pred.main[["labels"]][,1], col.name="rudensky_main")
        
        cat("Doing fine...\n")        
        pred.fine <- SingleR(method = "single", assay(sce,"logcounts"), assay(ref2,"logcounts"), ref2$labels.main, fine.tune = TRUE)
        so <- AddMetaData(so, pred.fine[["labels"]][,1], col.name="bulk_rna")

        return(so)
    }

    so <- annotations(so,ref)
    print("done")

    numColors = max(length(unique(so@meta.data$bulk_rna)),length(unique(so@meta.data$rudensky_main)))
    colpaired = colorRampPalette(brewer.pal(12,"Paired"))
    cols=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075",colpaired(numColors))
    
    imageWidth = 6000
    imageHeight = 3000
    dpi = 300

    png(
      filename="rudensky_annotation_reclustered.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p1 = DimPlot(so, reduction="umap", group.by="rudensky_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle("Rudensky DC Main Type Annotations")
    p2 = DimPlot(so, reduction="umap", group.by="bulk_rna")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle("Bulk RNA-seq Annotations")
    
    print(plot_grid(p1,p2,nrow=1))

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

print("template_function_rudensky_annotation_reclustered.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_reclusterd_so<-readRDS(paste0(rds_output,"/var_reclusterd_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_reclusterd_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_reclusterd_so<-as.data.frame(var_reclusterd_so)}else{var_reclusterd_so <- var_reclusterd_so}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Rudensky_se<-readRDS(paste0(rds_output,"/var_Rudensky_se.rds"))
Input_is_Seurat_count <- 0
for(item in var_Rudensky_se){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Rudensky_se<-as.data.frame(var_Rudensky_se)}else{var_Rudensky_se <- var_Rudensky_se}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_format_DC_RNA_seq<-readRDS(paste0(rds_output,"/var_format_DC_RNA_seq.rds"))
Input_is_Seurat_count <- 0
for(item in var_format_DC_RNA_seq){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_format_DC_RNA_seq<-as.data.frame(var_format_DC_RNA_seq)}else{var_format_DC_RNA_seq <- var_format_DC_RNA_seq}
invisible(graphics.off())
var_rudensky_annotation_reclustered<-rudensky_annotation_reclustered(var_reclusterd_so,var_Rudensky_se,var_format_DC_RNA_seq)
invisible(graphics.off())
saveRDS(var_rudensky_annotation_reclustered, paste0(rds_output,"/var_rudensky_annotation_reclustered.rds"))
