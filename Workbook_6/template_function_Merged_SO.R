# Combine & Renormalize [scRNA-seq] [CCBR] (7f919d9a-be13-43dd-967f-cbbc8c1cf813): v158
Merged_SO <- function(PCA_and_Elbow_Plots) {
    #image: png
    
    imageType = "png"
    
    suppressMessages(library(Seurat))
    suppressMessages(library(ggplot2))
    suppressMessages(library(gridExtra))
    suppressMessages(library(RColorBrewer))

    SO = PCA_and_Elbow_Plots$value

    #initialize Citeseq functionality as false, 
    #later the template will check for a Protein assay and run if it finds it
    doCiteSeq <- FALSE

    doMergeData <- !FALSE
    dat = vector()
    integratedata = FALSE

    if (length(SO) > 1) {
    for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
    SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
    allgenes <- rownames(SO_merge)
    SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
    } else {
    SO_merge <- SO[[1]]
    allgenes <- rownames(SO_merge)
    SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
    }

    if (!("orig.ident" %in% colnames(SO_merge@meta.data))) {
        SO_merge@meta.data$orig.ident <- SO_merge@meta.data$orig_ident
    }

    if ("Protein" %in% names(SO_merge@assays)){
        doCiteSeq <-TRUE
    }

    if(FALSE){
    SO_merge <- ScaleData(SO_merge, assay = "HTO")
    }
    
    npcs = 18
    Do_SCTransform = TRUE
    vars_to_regress = c("S.Score","G2M.Score")

if (Do_SCTransform){
     if(is.null(vars_to_regress)){
        SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, return.only.var.genes = FALSE)}
     else{       
        SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress,return.only.var.genes = FALSE) 
}
}
else{
    all.genes <- rownames(SO_merge)
    if(is.null(vars_to_regress)){
        SO_merge <- SO_merge
    }
    else{
        SO_merge <- ScaleData(SO_merge, features=all.genes, assay = "RNA", vars.to.regress=vars_to_regress) 
    }
 DefaultAssay(SO_merge) <- "RNA"   
}

if (length(SO)>1) {
        all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)
    if(integratedata==TRUE){
            integ_features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
            if(!is.null(SO[[1]]@assays$SCT)){
                SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ_features)
                k.filter <- min(200, min(sapply(SO, ncol)))
                integ_anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ_features)
                SO_merge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",features.to.integrate = all_features)
                SO_merge <- ScaleData(SO_merge,features=all_features)
            }
            else{
                k.filter <- min(200, min(sapply(SO, ncol)))
                integ_anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ_features)
                SO_merge <- IntegrateData(anchorset = integ_anchors,features.to.integrate = all_features)
                SO_merge <- ScaleData(SO_merge,features=all_features)  
            }}
    }

    SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
    SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
    SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
    SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
    SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)
    

    #check for CITE-seq data and if so, run reductions
     if(doCiteSeq) {
        # Moved below integration step. SO_merge is recreated and this information was lost
        SO_merge <- ScaleData(SO_merge, assay = "Protein")

        print("finding protein variable features...")
        VariableFeatures(SO_merge,assay="Protein") <- rownames(SO_merge$Protein)
        #Support for partial
        if(all(sapply(seq_along(SO),function(i) "Protein" %in% names(SO[[i]]@assays)))){
            print("running protein pca...")
            SO_merge <- RunPCA(object = SO_merge, assay="Protein",npcs = npcs,verbose = FALSE,reduction.name="protein_pca",seed.use = 42)
            SO_merge <- RunUMAP(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein), reduction.name="protein_umap",seed.use=42)
            SO_merge <- RunTSNE(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein),seed.use = 1,reduction.name="protein_tsne",check_duplicates=F)
            SO_merge <- FindNeighbors(SO_merge, assay="Protein",graph.name="Protein_snn",features=rownames(SO_merge$Protein))
        }else{
            doCiteSeq <- FALSE #set to false so we don't cluster protein
        }
        
    } else {
        doCiteSeq <- FALSE
    }

    for (i in seq(0.2,1.2,0.2)){
    SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
        if(doCiteSeq){
            SO_merge <- FindClusters(SO_merge, graph.name="Protein_snn",resolution = i, algorithm = 1)
        }
    }
    print("Clustering successful!")
    
    n <- 60
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    grobsList = list()
    if(TRUE){
    p1 <- DimPlot(SO_merge, reduction = "tsne", group.by = "orig.ident", repel = TRUE,          pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) +         ggtitle("RNA TSNE")
    grobsList[[length(grobsList)+1]] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))

    }
    if(TRUE){
    p2 <- DimPlot(SO_merge, reduction = "umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("RNA UMAP")
    grobsList[[length(grobsList)+1]] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))

    }
    if(TRUE & doCiteSeq){ 
    p3 <- DimPlot(SO_merge, reduction = "protein_tsne", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody TSNE")
    grobsList[[length(grobsList)+1]] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))

    }
    if(TRUE & doCiteSeq){ 
    p4 <- DimPlot(SO_merge, reduction = "protein_umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody UMAP")
    grobsList[[length(grobsList)+1]] <- p4
    print("Added Antibody UMAP")
    print(length(grobsList))

    }

    
    n = ceiling(length(grobsList)^0.5)
    m=ceiling(length(grobsList)/n)
    imageWidth = 1200*n
    imageHeight = 1200*m
    dpi = 300

    grobs=arrangeGrob(grobs=grobsList,ncol=n)

    if (imageType == 'png') {
    png(
      filename="Merged_SO.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=2,
      bg="white",
      res=dpi,
      type="cairo")
    } else {
        library(svglite)
        svglite::svglite(
        file="Merged_SO.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    plot(grobs)

    slot(SO_merge,"commands") <- list()
    cat("\nPCA Object Checksum:\n")
    print(digest::digest(SO_merge))
    return(list(value=SO_merge))
}

#################################################
## Global imports and functions included below ##
#################################################
# 
#   }
# 
# 

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Merged_SO.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_PCA_and_Elbow_Plots<-readRDS(paste0(rds_output,"/var_PCA_and_Elbow_Plots.rds"))
Input_is_Seurat_count <- 0
for(item in var_PCA_and_Elbow_Plots){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_PCA_and_Elbow_Plots<-as.data.frame(var_PCA_and_Elbow_Plots)}else{var_PCA_and_Elbow_Plots <- var_PCA_and_Elbow_Plots}
invisible(graphics.off())
var_Merged_SO<-Merged_SO(var_PCA_and_Elbow_Plots)
invisible(graphics.off())
saveRDS(var_Merged_SO, paste0(rds_output,"/var_Merged_SO.rds"))
