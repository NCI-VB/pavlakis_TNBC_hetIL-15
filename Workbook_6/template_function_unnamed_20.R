# [scRNA-Seq][CCBR] Color by Antibody (a2b6e4d0-bc67-4e21-aac0-38bbd9089dfa): v3
unnamed_20 <- function(Annotate_Cell_Types,Sample_Names) {
    #image: png
    
    imageType = "png"
    suppressMessages(library(Seurat))
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    suppressMessages(library(gridExtra))

    #Parameters
    seurat_object=Annotate_Cell_Types
    doCiteSeq=FALSE
    samples = c("Treated","Control")
    antibody = c("CD64")
    Reduction = "tsne"
    point_transparency =0.5
    point_shape =16
    point_size=1
    returnSO =FALSE
    Number_of_rows = 0
    point_color = "red"

    SO = seurat_object$value
    print(SO)

    # Fix for underscore
    colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))

    if("active.ident" %in% slotNames(SO)){
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
    } else {
    sample_name = as.factor(SO@meta.data$ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
    }

    print("Missing antibodies(s):")
    print(antibody[!antibody %in% rownames(SO.sub$Protein@scale.data)])
    antibody = antibody[antibody %in% rownames(SO.sub$Protein@scale.data)]

    plotgene <- function(antibody){
    antibody.mat=SO.sub$Protein@scale.data[antibody,]
    antibody.quant=quantile(antibody.mat[antibody.mat>1],probs=c(.1,.5,.90))
    antibody.mat[antibody.mat>antibody.quant[3]]=antibody.quant[3]
    antibody.mat[antibody.mat<antibody.quant[1]]=0

    if (!(doCiteSeq)) {
    if(Reduction == "tsne"){
    p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, antibody=antibody.mat)
    }
    else if(Reduction == "umap"){
    p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, antibody=antibody.mat)
    }
    else{
    p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, antibody=antibody.mat)
    } #if CITEseq is chosen then:
    } else {
        if(Reduction == "tsne"){
    p1 <- DimPlot(SO.sub, reduction = "protein_tsne", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, antibody=antibody.mat)
    }
    else if(Reduction == "umap"){
    p1 <- DimPlot(SO.sub, reduction = "protein_umap", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, antibody=antibody.mat)
    }
    else{
    p1 <- DimPlot(SO.sub, reduction = "protein_pca", group.by = "ident")
    clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, antibody=antibody.mat)
    }
    }
    
    clusmat %>% dplyr::arrange(antibody) -> clusmat
    g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    ggtitle(antibody) +
    geom_point(aes(colour=antibody),alpha=point_transparency,shape=point_shape, size=point_size) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.text=element_text(size=rel(0.5)) )+
    scale_color_gradient(limits = c(0, antibody.quant[3]),low = "lightgrey", high = point_color) +
    xlab(paste(Reduction,"-1",sep='')) + ylab(paste(Reduction,"-2",sep=''))  
    return(g)
    }

    if (Number_of_rows== 0) {
        n = ceiling(length(antibody)^0.5)
    } else {
        n = Number_of_rows  
    }

    m=ceiling(length(antibody)/n)
    imageWidth = 1200*m
    imageHeight = 1200*n
    dpi = 300

    if (imageType == 'png') {
    png(
      filename="unnamed_20.png",
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
        file="unnamed_20.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

grob <- lapply(seq_along(antibody),function(x) plotgene(antibody[x]))
gridExtra::grid.arrange(grobs=grob,nrow=n)

     if(returnSO){
    return(list(value=SO))
  else{
    antibody = as.data.frame(antibody) 
  return(antibody)
}

}

#################################################
## Global imports and functions included below ##
#################################################




    

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_unnamed_20.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Annotate_Cell_Types<-readRDS(paste0(rds_output,"/var_Annotate_Cell_Types.rds"))
Input_is_Seurat_count <- 0
for(item in var_Annotate_Cell_Types){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Annotate_Cell_Types<-as.data.frame(var_Annotate_Cell_Types)}else{var_Annotate_Cell_Types <- var_Annotate_Cell_Types}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Sample_Names<-readRDS(paste0(rds_output,"/var_Sample_Names.rds"))
Input_is_Seurat_count <- 0
for(item in var_Sample_Names){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Sample_Names<-as.data.frame(var_Sample_Names)}else{var_Sample_Names <- var_Sample_Names}
invisible(graphics.off())
var_unnamed_20<-unnamed_20(var_Annotate_Cell_Types,var_Sample_Names)
invisible(graphics.off())
saveRDS(var_unnamed_20, paste0(rds_output,"/var_unnamed_20.rds"))
