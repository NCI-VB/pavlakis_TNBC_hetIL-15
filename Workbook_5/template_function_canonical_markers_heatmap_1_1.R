canonical_markers_heatmap_1_1 <- function(Tmm_name_fix, edited_markers) {
    #image:png
  
  # list.of.packages <- c("pheatmap")
  # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # if(length(new.packages)) BiocManager::install(new.packages)
  install.packages("pheatmap")
    
    library(tidyverse)
    library(scales)
    library(Seurat) #Match color palette
    library(pheatmap)
    library(RColorBrewer)

    df <- Tmm_name_fix #%>% select(-CCS231_Macro)
    #df <- df %>% select(-one_of(grep("CCS231_Macro",colnames(df),value=TRUE)))
    marker.df <- edited_markers
    additional_markers <- c("Cd24a","Spn","Ly6c1","Ly6c2","Ly6g","Adgre1","Siglec1","Vcam1","Cd36")
    tmp.df <- as.data.frame(matrix(nrow = length(additional_markers), ncol = 4, dimnames = list(NULL,c("Genes","Origin","Species","Matches"))))
    tmp.df$Genes   <- additional_markers
    tmp.df$Matches <- additional_markers
    tmp.df$Origin <- "Canonical DC markers"
    tmp.df$Species <- "Mouse"
    marker.df <- rbind(marker.df,tmp.df)

    features_of_interest <- marker.df$Genes[marker.df$Genes %in% df$Gene]
    
    new_order <- c("Fcgr1","Siglec1","Ly6c2","Cx3cr1","Ly6c1","Adgre1","Fcgr3","Csf1r","Batf3","Itgax","Rbpj","Cd207","Cd24a","Cd209a","Xcr1","Irf4","Itgae","Itgam","Sirpa","Irf8","Flt3","Siglech","Clec9a","Irf7","Dpp4","Mgl2","Bcl11a","Tcf4","Tbx21","Spn","Bst2","Cd36","Vcam1","Ccr7","Cd14","Cd8a")

    features_of_interest <- features_of_interest[match(new_order,features_of_interest)]
    df <- df %>% filter(Gene %in% c(features_of_interest))

    row.names(df) <- df$Gene
    df$Gene <- NULL

    df <- df[features_of_interest, ]
    # Scale data
    tmean.scale = t(scale(t(df)))
    tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
    tmean.scale = na.omit(tmean.scale)

    # Setup annotation
    default.colors <- hue_pal()(7)
    marker.colors <- as.list(default.colors[1:3])
    sample.colors <- as.list(default.colors[4:7])
    sample.colors <- c("#eb3423","#003df5","#4eac5b","#7f7f7f")
    print(sample.colors)

    # 
    annotation_row <- marker.df %>% filter(Genes %in% features_of_interest) %>% select(Genes,Origin)
    row.names(annotation_row) <- annotation_row[["Genes"]]
    annotation_row[["Genes"]] <- NULL
    colnames(annotation_row) <- "Marker"
    
    annotation_row[["Marker"]] <- factor(as.character(annotation_row[["Marker"]]), levels = unique(as.character(annotation_row[["Marker"]])))

    annot_levels <- levels(annotation_row[["Marker"]])
    marker.colors <- marker.colors[seq_along(annot_levels)]
    names(marker.colors) <- annot_levels
    

    col.text <- c("CD103+ cDC1","CD11b+ cDC2","CD103+CD11b+ DC","Macrophage")
    annotation_col <- data.frame(Sample=factor(col.text, levels=col.text),row.names=colnames(df))
    names(sample.colors) <- col.text
    print(annotation_col)

    # Combine annotation colors
    annot_col <- c(sample.colors,marker.colors)
    annot_col <- list("Sample" = sample.colors)

    print(annot_col)
    
    # Set up colors
    col.pal <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #PurpleAndYellow(k = 100)
    breaks = seq(-2, 2, length=100)
    legbreaks = seq(-2, 2, length=5)   
    breaks = sapply(breaks, signif, 4)
    legbreaks = sapply(legbreaks, signif, 4)

    imageWidth = 3000
    imageHeight = 6000*0.54
    dpi = 300

    png(
      filename="canonical_markers_heatmap_1_1.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    
    p <- pheatmap(
            tmean.scale, 
            color=col.pal,
            #legend_breaks=legbreaks,
            #legend=FALSE,
            cellwidth=35, 
            # cellheight=10, 
            scale="none",
            treeheight_col=0,
            treeheight_row=0,
            kmeans_k=NA,
            breaks=breaks,
            # height=80,
            fontsize=14,
            #fontsize_row=4,
            #fontsize_col=8,
            show_rownames=TRUE, 
            show_colnames=FALSE,
            cluster_rows=FALSE, 
            cluster_cols=FALSE,
            cutree_rows=1,
            #clustering_distance_rows=drows1, 
            #clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            #annotation_row = annotation_row,
            annotation_colors = annot_col,
            gaps_col = c(3),
            annotation_names_col = FALSE#,
            #labels_col = labels_col
        )
    print(p)
# auto removed:     return(NULL)
}

print("template_function_canonical_markers_heatmap_1_1.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Tmm_name_fix<-readRDS(paste0(rds_output,"/var_Tmm_name_fix.rds"))
Input_is_Seurat_count <- 0
for(item in var_Tmm_name_fix){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Tmm_name_fix<-as.data.frame(var_Tmm_name_fix)}else{var_Tmm_name_fix <- var_Tmm_name_fix}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_edited_markers<-readRDS(paste0(rds_output,"/var_edited_markers.rds"))
var_edited_markers <- unique(var_edited_markers) #Itgae is duplicated in this list
Input_is_Seurat_count <- 0
for(item in var_edited_markers){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_edited_markers<-as.data.frame(var_edited_markers)}else{var_edited_markers <- var_edited_markers}
invisible(graphics.off())
var_canonical_markers_heatmap_1_1<-canonical_markers_heatmap_1_1(var_Tmm_name_fix,var_edited_markers)
invisible(graphics.off())
saveRDS(var_canonical_markers_heatmap_1_1, paste0(rds_output,"/var_canonical_markers_heatmap_1_1.rds"))
