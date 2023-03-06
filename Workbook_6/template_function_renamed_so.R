# Rename Seurat Clusters [scRNA-seq] [CCBR] (19334dc7-78d4-428a-8148-82e4a9bc7f2d): v10
renamed_so <- function(reclusterd_so,reclustered_metadata) {
    #image:png

    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("Seurat"))
    suppressMessages(library("cowplot"))
    suppressMessages(library("scales"))

    so <- reclusterd_so$value

    ident_to_rename = "SCT_snn_res_0_6"
    new_ident = "cell_identities"
    make_new_active = TRUE
    reduction = "umap"
    relevel_identities <- TRUE
    relevel_alphabetically <-TRUE
    relevel_order <- c("Cluster A","Cluster B","Cluster C")

    so <- suppressMessages(expr = StashIdent(object = so, save.name = 'old.ident'))
    #active ident is now stored as old.ident in metadata

    if(! ident_to_rename %in% colnames(so@meta.data)){
        ident_to_rename = gsub("_(\\d+)_(\\d+)$", ".\\1.\\2",ident_to_rename) #For resolution
        ident_to_rename = gsub("orig_ident","orig.ident",ident_to_rename)

        if(! ident_to_rename %in% colnames(so@meta.data)){
            stop("Unable to find ident in Seurat object metadata")
        }
    }

    # Set new idents
    Idents(so) <- ident_to_rename
    cat("Current Seurat Clusters:\n")
    print(levels(so))
    #active ident is now ident_to_rename

    #set new idents
    new.idents <- c("Monocyte 1","Monocyte 2","cDC2","CCR7 High DC","CD103intCD11b+ DC","Siglec-H DC","CCR7 High DC","Siglec-H DC","cDC1","CCR7 High DC","Other")
    
    if(length(new.idents) != length(levels(so))){
        stop("New identities must be the same length as current identities!")
    }

    names(new.idents) <- levels(so)
    so <- RenameIdents(object = so, new.idents)
    so <- suppressMessages(expr = StashIdent(object = so, save.name = new_ident))
    #active idents have been renamed, and their values have been stashed in a new column called new_ident

    #Relevel
    if(relevel_identities){

        #alphabetically or in order
        if(relevel_alphabetically){
            so@meta.data[[new_ident]] <- factor(as.character(so@meta.data[[new_ident]]))

        }else{
            if(length(relevel_order) != length(levels(so))){
                stop("Relevel order must be same length as SO levels!")
            }
            so@meta.data[[new_ident]] <- factor(as.character(so@meta.data[[new_ident]]), levels = relevel_order)
        }
        Idents(so) <- new_ident
        #fix these in meta.data then set them as active idents
    }

    # Make new ident
    if(make_new_active){
        so <- suppressMessages(expr = StashIdent(object = so, save.name = 'ident'))
    }else{
        Idents(so) <- "old.ident"
        so <- suppressMessages(expr = StashIdent(object = so, save.name = 'ident'))
    }

    #color palette
    max_colors <- length(levels(so@meta.data[[ident_to_rename]]))
    
    col.pal1 <- hue_pal()(max_colors)
    names(col.pal1) <- levels(so@meta.data[[ident_to_rename]])

    col.pal2 <- col.pal1
    names(col.pal2) <- new.idents
    filter <- duplicated(names(col.pal2))
    col.pal2 <- col.pal2[!filter]

    print(col.pal1)
    print(col.pal2)
    imageWidth = 2600*2
    imageHeight = 2000
    dpi = 300

    p1 = DimPlot(so, reduction=reduction, group.by=ident_to_rename) + scale_color_manual(values = col.pal1) + theme(legend.position="right",plot.title=element_text(hjust=0.5)) + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle(ident_to_rename)

    p2 = DimPlot(so, reduction=reduction, group.by=new_ident) + scale_color_manual(values = col.pal2) + theme(legend.position="right",plot.title=element_text(hjust=0.5)) + guides(override.aes = list(size=2),colour=guide_legend(ncol=4)) + ggtitle(new_ident)
    
    png(
      filename="renamed_so.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    print(plot_grid(p1,p2))
    dev.off()
    return(list(value=so))
}

print("template_function_renamed_so.R #########################################################################")
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
var_reclustered_metadata<-readRDS(paste0(rds_output,"/var_reclustered_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_reclustered_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_reclustered_metadata<-as.data.frame(var_reclustered_metadata)}else{var_reclustered_metadata <- var_reclustered_metadata}
invisible(graphics.off())
var_renamed_so<-renamed_so(var_reclusterd_so,var_reclustered_metadata)
invisible(graphics.off())
saveRDS(var_renamed_so, paste0(rds_output,"/var_renamed_so.rds"))
