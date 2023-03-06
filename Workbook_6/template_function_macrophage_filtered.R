# [scRNA-Seq][CCBR] Filter Seurat Object by Metadata (ec3f23f9-bcba-4f3a-8a08-8ba611fbb6c7): v114
macrophage_filtered <- function(Annotate_Cell_Types, Metadata_Table, Sample_Names) {
    
    ## Libraries
    suppressMessages(library(ggplot2))
    suppressMessages(library(Seurat))
    suppressWarnings(library(gridExtra))
    suppressMessages(library(grid))
    suppressMessages(library(gridBase))
    suppressMessages(library(cowplot))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(colorspace))
    suppressMessages(library(tidyverse))
    
    ## Inputs.
    so <- Annotate_Cell_Types$value
    reduction = "tsne"
    doCiteSeq <- FALSE
    #image: png
    imageType = "png"
    Keep <- FALSE
    plotinteractive = FALSE
    seed = 10
    metacols = "SCT_snn_res_0_6"
    cols1 <- c("aquamarine3","salmon1","lightskyblue3","plum3","darkolivegreen3","goldenrod1","burlywood2","gray70","firebrick2","steelblue","palegreen4","orchid4","darkorange1","yellow","sienna","palevioletred1","gray60","cyan4","darkorange3","mediumpurple3","violetred2","olivedrab","darkgoldenrod2","darkgoldenrod","gray40","palegreen3","thistle3","khaki1","deeppink2","chocolate3","paleturquoise3","wheat1","lightsteelblue","salmon","sandybrown","darkolivegreen2","thistle2","gray85","orchid3","darkseagreen1","lightgoldenrod1","lightskyblue2","dodgerblue3","darkseagreen3","forestgreen","lightpink2","mediumpurple4","lightpink1","thistle","navajowhite","lemonchiffon","bisque2","mistyrose","gray95","lightcyan3","peachpuff2","lightsteelblue2","lightyellow2","moccasin","gray80","antiquewhite2","lightgrey")

    samples = eval(parse(text=gsub('\\[\\]','c()','c("Control","Treated")')))

    if (length(samples) == 0) {
        samples = unique(SO@meta.data$sample_name)
    }

    ## Replace dots in metadata column names with underscores.
    colnames(so@meta.data) = gsub("\\.", "_", colnames(so@meta.data))
    
    ## If you have protien data, then ...
    if (doCiteSeq) {
        reduction =paste("protein_", reduction,sep='')
    }

    ## Set image dimensions within Vector.
    imageWidth = 2000 * 2
    imageHeight = 2000
    dpi = 300

    ## Set image format (png or svg) based on user input.
    if (imageType == 'png') {
        png(
            filename = "macrophage_filtered.png",
            width = imageWidth,
            height = imageHeight,
            units = "px",
            pointsize = 4,
            bg = "white",
            res = dpi,
            type = "cairo")
    } else {
        library(svglite)
        svglite::svglite(
            file = "macrophage_filtered.png",
            width = round(imageWidth/dpi,digits=2),
            height = round(imageHeight/dpi,digits=2),
            pointsize = 1,
            bg = "white")
    }

    ## Original color-picking code.
    n <- 2e3
    set.seed(seed)
    ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    ourColorSpace <- as(ourColorSpace, "LAB")
    distinctColorPalette <-function(k=1,seed) {
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(currentColorSpace, k, iter.max=20)
        colors <- unname(hex(LAB(km$centers)))
        return(colors)
    }

    ## User-selected metadata column is used to set idents.
    Filter.orig = so@meta.data[[metacols[1]]]
    colname <- metacols[1]

    ident_of_interest = as.factor(so@meta.data[[colname]])
    names(ident_of_interest)=names(so@active.ident)
    so@active.ident <- as.factor(vector())
    so@active.ident <- ident_of_interest
    
    ## Get colors from user parameter and add more if the default list is too short.
    if(class(so@meta.data[[metacols[1]]]) != "numeric"){
        q = length(levels(as.factor(Filter.orig)))
        if(length(cols1) < q) {
            r = q - length(cols1)
            more_cols = distinctColorPalette(r,10)
            cols1 <- c(cols1, more_cols)
        }
        names(cols1) <- levels(as.factor(Filter.orig))

    ## Keep or remove cells based on user input values.
        if (Keep) {
            subsetValue <- c("0","1","2","6","8","9","11","12")
            metaCol <- unique(so@meta.data[[metacols[1]]])
            print("Missing values:")
            print(setdiff(subsetValue,metaCol))
            subsetValue <- intersect(metaCol,subsetValue)
        } else {
            metaCol <- unique(so@meta.data[[colname]])
            valsToRemove <- c("0","1","2","6","8","9","11","12")
            subsetValue <- setdiff(metaCol, valsToRemove)
        }

        ## Subset Seurat object.
        #SO.sub <-SubsetData(so, ident.use=subsetValue)
        SO.sub <- subset(so, idents = subsetValue)

        ## Log output of tables of cell types by samples before and after filtes.
        print("Breakdown of filtered data:")
        print(table(so@meta.data[[metacols[1]]],so@meta.data$orig_ident))
        cat("\n")
        print("After Filtering:")
        print(table(SO.sub@meta.data[[metacols[1]]],SO.sub@meta.data$orig_ident))
    
        ## Set filter for the subsetted SO.
        SO.sub@meta.data[[colname]] <- as.factor(as.character(SO.sub@meta.data[[colname]])) #Relevel Factors

        Filter.sub = SO.sub@meta.data[[colname]]

        ## More color stuff.
        #Set colors for unfiltered and filtered data by sample name:
        n = length(levels(as.factor(Filter.sub)))
        idx = vector("list", n)
        names(idx) <- levels(as.factor(Filter.sub))
        for (i in 1:n) {
            id = Filter.orig %in% levels(as.factor(Filter.sub))[i]
            idx[[i]] <- rownames(so@meta.data)[id]
        }
        cols2 <- cols1[levels(as.factor(Filter.sub))]
        
        ## Make before and after plots.
        title <- paste0("filtered by ", metacols[1], " and split by ", metacols[2])
        p1 = DimPlot(so, reduction=reduction, 
                group.by=colname,
                pt.size=0.1) + 
            theme_classic() + 
            scale_color_manual(values=cols1) + 
            theme(legend.position="right") +
            guides(colour=guide_legend(ncol=1,override.aes = list(size = 2))) +
            ggtitle(colname)
        p2 = DimPlot(so, reduction=reduction, 
                cells.highlight = idx, 
                cols.highlight= rev(cols2[1:n]), 
                sizes.highlight = 0.5) + 
            theme_classic() + 
            theme(legend.position="right")+
            guides(colour=guide_legend(ncol=1,reverse=TRUE,override.aes = list(size = 2))) +
            ggtitle(title)

    ## Else, filter on numeric data with a user defined threshold and direction.
    } else {
        filterDirection <-"greater than"
        metaCol <- unique(so@meta.data[["SCT_snn_res_0_6"]])
        value <- 0.5
        if (filterDirection =="greater than") {
            SO.sub <- subset(so, subset = SCT_snn_res_0_6 > 0.5)
        } else {
            SO.sub <- subset(so, subset = SCT_snn_res_0_6 < 0.5)
        }

        drawtsne <- function(SO,reduction,scalecol,colgrad){
        
        SO.clus <- SO@meta.data[[m]]
    
        p1 <- DimPlot(SO, reduction = reduction, group.by = "ident")
        class(p1$data$ident) <- "numeric"

        if(reduction=="tsne"){
            clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO@meta.data[[m]]))
        } else if(reduction=="umap"){
            clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO@meta.data[[m]]))
        } else if (reduction=="pca"){
            clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO@meta.data[[m]]))
        } else if (reduction=="protein_tsne"){
            clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.numeric(SO@meta.data[[m]]))
        } else if (reduction=="protein_umap"){
            clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.numeric(SO@meta.data[[m]]))
        } else {
            clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.numeric(SO@meta.data[[m]]))
        }
        
        clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
        title=as.character(m)
        clusmat %>% dplyr::arrange(clusid) -> clusmat

        p2 <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
            theme_bw() +
            theme(legend.title=element_blank()) +
            geom_point(aes(colour=clusid),alpha=0.5,shape = 20,size=0.1) +
            #scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(min, midpt2,midpt,midpt3, max), limits = c(0, 1))) +
            #scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = c(min, midpt2,midpt,midpt3, max)) + 
            #scale_color_gradientn(colors=brewer.pal(n = 5, name = "RdBu"), values = scales::rescale(c(min, midpt2,midpt,midpt3, max))) +
            
            scale_color_gradientn(colors=brewer.pal(n = 5, name = colgrad), values = scalecol) +
            
            guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
            ggtitle(title) +
            xlab("umap-1") + ylab("umap-2")
        return(p2)
        }

    m = metacols
        clusid = so@meta.data[[m]]
        maxclus = max(clusid)
        clusmid = 0.01/maxclus        
        min = min(clusid)
        midpt1 = 0.99*value
        midpt = value
        midpt2 = 1.01*value
        max = max(clusid)
        colpoints <- c(min,midpt1,midpt,midpt2,max)
        colpoints <- scales::rescale(colpoints,c(0,1))

    p1 <- drawtsne(so,reduction,colpoints,"RdBu")

    clusid = scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
        clus.quant=quantile(clusid[clusid>0],probs=c(0,.25,.5,.75,1))
        min = clus.quant[1]
        midpt = clus.quant[3]
        midpt3 = clus.quant[2]
        midpt4 = clus.quant[4]
        max = clus.quant[5]  
        colpoints2 <- c(min,midpt3,midpt,midpt4,max)
     
    p2 <- drawtsne(SO.sub,reduction,colpoints2,"Blues")
    
    }

    ## If interactive plot requested, then ...
    if (plotinteractive == TRUE) {
        gp1 <- ggplotly(p1)
        gp2 <- ggplotly(p2)
        p <- subplot(gp1, gp2, nrows=2)
        print(p)
    }
    ## Else, print non-interactive plot.
    else {
        print(plot_grid(p1,p2,nrow=1))
    }
    
    ## Return the subsetted Seurat object.
    return(list(value=SO.sub))
}

## Commented out below on 7/28/21 because it had started failing in training. It is believed to be vestigial. -- Josh M.

#################################################
## Global imports and functions included below ##
#################################################

# 
#   }
# 
# 

# Functions defined here will be available to call in
# the code for any table.

print("template_function_macrophage_filtered.R #########################################################################")
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
var_Metadata_Table<-readRDS(paste0(rds_output,"/var_Metadata_Table.rds"))
Input_is_Seurat_count <- 0
for(item in var_Metadata_Table){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Metadata_Table<-as.data.frame(var_Metadata_Table)}else{var_Metadata_Table <- var_Metadata_Table}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Sample_Names<-readRDS(paste0(rds_output,"/var_Sample_Names.rds"))
Input_is_Seurat_count <- 0
for(item in var_Sample_Names){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Sample_Names<-as.data.frame(var_Sample_Names)}else{var_Sample_Names <- var_Sample_Names}
invisible(graphics.off())
var_macrophage_filtered<-macrophage_filtered(var_Annotate_Cell_Types,var_Metadata_Table,var_Sample_Names)
invisible(graphics.off())
saveRDS(var_macrophage_filtered, paste0(rds_output,"/var_macrophage_filtered.rds"))
