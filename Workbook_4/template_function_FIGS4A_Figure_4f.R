FIGS4A_Figure_4f <- function(DC_CPM_TMM_FILTERED, DC_metadata){
    # image: png
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(ggrepel))

    samples_to_include = c("CS234_cDC1s","CS234_cDC2s","CS234_NewDCs","CCS231_Macro")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]

    DC_CPM_TMM_FILTERED %>% dplyr::select(append("Gene",samples_to_include)) -> spark.df
    
    spark.df %>% dplyr::collect() -> edf.orig
    df.cols <- colnames(edf.orig)
    #print(df.cols)
    edf.orig$Gene <- apply(array(as.character(edf.orig$Gene)),1,function(z) unlist(strsplit(z, "_"))[2])
    edf.orig$sd <- apply(edf.orig[,-1], 1, sd)
    edf.orig <- edf.orig %>% arrange(desc(sd)) %>% select(one_of(df.cols))
    edf.orig <- edf.orig[!duplicated(edf.orig$Gene),]
    #spark_df <- df %>% dplyr::as.DataFrame()
    print(edf.orig)

    Sampinfo <- dplyr::collect(DC_metadata)
    # cell <- Sampinfo$sample_id
    #Sampinfo = Sampinfo[, colSums(is.na(Sampinfo)) == 0]
    rownames(Sampinfo) <- Sampinfo$sample_id

    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo$sample_id), ]
    #Sampinfo = na.omit(Sampinfo)
    Sampinfo = Sampinfo[complete.cases(Sampinfo[, "sample_id"]),]
    print(paste0("Total number of genes in input: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo$sample_id,colnames(edf.orig))]

    idx=!rowMeans(edf)==0
    edf=edf[idx,]
    edf.orig = edf.orig[idx,]
    print(paste0("Total number of genes included: ", nrow(edf.orig)))
    edf <- log2(edf+1)
    print(head(edf))
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- Sampinfo$cell
    pca.df$sample <- Sampinfo$cell
    
    # Manual changes to sample names
    replacements = c("")

    #plotcolors <- c("darkred","purple3","cadetblue","coral","deeppink","darkblue","darkgoldenrod","darkolivegreen3")
    plotcolors <- c("#eb3423","#003df5","#4eac5b","#7f7f7f")
    if (length(unique(Sampinfo$Group)) > length(plotcolors)) {
        plotcolors <- c(plotcolors, rep("black", length(unique(Sampinfo$Group)) - length(plotcolors)))
    }

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0("PC1 ", perc.var[1],"%")
    pc.y.lab <- paste0("PC2 ", perc.var[2],"%")
    
    labelpos <- pca.df
    labelpos$mean_y <- pca.df$PC2+2
    labelpos$mean_x <- pca.df$PC1+2
    print(labelpos$sample)

    labelpos$sample[3] <- "CD103^{int}~~CD11b^{+phantom(0)}~DC"

    if (TRUE){
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="none") +
      geom_point(aes(color=pca.df$celltype), size=14) +
    #   geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
        #   label=sample, color=celltype, vjust="inward", hjust="inward"), size=5, show.legend=FALSE) +
    #   geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
    #       label=sample, color=celltype, vjust=ifelse(PC2 > 0, 2, -0.5), hjust=ifelse(PC1 > 0, 1, -0.1)), size=13, show.legend=FALSE,parse=TRUE) +
      theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size=40),
                axis.title = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
    }

    else{
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }

    print(g)

    rownames(edf)=edf.orig[,1]
    edf.df = as.data.frame(edf)
    edf.df %>% rownames_to_column("Gene") -> edf.df
    return(edf.df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_FIGS4A_Figure_4f.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DC_CPM_TMM_FILTERED<-readRDS(paste0(rds_output,"/var_DC_CPM_TMM_FILTERED.rds"))
Input_is_Seurat_count <- 0
for(item in var_DC_CPM_TMM_FILTERED){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DC_CPM_TMM_FILTERED<-as.data.frame(var_DC_CPM_TMM_FILTERED)}else{var_DC_CPM_TMM_FILTERED <- var_DC_CPM_TMM_FILTERED}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DC_metadata<-readRDS(paste0(rds_output,"/var_DC_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_DC_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DC_metadata<-as.data.frame(var_DC_metadata)}else{var_DC_metadata <- var_DC_metadata}
invisible(graphics.off())
var_FIGS4A_Figure_4f<-FIGS4A_Figure_4f(var_DC_CPM_TMM_FILTERED,var_DC_metadata)
invisible(graphics.off())
saveRDS(var_FIGS4A_Figure_4f, paste0(rds_output,"/var_FIGS4A_Figure_4f.rds"))
