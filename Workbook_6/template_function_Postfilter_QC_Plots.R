# Post-filter QC Plots [scRNA-seq] [CCBR] (44f7a498-aa69-429d-a9e8-4b6c8650bb05): v37
Postfilter_QC_Plots <- function(Filter_and_QC_Samples) {
    #image: png
    imageType = "png"
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    suppressMessages(library(gtable))
    suppressMessages(library(gridExtra))
    suppressMessages(library(reshape2))
    suppressMessages(library(RColorBrewer))

    obj.list <- Filter_and_QC_Samples$value
    #obj.list <- so.list

    all.columns <- unique(unlist(sapply(seq_along(obj.list), function(i) colnames(obj.list[[i]]@meta.data))))
    qc.df <- array(0,dim=c(0,3))
    for (i in 1:length(obj.list)){
        so <- obj.list[[i]]
        #print(so)

        #Add missing columns to metadata
        missing.columns <- setdiff(all.columns,colnames(so@meta.data))
        
        for(i in missing.columns){
            so <- AddMetaData(so,rep(0,ncol(so)), i)
        }
        df.m <- melt(so@meta.data)
        qc.df <- rbind(qc.df,df.m)
    }

    qfilter <- function(x){
        library(dplyr)
        qc.df %>% dplyr::filter(variable == x)
    }

    col1=brewer.pal(8, "Set3")[-2] 
    col2=c(col1,brewer.pal(8,"Set2")[3:6])
    col3=c(col2,"#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9","#808080","#A9A9A9","#8B7355")

    plothist <- function(count.df){
    g=ggplot(count.df) + 
    theme_bw() +
    geom_density(aes(x = value, colour = orig.ident)) +
    labs(x = NULL) +
    theme(legend.position='right',legend.text=element_text(size=10),
    legend.title=element_blank()) + 
    ggtitle(count.df$variable[1]) +
    scale_x_continuous(trans='log10') + 
    scale_color_manual(values = col3) 
    #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
    return(g)
    #print(g)
# auto removed:     #return(NULL)
    }

    plotviolin <- function(count.df){
    axislab = unique(count.df$orig.ident)

    v=ggplot(count.df, aes(x=orig.ident, y=value)) +
    ggtitle(count.df$variable[1]) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.text=element_text(size=rel(1)),
          legend.title=element_blank(), axis.text=element_text(size=10),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          #axis.text.x=element_text(angle=45,hjust=1),
          plot.title = element_text(size = 20, face = "bold")) +
    geom_violin(aes(fill=as.factor(orig.ident))) +  
    scale_fill_manual(values = col3) +
    geom_boxplot(width=.1) +
    #labs(colour = n, y=m) +
    #geom_jitter(height = 0, width = 0.1, size = 0.1) +
    #scale_y_continuous(trans='log2') + 
    scale_x_discrete(limits = as.vector(axislab)) 
    return(v)
    #print(v)
# auto removed:     #return(NULL)
    }

    plotscatter <- function(count.df,counts){
            count.df %>% dplyr::mutate("value2"=counts) -> count.df 
            ylab = as.character(unique(count.df$variable))
            xlab = "RNA Count"
            name = paste(ylab,"vs.",xlab)          
            g <- ggplot(count.df, aes(x=value2, y=value,color = orig.ident)) +
                geom_point(size = 0.5) + 
                theme_classic() +
                theme(legend.position='right',legend.text=element_text(size=10),
    legend.title=element_blank()) + 
                guides(colour = guide_legend(override.aes = list(size=2))) +
                scale_color_manual(values = col3) +
                labs(title=name, x = xlab, y = ylab)
                ggtitle(name) 
            #print(g)
# auto removed:             #return(NULL)
            return(g)
        }
    
    
    
    useSpark <- FALSE

    if (useSpark) {
        qc.count <- spark.lapply(unique(qc.df$variable), function(x) {qfilter(x)})
        qc.count[[1]] %>% dplyr::filter(variable=="nCount_RNA") %>% pull(value) -> RNAcounts
    } else {
        qc.count <- lapply(unique(qc.df$variable), function(x) {qfilter(x)})
        qc.count[[1]] %>% dplyr::filter(variable=="nCount_RNA") %>% pull(value) -> RNAcounts
    }
    grobs <- lapply(seq_along(qc.count), function(x) {arrangeGrob(grobs = list(plotscatter(qc.count[[x]],RNAcounts),plothist(qc.count[[x]]),plotviolin(qc.count[[x]])),nrow=1,ncol=3)})
    
    imageWidth = 5000
    imageHeight = 1000*length(grobs)
    dpi = 300

    if (imageType == 'png') {
    png(
      filename="Postfilter_QC_Plots.png",
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
        file="Postfilter_QC_Plots.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    grid.arrange(grobs = grobs, nrow = length(grobs))

so@meta.data %>% rownames_to_column("Barcode") -> meta.df

cat("\nReturn objects checksum:\n")
print(digest::digest(obj.list))

return(list(value=obj.list))

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

print("template_function_Postfilter_QC_Plots.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Filter_and_QC_Samples<-readRDS(paste0(rds_output,"/var_Filter_and_QC_Samples.rds"))
Input_is_Seurat_count <- 0
# for(item in var_Filter_and_QC_Samples){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
# if(Input_is_Seurat_count == 0 ){
# var_Filter_and_QC_Samples<-as.data.frame(var_Filter_and_QC_Samples)}else{var_Filter_and_QC_Samples <- var_Filter_and_QC_Samples}
invisible(graphics.off())
var_Postfilter_QC_Plots<-Postfilter_QC_Plots(var_Filter_and_QC_Samples)
invisible(graphics.off())
saveRDS(var_Postfilter_QC_Plots, paste0(rds_output,"/var_Postfilter_QC_Plots.rds"))
