FIG5B_renamed_density <- function(renamed_so, unnamed_20) {
    #image: png

    library(Seurat)
    library(ggplot2)
    library(viridis)
    library(scales)
    library(tidyverse)
    library(cowplot)

    so <- renamed_so$value
    print(colnames(so[[]]))
    metadata_column <- "ident"
    clusters <- levels(so@meta.data[[metadata_column]])

    sample_name = as.factor(so@meta.data[[metadata_column]])
    names(sample_name)=names(so@active.ident)
    so@active.ident <- as.factor(vector())
    so@active.ident <- sample_name

    counts <- as.data.frame(as.data.frame.matrix(table(so@meta.data[[metadata_column]],so@meta.data$orig.ident)))

    counts$total <- apply(counts,1,sum)

    counts$control_frequency <- counts$Control / counts$total
    counts$treated_frequency <- counts$Treated / counts$total
    counts$score <- scales::rescale(-2*counts$control_frequency+1, from=c(-1,1), to=c(0,1))
    print(counts)
    embeddings <- as.data.frame(Embeddings(so, reduction = "umap"))

    embeddings$score <- apply(array(row.names(embeddings)), 1, function(z) counts[as.character(so@meta.data[z,metadata_column]), "score"])
    
    p <- ggplot(embeddings, aes(x=UMAP_1,y=UMAP_2,color=score))+
        theme_cowplot() +
        geom_point()+
        scale_color_viridis()+
        theme(legend.key.height = unit(2.5, "cm"),
                    axis.title = element_blank(),
                    #axis.text = element_blank(),
                    axis.line = element_line(colour = 'black', size = 1),
                    axis.ticks = element_line(colour = 'black', size = 1))

        #guides(override.aes = list(size=4),colour=guide_legend(ncol=4)) + 
        #guides(color = guide_legend(override.aes = list(size = 6)))
    print(p)
# auto removed:     return(NULL)
}

print("template_function_FIG5B_renamed_density.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_unnamed_20<-readRDS(paste0(rds_output,"/var_unnamed_20.rds"))
Input_is_Seurat_count <- 0
for(item in var_unnamed_20){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_unnamed_20<-as.data.frame(var_unnamed_20)}else{var_unnamed_20 <- var_unnamed_20}
invisible(graphics.off())
var_FIG5B_renamed_density<-FIG5B_renamed_density(var_renamed_so,var_unnamed_20)
invisible(graphics.off())
saveRDS(var_FIG5B_renamed_density, paste0(rds_output,"/var_FIG5B_renamed_density.rds"))
