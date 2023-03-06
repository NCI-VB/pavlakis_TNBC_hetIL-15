FIGS5A_CD24_expression <- function(renamed_so) {  # produces 2 color channels and the overlay

require(ggplot2)
require(gridExtra)
library(Seurat)
library(tidyverse)
library(ggnewscale)
    so = renamed_so$value

  
    gg.overlay <- function(so,df,marker1,marker2){
    
        df %>% dplyr::mutate(avg = mark2.scale+mark1.scale) %>% dplyr::arrange(avg) -> df
    
        gg <- ggplot(df, aes(dr1,dr2))+
            geom_point(aes(color=mark1.scale), size=0)+
            scale_color_gradient("SCT",low="black",high=rgb(1,0,0))+
            new_scale_color() +
            geom_point(aes(color=mark2.scale), size=0)+
            scale_color_gradient("Protein", low="black",high=rgb(0,1,0))+
            geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape=16,size=0.5, alpha=0.5)+
            theme_classic() +
            xlab("umap-1") + 
            ylab("umap-2") +
            ggtitle("CD24 Expression") +
            coord_fixed()

            return(gg)
    }

    t1 = 0.5
    t2 = 0.5
    addlines = TRUE

    marker1 <- "Cd24a"
    marker2 <- "CD24"
    scale1 <- TRUE
    scale2 <- TRUE

    mark1 = so@assays$SCT@scale.data[marker1,]
    if(scale1){
    q1 = quantile(mark1,0.80)
    q0 = quantile(mark1,0.5)
    mark1[mark1<q0]=q0
    mark1[mark1>q1]=q1
    }
    mark1.scale <- scales::rescale(mark1, to=c(0,1))

    mark2 = so@assays$Protein@scale.data[marker2,]
    if(scale2){
    q1 = quantile(mark2,0.90)
    q0 = quantile(mark2,0.5)
    print(q0)
    mark2[mark2<q0]=q0
    mark2[mark2>q1]=q1
    }
    mark2.scale <- scales::rescale(mark2, to=c(0,1))

    df <- data.frame(cbind(dr1=so@reductions$umap@cell.embeddings[,1],
                           dr2=so@reductions$umap@cell.embeddings[,2],
                            mark1.scale,mark2.scale))
    p <- gg.overlay(so,df,marker1,marker2)

    imageWidth = 1300
    imageHeight = 1200
    dpi = 300

    png(
      filename="FIGS5A_CD24_expression.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    print(p)
     
# auto removed:     return(NULL)
}

print("template_function_FIGS5A_CD24_expression.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
invisible(graphics.off())
var_FIGS5A_CD24_expression<-FIGS5A_CD24_expression(var_renamed_so)
invisible(graphics.off())
saveRDS(var_FIGS5A_CD24_expression, paste0(rds_output,"/var_FIGS5A_CD24_expression.rds"))
