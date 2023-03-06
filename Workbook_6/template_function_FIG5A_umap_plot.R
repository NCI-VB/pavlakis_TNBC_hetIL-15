FIG5A_umap_plot <- function(renamed_so) {
    
    library(Seurat)
    library(scales)
    library(colorspace)

    so <- renamed_so$value

    new_levels <- c("CCR7 High DC","CD103intCD11b+ DC","cDC1","cDC2","Siglec-H DC","Monocyte 1","Monocyte 2","Other")
    levels(so) <- new_levels
    
        #   Red: #eb3423 (cDC1)
    #   Blue: #003df5 (cDC2)
    #   Green #4eac5b (New)
    #   Gold: lightgoldenrod

    default.colors <- c(hue_pal()(length(x = levels(x = so))))
    default.colors <- qualitative_hcl(n=length(x = levels(x = so)), palette = "Dark 2")
    default.colors[1] <- "#e62ee2"
    default.colors[2] <- "#4eac5b"
    default.colors[3] <- "#eb3423"
    default.colors[4] <- "#003df5"
    default.colors[5] <- "#e68643"
    default.colors[6] <- default.colors[7]
    default.colors[7] <- "#34c3eb"
    p <- DimPlot(so, reduction = "umap", label = FALSE, pt.size = 1, cols = default.colors )
    p <- p + theme(axis.title = element_blank(),
                    #axis.text = element_blank(),
                    axis.line = element_line(colour = 'black', size = 1),
                    axis.ticks = element_line(colour = 'black', size = 1),)
    print(p)

# auto removed:     return(NULL)
}

print("template_function_FIG5A_umap_plot.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
invisible(graphics.off())
var_FIG5A_umap_plot<-FIG5A_umap_plot(var_renamed_so)
invisible(graphics.off())
saveRDS(var_FIG5A_umap_plot, paste0(rds_output,"/var_FIG5A_umap_plot.rds"))
