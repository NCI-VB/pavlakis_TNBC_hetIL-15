cd11_heatmap <- function(renamed_so) {
    
    suppressMessages(library("Seurat"))
    suppressMessages(library("tidyverse"))
    library(colorspace)
    so <- renamed_so$value
    
    new_levels <- c("CCR7 High DC","CD103intCD11b+ DC","cDC1","cDC2","Siglec-H DC","Monocyte 1","Monocyte 2","Other")
    levels(so) <- new_levels

    # so@active.ident <- factor(as.character(so@active.ident),levels = new_levels)
    # so[["ident"]] <- factor(as.character(so[["ident"]]), levels = new_levels)

    so <- subset(x = so, idents = "Other", invert = TRUE)

    markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    default.colors <- qualitative_hcl(n=length(x = levels(x = so)), palette = "Dark 2")
    default.colors[1] <- "#e62ee2"
    default.colors[2] <- "#4eac5b"
    default.colors[3] <- "#eb3423"
    default.colors[4] <- "#003df5"
    default.colors[5] <- "#e68643"
    default.colors[6] <- default.colors[7]
    default.colors[7] <- "#34c3eb"
    p <- DoHeatmap(so, features = top10$gene, angle=45, group.colors=default.colors) + NoLegend()+
        theme(axis.text.y = element_text(size=8))
    print(p)
    return(markers)

}

print("template_function_cd11_heatmap.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
invisible(graphics.off())
var_cd11_heatmap<-cd11_heatmap(var_renamed_so)
invisible(graphics.off())
saveRDS(var_cd11_heatmap, paste0(rds_output,"/var_cd11_heatmap.rds"))
