all_markers <- function(renamed_so) {
    
    suppressMessages(library("Seurat"))
    suppressMessages(library("tidyverse"))

    so <- renamed_so$value

    new_levels <- c("CCR7 High DC","CD103intCD11b+ DC","cDC1","cDC2","Siglec-H DC","Monocyte 1","Monocyte 2","Other")
    levels(so) <- new_levels

    so <- subset(x = so, idents = "Other", invert = TRUE)

    markers <- FindAllMarkers(so, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)

    return(markers)
}

print("template_function_all_markers.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
invisible(graphics.off())
var_all_markers<-all_markers(var_renamed_so)
invisible(graphics.off())
saveRDS(var_all_markers, paste0(rds_output,"/var_all_markers.rds"))
