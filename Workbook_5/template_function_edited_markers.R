edited_markers <- function(Canonical_dc_markers) {
    
    df <- Canonical_dc_markers

    tmp <- df[1,]
    tmp[1,c(1,4)] <- "Itgae"

    df <- rbind(df,tmp)
    return(df)
}

print("template_function_edited_markers.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Canonical_dc_markers<-readRDS(paste0(rds_output,"/var_Canonical_dc_markers.rds"))
Input_is_Seurat_count <- 0
for(item in var_Canonical_dc_markers){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Canonical_dc_markers<-as.data.frame(var_Canonical_dc_markers)}else{var_Canonical_dc_markers <- var_Canonical_dc_markers}
invisible(graphics.off())
var_edited_markers<-edited_markers(var_Canonical_dc_markers)
invisible(graphics.off())
saveRDS(var_edited_markers, paste0(rds_output,"/var_edited_markers.rds"))
