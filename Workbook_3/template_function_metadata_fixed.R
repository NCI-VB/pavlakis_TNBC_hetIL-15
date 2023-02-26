metadata_fixed <- function(metadata) {
  library(tidyverse)
  tab <- dplyr::collect(metadata)
  tab <- tab %>% mutate(sample_name = gsub("\\s+","_", sample_name))
  tab <- tab %>% mutate(combined_treatment = paste(treatment,timepoint,sep="_"))
  return(tab)  
}

print("template_function_metadata_fixed.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_metadata<-readRDS(paste0(rds_output,"/var_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_metadata<-as.data.frame(var_metadata)}else{var_metadata <- var_metadata}
invisible(graphics.off())
var_metadata_fixed<-metadata_fixed(var_metadata)
invisible(graphics.off())
saveRDS(var_metadata_fixed, paste0(rds_output,"/var_metadata_fixed.rds"))
