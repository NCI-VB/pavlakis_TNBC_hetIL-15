count_table_corrected <- function(count_table) {
  library(tidyverse)
  tab <- dplyr::collect(count_table)
  colnames(tab) <- gsub("\\s+","_", colnames(tab))
  colnames(tab) <- gsub("__","_", colnames(tab))
  return(tab)      
}

print("template_function_count_table_corrected.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_count_table<-readRDS(paste0(rds_output,"/var_count_table.rds"))
Input_is_Seurat_count <- 0
for(item in var_count_table){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_count_table<-as.data.frame(var_count_table)}else{var_count_table <- var_count_table}
invisible(graphics.off())
var_count_table_corrected<-count_table_corrected(var_count_table)
invisible(graphics.off())
saveRDS(var_count_table_corrected, paste0(rds_output,"/var_count_table_corrected.rds"))
