# [CCBR] Select Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v11
remove_ercc <- function(count_table_corrected) {
    library(tidyverse)
    tab <- dplyr::collect(count_table_corrected)
    if(FALSE){
        x=paste0("^", c("POS","NEG"), "$", collapse = "|")    
    }
    else{
        x=paste(c("POS","NEG"),collapse="|")
    }
    if(FALSE){
    tab <- dplyr::filter(tab,grepl(x,Name,ignore.case = TRUE) )
    }
    else{
    tab <- dplyr::filter(tab,!grepl(x,Name,ignore.case = TRUE) )    
    }
    return(tab)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_remove_ercc.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_count_table_corrected<-readRDS(paste0(rds_output,"/var_count_table_corrected.rds"))
Input_is_Seurat_count <- 0
for(item in var_count_table_corrected){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_count_table_corrected<-as.data.frame(var_count_table_corrected)}else{var_count_table_corrected <- var_count_table_corrected}
invisible(graphics.off())
var_remove_ercc<-remove_ercc(var_count_table_corrected)
invisible(graphics.off())
saveRDS(var_remove_ercc, paste0(rds_output,"/var_remove_ercc.rds"))
