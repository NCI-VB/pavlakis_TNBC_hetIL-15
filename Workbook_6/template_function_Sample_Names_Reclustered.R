# Sample Names [scRNA-seq] [CCBR] (51ea15ce-2be0-415f-902b-3c86175eb6cd): v37
Sample_Names_Reclustered <- function(reclusterd_so) {
    
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    if(FALSE){
        SO <- SO_all
        print(SO)
    } else {
        SO = reclusterd_so$value
        print(SO)
    }
    
    if("meta.data" %in% slotNames(SO)){
        if ("orig.ident" %in% colnames(SO@meta.data)) {
            dim.df <- as.data.frame(t(as.matrix(table(SO@meta.data$orig.ident))))
        } else {
            dim.df <- as.data.frame(t(as.matrix(table(SO@meta.data$orig_ident))))
        }
        colnames(dim.df) <- lapply(colnames(dim.df), function(x) gsub(".h5", "", x))   
    } else {
        if ("orig.ident" %in% colnames(SO$RNA@meta.data)) {
            dim.df <- as.data.frame(t(as.matrix(table(SO$RNA@meta.data$orig.ident))))
        } else {
            dim.df <- as.data.frame(t(as.matrix(table(SO$RNA@meta.data$orig_ident))))
        }
        dim.df <- as.data.frame(t(as.matrix(table(SO$RNA@meta.data$orig.ident))))
        colnames(dim.df) <- lapply(colnames(dim.df), function(x) gsub(".h5", "", x))
    }

    return(dim.df)
}

## Commenting out the below to unblock training. I believe this is vestigial code and it's started throwing errors.
## Today is 7/27/21. -- Josh

#################################################
## Global imports and functions included below ##
#################################################
# 
#   }
# 
# 

    

print("template_function_Sample_Names_Reclustered.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_reclusterd_so<-readRDS(paste0(rds_output,"/var_reclusterd_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_reclusterd_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_reclusterd_so<-as.data.frame(var_reclusterd_so)}else{var_reclusterd_so <- var_reclusterd_so}
invisible(graphics.off())
var_Sample_Names_Reclustered<-Sample_Names_Reclustered(var_reclusterd_so)
invisible(graphics.off())
saveRDS(var_Sample_Names_Reclustered, paste0(rds_output,"/var_Sample_Names_Reclustered.rds"))
