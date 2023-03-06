# Sample Names [scRNA-seq] [CCBR] (51ea15ce-2be0-415f-902b-3c86175eb6cd): v31
Sample_Names <- function(Annotate_Cell_Types) {
    
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    if(FALSE){
        SO <- SO_all
        print(SO)
    } else {
        SO = Annotate_Cell_Types$value
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

    

print("template_function_Sample_Names.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Annotate_Cell_Types<-readRDS(paste0(rds_output,"/var_Annotate_Cell_Types.rds"))
Input_is_Seurat_count <- 0
for(item in var_Annotate_Cell_Types){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Annotate_Cell_Types<-as.data.frame(var_Annotate_Cell_Types)}else{var_Annotate_Cell_Types <- var_Annotate_Cell_Types}
invisible(graphics.off())
var_Sample_Names<-Sample_Names(var_Annotate_Cell_Types)
invisible(graphics.off())
saveRDS(var_Sample_Names, paste0(rds_output,"/var_Sample_Names.rds"))
