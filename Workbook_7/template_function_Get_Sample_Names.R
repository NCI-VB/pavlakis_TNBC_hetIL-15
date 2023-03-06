Get_Sample_Names <- function(Labeled_so) {
    
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    if(FALSE){
        SO <- SO_all
        print(SO)
    } else {
        SO = Labeled_so
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

#################################################
## Global imports and functions included below ##
#################################################




    

print("template_function_Get_Sample_Names.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
nidap_output <- paste0(currentdir,'/nidap_downloads')
var_Labeled_so<-readRDS(paste0(nidap_output,'/RObjectdata.rds'))
Input_is_Seurat_count <- 0

if(Input_is_Seurat_count == 0 ){
}
invisible(graphics.off())
var_Get_Sample_Names<-Get_Sample_Names(var_Labeled_so)
invisible(graphics.off())
saveRDS(var_Get_Sample_Names, paste0(rds_output,"/var_Get_Sample_Names.rds"))
