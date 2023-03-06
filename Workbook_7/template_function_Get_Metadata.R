Get_Metadata <- function(Labeled_so) {
    
    ## Libraries.
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    ## Cached SO stuff.
    if(FALSE){
        SO <- SO_all  
        print(SO) 
    } else {
        SO = Labeled_so
        print(SO)
    } 
    
    ## Log output.
    print(SO)

    ## Cell embeddings stuff.
    if (FALSE) {
        reds = names(SO@reductions)
        for (i in seq_along(reds)){
            SO = Seurat::AddMetaData(SO,as.data.frame(SO@reductions[[i]]@cell.embeddings))
        }
    }

    ## Checking for where to pull meta.data.
    if("meta.data" %in% slotNames(SO)) {
        met.df <- SO@meta.data 
    } else {
        met.df <- SO$RNA@meta.data
    }
    
    ## If no barcode column, get them from the rownames.
    if (!("Barcode" %in% colnames(met.df))) {
        met.df %>% rownames_to_column("Barcode") -> met.df
    }
    
    ## Return metadata table.
    return(met.df)
}

#################################################
## Global imports and functions included below ##
#################################################




print("template_function_Get_Metadata.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
nidap_output <- paste0(currentdir,'/nidap_downloads')
var_Labeled_so<-readRDS(paste0(nidap_output,'/RObjectdata.rds'))
Input_is_Seurat_count <- 0

if(Input_is_Seurat_count == 0 ){
}
invisible(graphics.off())
var_Get_Metadata<-Get_Metadata(var_Labeled_so)
invisible(graphics.off())
saveRDS(var_Get_Metadata, paste0(rds_output,"/var_Get_Metadata.rds"))
