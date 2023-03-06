# Metadata Table [scRNA-seq] [CCBR] (6a8139b7-45b4-4c6a-8648-c8ec34e6fc60): v31
macrophage_filtered_metadata <- function(macrophage_filtered) {
    
    ## Libraries.
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    ## Cached SO stuff.
    if(FALSE){
        SO <- SO_all  
        #print(SO) 
    } else {
        SO = macrophage_filtered$value
        #print(SO)
    } 
    
    ## Log output.
    #print(SO)

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

## Commenting out the below to unblock training. I believe this is vestigial code and it's started throwing errors.
## Today is 7/27/21. -- Josh

#################################################
## Global imports and functions included below ##
#################################################
# 
#   }
# 
# 

print("template_function_macrophage_filtered_metadata.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_macrophage_filtered<-readRDS(paste0(rds_output,"/var_macrophage_filtered.rds"))
Input_is_Seurat_count <- 0
for(item in var_macrophage_filtered){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_macrophage_filtered<-as.data.frame(var_macrophage_filtered)}else{var_macrophage_filtered <- var_macrophage_filtered}
invisible(graphics.off())
var_macrophage_filtered_metadata<-macrophage_filtered_metadata(var_macrophage_filtered)
invisible(graphics.off())
saveRDS(var_macrophage_filtered_metadata, paste0(rds_output,"/var_macrophage_filtered_metadata.rds"))
