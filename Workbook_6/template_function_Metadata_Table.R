# Metadata Table [scRNA-seq] [CCBR] (6a8139b7-45b4-4c6a-8648-c8ec34e6fc60): v31
Metadata_Table <- function(Annotate_Cell_Types) {
    
    ## Libraries.
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    
    ## Cached SO stuff.
    if(FALSE){
        SO <- SO_all  
        #print(SO) 
    } else {
        SO = Annotate_Cell_Types$value
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

print("template_function_Metadata_Table.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Annotate_Cell_Types<-readRDS(paste0(rds_output,"/var_Annotate_Cell_Types.rds"))
Input_is_Seurat_count <- 0
for(item in var_Annotate_Cell_Types){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Annotate_Cell_Types<-as.data.frame(var_Annotate_Cell_Types)}else{var_Annotate_Cell_Types <- var_Annotate_Cell_Types}
invisible(graphics.off())
var_Metadata_Table<-Metadata_Table(var_Annotate_Cell_Types)
invisible(graphics.off())
saveRDS(var_Metadata_Table, paste0(rds_output,"/var_Metadata_Table.rds"))
