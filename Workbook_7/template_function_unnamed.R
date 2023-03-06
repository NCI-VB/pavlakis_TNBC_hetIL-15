unnamed <- function(Labeled_so, Get_Metadata) {
    
    ## Libraries.
    suppressMessages(library(Seurat))
    suppressMessages(library(tidyverse))
    suppressMessages(library(reshape2))

    ## Inputs.
    SO <- Labeled_so # Imports SO of interest
    meta.df <- dplyr::collect(Get_Metadata)
    selected_column <- "SCT_snn_res_0_8_merged"
    colnames(SO@meta.data) <- gsub("\\.","_",colnames(SO@meta.data))

    ## Find user-selected column with cell types in it.
    search_vec <- colnames(SO@meta.data)
    CellType <- SO@meta.data[grep(selected_column, search_vec, fixed = TRUE)]

    ## Unprocessed cell count table.
    cell_table <- as.data.frame(table(CellType[selected_column]))
    ## First convert to a diagonal matrix.
    count_matrix <- dcast(cell_table, Var1 ~ Var1)[,-1]
    ## ID for grouping.
    ID <- 1
    bound_matrix <- cbind(ID,count_matrix)

    # Final table to display
    Display <- bound_matrix %>% group_by(ID) %>% summarise_each(funs(sum(., na.rm = TRUE))) %>% subset(select = c(-1))}
    


#################################################
## Global imports and functions included below ##
#################################################




    

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.




print("template_function_unnamed.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
nidap_output <- paste0(currentdir,'/nidap_downloads')
var_Labeled_so<-readRDS(paste0(nidap_output,'/RObjectdata.rds'))
Input_is_Seurat_count <- 0

if(Input_is_Seurat_count == 0 ){
}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Get_Metadata<-readRDS(paste0(rds_output,"/var_Get_Metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_Get_Metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Get_Metadata<-as.data.frame(var_Get_Metadata)}else{var_Get_Metadata <- var_Get_Metadata}
invisible(graphics.off())
var_unnamed<-unnamed(var_Labeled_so,var_Get_Metadata)
invisible(graphics.off())
saveRDS(var_unnamed, paste0(rds_output,"/var_unnamed.rds"))
