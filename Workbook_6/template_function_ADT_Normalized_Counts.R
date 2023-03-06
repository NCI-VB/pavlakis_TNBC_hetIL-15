ADT_Normalized_Counts <- function(renamed_so) {
    
    library(Seurat)
    library(tidyverse)
    library(stringr)

    so <- renamed_so$value

    df <- so@assays$Protein@data

    df.t <- as.data.frame(x = t(x = as.matrix(x = df))) %>% rownames_to_column("Barcode") #Ab in colnames

    df.t <- df.t %>% mutate(Group = so@meta.data[ Barcode, "ident"],Sample = unlist(str_split(Barcode,"_"))[1]) %>% select(Barcode,Sample,Group,everything())

    return(df.t)
}

print("template_function_ADT_Normalized_Counts.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
invisible(graphics.off())
var_ADT_Normalized_Counts<-ADT_Normalized_Counts(var_renamed_so)
invisible(graphics.off())
saveRDS(var_ADT_Normalized_Counts, paste0(rds_output,"/var_ADT_Normalized_Counts.rds"))
