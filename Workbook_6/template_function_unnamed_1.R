unnamed_1 <- function(ADT_Raw_Counts) {
    
    df <- ADT_Raw_Counts

    df[,c(1,2)] <- NULL

    print(apply(df,2,max))
# auto removed:     return(NULL)

}

print("template_function_unnamed_1.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ADT_Raw_Counts<-readRDS(paste0(rds_output,"/var_ADT_Raw_Counts.rds"))
Input_is_Seurat_count <- 0
for(item in var_ADT_Raw_Counts){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ADT_Raw_Counts<-as.data.frame(var_ADT_Raw_Counts)}else{var_ADT_Raw_Counts <- var_ADT_Raw_Counts}
invisible(graphics.off())
var_unnamed_1<-unnamed_1(var_ADT_Raw_Counts)
invisible(graphics.off())
saveRDS(var_unnamed_1, paste0(rds_output,"/var_unnamed_1.rds"))
