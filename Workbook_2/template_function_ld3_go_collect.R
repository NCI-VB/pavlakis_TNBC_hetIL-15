# [CCBR] Collect L2P Pathway Gene Hits (1c6d4aa3-b246-461e-82f6-b82946fcb785): v8
ld3_go_collect <- function(ld3_go_pathway_up) {
    ## Set-up.
    suppressMessages(library(dplyr))
    df <- dplyr::collect(ld3_go_pathway_up)

    ## Make path_id, subset df, transpose, and move paths to colnames.
    df$path_ID <- paste(df$pathwayname, df$pathwayaccessionidentifier, sep = "_")
    df %>% dplyr::select(path_ID, genesinpathway) -> df2
    rownames(df2) <- df2$path_ID
    df2$path_ID <- NULL
    df3 = t(df2)
    pathways <- colnames(df3)

    ## Expand gene lists to create list of lists.
    pathways_gene_lists <- list()
    for(pathway in pathways){
        genes <- strsplit(df3[1,pathway], " ")
        pathways_gene_lists[pathway] <- genes
        }
    
    ## Create matrix big enough for final needs containing all NAs.
    colnum = length(pathways)
    rownum = max(unlist(lapply(pathways_gene_lists,FUN=length)))
    pathMat = matrix(nrow=rownum,ncol=colnum)
    pathMat[] = NA

    ## Replace NAs in each column with genes from one path, add path names as colnames, and convert to data frame.
    for (i in 1:colnum){
        pathLen = length(pathways_gene_lists[[pathways[i]]])
        pathMat[1:pathLen,i] = pathways_gene_lists[[pathways[i]]]
        }
    pathways = gsub("[ ,;:{}()/\n\t=]", "_", pathways)
    colnames(pathMat) = pathways
    path_df = as.data.frame(pathMat)

return(path_df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_ld3_go_collect.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ld3_go_pathway_up<-readRDS(paste0(rds_output,"/var_ld3_go_pathway_up.rds"))
Input_is_Seurat_count <- 0
for(item in var_ld3_go_pathway_up){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ld3_go_pathway_up<-as.data.frame(var_ld3_go_pathway_up)}else{var_ld3_go_pathway_up <- var_ld3_go_pathway_up}
invisible(graphics.off())
var_ld3_go_collect<-ld3_go_collect(var_ld3_go_pathway_up)
invisible(graphics.off())
saveRDS(var_ld3_go_collect, paste0(rds_output,"/var_ld3_go_collect.rds"))
