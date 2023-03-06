FIG5D_gsva_split_treatment <- function(renamed_so, msigdb_v6_2_with_orthologs) {
    
    suppressMessages(library(Seurat))
    suppressMessages(library(GSVA))
    suppressMessages(library(tidyverse))
    suppressMessages(library(pheatmap))
    
    object <- renamed_so$value
    db <- msigdb_v6_2_with_orthologs

    assay <- "RNA"
    slot  <- "counts"
    clusters_of_interest <- "cell_identities"
    species_of_interest <- "Mouse"
    collections_of_interest <- unique(db$collection) #"BP: GO biological process"# unique(db$collection) # "BP: GO biological process" #
    min.sz <- 15
    max.sz <- 500
    max_pathways <- 50

    # set idents
    Idents(object) <- clusters_of_interest
    object <- subset(object, idents = "Other", invert = TRUE)

    new_idents <- paste(object@meta.data$orig.ident,object@meta.data$cell_identities, sep="_")
    object <- AddMetaData(object, new_idents, "merged_idents")
    Idents(object) <- "merged_idents"
    # get the data
    raw_data <- Seurat::GetAssayData(object, assay = assay, slot = slot)

    # get the identis
    cell_ids <- as.character( Seurat::Idents(object) )

    av_counts <- apply(raw_data, 1, function(row_data) {
        by(row_data, cell_ids, mean)
    })

    # convert to a data.frame
    av_counts <- t( data.frame(av_counts) )

    # prepare gene sets
    db <- db %>% dplyr::filter(species == species_of_interest, collection %in% collections_of_interest)
    db <- db %>% group_by(collection, gene_set_name) %>% summarize(gene_symbol = as.list(strsplit(paste0(unique(gene_symbol), collapse = " "), " ")))
    Ldb <- db$gene_symbol
    names(Ldb) = db$gene_set_name

    filter <- duplicated(names(Ldb))
    Ldb <- Ldb[!filter]
    
    ## run gsea
    set.seed(42)
    gsva.es <- gsva(av_counts, Ldb, method = "ssgsea", min.sz = min.sz, max.sz = max.sz, verbose=FALSE)

    
    dendritic_cell_pathways <- grep("dendritic",row.names(gsva.es), value = TRUE, ignore.case = TRUE)
    print(dendritic_cell_pathways)
    dendritic_cell_pathways <- dendritic_cell_pathways[!grepl("spine",dendritic_cell_pathways,ignore.case=TRUE)]
    dendritic_cell_pathways <- dendritic_cell_pathways[!grepl("axo_dendritic",dendritic_cell_pathways,ignore.case=TRUE)]
    dendritic_cell_pathways <- dendritic_cell_pathways[!grepl("somatodendritic",dendritic_cell_pathways,ignore.case=TRUE)]
    dendritic_cell_pathways <- dendritic_cell_pathways[!grepl("shaft",dendritic_cell_pathways,ignore.case=TRUE)]

    gsva.es <- gsva.es[dendritic_cell_pathways, ]

    pathway_names <- row.names(gsva.es)
    # find the maximum differently expressed pathway

    max_difference <- do.call(rbind, apply(gsva.es, 1, function(row) {
        values <- as.numeric(row)
        return(data.frame(min = min(values), max = max(values)))
    }))
    max_difference$name <- pathway_names
    max_difference <- max_difference %>% select(name,min,max)

    max_difference$diff <- max_difference$max - max_difference$min

    # sort based on the difference
    max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

    gsva.es <- gsva.es[max_difference$name, ]
    gsva.t <- as.data.frame(t(gsva.es))

    max_pathways <- min(max_pathways,length(max_difference$name))
    pathway_ids <- max_difference$name[1:max_pathways]
    
    # get the expression values as a matrix
    expression_matrix <- as.matrix(gsva.es[pathway_ids, ])

    # merge default with user parameters
    heatmap_params <- list(
        mat = expression_matrix,
        #color = rev(RColorBrewer::brewer.pal(9, "RdYlBu")),
        kmeans_k = NA,
        cellwidth = 20,
        cellheight = 20,
        scale = "row",
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        show_rownames = T, show_colnames = T)
    
    p <- do.call(pheatmap, heatmap_params)

    print(p)

    gsva.t <- gsva.t %>% rownames_to_column("Barcode")
    return(gsva.t)
}

print("template_function_FIG5D_gsva_split_treatment.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_msigdb_v6_2_with_orthologs<-readRDS(paste0(rds_output,"/var_msigdb_v6_2_with_orthologs.rds"))
Input_is_Seurat_count <- 0
for(item in var_msigdb_v6_2_with_orthologs){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_msigdb_v6_2_with_orthologs<-as.data.frame(var_msigdb_v6_2_with_orthologs)}else{var_msigdb_v6_2_with_orthologs <- var_msigdb_v6_2_with_orthologs}
invisible(graphics.off())
var_FIG5D_gsva_split_treatment<-FIG5D_gsva_split_treatment(var_renamed_so,var_msigdb_v6_2_with_orthologs)
invisible(graphics.off())
saveRDS(var_FIG5D_gsva_split_treatment, paste0(rds_output,"/var_FIG5D_gsva_split_treatment.rds"))
