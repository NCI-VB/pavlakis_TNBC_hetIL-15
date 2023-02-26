FIGS2D_S4E_L_leukocyte_migration_new <- function(ld3_go_collect,DEG_lymph, metadata_fixed, Ortholog_Map_for_RNA_Seq) {
    # image: png
    suppressMessages(library(dplyr))
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))

    orthology_table = Ortholog_Map_for_RNA_Seq
    table_to_convert = ld3_go_collect

    gene_column_to_convert = "leukocyte_migration_GO_0050900"
    organism = "Mouse"
    if(organism != "Human"){
        organism_of_interest = "Mouse"
        orthology_table = dplyr::filter(orthology_table,           
        orthology_table$Organism==organism_of_interest)
        if ("from human"=="from human") {
            orthology_reference_column = "Human_gene_name"
            orthology_conversion_column = "Nonhuman_gene_name"
        } else {
            orthology_reference_column = "Nonhuman_gene_name"
            orthology_conversion_column = "Human_gene_name"
    }
        orthology_table %>% dplyr::rename("orthology_reference" = orthology_reference_column) %>%(orthology_reference_column,           "orthology_reference") %>%
dplyr::rename("orthology_conversion" = orthology_conversion_column ) %>% dplyr::select("orthology_reference", "orthology_conversion") -> orthology_table
        table_to_convert %>% dplyr::join(orthology_table, table_to_convert[[gene_column_to_convert]]==orthology_table$"orthology_reference") -> table_to_convert
    table_to_convert %>%
        dplyr::drop("orthology_reference") %>%
        dplyr::withColumn(gene_column_to_convert, table_to_convert$"orthology_conversion") %>%
        dplyr::drop("orthology_conversion") -> table_to_convert
        genes = table_to_convert %>% dplyr::select("leukocyte_migration_GO_0050900")
        genes = dplyr::distinct(genes)
    } else {
    genes = table_to_convert %>% dplyr::select("leukocyte_migration_GO_0050900") 
    genes = dplyr::distinct(genes)  
    }
    
    pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
        if (n < 1L) 
            return(character(0L))
        h <- rep(h, length.out = 2L)
        c <- c[1L]
        l <- rep(l, length.out = 2L)
        power <- rep(power, length.out = 2L)
        rval <- seq(1, -1, length = n)
        rval <- hex(
            polarLUV(
                L = l[2L] - diff(l) * abs(rval)^power[2L], 
                C = c * abs(rval)^power[1L],
                H = ifelse(rval > 0, h[1L], h[2L])
            ),
            fixup=fixup, ...
        )
        if (!missing(alpha)) {
            alpha <- pmax(pmin(alpha, 1), 0)
            alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                width = 2L, upper.case = TRUE)
            rval <- paste(rval, alpha, sep = "")
        }
        return(rval)
    }
    #Color selections for heatmap:
    np0 = pal(100)
    np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1)  #Blue to Red
    np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2))  #Red to Vanilla
    np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) #Violet to Pink
    np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
    np5 = colorRampPalette(c("steelblue","white", "red"))(100) #Steelblue to White to Red

    np = list(np0, np1, np2, np3, np4, np5)
    names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")

    doheatmap <- function(dat, clus, clus2, rn, cn, col) {
        require(pheatmap)
        require(dendsort)
        if (TRUE) {
            tmean.scale = t(scale(t(dat)))
            tmean.scale = tmean.scale[is.finite(rowSums(tmean.scale)),]
        } else {
            tmean.scale = dat
        }
        col.pal <- np[[col]]
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "euclidean"
        dcols1 <- "euclidean"
        minx = min(tmean.scale)
        maxx = max(tmean.scale)

        if (F) {
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            #absmax = ceiling(max(abs(c(minx, maxx))))
            #breaks = c(-1*absmax, seq(0, 1, length=98), absmax)
            #legbreaks = c(-1*absmax, 0, absmax)
            breaks = seq(-2, 2, length=100)
            legbreaks = seq(-2, 2, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        #Run cluster method using 
        hc = hclust(dist(t(tmean.scale)), method="complete")
        hcrow = hclust(dist(tmean.scale), method="complete")
        
        if (FALSE) {
            sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
        } else {
            sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        }
        if (clus2) {
            rowclus <- sort_hclust(hcrow)
        } else {
            rowclus = FALSE
        }
        #print('sorted the clusters')

        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

    pathname <- stringr::str_replace_all("Leukocyte_Migration", "_", " ") 
   pathname <- stringr::str_wrap(pathname,30)

        hm.parameters <- list(
            tmean.scale, 
            color=col.pal,
            legend_breaks=legbreaks,
            cellwidth=14, 
            #cellheight=16, 
            scale="none",
            treeheight_col=treeheight,
            treeheight_row=treeheight,
            kmeans_k=NA,
            breaks=breaks,
            # height=80,
            fontsize=14,
            fontsize_row=14,
            fontsize_col=10,
            show_rownames=rn, 
            show_colnames=cn,
            main=pathname,
            clustering_method="complete",
            cluster_rows=rowclus, 
            cluster_cols=clus,
            cutree_rows=1,
            clustering_distance_rows=drows1, 
            clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            annotation_colors = annot_col,
            labels_col = labels_col,
            annotation_names_col = FALSE
        )

        mat = t(tmean.scale)
        # print('calculated mat')
        
        callback = function(hc, mat) {
            # print('inside the callback')
            dend=rev(dendsort(as.dendrogram(hc)))
            # print ('reversed the dendsorted hc')
            dend %>% dendextend::rotate(c(1:length(dend))) -> dend
            as.hclust(dend)
        }
        do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
    }

    samples_to_include = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","LN_IL15_Treated_336")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]

    spark_df <- DEG_lymph

    if (FALSE == FALSE) {
        genes_to_include_df = genes
        genes_to_include_df <- dplyr::withColumnRenamed(genes_to_include_df, "leukocyte_migration_GO_0050900", "gene_name_column_temp")
        spark_df %>% 
            dplyr::join(genes_to_include_df, spark_df$Gene==genes_to_include_df$gene_name_column_temp, "inner") %>%
            dplyr::drop("gene_name_column_temp") -> spark_df
    }
    
    spark_df %>% dplyr::select(append("Gene", samples_to_include)) -> spark_df

    df=dplyr::collect(spark_df)
    df %>% group_by(Gene) %>% summarise_all(funs(sum)) -> df
    df %>% dplyr::filter(!is.na(Gene)) -> df
    df.mat = df[ , (colnames(df) != "Gene" )]
    df.mat <- as.data.frame(df.mat)
    row.names(df.mat) <- df$Gene 
    

    annot <- dplyr::collect(metadata_fixed)
    annot %>% dplyr::filter(sample_name %in% samples_to_include) -> annot
    timepoint <- as.character(annot$timepoint)
    timepoint[which(timepoint == "day_1")] <- "1st injection"
    timepoint[which(timepoint == "day_2")] <- "2nd injection"
    timepoint[which(timepoint == "day_3")] <- "3rd injection"
    annot$Timepoint <- factor(timepoint)

    
    Treatment <- as.character(annot$treatment)
    Treatment[which(Treatment == "control")] <- "Control"
    Treatment[which(Treatment == "IL15")] <- "hetIL-15"
    annot$Treatment <- factor(Treatment)
    
    if(TRUE){
      annot %>% arrange_(.dots=c("Treatment","Timepoint")) -> annot
      df.mat <- df.mat[,match(annot$sample_name,colnames(df.mat))] 
    }

    groups = c("Treatment","Timepoint")

    annot %>% dplyr::select(groups) -> annotation_col    
    annotation_col = as.data.frame(unclass(annotation_col))
    rownames(annotation_col) <- annot$sample_name
    annot_col = list()
    #colors <- c("#7E7D7D","#023CF3","greenyellow","darkviolet","darkorange","darkturquoise","darkblue","chocolate","deeppink","royalblue")
    colors <- c("darkred","greenyellow","darkviolet","black","darkorange","darkorchid","darkturquoise","darkblue","azure","cadetblue","chocolate","deeppink","lavender")

    b=1
    i=1
    while (i <= length(groups)){
        nam <- groups[i]
        grp <- as.factor(annotation_col[,i])
        c <- b+length(levels(grp))-1
        col = colors[b:c]
        names(col) <- levels(grp)
        assign(nam,col)
        annot_col = append(annot_col,mget(nam))
        b = b+c
        i=i+1
    }

    print(paste0("The total number of genes in heatmap: ", nrow(df.mat)))

    labels_col <- colnames(df.mat)

    manually_replace_sample_names = FALSE
    if (manually_replace_sample_names) {
        replacements = c("")
        old <- c()
        new <- c()
        for (i in 1:length(replacements)) {
            old[i] <- strsplit(replacements[i], ": ?")[[1]][1]
            new[i] <- strsplit(replacements[i], ": ?")[[1]][2]
        }
        df.relabel <- as.data.frame(cbind(old, new), stringsAsFactors=FALSE)
        labels_col %>% replace(match(df.relabel$old, labels_col), df.relabel$new) -> labels_col
    }

    png(
      filename="FIGS2D_S4E_L_leukocyte_migration_new.png",
      width=14,
      height=10,
      units="in",
      pointsize=4,
      bg="white",
      res=300,
      type="cairo")
    p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, rn=TRUE, cn=FALSE, col="Bu Yl Rd")
    print(p)

    return(DEG_lymph)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_FIGS2D_S4E_L_leukocyte_migration_new.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ld3_go_collect<-readRDS(paste0(rds_output,"/var_ld3_go_collect.rds"))
Input_is_Seurat_count <- 0
for(item in var_ld3_go_collect){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ld3_go_collect<-as.data.frame(var_ld3_go_collect)}else{var_ld3_go_collect <- var_ld3_go_collect}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_lymph<-readRDS(paste0(rds_output,"/var_DEG_lymph.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_lymph){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_lymph<-as.data.frame(var_DEG_lymph)}else{var_DEG_lymph <- var_DEG_lymph}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_metadata_fixed<-readRDS(paste0(rds_output,"/var_metadata_fixed.rds"))
Input_is_Seurat_count <- 0
for(item in var_metadata_fixed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_metadata_fixed<-as.data.frame(var_metadata_fixed)}else{var_metadata_fixed <- var_metadata_fixed}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Ortholog_Map_for_RNA_Seq<-readRDS(paste0(rds_output,"/var_Ortholog_Map_for_RNA_Seq.rds"))
Input_is_Seurat_count <- 0
for(item in var_Ortholog_Map_for_RNA_Seq){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Ortholog_Map_for_RNA_Seq<-as.data.frame(var_Ortholog_Map_for_RNA_Seq)}else{var_Ortholog_Map_for_RNA_Seq <- var_Ortholog_Map_for_RNA_Seq}
invisible(graphics.off())
var_FIGS2D_S4E_L_leukocyte_migration_new<-FIGS2D_S4E_L_leukocyte_migration_new(var_ld3_go_collect,var_DEG_lymph,var_metadata_fixed,var_Ortholog_Map_for_RNA_Seq)
invisible(graphics.off())
saveRDS(var_FIGS2D_S4E_L_leukocyte_migration_new, paste0(rds_output,"/var_FIGS2D_S4E_L_leukocyte_migration_new.rds"))
