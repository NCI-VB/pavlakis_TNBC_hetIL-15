# [CCBR] List to Pathway (l2p) Enrichment (e193c7d0-2c92-475b-b565-14e867487c58): v112
ld3_go_pathway_up <- function(DEG_lymph,Ortholog_Map_for_RNA_Seq) {
    # image: png
    suppressMessages(library(l2p))
    suppressMessages(library(magrittr))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    #suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(RColorBrewer))
    DEG_lymph %>% dplyr::select("Gene", "IL15_day_3-control_day_3_tstat") %>% dplyr::collect() -> genesmat
    genesmat %>% dplyr::arrange(desc(`IL15_day_3-control_day_3_tstat`)) -> genesmat
    
    sort_descending = TRUE
    if (sort_descending) {
        genes_to_include = head(genesmat["Gene"], 150)
    } else {
        genes_to_include = tail(genesmat["Gene"], 150)
    }
    genes_to_include <- as.vector(unique(unlist(genes_to_include)))
    genes_universe = as.vector(unique(unlist(genesmat["Gene"])))
    gene_set_sources_to_include = c("GO")
    categories_string <- paste(gene_set_sources_to_include, collapse=",")

    organism = "Mouse"
    if (organism != "Human") {
        organism_of_interest = "Mouse"
        orthology_table = dplyr::filter(Ortholog_Map_for_RNA_Seq,     
        Ortholog_Map_for_RNA_Seq$Organism==organism_of_interest)
        
        if ("from human"=="from human") {
            orthology_reference_column = "Human_gene_name"
            orthology_conversion_column = "Nonhuman_gene_name"
        } else {
            orthology_reference_column = "Nonhuman_gene_name"
            orthology_conversion_column = "Human_gene_name"
    }
        orthology_table %>% dplyr::rename("orthology_reference" = orthology_reference_column) %>%(orthology_reference_column,           "orthology_reference") %>%
dplyr::rename("orthology_conversion" = orthology_conversion_column ) %>% dplyr::select("orthology_reference", "orthology_conversion") -> orthology_table
        orthology_table <- dplyr::collect(orthology_table)
        orthology_table %>% dplyr::filter(orthology_conversion %in% genes_to_include) %>% dplyr::select(orthology_reference) -> genes_to_include
        genes_to_include <- as.vector(genes_to_include$orthology_reference)
        genes_to_include <- unique(genes_to_include)
        orthology_table %>% dplyr::filter(orthology_conversion %in% genes_universe) %>% dplyr::select(orthology_reference) -> genes_universe
        genes_universe <- as.vector(genes_universe$orthology_reference)
        genes_universe <- unique(genes_universe)

    }
    else{
    genes_to_include = genes_to_include
    genes_universe = genes_universe
    }

    use_built_in_gene_universe = FALSE
    if (use_built_in_gene_universe) {
        x <- l2pwcats(genes_to_include, categories_string)
        print("Using built-in gene universe.")
    } else {
        x <- l2puwcats(genes_to_include, genes_universe, categories_string)
        print("Using all genes in differential expression analysis as gene universe.")
    }
    x %>%
        dplyr::arrange(pval) %>% dplyr::mutate(hitsPerc=(pwhitcount/(pwnohitcount+pwhitcount))*100) %>% dplyr::mutate(pathtotal=pwnohitcount+pwhitcount) %>% dplyr::filter(ratio >= 0) %>%
        dplyr::select("pathwayname", "category", "pathwayaccessionidentifier", "pval", "fdr", "pwhitcount", "genesinpathway", "pwnohitcount","pathtotal","hitsPerc", "inputcount", "pwuniverseminuslist","ratio")  %>% dplyr::rename(diff_ratio = ratio) %>% dplyr::filter(pval < 0.05)%>% dplyr::filter(pwhitcount >= 5) -> x
  #allgenes <- lapply(x$pathwayaccessionidentifier,function(x) l2pgetgenes4acc(x)) 
  #x %>% mutate("allgenes" = paste(allgenes,sep="",collapse=" ")) -> x
    print(paste0("Total number of pathways: ", nrow(x)))
    goResults <- x
    goResults %>% top_n(10, wt=-log(pval)) %>%
        dplyr::arrange(-log(pval)) -> goResults
    minp = min(goResults$pval) - 0.1*min(goResults$pval)
    maxp = max(goResults$pval) + 0.1*max(goResults$pval)
    #print(goResults$pval)
    sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
    goResults %>% dplyr::mutate(pathwayname2 = stringr::str_replace_all(pathwayname, "_", " ")) -> goResults
    goResults %>% dplyr::mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults
    
    if (FALSE){
    goResults %>% dplyr::mutate(percorder = order(goResults$pval)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$pval))
    xmax = max(goResults$pval) 
    gplot <- goResults %>% 
               ggplot(aes(x=pval,
               y=pathwayname2, 
               colour=hitsPerc, 
               size=pwhitcount)) +
        geom_point() +
        theme(text = element_text(size=8)) +
        xlim(xmin,xmax) +
        expand_limits(colour = seq(minp, maxp, by = 10),
                size = seq(0, sizemax,by=10)) +
        labs(x="p value", y="GO term", colour="Hits (%)", size="Count") 
    print(gplot)
    }
    else{
    goResults %>% dplyr::mutate(percorder = order(goResults$hitsPerc)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$hitsPerc)-5)
    xmax = ceiling(max(goResults$hitsPerc)+5) 
    gplot <- goResults %>% 
               ggplot(aes(x=hitsPerc,
               y=pathwayname2, 
               colour=pval, 
               size=pwhitcount)) +
        geom_point() +
        theme(text = element_text(size=8)) +
        xlim(xmin,xmax) +
        expand_limits(colour = seq(minp, maxp, by = 10),
                size = seq(0, sizemax,by=10)) +
        labs(x="Hits (%)", y="GO term", colour="p value", size="Count") 
    print(gplot)
    }
    return(x)
}

#################################################
## Global imports and functions included below ##
#################################################
install_bioconductor_package <- function(pkg) {
  }
install_bioconductor_package("GenomeInfoDbData_1.2.1_r351")
suppressMessages(library(GenomeInfoDbData))

# Functions defined here will be available to call in
# the code for any table.

print("template_function_ld3_go_pathway_up.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_lymph<-readRDS(paste0(rds_output,"/var_DEG_lymph.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_lymph){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_lymph<-as.data.frame(var_DEG_lymph)}else{var_DEG_lymph <- var_DEG_lymph}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Ortholog_Map_for_RNA_Seq<-readRDS(paste0(rds_output,"/var_Ortholog_Map_for_RNA_Seq.rds"))
Input_is_Seurat_count <- 0
for(item in var_Ortholog_Map_for_RNA_Seq){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Ortholog_Map_for_RNA_Seq<-as.data.frame(var_Ortholog_Map_for_RNA_Seq)}else{var_Ortholog_Map_for_RNA_Seq <- var_Ortholog_Map_for_RNA_Seq}
invisible(graphics.off())
var_ld3_go_pathway_up<-ld3_go_pathway_up(var_DEG_lymph,var_Ortholog_Map_for_RNA_Seq)
invisible(graphics.off())
saveRDS(var_ld3_go_pathway_up, paste0(rds_output,"/var_ld3_go_pathway_up.rds"))
