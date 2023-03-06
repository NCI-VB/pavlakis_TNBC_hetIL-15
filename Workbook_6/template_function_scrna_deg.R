# DEG (Gene Expression Markers) [scRNA-seq] [CCBR] (54b6dd44-e233-4fb8-9bae-6ab0cb46e399): v58
scrna_deg <- function(reclusterd_so,Sample_Names_Reclustered,deg_metadata) {
    
    # Faking package "GenomeInfoDbData". See bottom of code for more details.
    if(!"GenomeInfoDbData" %in% loadedNamespaces()){
        fakepackage("GenomeInfoDbData")
    }
    
    suppressMessages(library(Seurat))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(scales))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggrepel))
    suppressMessages(library(gdata))
    suppressMessages(library(reshape2))
    suppressMessages(library(tools))
    suppressMessages(library(grid))
    suppressMessages(library(gridBase))
    suppressMessages(library(gridExtra))
    suppressMessages(library(parallel))
    suppressMessages(library(MAST))

    SO = reclusterd_so$value

     #collect parameters here:
    contrasts <- c("0-all","1-2")
    useSpark <- FALSE
    
    samples = eval(parse(text=gsub('\\[\\]','c()','c("Control","Treated")')))
    
    if (length(samples) == 0) {
        samples = unique(SO@meta.data$sample_name)
    }

    colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))
    if("active.ident" %in% slotNames(SO)){
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
    } else {
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
    }
    
    print("selected samples:")
    print(SO.sub)

    meta.df <- dplyr::collect(deg_metadata)
    colnames(SO.sub@meta.data) = gsub("\\.","_",colnames(SO.sub@meta.data))

    #define contrasts
    newcont <- list()
    for (i in 1:length(contrasts)){
       newcont[[i]] <- c(paste(unlist(strsplit(contrasts[i],"-"))))
    }
    contrasts <- newcont

    #ERROR CATCHING
    #collect valid names of valid columns
    validColumns <- character()
    for (i in colnames(meta.df)) {
        if (!any(is.na(meta.df[[i]]))) {
            validColumns <-c(validColumns,i)
        }
    }

    param2test <- "SCT_snn_res_0_6"

    if (param2test =="") {
        mcols = colnames(SO.sub@meta.data)
        param2test <-mcols[grepl("RNA_snn",mcols)][[1]]
        print(paste("No parameter selected, defaulting to",param2test))
    }

    contrastTarget <- SO.sub@meta.data[[param2test]]
    contrastType <- param2test
    contrastCounts = as.data.frame(table(contrastTarget))
    validContrasts = subset(contrastCounts, Freq>2)[[1]]

    #catch malformed contrasts
    for (i in contrasts) {
        if (!(i[[1]] %in% contrastTarget)) {
            print(paste(i[[1]],"is not a valid contrast for contrast type:", contrastType))
            print("Please see below for an example of valid contrasts for your selected contrast type.")
            print(validContrasts)
            stop("You have entered an invalid group to contrast against.")
        } else if (!(i[[2]] %in% contrastTarget) & (i[[2]] != "all")) {
            print(paste(i[[2]],"is not a valid contrast for contrast type:", contrastType))
            print("Please see below for an example of valid contrasts for your selected contrast type.")
            print(validContrasts)
            stop("You have entered an invalid group to contrast against.")
        } else if (length(i)>2) {
            print("Contrasts are as follows..")
            print(i)
            stop("The console says there are too many inputs in your contrasts. A contrast should only contain Group1-Group2, but the console thinks you have inputed Group1-Group2-Group3")
        } else if (!(i[[2]] %in% validContrasts) & (i[[2]] != "all")) {
            print(paste(i[[2]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
            stop("You have entered an invalid group to contrast against.")
        } else if (!(i[[1]] %in% validContrasts)) {
            print(paste(i[[1]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
            stop("You have entered an invalid group to contrast against.")
        }
    }

    #print out contrast cell contrastCounts
    for (i in seq_along(contrasts)) {
        firstGroup <- contrasts[[i]][[1]]
        firstGroupCount <- subset(contrastCounts, contrastTarget == firstGroup)$Freq
        if  (contrasts[[i]][[2]]!= "all") {
            secondGroup <-contrasts[[i]][[2]]
            secondGroupCount <-subset(contrastCounts, contrastTarget == secondGroup)$Freq      
        print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. cluster",secondGroup,"with",secondGroupCount,"cells."))
        } else {
            secondGroupCount <-ncol(SO.sub)-firstGroupCount
            print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. all other clusters, totalling",secondGroupCount,"cells."))
        } 
    }

    #define and call function for running DEG
    get_deg_table <- function(n) {
        library(Seurat)

        firstCluster <-n[[1]]
        secondCluster <- n[[2]]

        if (n[[2]]=='all') {
            secondCluster <- NULL
        }
        
        Idents(SO.sub) <- param2test
        markers = FindMarkers(SO.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = "MAST", logfc.threshold = 0.25, verbose=FALSE, assay = "RNA")
        colnames(markers) <- chartr(old=" ",new="_",paste(colnames(markers), n[[1]],"vs",n[[2]],sep = "_"))
        return(markers)
    }

    if (useSpark) {
        deg_tables <- spark.lapply(contrasts, get_deg_table) 
    } else {
        deg_tables <- lapply(contrasts, get_deg_table) 
    }

    for(i in seq_along(deg_tables)){
        degtab <- deg_tables[[i]]
        degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos 
        degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
        print(paste0("The number of upregulated genes at p<0.05 in contrast number ", i, " is:"))
        print(pos[1])
        print(paste0("The number of downregulated genes at p<0.05 in contrast number ", i, " is:"))
        print(neg[1]) 
    }

    #Merge the deg tables together
    out_df <- NULL
    for (i in deg_tables) {
        if (is.null(out_df)) {
            out_df <- deg_tables[1]
            out_df <- as.data.frame(out_df)
        } else {
            out_df <- merge(out_df, i, by="row.names", all=TRUE)
            rownames(out_df) <- out_df$Row.names #set the rownames
            out_df$Row.names <- NULL #drop the row.names columns which we no longer need
        }
    }
    
    out_df$Gene <- rownames(out_df)
    out_df$Row.names <- NULL
    out_df <- out_df %>% dplyr::select(Gene, everything())
    return(out_df)

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

# Functions defined here will be available to call in
# the code for any table.

# Using "fakepackage" to fake the existance of GenomeInfoDbData in the namespace.
# This is needed because, while none of its functions are used, it is a dependency of one of the dependencies of MAST.
# For reasons unclear, we can't get GenomeInfoDbData to load in normal ways.
fakepackage <- function(name, exported = NULL, unexported = NULL, attach = TRUE) {
  # fetch and eval call to create `makeNamespace`
  eval(body(loadNamespace)[[c(8, 4, 4)]])
  # create an empty namespace
  ns <- makeNamespace(name)
  # makethis namespace the closure env of our input functions
  exported <- lapply(exported, `environment<-`, ns)
  unexported <- lapply(unexported, `environment<-`, ns)
  # place these in the namespace
  list2env(exported, ns)
  list2env(unexported, ns)
  # export relevant functions
  namespaceExport(ns, names(exported))
  if(attach) {
    # copy exported funs to "package:pkg" envir of the search path
    attach(exported, name = paste0("package:", name))
  invisible()
}

print("template_function_scrna_deg.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_reclusterd_so<-readRDS(paste0(rds_output,"/var_reclusterd_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_reclusterd_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_reclusterd_so<-as.data.frame(var_reclusterd_so)}else{var_reclusterd_so <- var_reclusterd_so}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Sample_Names_Reclustered<-readRDS(paste0(rds_output,"/var_Sample_Names_Reclustered.rds"))
Input_is_Seurat_count <- 0
for(item in var_Sample_Names_Reclustered){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Sample_Names_Reclustered<-as.data.frame(var_Sample_Names_Reclustered)}else{var_Sample_Names_Reclustered <- var_Sample_Names_Reclustered}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_deg_metadata<-readRDS(paste0(rds_output,"/var_deg_metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_deg_metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_deg_metadata<-as.data.frame(var_deg_metadata)}else{var_deg_metadata <- var_deg_metadata}
invisible(graphics.off())
var_scrna_deg<-scrna_deg(var_reclusterd_so,var_Sample_Names_Reclustered,var_deg_metadata)
invisible(graphics.off())
saveRDS(var_scrna_deg, paste0(rds_output,"/var_scrna_deg.rds"))
