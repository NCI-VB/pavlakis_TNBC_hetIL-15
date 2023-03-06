cytokine_violins_blank <- function(Labeled_so, Get_Metadata, unnamed, Gene_list) {

    library(Seurat)
    library(reshape2)
    library(cowplot)
    library(rlang)
    
    so <- Labeled_so
    so <- suppressMessages(expr = StashIdent(object = so, save.name = 'ident'))

    assay <- "SCT" #Options RNA/SCT/Protein
    slot  <- "data" #Options: counts/data/scale.data

    #group to plot. This is a metadata column.
    ident_of_interest  <- "ident"

    #groups of interest. From. Metadata column count table
    groups_of_interest <- c("CD103intCD11b+ DC","cDC1","cDC2")
    genes_of_interest <- c("Ermap","Mki67","Klf1","Klf2")

    gene.df <- Gene_list
    genes_of_interest <- gene.df[["A"]]

    scale_data = FALSE
    filter_outliers = FALSE
    reorder_ident = TRUE
    rename_ident = ""
    ylimit = 0
    plot_style = "grid"
    jitter_points <- FALSE
    log_scale_data <- FALSE

    ######################
    ### Error checking ###
    ######################

    gene_filter <- genes_of_interest %in% rownames(x = GetAssayData(object = so, slot = slot, assay = assay))
    genes_of_interest <- genes_of_interest[gene_filter]
    missing_genes <- genes_of_interest[!gene_filter]
    if(length(missing_genes) > 0){
        cat("The following genes are missing from the dataset:\n")
        print(missing_genes)
    }

    if(length(genes_of_interest) == 0){
        stop("No query genes were found in the dataset.")
    }

    if(!ident_of_interest %in% colnames(so@meta.data)){
        colnames(so@meta.data) <- gsub("\\.","_",colnames(so@meta.data))
        if(!ident_of_interest %in% colnames(so@meta.data)){
            stop("Unable to find ident of interest in metadata.")
        }
    }

    group_filter <- groups_of_interest %in% so@meta.data[[ident_of_interest]]
    groups_of_interest <- groups_of_interest[group_filter]
    missing_groups <- groups_of_interest[!group_filter]
    if(length(missing_groups) > 0){
        cat("The following groups are missing from the selected ident:\n")
        print(missing_groups)
    }

    if(length(groups_of_interest) == 0){
        stop("No groups were found in the selected ident.")
    }

    if(rename_ident %in% c("Gene","Expression","scaled")){
        stop("New ident name cannot be one of Gene, Expression, or scaled.")
    }

    ##########################
    ### End error checking ###
    ##########################

    # deal with limits
    if(ylimit == 0){
        ylimit <- NULL
    }

    #set idents
    if("active.ident" %in% slotNames(so)){
        new_idents = as.factor(so@meta.data[[ident_of_interest]])
        names(new_idents) = names(so@active.ident)
        so@active.ident <- as.factor(vector())
        so@active.ident <- new_idents
        so.sub = subset(so, ident = groups_of_interest)
    } else {
        new_idents = as.factor(so@meta.data[[ident_of_interest]])
        names(sample_name)=so@meta.data[["Barcode"]]
        so@active.ident <- as.factor(vector())
        so@active.ident <- new_idents
        so.sub = subset(so, ident = groups_of_interest)
    }

    DefaultAssay(object = so) <- assay
    data <- FetchData(object = so.sub, vars = genes_of_interest, slot = slot)

    data[[ident_of_interest]] <- so.sub@meta.data[row.names(data),ident_of_interest]

    df.melt <- melt(data)
    
    if(!is.empty(rename_ident)){
        ident_of_interest <- rename_ident
    }
    colnames(df.melt) <- c(ident_of_interest,"Gene","Expression")

    #check to see if ident of interest looks numeric
    if(suppressWarnings(all(!is.na(as.numeric(as.character(df.melt[[ident_of_interest]])))))){
        ident.values <- strtoi(df.melt[[ident_of_interest]])
        ident.levels <- unique(ident.values)[order(unique(ident.values))]
        df.melt[[ident_of_interest]] <- factor(ident.values, levels = ident.levels)
    }else if(reorder_ident){
        # if non-numeric, place in order of groups of interests
        df.melt[[ident_of_interest]] <- factor(df.melt[[ident_of_interest]], levels = groups_of_interest)        
    }

    # Filter outliers
    if(filter_outliers){
        for(gene in genes_of_interest){
            for(group in groups_of_interest){
                current.ind <- which(df.melt[["Gene"]] == gene & df.melt[[ident_of_interest]] == group)
                df.melt[current.ind,"Expression"] <- remove_outliers_func(df.melt[current.ind,"Expression", drop = TRUE])
            }
        }
    }

    # Scale data to y limit
    if(scale_data){
        expression_data = "scaled"
        axis.title.y = "Expression (scaled)"
        ylimit <- ylimit %||% 1
        df.melt <- df.melt %>% group_by(Gene) %>% mutate(scaled = scales::rescale(Expression, to=c(0,ylimit)))
    }else{
        expression_data <- axis.title.y <- "Expression"
    }

    imageWidth = 2000*2
    imageHeight = 2000*4
    png(
      filename="cytokine_violins_blank.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=300,
      type="cairo")

    # #Reorder
    # #Change color
    #   Red: #eb3423 (cDC1)
    #   Blue: #003df5 (cDC2)
    #   Green #4eac5b (New)
    #   Gold: lightgoldenrod

    # #Axis thicker
    # #Remove text
    #groups_of_interest <- c("Monocytes","CD103intCD11b+ DC","cDC1","cDC2")
    df.melt$Gene <- factor(as.character(df.melt$Gene), levels = c("Ccl6","Ccl9","Cxcl2","Ccl17","Ccl2","Ccl4","Ccl22","Ccl24"))
    
    g <- ggplot(df.melt, aes_string(x=ident_of_interest, y=expression_data, fill=ident_of_interest)) +
        geom_violin(aes_string(fill = ident_of_interest), scale="width", trim = FALSE, show.legend = FALSE) +
        scale_fill_manual(values=c("#4eac5b", "#eb3423","#003df5"))+
        #geom_jitter(height = 0, width = 0.05, size=0.1) +
        theme_classic() + 
        #scale_fill_manual(values=cols) + 
        labs(y=axis.title.y) +
        theme(  strip.text.x = element_text( color="white", face="plain", size=32),
                axis.line = element_line(colour = 'black', size = 2),
                axis.ticks = element_line(colour = 'black', size = 2),
                axis.ticks.length = unit(12, "points"),
                axis.title = element_blank(),
                axis.text = element_blank())

    if(!is.null(ylimit)){
        g <- g + ylim(0,ylimit)
    }

    if(jitter_points){
        g <- g + geom_jitter(height = 0)
    }

    if(log_scale_data){
        g <- g + scale_y_log10()
    }

    # Plot after jitter if wanted
    g <- g + geom_boxplot(width=0.1, fill="white")
    
    # Plot styles
    ncol = 2 #ceiling(length(unique(df.melt$Gene))^0.5)
    nrow = 4 #ceiling(length(unique(df.melt$Gene)) / ncol)
    if(plot_style == "rows"){
        g <- g + facet_grid(rows = vars(Gene))
    }else{
        g <- g + facet_wrap(~Gene, nrow = nrow, ncol = ncol)
        if(plot_style == "labeled"){
            plots <- splitFacet(g)
            plots <- lapply(seq_along(plots), function(i) plots[[i]] + ggtitle(genes_of_interest[i]) + theme(plot.title = element_text(hjust = 0.5)) )
            g <- plot_grid(plotlist = plots, nrow = nrow, ncol = ncol, labels = LETTERS[seq( from = 1, to = length(plots) )])
        }
    }
    print(g)
# auto removed:     return(NULL)
}

splitFacet <- function(x){
    # https://stackoverflow.com/questions/30510898/split-facet-plot-into-list-of-plots/52225543
    facet_vars <- names(x$facet$params$facets)         # 1
    x$facet    <- ggplot2::ggplot()$facet              # 2
    datasets   <- split(x$data, x$data[facet_vars])    # 3
    new_plots  <- lapply(datasets,function(new_data) { # 4
        x$data <- new_data
        x})
}

# from rapportools
is.empty <- function(x, trim = TRUE, ...) {
    if (length(x) <= 1) {
        if (is.null(x))
            return (TRUE)
        if (length(x) == 0)
            return (TRUE)
        if (is.na(x) || is.nan(x))
            return (TRUE)
        if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
            return (TRUE)
        if (is.logical(x) && !isTRUE(x))
            return (TRUE)
        if (is.numeric(x) && x == 0)
            return (TRUE)
        return (FALSE)
    } else
        sapply(x, is.empty, trim = trim, ...)
}

# from rapportools
trim.space <- function(x, what = c('both', 'leading', 'trailing', 'none'), space.regex = '[:space:]', ...){
    if (missing(x))
        stop('nothing to trim spaces to =(')
    re <- switch(match.arg(what),
                 both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
                 leading  = sprintf('^[%s]+', space.regex),
                 trailing = sprintf('[%s]+$', space.regex),
                 none     = {
                     return (x)
                 })
    vgsub(re, '', x, ...)
}

# from rapportools
vgsub <- function(pattern, replacement, x, ...){
    for(i in 1:length(pattern))
        x <- gsub(pattern[i], replacement[i], x, ...)
    x
}

remove_outliers_func <- function(x, na.rm = TRUE){
    qnt <- quantile(x, probs=c(0.25,0.75), na.rm = na.rm)
    H <- 3*IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
}

print("template_function_cytokine_violins_blank.R #########################################################################")
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
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_unnamed<-readRDS(paste0(rds_output,"/var_unnamed.rds"))
Input_is_Seurat_count <- 0
for(item in var_unnamed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_unnamed<-as.data.frame(var_unnamed)}else{var_unnamed <- var_unnamed}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Gene_list<-readRDS(paste0(rds_output,"/var_Gene_list.rds"))
Input_is_Seurat_count <- 0
for(item in var_Gene_list){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Gene_list<-as.data.frame(var_Gene_list)}else{var_Gene_list <- var_Gene_list}
invisible(graphics.off())
var_cytokine_violins_blank<-cytokine_violins_blank(var_Labeled_so,var_Get_Metadata,var_unnamed,var_Gene_list)
invisible(graphics.off())
saveRDS(var_cytokine_violins_blank, paste0(rds_output,"/var_cytokine_violins_blank.rds"))
