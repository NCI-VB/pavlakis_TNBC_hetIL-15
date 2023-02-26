# [CCBR] Filter Low Count Genes For Bulk RNAseq Data (59e7ee62-1715-4276-b370-e8395168f9d8): v96
lymph_node_filtered <- function(remove_ercc, metadata_fixed) {
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(reshape))
    suppressMessages(library(gridExtra))

    metadata <- dplyr::collect(metadata_fixed)
    spark.df = remove_ercc
    samples_to_include = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","LN_IL15_Treated_336")
    samples_to_include <- samples_to_include[samples_to_include %in% metadata$sample_name]
    samples_to_include <- samples_to_include[samples_to_include != "Name"]

    spark.df %>% dplyr::select(append("Name", samples_to_include)) -> spark.df

    df <- dplyr::dropna(spark.df)
    df <- dplyr::collect(df)
    print(paste("filtered for samples =", ncol(df)-1))

    df %>% dplyr::group_by(Name) %>% summarize_all(sum)  %>% data.frame -> df
    print(paste0("Number of genes before filtering: ", nrow(df)))
    
    if (TRUE) {
        rownames(df) <- df$Name
        df$Name <- NULL
        if(TRUE){
        counts <- edgeR::cpm(as.matrix(df)) > 5 }# boolean matrix}
        else{
        counts <- as.matrix(df) > 5
        }
        tcounts <- as.data.frame(t(counts))
        colnum <- dim(counts)[1] # number of genes
        
        rownames(metadata) <- metadata$sample_name
        tcounts <- merge(metadata["treatment"], tcounts, by="row.names")
        tcounts$Row.names <- NULL
        melted <- melt(tcounts, id.vars="treatment")
        tcounts.tot <- dplyr::summarise(dplyr::group_by(melted, treatment, variable), sum=sum(value))
        tcounts.tot %>% tidyr::spread(variable, sum) -> tcounts.group
        colSums(tcounts.group[(1:colnum+1)]>=5) >= 1 -> tcounts.keep
        df <- df[tcounts.keep, ]
        df %>% rownames_to_column("Name") -> df
    } else {
        if (TRUE){
df$isexpr1 <- rowSums(edgeR::cpm(as.matrix(df[, -1])) > 5) >= 1} else {
df$isexpr1 <- rowSums(as.matrix(df[, -1]) > 5) >= 1}

    df <- as.data.frame(df[df$isexpr1, ])
    df$isexpr1 <- NULL
    }

    colnames(df)[colnames(df)=="Name"] <- "Gene"

    print(paste0("Number of genes after filtering: ", nrow(df)))
    
    df -> edf.orig

    metadata[is.na(metadata)]<- "NA"
    metadata = metadata[, colSums(is.na(metadata)) == 0]
    rownames(metadata) <- metadata$sample_name
    metadata <- metadata[match(colnames(edf.orig[,-1]), metadata$sample_name), ]
    metadata = metadata[!is.na(metadata["sample_name"]),]
    print(paste0("Total number of genes included: ", nrow(edf.orig)))

    edf <- edf.orig[,match(metadata$sample_name,colnames(edf.orig))]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- metadata$treatment
    pca.df$sample <- metadata$combined_treatment
    
    # Manual changes to sample names
    replacements = c("")

    plotcolors <- c("darkred","yellowgreen","purple3","black","darkolivegreen","lightgreen","lightgoldenrodyellow","cadetblue","azure","chocolate","coral","deeppink")
    if (length(unique(metadata$Group)) > length(plotcolors)) {
        plotcolors <- c(plotcolors, rep("black", length(unique(metadata$Group)) - length(plotcolors)))
    }

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0("PC1 ", perc.var[1],"%")
    pc.y.lab <- paste0("PC2 ", perc.var[2],"%")
    
    labelpos <- pca.df
    labelpos$mean_y <- pca.df$PC2+2
    labelpos$mean_x <- pca.df$PC1+2
    
    if (TRUE){
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=1) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
          label=sample, color=celltype, vjust="inward", hjust="inward"), size=3, show.legend=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
    }

    else{
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }
    par(mfrow = c(2,1))
   
   df.filt <- df[,match(metadata$sample_name,colnames(df))]
    print(paste0("Total number of samples included: ", ncol(df.filt)))
    rownames(df.filt)=df[,1]
    df.filt %>% rownames_to_column("Name") -> df.filt
    df.m <- melt(df.filt,id.vars=c("Name"))
    if(TRUE){
    df.m %>% dplyr::mutate(value = log2(value+0.5)) -> df.m
    }
    n <- 40
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    if (length(unique(metadata$sample_name)) > length(cols)) {
        cols <- c(cols, rep("black", length(unique(metadata$sample_name)) - length(cols)))
    }

    if(FALSE){
        xmin = -1
        xmax = 1
    }
    else{
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

    if(FALSE){
    df.m %>% mutate(colgroup = metadata[variable,"treatment"]) -> df.m
    g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = colgroup, linetype = colgroup))+
        xlab("filtered counts") + ylab("density") +
         theme_bw() +
        theme(legend.position='top', legend.text = element_text(size = 10),legend.title=element_blank()) + #scale_x_log10() + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) + 
        scale_colour_manual(values=cols) 
    }
    else{
        df.m$variable = metadata[df.m$variable,"sample_name"]
     n=length(unique(df.m$variable))
     if(n>160){
     m=ceiling(n/4)
     n=m*4
     color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
     cols=sample(color, n)}
     else{
         m=40
     }
        g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = variable, linetype = variable) ) + 
        xlab("filterd counts") + ylab("density") +
        theme_bw() +
        theme(legend.position='top', legend.title=element_blank(), legend.text = element_text(size = 10)) + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),m)) + 
        scale_colour_manual(values=cols) +
       guides(linetype = guide_legend(ncol = 6))
    }

    imageWidth = 3000
    imageHeight = 1500*1
    dpi = 300

    png(
      filename="lymph_node_filtered.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    require(gridExtra)
    plots = grid.arrange(g,g2, nrow=1)
    print(plots)

    return(df)
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

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_lymph_node_filtered.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_remove_ercc<-readRDS(paste0(rds_output,"/var_remove_ercc.rds"))
Input_is_Seurat_count <- 0
for(item in var_remove_ercc){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_remove_ercc<-as.data.frame(var_remove_ercc)}else{var_remove_ercc <- var_remove_ercc}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_metadata_fixed<-readRDS(paste0(rds_output,"/var_metadata_fixed.rds"))
Input_is_Seurat_count <- 0
for(item in var_metadata_fixed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_metadata_fixed<-as.data.frame(var_metadata_fixed)}else{var_metadata_fixed <- var_metadata_fixed}
invisible(graphics.off())
var_lymph_node_filtered<-lymph_node_filtered(var_remove_ercc,var_metadata_fixed)
invisible(graphics.off())
saveRDS(var_lymph_node_filtered, paste0(rds_output,"/var_lymph_node_filtered.rds"))
