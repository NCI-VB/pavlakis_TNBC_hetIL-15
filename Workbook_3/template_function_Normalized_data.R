# [CCBR] Normalization of Bulk RNA-seq Data (afc2524c-9bae-4873-98c1-e06ef8f4632b): v121
Normalized_data <- function(unnamed, metadata_fixed) {
    #image: png
    imageType = "png"
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(reshape))
    suppressMessages(library(gridExtra))

    df <- dplyr::collect(unnamed)
    
    samples_for_deg_analysis = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","Tumor_CTRL_321","Tumor_CTRL_327","Tumor_CTRL_331","Tumor_IL15_320","Tumor_IL15_338","Tumor_IL15_339","Tumor_CTRL_203","Tumor_CTRL_205","Tumor_CTRL_206","Tumor_IL15_209","Tumor_IL15_215","Tumor_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","Tumor_Control_301","Tumor_Control_304","Tumor_Control_324","Tumor_IL15_Treated_684","Tumor_IL15_Treated_336","Tumor_IL15_Treated_340")
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis !=            "Gene"]
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis !=            "GeneName"]
    df.m <- df[,samples_for_deg_analysis]
    gene_names <- NULL
    gene_names$GeneID <- df[,1]
    targetfile <- dplyr::collect(metadata_fixed)
    targetfile <- targetfile[match(colnames(df.m),targetfile$sample_name),]
    targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
    df.m <- df.m[,match(targetfile$sample_name,colnames(df.m))]
    if(FALSE){
    x <- DGEList(counts=2^df.m, genes=gene_names)
    } else {
    x <- DGEList(counts=df.m, genes=gene_names)     
    }
    dm.formula <- as.formula(paste("~0 +", paste(c("tissue","timepoint"), sep="+", collapse="+")))
    design=model.matrix(dm.formula, targetfile)
    colnames(design) <- str_replace_all(colnames(design), c("tissue","timepoint"), "")
    v <- voom(x,design=design,normalize="quantile")
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    print(paste0("Total number of genes included: ", nrow(df.voom)))

    samples_to_include = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","Tumor_CTRL_321","Tumor_CTRL_327","Tumor_CTRL_331","Tumor_IL15_320","Tumor_IL15_338","Tumor_IL15_339","Tumor_CTRL_203","Tumor_CTRL_205","Tumor_CTRL_206","Tumor_IL15_209","Tumor_IL15_215","Tumor_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","Tumor_Control_301","Tumor_Control_304","Tumor_Control_324","Tumor_IL15_Treated_684","Tumor_IL15_Treated_336","Tumor_IL15_Treated_340")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    #df %>% dplyr::select(append("Gene", samples_to_include)) -> spark.df
    df.voom -> df
    df -> edf.orig

    Sampinfo <- dplyr::collect(metadata_fixed)
    # cell <- Sampinfo$sample_name
    Sampinfo = Sampinfo[, colSums(is.na(Sampinfo)) == 0]
    rownames(Sampinfo) <- Sampinfo$sample_name
    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo$sample_name), ]
    Sampinfo = na.omit(Sampinfo)
    print(paste0("Total number of genes included: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo$sample_name,colnames(edf.orig))]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- Sampinfo$timepoint
    pca.df$sample <- Sampinfo$tissue
    
    # Manual changes to sample names
    replacements = c("")

    plotcolors <- c("darkred","yellowgreen","purple3","black","darkolivegreen","lightgreen","lightgoldenrodyellow","cadetblue","azure","chocolate","coral","deeppink")
    if (length(unique(Sampinfo$timepoint)) > length(plotcolors)) {
        plotcolors <- c(plotcolors, rep("black", length(unique(Sampinfo$timepoint)) - length(plotcolors)))
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
    } else {
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
   
    #samples_to_include = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","Tumor_CTRL_321","Tumor_CTRL_327","Tumor_CTRL_331","Tumor_IL15_320","Tumor_IL15_338","Tumor_IL15_339","Tumor_CTRL_203","Tumor_CTRL_205","Tumor_CTRL_206","Tumor_IL15_209","Tumor_IL15_215","Tumor_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","Tumor_Control_301","Tumor_Control_304","Tumor_Control_324","Tumor_IL15_Treated_684","Tumor_IL15_Treated_336","Tumor_IL15_Treated_340")
   # samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    df.filt <- df %>% dplyr::select(samples_to_include)
   # Sampinfo <- dplyr::collect(metadata_fixed)
        # cell <- Sampinfo$SampleName
        rownames(Sampinfo) <- Sampinfo$sample_name
        Sampinfo = Sampinfo[complete.cases(Sampinfo[, "sample_name"]),]
        print(paste0("Total number of samples included: ", nrow(Sampinfo)))
    df.filt <- df.filt[,match(Sampinfo$sample_name,colnames(df.filt))]
    rownames(df.filt)=df[,1]
    df.filt %>% rownames_to_column("Gene") -> df.filt
    df.m <- melt(df.filt,id.vars=c("Gene"))
    if(FALSE){
    df.m %>% mutate(value = log2(value+0.5)) -> df.m
    }
    n <- 40
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    if (length(unique(Sampinfo$sample_name)) > length(cols)) {
        cols <- c(cols, rep("black", length(unique(Sampinfo$sample_name)) - length(cols)))
    }

    if(FALSE){
        xmin = -1
        xmax = 1
    } else {
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

  if(FALSE){
    df.m %>% mutate(colgroup = Sampinfo[variable,"tissue"]) -> df.m
    g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = colgroup, linetype = colgroup))+
        xlab("normalized counts") + ylab("density") +
         theme_bw() +
        theme(legend.position='top', legend.text = element_text(size = 10),legend.title=element_blank()) + #scale_x_log10() + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) + 
        scale_colour_manual(values=cols) 
    } else {
        df.m$variable = Sampinfo[df.m$variable,"treatment"]
     n=length(unique(df.m$variable))
     if(n>160){
     m=ceiling(n/4)
     n=m*4
     color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
     cols=sample(color, n)
    } else {
         m=40
    }
        g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = variable, linetype = variable) ) + 
        xlab("normalized counts") + ylab("density") +
        theme_bw() +
        theme(legend.position='top', legend.title=element_blank(), legend.text = element_text(size = 10)) + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),m)) + 
        scale_colour_manual(values=cols) + guides(linetype = guide_legend(ncol = 6))
    }

    imageWidth = 3000
    imageHeight = 1500*1
    dpi = 300

    ## Choice of image output format: PNG or SVG.
    if (imageType == 'png') {
    png(
      filename="Normalized_data.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    } else {
        library(svglite)
        svglite::svglite(
        file="Normalized_data.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    require(gridExtra)
    plots = grid.arrange(g,g2,nrow=1)
    print(plots)

    return(df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Normalized_data.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_unnamed<-readRDS(paste0(rds_output,"/var_unnamed.rds"))
Input_is_Seurat_count <- 0
for(item in var_unnamed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_unnamed<-as.data.frame(var_unnamed)}else{var_unnamed <- var_unnamed}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_metadata_fixed<-readRDS(paste0(rds_output,"/var_metadata_fixed.rds"))
Input_is_Seurat_count <- 0
for(item in var_metadata_fixed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_metadata_fixed<-as.data.frame(var_metadata_fixed)}else{var_metadata_fixed <- var_metadata_fixed}
invisible(graphics.off())
var_Normalized_data<-Normalized_data(var_unnamed,var_metadata_fixed)
invisible(graphics.off())
saveRDS(var_Normalized_data, paste0(rds_output,"/var_Normalized_data.rds"))
