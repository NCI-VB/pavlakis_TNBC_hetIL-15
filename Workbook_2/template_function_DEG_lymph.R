# [CCBR] Voom Normalization and DEG Analysis (de953b9c-b7f3-49ac-b883-55070d759e5d): v79
DEG_lymph <- function(lymph_node_filtered, metadata_fixed) {
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(stringr))
    df <- dplyr::collect(lymph_node_filtered)

    samples_for_deg_analysis = c("Lymph_CTRL_321","Lymph_CTRL_327","Lymph_CTRL_331","Lymph_IL15_320","Lymph_IL15_338","Lymph_IL15_339","Lymph_CTRL_203","Lymph_CTRL_205","Lymph_CTRL_206","Lymph_IL15_209","Lymph_IL15_215","Lymph_IL15_216","LN_Control_301","LN_Control_304","LN_Control_324","LN_IL15_Treated_684","LN_IL15_Treated_340","LN_IL15_Treated_336")
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis != "Gene"]
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis != "GeneName"]

    df.m <- df[,samples_for_deg_analysis]

    gene_names <- NULL
    gene_names$GeneID <- df[,1]
    
    targetfile <- dplyr::collect(metadata_fixed)
    targetfile <- targetfile[match(colnames(df.m),targetfile$sample_name),]
    targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
    df.m <- df.m[,match(targetfile$sample_name,colnames(df.m))]
    if(FALSE){
    x <- DGEList(counts=2^df.m, genes=gene_names)
    }
    else{
    x <- DGEList(counts=df.m, genes=gene_names)     
    }
    
    dm.formula <- as.formula(paste("~0 +", paste(c("combined_treatment"), sep="+", collapse="+")))
    design=model.matrix(dm.formula, targetfile)
    print(design)
    colnames(design) <- str_replace_all(colnames(design), "combined_treatment", "")
    
    if ("quantile" %in% c("TMM","TMMwzp","RLE","upperquartile")){
    x <- calcNormFactors(x, method = "quantile") 
    rownames(x) <- x$genes$GeneID
    v <- voom_v2(x,design=design,normalize="none")
    }
    else {
    v <- voom_v2(x,design=design,normalize="quantile",plot=TRUE)
    }
    
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    fit <- lmFit(v, design)
    cm <- makeContrasts(contrasts = c("IL15_day_1-control_day_1","IL15_day_2-control_day_2","IL15_day_3-control_day_3"), levels=design)
    print(colnames(design))
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    logFC = fit2$coefficients
    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
    tstat = fit2$t
    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
    FC = 2^fit2$coefficients
    FC = ifelse(FC<1,-1/FC,FC)
    colnames(FC)=paste(colnames(FC),"FC",sep="_")
    pvalall=fit2$p.value
    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
    
    
    if(TRUE){
    tve <- t(v$E)
    mean.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=targetfile[targetfile$sample_name==Sample,"combined_treatment"]) %>% group_by(group) %>% summarise_all(funs(mean)) %>% as.data.frame()
    mean.df[,-c(1,2)] %>% as.matrix() %>% t() -> mean
    colnames(mean) <- mean.df[,1]
    colnames(mean)=paste("Mean",colnames(mean),sep="_")
    colnames(mean) = gsub("\\.", "_", colnames(mean))
    sd.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=targetfile[targetfile$sample_name==Sample,"combined_treatment"]) %>% group_by(group) %>% summarise_all(funs(sd)) %>% as.data.frame()
    sd.df[,-c(1,2)] %>% as.matrix() %>% t() -> sd
    colnames(sd) <- sd.df[,1]
    colnames(sd)=paste("SD",colnames(sd),sep="_")
    colnames(sd) = gsub("\\.", "_", colnames(sd))
   finalres=as.data.frame(cbind(mean,sd,v$E,FC, logFC, tstat, pvalall, pvaladjall)) 
}
   else{
    finalres=as.data.frame(cbind(v$E,FC, logFC, tstat, pvalall, pvaladjall))
    }
    finalres %>% rownames_to_column("Gene") -> finalres
    print(paste0("Total number of genes included: ", nrow(finalres)))
    return(finalres) 
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#voom function was modified to incrase dot size,line-width, x and y Label size, title label size, margin between tick, tick-values, and labels.

voom_v2 <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
                     span = 0.5, plot = FALSE, save.plot = FALSE, ...) 
{
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    par(mgp=c(7,2.5,0),mar=c(9,11,6,4)+.1)
    #The 'mar' argument of 'par' sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. 
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
         pch = 16, cex = 1.25, cex.lab=4.5,cex.axis=3)
    title(main="voom: Mean-variance trend", cex.main=4.5)
    lines(l, col = "red",lwd=3)
    #mtext("Sqrt( standard deviation )", side=2, line=6, cex=4)
    #mtext("log2( count size + 0.5 )", side=1, line=6, cex=4)
    
  }
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets)) 
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}

print("template_function_DEG_lymph.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_lymph_node_filtered<-readRDS(paste0(rds_output,"/var_lymph_node_filtered.rds"))
Input_is_Seurat_count <- 0
for(item in var_lymph_node_filtered){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_lymph_node_filtered<-as.data.frame(var_lymph_node_filtered)}else{var_lymph_node_filtered <- var_lymph_node_filtered}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_metadata_fixed<-readRDS(paste0(rds_output,"/var_metadata_fixed.rds"))
Input_is_Seurat_count <- 0
for(item in var_metadata_fixed){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_metadata_fixed<-as.data.frame(var_metadata_fixed)}else{var_metadata_fixed <- var_metadata_fixed}
invisible(graphics.off())
var_DEG_lymph<-DEG_lymph(var_lymph_node_filtered,var_metadata_fixed)
invisible(graphics.off())
saveRDS(var_DEG_lymph, paste0(rds_output,"/var_DEG_lymph.rds"))
