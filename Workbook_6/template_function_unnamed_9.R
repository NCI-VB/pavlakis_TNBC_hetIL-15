# PCA & Normalization [scRNA-seq] [CCBR] (4459a0f4-c07b-4d42-98e6-22f2e03230ef): v81
unnamed_9 <- function(Postfilter_QC_Plots) {
    #image: png
    imageType = "png"
    suppressMessages(library(Seurat))
    suppressMessages(library(gridExtra))
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))

    SO = Postfilter_QC_Plots$value

    #in case you want to redo this on a merged SO
    if (class(SO) =="Seurat") {
        x =list()
        x[[1]] <- SO
        SO <- x
    }
    vars_to_regress <- c("S.Score","G2M.Score")
    doJackStraw <- TRUE
    npcs = 25
    linearScale = FALSE
    
    # Linearly scale data without regressing anything.
    scale_so <- function(so){
        so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
        so$CC.Difference <- so$S.Score - so$G2M.Score
        so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst")
        all.genes <- rownames(so)
        so <- ScaleData(so,features=all.genes)
        return(so)
    }

    # Make PCA without regressing anything, and using only SCTransform().
    pca_noregress <- function(so) {
        so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE)
        so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
        return(so)
    }

    # Make PCA with SCTransform() and optional ScaleData, and do so with
    # both regression (if user requests) and on all genes.
    pca <- function(so) {
        # If user sets Linear Scaling toggle TRUE, also run ScaleData().
        # Use case: user has legacy project from Seurat 2 and wants to keep
        # methods consistent with pre-SCT Seurat.
        if(linearScale == TRUE) {
            all.genes <- rownames(so)
            if(is.null(vars_to_regress)){     
                so <- so
            }
            else{
                so <- ScaleData(so, features=all.genes, vars.to.regress = vars_to_regress) 
            }
        }
        # Run SCTransform().
        if(is.null(vars_to_regress)){
            so <- so
        }
        else { 
            so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
        }
        # Make PCA using last transform run, which will always be that from
        # SCTransform().
        so <- RunPCA(object = so, npcs = npcs)
        slot(so,"commands") <- list()
        return(so)
    }

    # Do transformation with and without regression using SCTransform()
    # and ScaleData().
    so_scale <- lapply(SO, scale_so) 
    so_orig <- lapply(so_scale, pca_noregress)
    so_list <- lapply(so_scale, pca) 

    vars_to_plot = c("S.Score","G2M.Score")
    imageCols = 2
    if (length(vars_to_plot) > 0) {
        imageCols <- imageCols + length(vars_to_plot)
    }
    
    if (doJackStraw) {
        imageCols <- imageCols + 2
    }

    plotPCA <- function(so,m){
    p1 <- DimPlot(so, reduction = "pca")
    clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=so@meta.data[[m]])

    sumpcsd = sum(so@reductions$pca@stdev)
    pcvar = (so@reductions$pca@stdev/sumpcsd)*100
    pcvar=formatC(pcvar,format = "g",digits=3)
    pcvar1 = pcvar[1] 
    pcvar2 = pcvar[2]
    pcvar 

    run.categ <- function(mat){
    colors=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9")

    g <- ggplot(mat, aes(x=umap1, y=umap2)) +
                theme_bw() +
                theme(legend.title=element_blank()) +
                geom_point(aes(colour=clusid),size=0.5) +
                scale_color_manual(values=colors) +
                xlab(paste0("PC-1 ",pcvar[1],"%")) + ylab(paste0("PC-2 ",pcvar[2],"%")) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(),legend.text=element_text(size=rel(1.5))) +
                guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
                ggtitle(so@project.name)  
                return(g)
    }

    run.cont <- function(mat,midpt,maxpt){
        g <- ggplot(mat, aes(x=umap1, y=umap2)) +
                theme_bw() +
                theme(legend.title=element_blank()) +
                geom_point(aes(colour=clusid),size=0.5) +
                scale_colour_gradient2(low = "blue",mid="lightgrey",high = "red",limits=c(0, maxpt),midpoint = midpt) +
                xlab(paste0("PC-1 ",pcvar[1],"%")) + ylab(paste0("PC-2 ",pcvar[2],"%")) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(),legend.text=element_text(size=rel(1.5))) +
                guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
                ggtitle(paste0(so@project.name,"_",m))  
                return(g)

    }
    if(class(clusmat$clusid) == "factor"){
        g = run.categ(clusmat)
        return(g)
    }
    else{
        clusmat %>% arrange(clusid) -> clusmat
        if(m=="percent.mt"){
          #mid = quantile(clusmat$clusid)[4]
          mid=5
          max = 10
          clusmat$clusid[clusmat$clusid > max] <- max
        }
        else{
          mid = quantile(clusmat$clusid)[3]
          max = quantile(clusmat$clusid,probs=0.95)
          clusmat$clusid[clusmat$clusid > max] <- max
        }
        
        g = run.cont(clusmat,mid,max)
        return(g)
    }

    }

    methods = c("Elbow")
    
    plotElbow <- function(so){

    if("Elbow" %in% methods){
        #Find Elbow:
        sumpcsd = sum(so@reductions$pca@stdev)
        pcvar = (so@reductions$pca@stdev/sumpcsd)*100
        cumu <- cumsum(pcvar)
        co1 <- which(cumu > 80 & pcvar < 5)[1]
        co2 <- sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 0.1), decreasing = T)[1] + 1
        pcs = min(co1,co2)
        lab = paste0("Elbow = ", pcs)
        xpos = pcs + 4
    }

    if("Marchenko-Pastur" %in% methods){
        #Using code from URD (https://rdrr.io/github/farrellja/URD/src/R/pca.R)
        pcaMarchenkoPastur <- function(M, N, pca.sdev, factor=1, do.print=T) {
        pca.eigenvalue <- (pca.sdev)^2
        marchenko.pastur.max <- (1+sqrt(M/N))^2
        pca.sig <- pca.eigenvalue > (marchenko.pastur.max * factor)
    if (do.print) {
        print(paste("Marchenko-Pastur eigenvalue null upper bound:", marchenko.pastur.max))
    if (factor != 1) {
      print(paste(length(which(pca.sig)), "PCs have eigenvalues larger than", factor, "times null upper bound."))
    } else {
      print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), "PCs have larger eigenvalues."))
    }}
  pca.sig.length = length(pca.sig[pca.sig==TRUE])
  return(pca.sig.length)
}

    M = dim(so$RNA@data)[1]
    N = dim(so$RNA@data)[2]
    pca.sdev = so@reductions$pca@stdev
    pca.sig.num = pcaMarchenkoPastur(M=M,N=N,pca.sdev = pca.sdev)
    lab2 = paste0("MP = ", pca.sig.num)
    xpos2 = pca.sig.num+4
    }
        ep <- ElbowPlot(so,ndims=30) + theme_bw() + 
        ggtitle(paste0(so@project.name," Elbow Plot")) 

        if(exists("lab")){
        ep <- ep + 
        geom_vline(xintercept = pcs, color="red") +
        annotate("text",  x=xpos, y = 4, label = lab, color="red",size=4) 
        }
        if(exists("lab2")){
         ep <- ep + 
         geom_vline(xintercept = pca.sig.num, color="blue") +
         annotate("text",  x=xpos2, y = 6, label = lab2, color="blue",size=4)
        }
        return(ep)
        #print(ep)
# auto removed:         #return(NULL)
    }

    if(is.null(vars_to_plot)){
        vars_to_plot = "nCount_RNA"
    }
    len <- length(vars_to_plot)*2
    grobsList <- vector(mode = "list", length = len)

    k=1
    for (i in 1:length(vars_to_plot)){ 
        grob <- lapply(so_orig, function(x) plotPCA(x,vars_to_plot[i]))
        grob=grid.arrange(grobs=grob,nrow=length(grob))
        grobsList[[k]] <- grob
        grob2 <- lapply(so_list, function(x) plotPCA(x,vars_to_plot[i]))
        grob2=grid.arrange(grobs=grob2,nrow=length(grob2))
        l=k+1
        grobsList[[l]] <- grob2
        k=k+2   
    }
    
    grob3 <- lapply(so_list, function(x) plotElbow(x))
    grob3=grid.arrange(grobs=grob3,nrow=length(grob3))

    grobsList[[length(grobsList)+1]] <- grob3
    
    if (doJackStraw) {
        grob4 <-lapply(so_list, function(x) JackStraw(x, reduction = "pca", dims = 5,
  num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
  maxit = 1000))
        grob4 <- lapply(grob4, function(x) ScoreJackStraw(x, dims = 1:5))
        grob4 <- lapply(grob4, function(x) JackStrawPlot(x, dims = 1:5))
        grob4=grid.arrange(grobs=grob4,nrow=length(grob4))
        grobsList[[length(grobsList)+1]] <- grob4
    }

    grobs <- arrangeGrob(grobs=grobsList,ncol=length(grobsList),newpage=F)

    imageWidth = 1000*2*imageCols
    imageHeight = 1000*length(so_list)
    dpi = 300

   if (imageType == "png") {
    png(
      filename="unnamed_9.png",
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
        file="unnamed_9.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }
    
    plot(grobs)

    cat("\nPCA Object Checksum:\n")
    print(digest::digest(so_list))
    return(list(value=so_list))
}

#################################################
## Global imports and functions included below ##
#################################################
# 
#   }
# 
# 

# Functions defined here will be available to call in
# the code for any table.

print("template_function_unnamed_9.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Postfilter_QC_Plots<-readRDS(paste0(rds_output,"/var_Postfilter_QC_Plots.rds"))
# Input_is_Seurat_count <- 0
# for(item in var_Postfilter_QC_Plots){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
# if(Input_is_Seurat_count == 0 ){
# var_Postfilter_QC_Plots<-as.data.frame(var_Postfilter_QC_Plots)}else{var_Postfilter_QC_Plots <- var_Postfilter_QC_Plots}
invisible(graphics.off())
var_unnamed_9<-unnamed_9(var_Postfilter_QC_Plots)
invisible(graphics.off())
saveRDS(var_unnamed_9, paste0(rds_output,"/var_unnamed_9.rds"))
