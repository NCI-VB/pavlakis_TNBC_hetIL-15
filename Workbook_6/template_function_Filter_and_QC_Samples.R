# Filter & QC Samples [scRNA-seq] [CCBR] (c4dccc54-4b4e-4860-af86-6398a6562013): v147
Filter_and_QC_Samples <- function(hetIL15_citeseq) {
    #image:png
    imageType = "png"
     suppressMessages(library(Seurat))
     suppressMessages(library(reshape2))
     suppressMessages(library(tidyverse))
     suppressMessages(library(gridExtra))
     suppressMessages(library(RColorBrewer))
     suppressMessages(library(stringr))

localFilePaths <- readRDS("./rds_output/var_hetIL15_citeseq.rds")
localFilePaths <- readRDS("./rds_output/var_hetIL15_citeseq.rds")

    subsetRegex = eval(parse(text=gsub('\\[\\]','c()','c("CD11.protein.h5")')))
    Keep <- TRUE
    if (length(subsetRegex) > 0) {
        if (Keep == TRUE){
        for (i in length(subsetRegex)) {
localFilePaths <- readRDS("./rds_output/var_hetIL15_citeseq.rds")
        }
        }
        else{
            for (i in length(subsetRegex)) {
localFilePaths <- readRDS("./rds_output/var_hetIL15_citeseq.rds")
        }
        }
    }
    
    ##force path to correct file
    localFilePaths <- c("nidap_downloads/C_TUMOR_CD11.protein.h5","nidap_downloads/IL15_TUMOR_CD11.protein.h5") 
    
    #obj.list <- lapply(localFilePaths, function(x) { return(c(Read10X_h5(x, use.names=FALSE), Read10X_h5(x, use.names=TRUE))) })
    obj.list <- lapply(localFilePaths, function(x) { return(Read10X_h5(x, use.names=TRUE)) })
    rename = TRUE
    if (rename == FALSE){
    names(obj.list) <- lapply(localFilePaths, basename)
    names(obj.list) <- sapply(names(obj.list), function(x) gsub("_filtered(\\w+)?.h5","", x))
    names(obj.list) <- sapply(names(obj.list), function(x) gsub("\\.(\\w+)?","", x))
    }
    else{
    names(obj.list) <- c("Control","Treated")
    obj.list <- obj.list[sort(names(obj.list))]
    }
     #obj.list <- hetIL15_citeseq$value
     
     mincells = 3
     mingenes = 200
     organism = "Mouse"

     if (organism == "Human"){
         mitoch = "^MT-"
     } else{
         mitoch = "^mt-"
         cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
         cc.genes$s.genes = str_to_title(cc.genes$s.genes)
     }

    seurat_object <- function(i) {
         if (class(obj.list[[i]]) == "dgCMatrix"){
        so.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
        } else {
        so.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
        }
        so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
        so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
        so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
        
        if (TRUE & "Antibody Capture" %in% names(obj.list[[i]])){
        antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
        so.nf[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so.nf)])
        so.nf <- NormalizeData(so.nf, assay = "Protein", normalization.method = "CLR")
        }
        if (FALSE){
        antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
        HTO = rownames(obj.list[[i]][2]$`Antibody Capture`)[grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
        so.nf[['HTO']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[HTO, colnames(x = so.nf)])
        so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
        #so.nf <- MULTIseqDemux(so.nf, assay="HTO")
        }

    #Filtered Seurat Object:
        if (class(obj.list[[i]]) == "dgCMatrix"){
        so <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
        } else {
        so <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
        }
        so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
        so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = mitoch)
       # so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
        so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
        
        if (TRUE & "Antibody Capture" %in% names(obj.list[[i]])){
        rownames(obj.list[[i]][2]$`Antibody Capture`) <- gsub(pattern = "_TotalSeqC", replacement = "", rownames(obj.list[[i]][2]$`Antibody Capture`))
        antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
        so[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so)])
        so <- NormalizeData(so, assay = "Protein", normalization.method = "CLR")
        }
        if (FALSE){
        rownames(obj.list[[i]][2]$`Antibody Capture`) <- gsub(pattern = "_TotalSeqC", replacement = "", rownames(obj.list[[i]][2]$`Antibody Capture`))
        HTO = rownames(obj.list[[i]][2]$`Antibody Capture`)[grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
        so[['HTO']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[HTO, colnames(x = so)])
        so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")
        #so <- MULTIseqDemux(so, assay="HTO")
        }
        if (FALSE) {
            allGenes = rownames(so)
            VDJgenes = c("TRBV","TRAV","TRBD","TRAJ","TRBJ")
            print("Removing VDJ genes. Genes removed...")
            for (j in VDJgenes) {
                print(allGenes[grepl(j, allGenes)])
                allGenes = allGenes[!grepl(j, allGenes)]  
            }
            so <- SubsetData(so,features = allGenes,assay="RNA")
        }

    cat("\n\n")
    cat(names(obj.list)[i],":\n")
    so.origcount = dim(so.nf)[2]
    cat(paste0("Original Cell Count=", so.origcount),"\n")

    #Start with filtering here:
        maxgenes = 2500
        complexity = 0.5
        MAD_gene <- TRUE
        ngenestdev <- mad(so@meta.data$nFeature_RNA)
        ngenemed <- median(so@meta.data$nFeature_RNA)
        ngenemaxlim <- ngenemed+(3*ngenestdev)
        gl = format(round(ngenemaxlim,0),nsmall=0)

        maxmitoch = 10

        MAD_mitoch <- TRUE
        mitostdev <- mad(so@meta.data$percent.mt)
        mitomed <- median(so@meta.data$percent.mt)
        mitomaxlim <- mitomed+(3*mitostdev)
        ml = format(round(mitomaxlim,2),nsmall=2)

      if (MAD_gene == TRUE & MAD_mitoch == TRUE)       {
        cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
        cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
        cat(paste0("Complexity Filter =",complexity,"\n"))
        so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
        perc.remain = (dim(so)[2]/so.origcount)*100
        perc.remain=formatC(perc.remain,format = "g",digits=3)
        cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
        cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
      }
       else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
        cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
        cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
        so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
        perc.remain = (dim(so)[2]/so.origcount)*100
        perc.remain=formatC(perc.remain,format = "g",digits=3)
        cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
        cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
        }
        else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
        cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
        cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
        so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
         perc.remain = (dim(so)[2]/so.origcount)*100
        perc.remain=formatC(perc.remain,format = "g",digits=3)
        cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
        cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
        }
        else {
        cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
        cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
        so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
        perc.remain = (dim(so)[2]/so.origcount)*100
        perc.remain=formatC(perc.remain,format = "g",digits=3)
        cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
        cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
        }

        plothist <- function(count.df,name){
            g=ggplot(count.df,aes(x=value,fill=filt)) + 
            theme_bw() +
            geom_histogram(binwidth=.05, alpha = 0.7, position="identity") +
            #geom_density(aes(x = value, colour = filt, linetype = filt)) +
            scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
            scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
            labs(x = NULL) +
            theme(plot.title = element_text(size=6),legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
            ggtitle(paste(name,count.df$variable[1])) +
            scale_x_continuous(trans='log10') + 
            scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
            return(g)
            #print(g)
# auto removed:             #return(NULL)
        }

        plotviolin <- function(count.df,name){
            axislab = unique(count.df$filt)
            col1=brewer.pal(8, "Set3")[-2] 
            col2=c(col1,brewer.pal(8,"Set2")[3:6])

            v = ggplot(count.df, aes(x=filt, y=value)) +
            ggtitle(paste(name,count.df$variable[1])) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 12, face = "bold")) +
            geom_violin(aes(fill=as.factor(filt))) +  
            scale_fill_manual(values = c("#00AFBB", "#FC4E07")) + 
            geom_boxplot(width=.1) +
            #labs(colour = n, y=m) +
            #geom_jitter(height = 0, width = 0.1, size = 0.1) +
            #scale_y_continuous(trans='log2') + 
            scale_x_discrete(limits = as.vector(axislab)) 
            return(v)
            #print(v)
# auto removed:             #return(NULL)
        }

        Runplots <- function(x,name){
                df.m %>% dplyr::filter(variable == x) -> count.df
                df2.m %>% dplyr::filter(variable == x) -> count2.df
                qc.df <- array(0,dim=c(0,4))
                qc.df <- rbind(qc.df,count2.df,count.df)
                if(FALSE){
                gg <- plothist(qc.df,name)}
                else{
                gg <- plotviolin(qc.df,name)
        }
        }

        RunScatter <- function(x,name){
            x <- as.character(x)
            scplot.m = so@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "filt")
            scplot2.m = so.nf@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "raw") 
            sc.plot.all = rbind(scplot2.m,scplot.m)
            g=ggplot(sc.plot.all,aes_string(x="nCount_RNA",y=x,color="filt")) + 
                geom_point(size = 0.5) + 
                theme_classic() +
                ggtitle(paste(name)) 
            #print(g)
# auto removed:             #return(NULL)
            return(g)
        }
      
        df.m <- melt(so@meta.data)
        df.m$filt <- "filt"
        df.m$filt <- as.factor(df.m$filt)
        df2.m <- melt(so.nf@meta.data)
        df2.m$filt <- "raw"
        df2.m$filt <- as.factor(df2.m$filt)

        v <- unique(df.m$variable)
        grob.list <- lapply(v,function(x){Runplots(x,so@project.name)})
        grob2.list <- lapply(v,function(x){RunScatter(x, so@project.name)})
        grob.all <- arrangeGrob(grobs = grob.list, ncol = length(grob.list))
        grob2.all <- arrangeGrob(grobs = grob2.list, ncol = length(grob2.list))
        
        slot(so,"commands") <- list() #clear command with timestamp for consistent checksum
        so2.list <- list(so,so.nf,grob.all,grob2.all)
        
        return(so2.list)
    }

     so.list <- lapply(seq_along(obj.list), seurat_object)
     
    so.f.list <- lapply(so.list,function(x) x[[1]])
    names(so.f.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))
   

    so.nf.list <- lapply(so.list,function(x) x[[2]])
    names(so.nf.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))
    
     so.final.list <<- so.f.list

    so.grobs.list <- lapply(so.list,function(x) x[[3]])
    so.grobs2.list <- lapply(so.list,function(x) x[[4]])

    cat("Final filtered samples:\n")
    print(so.f.list)
    cat("Final filtered sample names:\n")
    print(names(so.f.list))   

    imageWidth = min(1000*length(so.list[[1]][[3]]),15000)
    imageHeight = min(1000*length(so.grobs.list)*2,24000)
    dpi = 300

    if (imageType == 'png') {
    png(
      filename="Filter_and_QC_Samples.png",
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
        file="Filter_and_QC_Samples.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    #grid.arrange(grobs = so.grobs.list, nrow = length(so.grobs.list))
grobdat = list()
    for(i in 1:length(so.grobs.list)){grobdat=append(grobdat,list(so.grobs.list[[i]])) }
    for(i in 1:length(so.grobs2.list)){grobdat=append(grobdat,list(so.grobs2.list[[i]])) }
    
    grid.arrange((arrangeGrob(grobs=grobdat,nrow=length(grobdat))),nrow=1)

    cat("Return objects checksum:\n")
    #print(digest::digest(so.final.list))
    return(list(value=so.final.list))  
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

print("template_function_Filter_and_QC_Samples.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_hetIL15_citeseq<-readRDS(paste0(rds_output,"/var_hetIL15_citeseq.rds"))
Input_is_Seurat_count <- 0
for(item in var_hetIL15_citeseq){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_hetIL15_citeseq<-as.data.frame(var_hetIL15_citeseq)}else{var_hetIL15_citeseq <- var_hetIL15_citeseq}
invisible(graphics.off())
var_Filter_and_QC_Samples<-Filter_and_QC_Samples(var_hetIL15_citeseq)
invisible(graphics.off())
saveRDS(var_Filter_and_QC_Samples, paste0(rds_output,"/var_Filter_and_QC_Samples.rds"))
