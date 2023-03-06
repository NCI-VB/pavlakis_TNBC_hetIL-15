Plot_Metadata <- function(Labeled_so,Get_Sample_Names,Get_Metadata) {
    #image: png
    
    imageType="png"

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

    SO = Labeled_so
    doCiteSeq <- FALSE
    
    summarizeCutOff <-min(5,20) 

    samples = eval(parse(text=gsub('\\[\\]','c()','c()')))
    if (length(samples) == 0) {
        print("No samples specified. Using all samples...")
        samples = unique(SO@meta.data$sample_name)
    }

    ## Goal is to have column 1 of the new metadata be named "orig.ident" for downstream compatibility.
    ## Check new metadata for "orig.ident" column, else fix the "orig_ident" column name, else print an error message.
    if ("orig.ident" %in% colnames(SO@meta.data)) { ## If orig.ident already is the first column ...
        print("Found orig.ident in column 1 of SO metadata.")
    } else if ("orig_ident" %in% colnames(SO@meta.data)) { ## Else if "orig_ident" is the first column ...
        colnames(SO@meta.data)[colnames(SO@meta.data) == "orig_ident"] <- "orig.ident"
        print("Found orig_ident in column 1 of new metadata table. Changed to orig.ident for downstream compatibility.")
    } else { ## Else print an error message explaining we expect one of the two above as the first column in the new metadata.
        print("ERROR: Found neither orig.ident nor orig_ident in column 1 of new metadata table. Please try again with a new metadata table with one of these as the column name of the first column in the dataframe.")
    }

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
    
    print("selected SO:")
    print(SO)

    meta.df <- dplyr::collect(Get_Metadata)
    colnames(SO.sub@meta.data) = gsub("\\.","_",colnames(SO.sub@meta.data))

    m = eval(parse(text=gsub('\\[\\]','c()','c()')))
    m = m[!grepl("Barcode",m)]
    if (length(m) == 0) {
        print("No metadata columns specified. Plotting sample names and RNA clusters...")
        x = colnames(SO.sub@meta.data)
        x = x[grepl("RNA",x)]
        m = c("sample_name",x)
    }

    #ERROR CATCHING
    #collect valid names of valid columns
    validColumns <- character()
    for (i in colnames(meta.df)) {
        if (!any(is.na(meta.df[[i]]))) {
            validColumns <-c(validColumns,i)
        }
    }

    colsToSummarize <- eval(parse(text=gsub('\\[\\]','c()','c()')))
    m = unique(c(m,colsToSummarize))

    if (length(colsToSummarize)>0) {
    #Optional Summarization of Metadata
        for (i in colsToSummarize) {
            col <- meta.df[[i]]
            valCount <- length(unique(col))

            if ((valCount >=summarizeCutOff) & (i != 'Barcode') & (!is.element(class(meta.df[[i]][1]),c("numeric","integer")))) {
                freqVals <- as.data.frame(-sort(-table(col)))$col[1:summarizeCutOff]
                print(freqVals)
                summarized_col = list()
                count <- 0
                for (j in col) {

                    print(j)
                    print(paste("count is",count))
                    
                    if (is.na(j) || is.null(j) || (j =="None")) {
                        count <- count + 1
                        summarized_col[count] <- "NULLorNA"
                        print("NULLorNA")
                    } else if (j %in% freqVals){
                        count <- count + 1
                        summarized_col[count] <- j
                        print("valid")
                    } else {
                        count <- count + 1
                        summarized_col[count] <- "Other"
                        print("Other")
                    }
                }
                meta.df[[i]] <- summarized_col
            }
        }
        #assign new metadata
        colnames(meta.df) = gsub("\\.","_",colnames(meta.df))
        SO.sub@meta.data <- meta.df
        colnames(SO.sub@meta.data) = gsub("\\.","_",colnames(SO.sub@meta.data))
    }

    drawMetadata <- function(m){
    #check if there are NaNs in metadata, if there are, catch 
    if (any(is.na(meta.df[[m]]))) {
        print("ERROR: Metadata column appears to contain NA values. This is not recommended for clustering plots.")
        print("Please review your selected metadata column")
        print(head(meta.df[[m]]))
        print("Below are valid metadata to select for this plot:")
        print(validColumns)
        stop("End of error message.")
    }

    reduction = "umap"

    if (!(doCiteSeq)) {
        if(reduction=="tsne"){
            p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident")
        } else if(reduction=="umap"){
            p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
        } else { 
            p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
        }
    } else {
        if(reduction=="tsne"){
            p1 <- DimPlot(SO.sub, reduction = "protein_tsne", group.by = "ident")
        } else if(reduction=="umap"){
            p1 <- DimPlot(SO.sub, reduction = "protein_umap", group.by = "ident")
        } else { 
            p1 <- DimPlot(SO.sub, reduction = "protein_pca", group.by = "ident")
        }
    }

    #Categorical/Qualitative Variables
    if (!is.element(class(meta.df[[m]][1]),c("numeric","integer"))){
        if (!(doCiteSeq)) {
            #plot RNA clusters
        if(reduction=="tsne"){
            clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.character(SO.sub@meta.data[[m]]))
        } else if(reduction=="umap"){
            clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.character(SO.sub@meta.data[[m]]))
        } else { 
            clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.character(SO.sub@meta.data[[m]]))
        }

        } else {
            #else plot Antibody clusters
            if(reduction=="tsne"){
            clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.character(SO.sub@meta.data[[m]]))
        } else if(reduction=="umap"){
            clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.character(SO.sub@meta.data[[m]]))
        } else { 
            clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.character(SO.sub@meta.data[[m]]))
        }
        }

        clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
        title=as.character(m)
        cols=list()

        #Lab.palette <- colorRampPalette(brewer.pal(12,"Set3"))
        Lab.palette <- colorRampPalette(brewer.pal(12,"Paired"))
        n=length(unique((SO.sub@meta.data[[m]])))
        cols[[1]]=brewer.pal(8, "Set3")[-2]  #Alternative
        cols[[2]]=brewer.pal(8, "Set1")
        cols[[3]]=c(cols[[1]],brewer.pal(8,"Set2")[3:6])
        cols[[4]]=c("#F8766D","#FF9912","#a100ff","#00BA38","#619CFF","#FF1493","#010407")
        cols[[5]]=c("blue","red","grey")
        cols[[6]]=Lab.palette(n)
        cols[[7]]=c("red","green","blue","orange","cyan","purple")
        cols[[8]]=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9","#808080","#A9A9A9","#8B7355")
        colnum = 8

        n = length(unique(clusmat$clusid))
        #col=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075",Lab.palette(max(0,n-15)))[1:n]
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
        col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

        #Select to add labels to plot or not:
        if(FALSE){
            g <- ggplot(clusmat) +
                theme_bw() +
                theme(legend.title=element_blank()) +
                geom_point(aes(x=umap1, y=umap2,colour=clusid),size=0.01) +
                scale_color_manual(values=col) +
                xlab(paste(reduction,"-1")) + ylab(paste(reduction,"-2")) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="top",
                panel.background = element_blank(), legend.text=element_text(size=rel(1))) +
                guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
                ggtitle(title) +
                geom_label_repel(data=umap.pos,aes(x=umap1.mean,y=umap2.mean,label=umap.pos$clusid),size=4)
            }
        else{
            g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
                theme_bw() +
                theme(legend.title=element_blank()) +
                geom_point(aes(colour=clusid),size=0.01) +
                scale_color_manual(values=col) +
                xlab(paste(reduction,"-1")) + ylab(paste(reduction,"-2")) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top",
                panel.background = element_blank(),legend.text=element_text(size=rel(1))) +
                guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
                ggtitle(title)    
            } 

    } else {
        ##THIS IS THE PART WE PLOT QUANTITATIVE DATA
        m = as.character(m)
        clusid = SO.sub@meta.data[[m]]
        clusid = scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
        clus.quant=quantile(clusid[clusid>0],probs=c(.1,.5,.9))
        midpt = clus.quant[2]
        midpt2 = clus.quant[1]
        midpt3 = clus.quant[3]
        #hist(clusid[!is.na(clusid)], breaks=100, main=m)
        #abline(v=midpt,col="red",lwd=2)

        if (!(doCiteSeq)) {
            #plot RNA clusters
        if(reduction=="tsne"){
            clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        } else if(reduction=="umap"){
            clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        } else { 
            clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        }

        } else {
            #else plot Antibody clusters
            if(reduction=="tsne"){
            clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        } else if(reduction=="umap"){
            clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        } else { 
            clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
        }
        }

            clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
            title=as.character(m)
        print(environmentName(environment(arrange))) 
        clusmat %>% dplyr::arrange(clusid) -> clusmat
        print(environmentName(environment(arrange))) 
        g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
            theme_bw() +
            theme(legend.title=element_blank()) +
            geom_point(aes(colour=clusid),size=1) +
            #scale_color_gradient2(low = "blue4", mid = "white", high = "red",
            #          midpoint = midpt[[p]], na.value="grey",limits = c(0, 1)) + 
            scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(0, midpt2,midpt,midpt3, 1), limits = c(0, 1))) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
            ggtitle(title) +
            xlab("umap-1") + ylab("umap-2")
    }

    return(g)
    }

grobs <- lapply(m, function(x) drawMetadata(x))

    if (0 == 0) {
        n = ceiling(length(m)^0.5)
    } else {
        n = 0
    }

    
    y=ceiling(length(m)/n)
    imageWidth = 2000*n
    imageHeight = 2000*y
    dpi = 300

    if (imageType == 'png') {
    png(
      filename="Plot_Metadata.png",
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
        file="Plot_Metadata.png",
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

grid.arrange(grobs = grobs,ncol = n,newpage=F)

if (FALSE){
return(list(value=SO))   
    } else{
 return(SO.sub@meta.data)
}
}

#################################################
## Global imports and functions included below ##
#################################################




    

print("template_function_Plot_Metadata.R #########################################################################")
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
var_Get_Sample_Names<-readRDS(paste0(rds_output,"/var_Get_Sample_Names.rds"))
Input_is_Seurat_count <- 0
for(item in var_Get_Sample_Names){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Get_Sample_Names<-as.data.frame(var_Get_Sample_Names)}else{var_Get_Sample_Names <- var_Get_Sample_Names}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Get_Metadata<-readRDS(paste0(rds_output,"/var_Get_Metadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_Get_Metadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Get_Metadata<-as.data.frame(var_Get_Metadata)}else{var_Get_Metadata <- var_Get_Metadata}
invisible(graphics.off())
var_Plot_Metadata<-Plot_Metadata(var_Labeled_so,var_Get_Sample_Names,var_Get_Metadata)
invisible(graphics.off())
saveRDS(var_Plot_Metadata, paste0(rds_output,"/var_Plot_Metadata.rds"))
