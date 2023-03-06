cd24_separate <- function(renamed_so,Sample_Names_Reclustered) {  # produces 2 color channels and the overlay
    require(ggplot2)
    require(gridExtra)
    library(Seurat)
    library(tidyverse)

  so = renamed_so$value
  #print(rownames(so[['Protein']]))

  samples = c("Control","Treated")

    ## Goal is to have column 1 of the new metadata be named "orig.ident" for downstream compatibility.
    ## Check new metadata for "orig.ident" column, else fix the "orig_ident" column name, else print an error message.
    if ("orig.ident" %in% colnames(so@meta.data)) { ## If orig.ident already is the first column ...
        print("Found orig.ident in column 1 of SO metadata.")
    } else if ("orig_ident" %in% colnames(so@meta.data)) { ## Else if "orig_ident" is the first column ...
        colnames(so@meta.data)[colnames(so@meta.data) == "orig_ident"] <- "orig.ident"
        print("Found orig_ident in column 1 of new metadata table. Changed to orig.ident for downstream compatibility.")
    } else { ## Else print an error message explaining we expect one of the two above as the first column in the new metadata.
        print("ERROR: Found neither orig.ident nor orig_ident in column 1 of new metadata table. Please try again with a new metadata table with one of these as the column name of the first column in the dataframe.")
    }

    if("active.ident" %in% slotNames(so)){
        sample_name = as.factor(so@meta.data$orig.ident)
        names(sample_name)=names(so@active.ident)
        so@active.ident <- as.factor(vector())
        so@active.ident <- sample_name
        so.sub = subset(so, ident = samples)
    } else {
        sample_name = as.factor(so@meta.data$orig.ident)
        names(sample_name)=names(so@active.ident)
        so@active.ident <- as.factor(vector())
        so@active.ident <- sample_name
        so.sub = subset(so, ident = samples)
    }
  
    gg.overlay <- function(so.sub,df,marker1,marker2){
    
    df %>% dplyr::arrange(mark1.scale) -> df    
    xmin = min(df$dr1) - 0.1*min(df$dr1)
    xmax = max(df$dr1) + 0.1*min(df$dr1)

    gg.z1 <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape=16,size=0.5, alpha=0.5)+
    theme_classic() +
    xlab("umap-1") + 
    ylab("umap-2") +
    ggtitle(marker1) +
    coord_fixed()
  
    df %>% dplyr::arrange(mark2.scale) -> df
    
    gg.z2 <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape=16,size=0.5, alpha=0.5)+
    theme_classic() +
    xlab("umap-1") + 
    ylab("umap-2") +
    ggtitle(marker2) +
    coord_fixed()
  
    df %>% dplyr::mutate(avg = mark2.scale+mark1.scale) %>% dplyr::arrange(avg) -> df

    gg <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape=16,size=0.5, alpha=0.5)+
    theme_classic() +
    xlab("umap-1") + 
    ylab("umap-2") +
    ggtitle("Combined") +
    coord_fixed()

    return(list(gg.z1,gg.z2,gg))
    }

    t1 = 0.5
    t2 = 0.5
    addlines = TRUE

    gg.overlay2 <- function(so.sub,df,marker1,marker2){
    
    df %>% dplyr::arrange(mark1.scale) -> df    

# Create unscaled axis labels
display_unscaled <- FALSE

if(display_unscaled == TRUE){
label1_min <- paste("unscaled min:", round(min(mark1),digits = 2))
label1_max <- paste("unscaled max:", round(max(mark1),digits = 2))
label1 <- paste(as.character(marker1), label1_min, label1_max, sep = "\n")

label2_min <- paste("unscaled min:", round(min(mark2),digits = 2))
label2_max <- paste("unscaled max:", round(max(mark2),digits = 2))
label2 <- paste(as.character(marker2), label2_min, label2_max, sep = "\n")} else {

    label1 <- as.character(marker1)
    label2 <- as.character(marker2)
    
}

    gg.z1 <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape = 20,size=0.5)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()
  
    df %>% dplyr::arrange(mark2.scale) -> df
    
    gg.z2 <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape = 20,size=0.5)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()
  
    df %>% dplyr::mutate(avg = mark2.scale+mark1.scale) %>% dplyr::arrange(avg) -> df
    
    gg <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape = 20,size=0.5)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()

    if(addlines==TRUE){
        gg.z1 <- gg.z1 + 
            geom_vline(xintercept=t1,linetype="dashed") +
            geom_hline(yintercept=t2,linetype="dashed") 
        gg.z2 <- gg.z2 + 
            geom_vline(xintercept=t1,linetype="dashed") +
            geom_hline(yintercept=t2,linetype="dashed") 
        gg <- gg +
            geom_vline(xintercept=t1,linetype="dashed") +
            geom_hline(yintercept=t2,linetype="dashed") 
    }

    return(list(gg.z1,gg.z2,gg))
    }

    marker1 <- "Cd24a"
    marker2 <- "CD24"
    scale1 <- TRUE
    scale2 <- TRUE

    mark1 = so.sub@assays$SCT@scale.data[marker1,]
    if(scale1){
    q1 = quantile(mark1,0.80)
    q0 = quantile(mark1,0.5)
    mark1[mark1<q0]=q0
    mark1[mark1>q1]=q1
    }
    mark1.scale <- scales::rescale(mark1, to=c(0,1))

    mark2 = so.sub@assays$Protein@scale.data[marker2,]
    if(scale2){
    q1 = quantile(mark2,0.90)
    q0 = quantile(mark2,0.5)
    print(q0)
    mark2[mark2<q0]=q0
    mark2[mark2>q1]=q1
    }
    mark2.scale <- scales::rescale(mark2, to=c(0,1))

    #mark1.scale <- scales::rescale(so.sub@assays$SCT@data[marker1,], to=c(0,1))
    #mark2.scale <- scales::rescale(so.sub@assays$Protein@data[marker2,], to=c(0,1))
    df <- data.frame(cbind(dr1=so.sub@reductions$umap@cell.embeddings[,1],
                           dr2=so.sub@reductions$umap@cell.embeddings[,2],
                            mark1.scale,mark2.scale))
    gg.list <- gg.overlay(so.sub,df,marker1,marker2)
    gg.list2 <- gg.overlay2(so.sub,df,marker1,marker2)

    n=3
    imageWidth = 1000*n
    imageHeight = 1000*2
    dpi = 300

    png(
      filename="cd24_separate.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

if (FALSE){
    x <- df$mark1.scale
    y <- df$mark2.scale

    df_heatmap <- data.frame(x = x, y = y,
                 d = densCols(x, y, colramp = colorRampPalette(rev                      (rainbow(10, end = 4/6)))))

    p <- ggplot(df_heatmap) +
        geom_point(aes(x, y, col = d), size = 1) +
        scale_color_identity() + xlab(marker1) + ylab(marker2) +
        theme_bw() + geom_vline(xintercept=t1,linetype="dashed") +
        geom_hline(yintercept=t2,linetype="dashed") 

grid.arrange(gg.list[[1]], gg.list[[2]], gg.list[[3]], gg.list2[[1]],gg.list2[[2]],gg.list2[[3]],p,ncol=3)} else {
    grid.arrange(gg.list[[1]], gg.list[[2]], gg.list[[3]], gg.list2[[1]],gg.list2[[2]],gg.list2[[3]],ncol=3)
}

     
    filterdata = FALSE
    if(filterdata){
        M1_filt_direction = "greater_than"
        M2_filt_direction = "greater_than"
    
    df = df %>% 
        mutate(sample=so.sub@meta.data$sample) %>% 
        mutate(cellbarcode=rownames(so.sub@meta.data))
    
    if(M1_filt_direction == "greater_than"){
            ind1 <- df$`mark1.scale` > t1
        }else{
            ind1 <- df$`mark1.scale` < t1
        }

    cat("\n")
    print("Marker 1 filter:")
    print(sum(ind1))

    if(M2_filt_direction == "greater_than"){
            ind2 <- df$`mark2.scale` > t2
        }else{
            ind2 <- df$`mark2.scale` < t2
        }
    
    cat("\n")
    print("Marker 2 filter:")
    print(sum(ind2))

    applyfilt1 <- TRUE
    applyfilt2 <- FALSE

    if(applyfilt1){
        if(applyfilt2){
          if(TRUE){
              df <- df[c(ind1&ind2),]
          }
          else{
              df <- df[c(ind1|ind2),]  
          }  
        }
        else{
          df <- df[ind1,]   
        }
    }
    else{
        if(applyfilt2){
            df <- df[ind2,]    
    }
    }

    colnames(df)[3:4]= c(marker1,marker2)
    so.sub@meta.data %>% 
            mutate(Marker=case_when(rownames(so.sub@meta.data) %in% df$cellbarcode ~ TRUE, TRUE ~ FALSE)) -> so.sub.df
    
    cat("\n")
    print("Final Breakdown:")
    print(addmargins(table(so.sub.df$Marker,so.sub.df$sample_name)))
    rownames(so.sub.df) <- rownames(so.sub@meta.data)
    so.sub@meta.data <- so.sub.df
    
    cat("\n")
    print("After Filter applied:")
    print(dim(df)[1])
    }
return(list(value=so.sub))
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_cd24_separate.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_renamed_so<-readRDS(paste0(rds_output,"/var_renamed_so.rds"))
Input_is_Seurat_count <- 0
for(item in var_renamed_so){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_renamed_so<-as.data.frame(var_renamed_so)}else{var_renamed_so <- var_renamed_so}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Sample_Names_Reclustered<-readRDS(paste0(rds_output,"/var_Sample_Names_Reclustered.rds"))
Input_is_Seurat_count <- 0
for(item in var_Sample_Names_Reclustered){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Sample_Names_Reclustered<-as.data.frame(var_Sample_Names_Reclustered)}else{var_Sample_Names_Reclustered <- var_Sample_Names_Reclustered}
invisible(graphics.off())
var_cd24_separate<-cd24_separate(var_renamed_so,var_Sample_Names_Reclustered)
invisible(graphics.off())
saveRDS(var_cd24_separate, paste0(rds_output,"/var_cd24_separate.rds"))
