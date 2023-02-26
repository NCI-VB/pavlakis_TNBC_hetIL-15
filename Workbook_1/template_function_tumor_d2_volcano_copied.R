tumor_d2_volcano_copied <- function(Deg_tumor) {

    suppressMessages(library(ggplot2))
    suppressMessages(library(dplyr))
    suppressMessages(library(ggrepel))
    genesmat <- dplyr::collect(Deg_tumor)

    genesmat %>% dplyr::arrange(`IL15_day_2-control_day_2_pval`) -> genesmat

   # print(genesmat)

    print(paste0("Total number of genes included in volcano plot: ", nrow(genesmat)))

    negative_log10_p_values <- -log10(genesmat$`IL15_day_2-control_day_2_pval`)
    ymax <- ceiling(max(negative_log10_p_values[is.finite(negative_log10_p_values)]))

    xmax = ceiling(max(genesmat$`IL15_day_2-control_day_2_logFC`))
    xmin = floor(min(genesmat$`IL15_day_2-control_day_2_logFC`))

    nudge_x_all <- 0.8
    nudge_y_all <- 0.8

   
    new_contrast_label <- "IL15_day_2-control_day_2_logFC"
    
    pointsize = 4

        png(
      filename="tumor_d2_volcano_copied.png",
      width=7,
      height=5,
      units="in",
      pointsize=4,
      bg="white",
      res=300,
      type="cairo")
      
    p <- ggplot(genesmat,
        aes_(x = as.name(new_contrast_label), y = quote(-log10(`IL15_day_2-control_day_2_pval`)))) + # modified by RAS
        theme_classic() +
        geom_point(
            color='gray',
            size = pointsize-3) +
        geom_vline(xintercept=c(-1,1), color='black', alpha=1.0,linetype="dotted") + 
        geom_hline(yintercept=1.5, color='black', alpha=1.0, linetype="dotdash") +  
        geom_point(
            data = subset(genesmat, Gene %in% c("Ccl9","Cxcr3","Ccl19")),
            color = 'blue',
            size = pointsize) +
                geom_point(
            data = subset(genesmat, Gene %in% c("Cd247","Zap70","Cd3d")), 
            color = 'green',
            size = pointsize) +
        geom_point(
            data = subset(genesmat, Gene %in% c("Gzma","Gzmb","Ctsw","Klrg1","Prf1")), 
            color = 'red',
            size = pointsize) +
        geom_text_repel(
            data = genesmat[genesmat$Gene %in% c("Cd247","Cd3d","Ccl19","Ccl9","Cxcr3","Zap70","Gzma","Gzmb","Ctsw","Klrg1","Prf1","Ifng"),], 
            label = genesmat[genesmat$Gene %in% c("Cd247","Cd3d","Ccl19","Ccl9","Cxcr3","Zap70","Gzma","Gzmb","Ctsw","Klrg1","Prf1","Ifng"),"Gene"], 
            color = "black",
            fontface = 1,
            segment.size = 0.5,
            nudge_x = nudge_x_all,
            nudge_y = nudge_y_all,
            size=7)+#,
            # nudge_y = nudge_y_all,
            # size = 3) +
        #xlim(-4.5,6) +
        ylim(0,ymax) +
        labs(x="Log2 Fold-change", y="-log10(p-value)") +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank())

    print(p)
# auto removed:     return(NULL)
    save_img=TRUE
    if(save_img){
        return(list(value=p))
    }else{
# auto removed:         return(NULL)
    }
}

#################################################
## Global imports and functions included below ##
#################################################
#BiocManager::install("GenomeInfoDbData", destdir="/scratch")

# Functions defined here will be available to call in
# the code for any table.

print("template_function_tumor_d2_volcano_copied.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Deg_tumor<-readRDS(paste0(rds_output,"/var_Deg_tumor.rds"))
Input_is_Seurat_count <- 0
for(item in var_Deg_tumor){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Deg_tumor<-as.data.frame(var_Deg_tumor)}else{var_Deg_tumor <- var_Deg_tumor}
invisible(graphics.off())
var_tumor_d2_volcano_copied<-tumor_d2_volcano_copied(var_Deg_tumor)
invisible(graphics.off())
saveRDS(var_tumor_d2_volcano_copied, paste0(rds_output,"/var_tumor_d2_volcano_copied.rds"))
