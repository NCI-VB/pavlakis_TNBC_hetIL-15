FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq <- function(renamed_so, Canonical_DC_Markers) {
    
    library(Seurat)
    library(rlang)
    library(RColorBrewer)
    library(tidyverse)
    library(cowplot)
    library("scales")
    library(colorspace)
    
    #Add other genes e
    
    so <-renamed_so$value

    so <- suppressMessages(expr = StashIdent(object = so, save.name = 'ident'))
    new_levels <- c("CCR7 High DC","CD103intCD11b+ DC","cDC1","cDC2","Siglec-H DC","Monocyte 1","Monocyte 2","Other")
    levels(so) <- new_levels
        
    so <- subset(so, idents = c("cDC1","Other"), invert = TRUE)
    
    assay = "SCT"
    scale.by = 'radius'
    split.by = "cluster"
    clustering_distance = "euclidean"
    clustering_method = "complete"
    dot.min = 0
    dot.scale = 6
    scale = TRUE
    col.min = -2.5
    col.max = 2.5

    features <- rownames(x = so)
    markers <- Canonical_DC_Markers
    markers$Origin <- factor(markers$Origin, levels=unique(markers$Origin))

    features_of_interest <- unique(markers$Matches)
    bad <- c("Ly6g","Itgad")

    features_of_interest <- features_of_interest[!features_of_interest %in% bad]
    if(!all(features_of_interest %in% features)){
        print(features_of_interest[!features_of_interest %in% features])
        stop("Could not find all features.")
    }

    default.colors <- c(hue_pal()(nlevels(markers$Origin)))
    names(default.colors) <- levels(markers$Origin)

    group.use2 <- markers$Origin
    names(group.use2) <- markers$Genes

    annot_cols <- t(default.colors)
    
    #define colors
    # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    # qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    # colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # colors = diverging_hcl(100,palette="Blue-Red 2")
    # colors = colors[c(1,length(colors))]
    # new.cluster.ids <- c("Monocytes", "CD103intCD11b+ DC", "cDC1","cDC2", "CCR7Hi", "Siglec-H DC")
    # so@meta.data$cluster <- factor(as.character(so@meta.data$cluster), levels = new.cluster.ids)
    # sample_name = so@meta.data$cluster
    # names(sample_name)=names(so@active.ident)
    # #so@active.ident <- as.factor(vector())
    # so@active.ident <- sample_name
    
    #so@meta.data$cluster <- factor(as.character(so@meta.data$cluster), levels=c("CCR7hi DC", "Siglec-H DC", "cDC1", "cDC2", "DC Cluster 1", "DC Cluster 2", "Monocyte"))    

    #so@meta.data$cluster <- factor(as.character(so@meta.data$cluster), levels = new.cluster.ids)
    #so@active.ident <- factor(as.character(so@active.ident), levels = new.cluster.ids)
    plot.data <- .DotPlot(so, features = features_of_interest, dot.scale = 8, cluster.idents=FALSE)

    df <- plot.data[[2]]
    genes <- unique(as.character(df$features.plot))
    samples <- unique(as.character(df$id))

    dist.df <- as.data.frame(matrix(nrow=length(genes), ncol=length(samples), dimnames=list(genes,samples)))
    dist.df[] <- 0

    for(s in samples){
        df.sub <- df %>% filter(id == s)
        for(g in genes){
            if(g %in% df.sub$features.plot){
                gene_idx <- which(as.character(df.sub$features.plot) == g)
                dist.df[g,s] <- df.sub[gene_idx,"avg.exp"]
            }
        }
    }
    tree_row = cluster_mat(dist.df, distance = clustering_distance, method = clustering_method)
    tree_row = identity2(tree_row, dist.df)
    features_of_interest <- features_of_interest[tree_row$order]

    imageWidth = 3000
    imageHeight = 6000*0.54
    dpi = 300

    png(
      filename="FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    plot.data <- .DotPlot(so, features = features_of_interest, dot.scale = 8, cluster.idents=TRUE)

    #do annotations
    # groups <- markers %>% group_by(Origin) %>% summarize(n()) %>% pull(`n()`)
    # annot_locations <- list()
    # i <- 1
    # for(n in seq_along(groups)){
    #     print(groups(n))
    #     annot_locations[[n]] <- c(i,(i+groups[n]))
    #     i <- i+groups[n]+1
    # }
    
    p <- plot.data[[1]] + coord_flip() + RotatedAxis()

    # #p <- p + ggplot2::annotate(geom = "text", x = -1, y = -1, label = "Helpful annotation", color = "red", angle = 90)
    # p <- p + annotation_raster(
    #       raster = annot_cols[rev(group.use2)],
    #       xmin = -Inf,
    #       xmax = Inf,
    #       ymin = 7.3,
    #       ymax = 8
    #     )

    # y.max <- length(group.use2)
    
    # pbuild <- ggplot_build(plot = p)
    # y <- markers %>% select(Matches,Origin)
    # row.names(y) <- gsub("Lamp3","Lamp3+",y$Matches)
    # y$Matches <- NULL
    # colnames(y) <- gsub("Origin","group",colnames(y))
    # y$group <- gsub("Lamp3","Lamp3+",y$group)
    # y$group <- gsub("Intestine","Intestinal CD103+CD11b+ DC",y$group)
    # y.max <- max(pbuild$layout$panel_params[[1]]$y.range)
    # # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot

    # y$y <- pbuild$layout$panel_params[[1]]$y.major %||% pbuild$layout$panel_params[[1]]$y$break_positions()
    # label.y.pos <- tapply(X = y$y, INDEX = y$group, FUN = median) * y.max
    # label.y.pos <- data.frame(group = names(x = label.y.pos), label.y.pos) 

    # x.range <- diff(x = pbuild$layout$panel_params[[1]]$x.range)
    # x.pos <- max(pbuild$layout$panel_params[[1]]$x.range) + x.range * 0.015
    # x.max <- x.pos + 0.02 * x.range
    # print(x.max * 0.03 * 0.5)
    # p <- p + geom_text(
    #             stat = "identity",
    #             data = label.y.pos,
    #             aes_string(label = 'group', x = 'label.y.pos'),
    #             y = x.max - 0.4,
    #             angle = 90,
    #             #hjust = hjust,
    #             size = 6,
    #             color = "#FFFFFF")

    # for(i in seq_along(annot_locations)){
    #     p <- p + ggplot2::annotate(geom = "text", y=7.5, x=((annot_locations[[i]][1]+annot_locations[[i]][2])/2)-1, label = levels(markers$Origin)[i], color = "black", angle = 90, size=5)
    #     #p <- p + geom_vline(xintercept=annot_locations[[i]][2]-(length(annot_locations)-i+0.5))
    # }
    print(p)
# auto removed:     return(NULL)

}

RotatedAxis <- function(...) {
  rotated.theme <- theme(
    # Rotate X axis text
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(rotated.theme)
}

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

.DotPlot <- function(
  object,
  assay = NULL,
  features,
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  id.levels <- levels(x = data.features$id)

  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )

  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
#   new.cluster.ids <- c("Monocytes", "CD103intCD11b+ DC", "cDC1","cDC2", "CCR7Hi", "Siglec-H DC")
#   data.plot$id <- factor(as.character(data.plot$id), levels = new.cluster.ids)
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot(font_size = 16)+
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    #plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    library(colorspace)
    #plot  <- plot + scale_color_continuous_diverging(palette = "Blue-Red 3") #scale_colour_gradientn(divergent_hcl(100, palette = "Blue-Red 3"))
    plot  <- plot + scale_colour_gradient2(  low = "blue",
                                                mid = "gray",
                                                high = "red", limits=c(-2.5,2.5)) #scale_colour_gradientn(divergent_hcl(100, palette = "Blue-Red 3"))

  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  return(list(plot,data.plot))
}

cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
        stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
        stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
        d = as.dist(1 - cor(t(mat)))
    }
    else{
        if(class(distance) == "dist"){
            d = distance
        }
        else{
            d = dist(mat, method = distance)
        }
    }
    
    return(hclust(d, method = method))
}

identity2 = function(x, ...){
    return(x)
}

print("template_function_FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq.R #########################################################################")
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
var_Canonical_DC_Markers<-readRDS(paste0(rds_output,"/var_Canonical_DC_Markers.rds"))
Input_is_Seurat_count <- 0
for(item in var_Canonical_DC_Markers){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Canonical_DC_Markers<-as.data.frame(var_Canonical_DC_Markers)}else{var_Canonical_DC_Markers <- var_Canonical_DC_Markers}
invisible(graphics.off())
var_FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq<-FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq(var_renamed_so,var_Canonical_DC_Markers)
invisible(graphics.off())
saveRDS(var_FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq, paste0(rds_output,"/var_FIGS5B_canonical_dc_bubbleplot_no_cdc1_citeseq.rds"))
