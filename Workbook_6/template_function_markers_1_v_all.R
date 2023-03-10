# Enhanced Volcano (0c91aa57-0f76-4513-a063-5f9263d65727): v29
markers_1_v_all <- function(deg_cluster_vs_all) {
    # image: png

    # Changelog
    # 2020-10-29 Add support for pval == 0

    suppressMessages(library(stringr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggrepel))
    
    #Class handling
    df <- deg_cluster_vs_all
    if(class(df) == "SparkDataFrame"){
        df <- dplyr::collect(df)
    }
  
    # Log fold-change column
    # p-value column
    # Position of p-value threshold line #default 0.05
    # Position of log fold-change threshold lines #default 1 lfc
    # Label column #default gene names
    # Specific genes to label
    # Number of genes to label
    
    #label topN genes or specific genes
    
    label.col <- "Gene"
    sig.col <- "p_val_adj_1_vs_all"
    lfc.col <- "avg_logFC_1_vs_all"
    columns_of_interest <- c(label.col,lfc.col,sig.col)

    df <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs

    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df) <- c(label.col,change_lfc_name,change_sig_name)

    cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        df %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df
    } else if (value_to_sort_the_output_dataset=="p-value") {
        df %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df
    }
    genes_to_label <- as.character(df[1:no_genes_to_label,label.col])

    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df[,label.col]
    additional_labels <- additional_labels[filter]
    missing_labels <- additional_labels[!filter]

    if(length(missing_labels) > 0){
        cat("Could not find:\n")
        print(missing_labels)
    }
    use_only_addition_labels <- FALSE
    if(use_only_addition_labels){
        genes_to_label <- additional_labels
    }else{
        genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    # Plot variables
    title <- "Volcano Plot"
    subtitle <- "Contrast 1 vs all"

    pCutoff  = 10e-6
    FCcutoff = 0.69

    significant=as.vector(table(abs( df[,change_lfc_name] ) > FCcutoff &
                                     df[,change_sig_name]   < pCutoff))[2]

    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df))
    maxy <-  max(-log10(df[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
        maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
        # Sometimes, pvalues == 0
        keep <- df[[change_sig_name]] > 0
        df[[change_sig_name]][!keep] <- min(df[[change_sig_name]][keep])
        shapeCustom[!keep] <- 17

        maxy <- -log10(min(df[[change_sig_name]][keep]))
        cat("Some p-values equal zero. Adjusting y-limits.\n")
        cat(paste0("Maxy adjusted: ",maxy,"\n"))

    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[change_sig_name]]) <= maxy
    df[[change_sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom)<- rep("Exact",length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
    
    #Remove if nothin' doin'
    if(all(shapeCustom == 19)){
        shapeCustom <- NULL
    }
    
    maxy <- ceiling(maxy)

    use_custom_xlab <- FALSE
    if(use_custom_xlab){
        xlab <- gsub("_"," ",change_lfc_name)
    }else{
        xlab <- bquote(~Log[2]~ "fold change")
    }

    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    labSize <- 4
    pointSize <- 2

    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

    png(
      filename="markers_1_v_all.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p <- EnhancedVolcano(   df,x=change_lfc_name,y=change_sig_name,
                            lab=df[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df[,change_lfc_name])) - xlim_additional,ceiling(max(df[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
                            )
    print(p)

    df$rank <- -log10(df[[change_sig_name]]) * sign(df[[change_lfc_name]]) #sig already -log10
    df <- df[order(df$rank, decreasing=TRUE), ]
    #df <- df %>% arrange(desc(rank))
    return(df)
}

# This is directly copied from https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html v1.6.0
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1,
    max(toptable[[x]], na.rm=TRUE) + 1),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = 'EnhancedVolcano',
  caption = paste0('Total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 3.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  .legend = c('NS','Log2 FC','P','P & Log2 FC'),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 4.0,
  shade = NULL,
  shadeLabel = NULL,
  shadeAlpha = 1/2,
  shadeFill = "grey",
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black')
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))

  if(!is.numeric(toptable[[y]])) {
    stop(paste(y, ' is not numeric!', sep=''))

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'
  toptable$Sig[(toptable[[y]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[y]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values=colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(Sig)),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes=list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes=list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = factor(Sig),
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values=c(
            NS=col[1],
            FC=col[2],
            P=col[3],
            FC_P=col[4]),
          labels=c(
            NS=legendLabels[1],
            FC=legendLabels[2],
            P=legendLabels[3],
            FC_P=legendLabels[4])) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)
    } else {
        plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes=list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = factor(Sig),
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)
    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes=list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = factor(Sig)),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]))
    } else {
      plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }
  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = factor(Sig),
            shape = factor(Sig)),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4])) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE)
    } else {
      plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = factor(Sig)),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE)
    }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")

  # Gridlines
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())

  # user has specified to draw with geom_text or geom_label?
  if (boxedLabels == FALSE) {
    # For labeling with geom_text/label_repel (connectors) and
    # geom_text/label (.., check_overlap = TRUE), 4 possible
    # scenarios can arise
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    } else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label=subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
     } else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label=subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.size = 0,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.size = 0,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }
  } else {
    # For labeling with geom_text/label_repel (connectors) and
    # geom_text/label (.., check_overlap = TRUE), 4 possible
    # scenarios can arise
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(
        data=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    } else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(
        data=subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label=subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(
        data=subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label=subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(
        data=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }

  # shading
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE) +

      scale_fill_identity(name = shadeLabel,
        labels = shadeLabel,
        guide = 'legend')

  return(plot)
}

#################################################
## Global imports and functions included below ##
#################################################

print("template_function_markers_1_v_all.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_deg_cluster_vs_all<-readRDS(paste0(rds_output,"/var_deg_cluster_vs_all.rds"))
Input_is_Seurat_count <- 0
for(item in var_deg_cluster_vs_all){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_deg_cluster_vs_all<-as.data.frame(var_deg_cluster_vs_all)}else{var_deg_cluster_vs_all <- var_deg_cluster_vs_all}
invisible(graphics.off())
var_markers_1_v_all<-markers_1_v_all(var_deg_cluster_vs_all)
invisible(graphics.off())
saveRDS(var_markers_1_v_all, paste0(rds_output,"/var_markers_1_v_all.rds"))
