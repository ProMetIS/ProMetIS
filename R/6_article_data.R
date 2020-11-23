# Score plot visualization for PCA and (O)PLS(-DA) models
# 
# @param ropls.model opls: object obtained with the 'ropls::opls' applied on an ExpressionSet
# @param label.c character(1): name of the pData column to be used for the labels
# @param color.c character(1): name of the pData column to be used for the colors
# @param info.vc character(): names of the pData columns to be used for the plotly info
# @param title.c character(1): plot title
# @param components.vi integer(2): number of the components to display as x and y axis
# @param palette.c character(1): name of the RColorBrewer palette (for qualitative factor)
# @param ellipse.l logical(1): should ellipses be drawn (for qualitative factor)
# @param plot.l logical(1): set to FALSE in the interactive(_plotly) mode to return
# the plot object only without displaying it (default: TRUE)
# @param size.ls list: sizes for axis labels (default: 16), axis text (default: 14),
# points (default: 3), labels (default = 5), title (default = 20), legend title (default: 15),
# legend text (default: 15)
# @param figure.c Character: either 'interactive' (respectively,
# 'interactive_plotly') for interactive display with ggplot2 (respectively,
# with plotly::ggplotly [default]), or 'my_scoreplot.pdf' (respectively
# 'my_scoreplot.html') for figure saving (only the extension matters) with
# ggplot2 (respectively, with plotly::ggplotly)
# @return invisible ggplot2 object
# @examples
# # loading the 'sacurine' dataset from the 'ropls' package
# data(sacurine, package = "ropls")
# # building the ExpressionSet object
# sacurine.eset <- Biobase::ExpressionSet(assayData = t(sacurine[["dataMatrix"]]),
#                                         phenoData = new("AnnotatedDataFrame",
#                                                         data = sacurine[["sampleMetadata"]]),
#                                         featureData = new("AnnotatedDataFrame",
#                                                           data = sacurine[["variableMetadata"]]),
#                                         experimentData = new("MIAME",
#                                                              title = "sacurine"))
# # computing the PCA
# sacurine.pca <- ropls::opls(sacurine.eset)
# # score plot
# score_plotly(sacurine.pca)
# score_plotly(sacurine.pca, color.c = "age")
# score_plotly(sacurine.pca, color.c = "gender")
score_plotly <- function(ropls.model,
                         label.c = "sampleNames",
                         color.c = "",
                         info.vc = "all",
                         title.c = "",
                         components.vi = c(1, 2),
                         palette.c = "Set1",
                         ellipse.l = TRUE,
                         plot.l = TRUE,
                         size.ls = list(axis_lab.i = 16,
                                        axis_text.i = 14,
                                        point.i = 3,
                                        label.i = 5,
                                        title.i = 20,
                                        legend_title.i = 15,
                                        legend_text.i = 15),
                         figure.c = c("interactive",
                                      "interactive_plotly",
                                      "my_scoreplot.pdf",
                                      "my_scoreplot.html")[2]) {
  
  
  # checking arguments and preparing data
  
  ## dataset [ExpressionSet]
  
  eset <- ropls::getEset(ropls.model)
  
  if (any(dim(eset) < 1))
    stop("ExpressionSet object could not be extracted from the 'ropls.model'")
  
  pdata.df <- Biobase::pData(eset)
  
  data.df <- cbind.data.frame(sampleNames = Biobase::sampleNames(eset),
                              pdata.df,
                              stringsAsFactors = FALSE)
  
  ## plotly info
  
  if (length(info.vc) == 1 && info.vc == "all") {
    text.df <- data.df
  } else {
    info_in_metadata.vl <- info.vc %in% colnames(data.df)
    if (any(!info_in_metadata.vl)) {
      if (all(!info_in_metadata.vl)) {
        stop("None of the selected info.vc names was found in the sampleMetadata:\n",
             paste(info.vc, collapse = ", "))
      } else {
        warning("The following columns were not found in the sampleMetadata:\n",
                paste(info.vc[!info_in_metadata.vl], collapse = ", "))
      }
    }
    text.df <- data.df[, info.vc[info_in_metadata.vl], drop = FALSE]
  }
  
  text.vc <- apply(text.df, 1,
                   function(row.vc)
                     paste(paste0(colnames(text.df), " = ", row.vc), collapse = "\n"))
  
  ## score vectors
  
  score.mn <- ropls::getScoreMN(ropls.model)
  data.df <- cbind.data.frame(data.df,
                              text = text.vc,
                              score.mn)
  
  ## labels
  
  stopifnot(length(label.c) == 1)
  
  if (label.c != "")
    stopifnot(label.c %in% colnames(data.df))
  
  ## colors
  
  stopifnot(length(color.c) == 1)
  
  if (color.c != "") {
    stopifnot(color.c %in% colnames(data.df))
    if (is.factor(data.df[, color.c]) || is.character(data.df[, color.c])) {
      color_type.c <- "qualitative"
    } else
      color_type.c <- "quantitative"
  }
  
  ## components
  
  stopifnot(length(components.vi) == 2)
  stopifnot(max(components.vi) <= ncol(score.mn))
  
  ## sizes
  
  size_default.vi <- c(axis_lab.i = 16,
                       axis_text.i = 14,
                       point.i = 3,
                       label.i = 5,
                       title.i = 20,                      
                       legend_title.i = 15,
                       legend_text.i = 15)
  
  for (size.c in names(size_default.vi)) {
    if (!(size.c %in% names(size.ls)))
      size.ls[[size.c]] <- size_default.vi[size.c]
  }
  
  ## plot
  
  if (!plot.l && !grepl("interactive", figure.c))
    stop("'plot.l' can be set to 'FALSE' only for the interactive(_plotly) figures.")
  
  ## filename extention
  
  filename_ext.c <-  utils::tail(unlist(strsplit(basename(figure.c), ".", fixed = TRUE)), 1)
  
  
  # starting the plot [ggplot]
  
  p <- eval(parse(text = paste0("ggplot2::ggplot(data.df, ggplot2::aes(x = ",
                                paste0('p', components.vi[1]),
                                ", y = ",
                                paste0('p', components.vi[2]),
                                ifelse(color.c != "",
                                       paste0(", color = ", color.c),
                                       ""),
                                ifelse(label.c != "",
                                       paste0(", label = ", label.c),
                                       ""),
                                ", text = text))")))
  
  ## text/points [geom_text/geom_point]
  
  if (label.c != "") {
    p <- p + ggplot2::geom_text(size = size.ls[["label.i"]], ggplot2::aes(fontface = "bold"))
  } else {
    p <- p + ggplot2::geom_point(size = size.ls[["point.i"]])
  }
  
  ## horizontal and vertical lines [geom_hline, geom_vline]
  
  p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0))
  
  ## ellipses [stat_ellipse]
  
  p <- p + eval(parse(text = paste0("ggplot2::stat_ellipse(ggplot2::aes(x = ",
                                    paste0('p', components.vi[1]),
                                    ", y = ",
                                    paste0('p', components.vi[2]),
                                    ", group = 1), type = 'norm')")))
  
  if (ellipse.l && color.c != "" && color_type.c == "qualitative")
    p <- p + eval(parse(text = paste0("ggplot2::stat_ellipse(ggplot2::aes(x = ",
                                      paste0('p', components.vi[1]),
                                      ", y = ",
                                      paste0('p', components.vi[2]),
                                      ", group = ",
                                      ifelse(color.c != "", color.c, 1),
                                      "), type = 'norm')")))
  
  # title and axis labels [labs]
  
  p <- p + ggplot2::labs(title = title.c,
                         x = paste0("t", components.vi[1],
                                    " (",
                                    round(ropls.model@modelDF[components.vi[1], "R2X"] * 100),
                                    "%)"),
                         y = paste0("t", components.vi[2],
                                    " (",
                                    round(ropls.model@modelDF[components.vi[2], "R2X"] * 100),
                                    "%)"))
  
  # theme [them_bw, theme]
  
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = size.ls[["title.i"]], face = "bold"),
                   axis.title.x = ggplot2::element_text(size = size.ls[["axis_lab.i"]], face = "bold"),
                   axis.title.y = ggplot2::element_text(size = size.ls[["axis_lab.i"]], face = "bold"),
                   axis.text = ggplot2::element_text(size = size.ls[["axis_text.i"]]),
                   legend.title = ggplot2::element_text(face = "bold", size = size.ls[["legend_title.i"]]),
                   legend.text = ggplot2::element_text(face = "bold", size = size.ls[["legend_text.i"]]))
  
  # palette [scale_colour_brewer, scale_colour_gradientn]
  
  if (color.c != "") {
    if (color_type.c == "qualitative") {
      if (palette.c != "")
        p <- p + ggplot2::scale_colour_brewer(palette = palette.c)
    } else
      p <- p + ggplot2::scale_colour_gradientn(colours = rev(rainbow(100, end = 4/6)))
  }
  
  # display/saving [plotly::ggplotly, plotly::layout, htmlwidgets::saveWidget, plotly::as_widget]
  
  if (figure.c == "interactive_plotly" || filename_ext.c == "html") {
    
    p <- plotly::ggplotly(p, tooltip = "text")
    
    p <- plotly::layout(p, hoverlabel = list(font = list(size = 20)))
    
    if (filename_ext.c == "html") {
      
      htmlwidgets::saveWidget(plotly::as_widget(p), figure.c)
      
      return(invisible(p))
      
      
    } else {
      
      if (plot.l)
        print(p)
      
      return(invisible(p))
      
    }
    
  } else if (figure.c == "interactive" || filename_ext.c == "pdf") {
    
    if (filename_ext.c == "pdf")
      grDevices::pdf(figure.c)

    if (plot.l)
      print(p)
    
    if (filename_ext.c == "pdf")
      grDevices::dev.off()
    
    return(invisible(p))
    
  }
  
}

.metabo_qcmetrics <- function(metabo.mset) {
  
  metrics.vc <- c("cor_univ", "cor_pca", "qc_spread")
  metrics.mn <- matrix(0, nrow = length(metabo.mset), ncol = length(metrics.vc),
                       dimnames = list(names(metabo.mset), metrics.vc))
  
  # testing if 'injectionOrder' column is present with no NA
  
  injOrd_test.vl <- sapply(names(metabo.mset),
                           function(set.c) {
                             pdata.df <- Biobase::pData(metabo.mset[[set.c]])
                             "injectionOrder" %in% colnames(pdata.df) &&
                               !any(is.na(pdata.df[, "injectionOrder"]))
                           })
  if (any(!injOrd_test.vl))
    stop("The 'injectionOrder' column from the sampleMetadata is either missing or contains missing values in the following dataset(s):\n",
        paste(names(injOrd_test.vl)[!injOrd_test.vl], collapse = ",\n"), sep = "")
  
  # significant correlation with injection order
  
  metabo.mset <- phenomis::hypotesting(metabo.mset, "spearman", "injectionOrder",
                                       figure.c = "none",
                                       report.c = "none")
  
  for (set.c in names(metabo.mset)) {
    
    # set.c <- "metabolomics_plasma_c18acqui_pos"
    # set.c <- "metabolomics_liver_hilic_neg"
 
    eset <- metabo.mset[[set.c]]
    
    # eset <- eset[, order(Biobase::pData(eset)[, "injectionOrder"])]
    
    eset.pca <- ropls::opls(eset, predI = 2, fig.pdfC = "none", info.txtC = "none")
    
    # significant correlation with injection order

    metrics.mn[set.c, "cor_univ"] <- sum(Biobase::fData(eset)[, "spearman_injectionOrder_signif"],
                                         na.rm = T) / dim(eset)["Features"]
    
    # correlation between scores and injection order
    
    scores.mn <- ropls::getScoreMN(eset.pca)
    
    metrics.mn[set.c, "cor_pca"] <- sqrt(tcrossprod(cor(Biobase::pData(eset)[, "injectionOrder"], scores.mn)))
    
    # QC spread
    
    scores_invcov.mn <- solve(stats::cov(scores.mn))
    
    scores_dist.vn <- apply(scores.mn,
                            1,
                            function(x)
                              t(as.matrix(x)) %*% scores_invcov.mn %*% as.matrix(x))
    
    stopifnot(identical(names(scores_dist.vn), Biobase::sampleNames(eset)))
    
    metrics.mn[set.c, "qc_spread"] <- max(scores_dist.vn[Biobase::pData(eset)[, "sampleType"] == "pool"]) / max(scores_dist.vn)
    
  }
  
  # return

  return(metrics.mn)
  
}

quality <- function(metabo.mset, figure.c = c("none", "interactive")[2]) {
  
  quality.vc <- c("compounds",
                  "correl_inj_pca",
                  "pool_spread_pca")
  quality.mn <- matrix(0,
                       nrow = length(metabo.mset),
                       ncol = length(quality.vc),
                       dimnames = list(names(metabo.mset),
                                       quality.vc))
  
  for (set.c in names(metabo.mset)) {
    # set.c <- "metabolomics_liver_hilic_neg"
    
    eset <- metabo.mset[[set.c]]
    
    ## compounds
    
    quality.mn[set.c, "compounds"] <- dim(eset)["Features"]
    
    ## correlation of the 3 first PCA components with the injection order (%)
    
    eset.pca <- ropls::opls(eset, predI = 3, fig.pdfC = "none", info.txtC = "none")
    
    scores.mn <- ropls::getScoreMN(eset.pca)
    
    inj_cor.vn <- drop(cor(Biobase::pData(eset)[, "injectionOrder"], scores.mn))
    
    quality.mn[set.c, "correl_inj_pca"] <- sqrt(sum(inj_cor.vn^2)) * 100
    
    ## ratio between the maximum QC spread and the maximum sample spread in the 3 first PCA component space
    ## assessed by the Mahalanobis distance (%)
    
    scores_invcov.mn <- solve(stats::cov(scores.mn))
    
    scores_dist.vn <- apply(scores.mn,
                            1,
                            function(x)
                              t(as.matrix(x)) %*% scores_invcov.mn %*% as.matrix(x))
    
    quality.mn[set.c, "pool_spread_pca"] <- max(scores_dist.vn[Biobase::pData(eset)[, "sampleType"] == "pool"]) / max(scores_dist.vn) * 100
    

  }
 
 return(quality.mn)
  
}
