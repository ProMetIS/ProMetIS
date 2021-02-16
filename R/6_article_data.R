.cv_ggplot <- function(data.tb,
                       title.c = "") {
  
  color.c <- x.c <- "type"
  y.c <- "value"
  
  p <- ggplot2::ggplot(data.tb, ggplot2::aes(x = type, y = value, color = type)) +
    ggplot2::geom_boxplot(outlier.size = 1) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(title = title.c, x = "", y = "CV (%)") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 11, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 11, face = "bold"),
                   axis.text = ggplot2::element_text(size = 11, face = "bold"),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.2)
  
  return(invisible(p))
  
}


.metabo_quality_metrics <- function(metabo.mset,
                                    figure.c = c("none", "interactive")[2]) {
  
  .hist <- function(metric.vn, bin.vn) {
    
    stopifnot(identical(metric.vn, sort(metric.vn)))
    
    ind.vi <- numeric(length(metric.vn))
    
    for (k in 1:(length(bin.vn) - 1))
      ind.vi <- ind.vi + as.numeric(bin.vn[k] <= metric.vn)
    
    names(ind.vi) <- names(metric.vn)
    
    return(ind.vi)
    
  }
  
  # Zhang, X., Dong, J., & Raftery, D. (2020).
  # Five easy metrics of data quality for LC-MS based global metabolomics.
  # Analytical Chemistry, 92(19), 12925–12933. https://doi.org/10.1021/acs.analchem.0c01493
  
  quality.vc <- c("features",
                  "groups",
                  "NA_0.2",
                  "correl_inj_test",
                  "correl_inj_pca",
                  "pool_spread_pca",
                  "poolCV_0.3",
                  "ICCpool")
  quality.mn <- matrix(0,
                       nrow = length(quality.vc),
                       ncol = length(metabo.mset),
                       dimnames = list(quality.vc,
                                       names(metabo.mset)))
  
  cv_bin.i <- 20
  cv_bin.vn <- seq(0, 1, length.out = cv_bin.i + 1)
  
  cpd.mn <- matrix(0,
                   nrow = cv_bin.i,
                   ncol = length(metabo.mset),
                   dimnames = list(1:cv_bin.i, names(metabo.mset)))
  
  int_bin.i <- 20
  
  cpd.ls <- vector(mode = "list", length = length(metabo.mset))
  names(cpd.ls) <- names(metabo.mset)
  cpd_temp.mn <- matrix(0,
                        nrow = cv_bin.i,
                        ncol = int_bin.i,
                        dimnames = list(cv_bin.vn[-1], 1:int_bin.i))
  
  icc.ls <- vector(mode = "list", length = length(metabo.mset))
  names(icc.ls) <- names(metabo.mset)
  
  for (set.c in names(metabo.mset)) {
    # set.c <- "metabolomics_liver_hilic_neg"
    
    eset <- metabo.mset[[set.c]]
    eset <- eset[, Biobase::pData(eset)[, "sampleType"] %in% c("pool", "sample")]
    
    ## Total number of features
    
    quality.mn["features", set.c] <- dim(eset)["Features"]
    
    ## Features with <= 20% NA (%)
    
    isna.vn <- apply(Biobase::exprs(eset), 1, function(feat.vn) sum(is.na(feat.vn)))
    
    quality.mn["NA_0.2", set.c] <- sum(isna.vn / dim(eset)["Samples"] <= 0.2) / dim(eset)["Features"] * 100
    
    ## Number of chemical groups
    
    eset <- phenomis::reducing(eset)
    
    quality.mn["groups", set.c] <- length(table(Biobase::fData(eset)[, "redund_group"]))
    
    ## Features correlated with the injection order: univariate (%)
    
    eset <- phenomis::hypotesting(eset, "spearman", "injectionOrder",
                                  figure.c = "none",
                                  report.c = "none")
    
    quality.mn["correl_inj_test", set.c] <- sum(Biobase::fData(eset)[, "spearman_injectionOrder_signif"],
                                                na.rm = TRUE) / dim(eset)["Features"] * 100
    
    ## Features correlated with the injection order: multivariate (%)
    
    eset.pca <- ropls::opls(eset, predI = 3, fig.pdfC = "none", info.txtC = "none")

    scores.mn <- ropls::getScoreMN(eset.pca)

    inj_cor.vn <- drop(cor(Biobase::pData(eset)[, "injectionOrder"], scores.mn))

    quality.mn["correl_inj_pca", set.c] <- sqrt(sum(inj_cor.vn^2)) * 100

    ## Spread of QC (%)

    scores_invcov.mn <- solve(stats::cov(scores.mn))

    scores_dist.vn <- apply(scores.mn,
                            1,
                            function(x)
                              t(as.matrix(x)) %*% scores_invcov.mn %*% as.matrix(x))

    quality.mn["pool_spread_pca", set.c] <- max(scores_dist.vn[Biobase::pData(eset)[, "sampleType"] == "pool"]) / max(scores_dist.vn) * 100
    
    ## CV pool <= 30% (%)
    
    qc.eset <- eset[, Biobase::pData(eset)[, "sampleType"] == "pool"]
    
    # qc.mn: matrix of QC intensities
    
    qc.mn <- Biobase::exprs(qc.eset)
    
    cv.vn <- apply(qc.mn, 1, function(row.vn) sd(row.vn, na.rm = TRUE) / mean(row.vn, na.rm = TRUE))
    
    quality.mn["poolCV_0.3", set.c] <- signif(sum(cv.vn <= 0.3, na.rm = TRUE) / length(cv.vn) * 100, 2)
    
    cv.vn <- cv.vn[order(cv.vn)]
    
    cv_ind.vi <- .hist(cv.vn, cv_bin.vn)
    
    # cpd.mn: cumulative percentage of compounds as a function of QC CV
    
    cpd.mn[names(table(cv_ind.vi)), set.c] <- table(cv_ind.vi)
    
    cpd.mn[, set.c] <- cumsum(cpd.mn[, set.c]) / length(cv_ind.vi) * 100
    
    # qcl: log10 transformation of QC intensities
    
    qcl.mn <- qc.mn
    qcl.mn[qcl.mn < .Machine$double.eps] <- NA
    qcl.mn <- log10(qcl.mn)
    
    int.vn <- apply(qcl.mn, 1, function(feat.vn) median(feat.vn, na.rm = TRUE))
    int.vn <- int.vn[order(int.vn)]
    
    int_bin.vn <- seq(min(int.vn), max(int.vn), length.out = int_bin.i + 1)
    
    cpd_set.mn <- cpd_temp.mn
    
    int_ind.vi <- .hist(int.vn, int_bin.vn)
    
    colnames(cpd_set.mn) <- int_bin.vn[-1]
    # colnames(cpd_set.mn) <- as.numeric(colnames(cpd_set.mn)) - as.numeric(colnames(cpd_set.mn)[which(table(int_ind.vi) == max(table(int_ind.vi)))])
    
    for (i in 1:nrow(cpd_set.mn)) {
      for (j in 1:ncol(cpd_set.mn)) {
        cpd_set.mn[i, j] <- length(intersect(names(cv_ind.vi)[cv_ind.vi == i],
                                             names(int_ind.vi)[int_ind.vi == j]))
      }
    }
    
    cpd_set.mn <- cpd_set.mn / sum(cpd_set.mn) * 100
    cpd_set.mn <- t(cpd_set.mn)
    cpd_set.mn <- cpd_set.mn[, ncol(cpd_set.mn):1]
    cpd_set.mn <- cpd_set.mn[nrow(cpd_set.mn):1, ]
    
    cpd.ls[[set.c]] <- cpd_set.mn
    
    
    ## ICC at most probable abundance
    
    int.vn <- apply(qcl.mn, 1, function(feat.vn) median(feat.vn, na.rm = TRUE))
    qcl.mn <- qcl.mn[order(int.vn), ]
    int.vn <- int.vn[order(int.vn)]
    
    int_bin.vn <- seq(min(int.vn), max(int.vn), length.out = int_bin.i + 1)
    
    int_ind.vi <- .hist(int.vn, int_bin.vn)
    
    icc_set.vn <- sapply(seq_len(int_bin.i), function(k) {
      irr::icc(qcl.mn[int_ind.vi <= k, , drop = FALSE])[["value"]]
    })
    names(icc_set.vn) <- signif(int_bin.vn[-1], 2)
    # signif(caTools::runmean(int_bin.vn, 2)[-1], 2)
    
    int_bin.tab <- table(int_ind.vi)
    names(icc_set.vn) <- as.numeric(names(icc_set.vn)) - as.numeric(names(icc_set.vn))[which(int_bin.tab == max(int_bin.tab))[1]]
    
    icc.ls[[set.c]] <- icc_set.vn
    
  }
  rownames(cpd.mn) <- cv_bin.vn[-1]
  
  
  quality.mn["ICCpool", ] <- sapply(icc.ls,
                                    function(icc_set.vn)
                                      as.numeric(icc_set.vn[which.min(abs(as.numeric(names(icc_set.vn))))]),
                                    USE.NAMES = FALSE) * 100
  
  names(icc.ls) <- names(cpd.ls) <- colnames(cpd.mn) <- colnames(quality.mn) <- gsub("metabolomics_", "",
                                                                                     colnames(quality.mn))
  

  return(list(quality.mn = quality.mn,
              cv_bin.vn = cv_bin.vn,
              cpd.mn = cpd.mn,
              cpd.ls = cpd.ls,
              icc.ls = icc.ls))
  
}


.metabo_quality_plot <- function(quality_metrics.ls,
                                 types.vc = c("all",
                                              "cv_percent",
                                              "cv_percent_intens",
                                              "icc_intens"),
                                 figure.c = "interactive") {
  
  quality.mn <- quality_metrics.ls[["quality.mn"]]
  cv_bin.vn <- quality_metrics.ls[["cv_bin.vn"]]
  cpd.mn <- quality_metrics.ls[["cpd.mn"]]
  icc.ls <- quality_metrics.ls[["icc.ls"]]
  
  print(colnames(cpd.mn))
  
  all_types.vc <- c("cv_percent",
                    "cv_percent_intens",
                    "icc_intens")
  
  if (types.vc == "all")
    types.vc <- all_types.vc
  
  stopifnot(all(types.vc %in% all_types.vc))
  
  if (figure.c != "interactive")
    pdf(figure.c)
  
  if ("cv_percent" %in% types.vc) {
    
    # Fig. 4a: CVs of the missing value-free compounds versus the accumulated percentage of compounds
    
    par(font = 2, font.axis = 2, font.lab = 2, lwd = 2)
    
    pal.vc <- RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7)]
    
    plot(c(0, 1), c(0, 100), type = "n",
         xlab = "Coefficient of variation (CV)",
         ylab = "Percentage of compounds (%)",
         las = 1,
         main = "Cumulative percent of compounds as a function of QC CVs",
         xaxs = "i",
         yaxs = "i")
    
    for (set.c in colnames(quality.mn)) {
      
      lines(cv_bin.vn, c(0, cpd.mn[, set.c]),
            col = pal.vc[which(colnames(quality.mn) == set.c)],
            lwd = 2)
      
      points(cv_bin.vn, c(0, cpd.mn[, set.c]),
             col = pal.vc[which(colnames(quality.mn) == set.c)],
             pch = 16)
      
    }
    
    abline(v = 0.3, col = "red", lty = "dashed")
    
    legend("bottomright", legend = paste0(colnames(quality.mn),
                                          " (", signif(cpd.mn[which(abs(cv_bin.vn[-1] - 0.3) < 0.04), ], 2),
                                          "%)"),
           bty = "n", text.col = pal.vc, text.font = 2)
    
  }
  
  if ("cv_percent_intens" %in% types.vc) {
    
    # Fig. 4b: CVs of the missing value-free compounds versus the normalized log10 abundance and percentage of compunds
    
    for (set.c in colnames(quality.mn)) {
      cpd_set.mn <- quality_metrics.ls[["cpd.ls"]][[set.c]]
      
      # int_prop_max.n <- as.numeric(rownames(which(cpd_set.mn == max(cpd_set.mn, na.rm = TRUE), arr.ind = TRUE))[1])
      # rownames(cpd_set.mn) <- as.numeric(rownames(cpd_set.mn)) - int_prop_max.n
      
      cpd_set.mn[cpd_set.mn == 0] <- NA
      
      ropls::view(cpd_set.mn, mainC = set.c,
                  rowLabC = "Log10(Intensity)",
                  rowAllL = TRUE,
                  rowMarN = 3.1,
                  colLabC = "Coefficient of variation",
                  colAllL = TRUE)
    }
    
  }
  
  if ("icc_intens" %in% types.vc) {
    
    # Fig. 5a: ICC values versus percentages of compounds sorted by abundance
    
    par(font = 2, font.axis = 2, font.lab = 2, lwd = 2)
    
    pal.vc <- RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7)]
    
    plot(range(as.numeric(c(sapply(icc.ls, names)))), c(0, 1), type = "n",
         xlab = "Log10(Intensity)",
         ylab = "ICC",
         las = 1,
         main = "Intra correlation coefficient in QCs as a function of cumulative intensity",
         xaxs = "i",
         yaxs = "i")
    
    for (set.c in colnames(quality.mn)) {
      
      lines(as.numeric(names(icc.ls[[set.c]])),
            icc.ls[[set.c]],
            col = pal.vc[which(colnames(quality.mn) == set.c)],
            lwd = 2)
      
      points(as.numeric(names(icc.ls[[set.c]])),
             icc.ls[[set.c]],
             col = pal.vc[which(colnames(quality.mn) == set.c)],
             pch = 16)
      
    }
    
    abline(v = 0, col = "red", lty = "dashed")
    
    legend("bottomright", paste0(colnames(quality.mn),
                                 " (", round(quality.mn["ICCpool", ]), "%)"),
           bty = "n", text.col = pal.vc, text.font = 2)
    
  }
  
  if (figure.c != "interactive")
    dev.off()
  
}


.metabo_drift_plot <- function(eset,
                               span.n = 1,
                               sample_intensity.c = "mean",
                               title.c = "",
                               boxplot.l = FALSE,
                               col_batch.c = "batch",
                               col_injectionOrder.c = "injectionOrder",
                               col_sampleType.c = "sampleType",
                               mar.vn = c(3.5, 3.6, 1.1, 0.6)) {
  
  sample_color.vc <- phenomis:::.sample_color_eset(eset = eset,
                                        col_sampleType.c = col_sampleType.c)
  
  graphics::par(mar = mar.vn)
  
  ## ordering
  
  texprs.mn <- t(Biobase::exprs(eset))
  pdata.df <- Biobase::pData(eset)
  
  pdata.df[, "ordIniVi"] <- 1:nrow(texprs.mn)
  
  if (col_injectionOrder.c %in% colnames(pdata.df)) {
    ordNamC <- "Injection Order"
    if (col_batch.c %in% colnames(pdata.df)) {
      ordVi <- order(pdata.df[, col_batch.c],
                     pdata.df[, col_injectionOrder.c])
    } else
      ordVi <- order(pdata.df[, col_injectionOrder.c])
  } else {
    ordNamC <- "Samples"
    ordVi <- 1:nrow(texprs.mn)
  }
  
  texprs.mn <- texprs.mn[ordVi, ]
  pdata.df <- pdata.df[ordVi, ]
  sample_color_ordered.vc <- sample_color.vc[ordVi]
  
  if (col_batch.c %in% colnames(pdata.df))
    batch.table <- table(pdata.df[, col_batch.c])
  
  sample_means.vn <- eval(parse(text = paste0("apply(texprs.mn, 1, function(obsVn) ",
                                              sample_intensity.c, "(obsVn, na.rm = TRUE))")))
  
  graphics::plot(sample_means.vn,
                 col = sample_color_ordered.vc,
                 pch = 18,
                 main = title.c,
                 type = "n",
                 xaxs = "i",
                 xlab = "",
                 ylab = "")
  
  if (boxplot.l) {
    
    for (samI in 1:nrow(texprs.mn))
      graphics::boxplot(texprs.mn[samI, ],
                        at = samI,
                        add = TRUE)
    
  }
  
  graphics::points(sample_means.vn,
                   col = sample_color_ordered.vc,
                   pch = 18)
  
  graphics::mtext(ordNamC,
                  cex = 0.8,
                  line = 2,
                  side = 1)
  
  graphics::mtext(paste0(toupper(substr(sample_intensity.c, 1, 1)), substr(sample_intensity.c, 2, nchar(sample_intensity.c)), " of variable intensities"),
                  cex = 0.8,
                  line = 2,
                  side = 2)
  
  batch_sample.vi <- intersect(1:nrow(pdata.df),
                               grep("sample", pdata.df[, col_sampleType.c]))

  # graphics::lines(batch_sample.vi,
  #                 phenomis:::.loess(sample_means.vn, batch_sample.vi, batch_sample.vi, span.n),
  #                 col = phenomis:::.sample_color_vector("sample"),
  #                 lwd = 2)
  
  batch_pool.vi <- intersect(1:nrow(pdata.df),
                             grep("^pool$", pdata.df[, col_sampleType.c]))
  
  graphics::lines(batch_sample.vi,
                  phenomis:::.loess(sample_means.vn, batch_pool.vi, batch_sample.vi, span.n),
                  col = phenomis:::.sample_color_vector("pool"),
                  lwd = 2)
  
  # legend
  
  if (col_sampleType.c %in% colnames(pdata.df)) {
    obsColVuc <- sample_color_ordered.vc[sort(unique(names(sample_color_ordered.vc)))]
    legOrdVc <- c("blank", paste0("pool", 8:1), "pool", "other", "sample")
    obsColVuc <- obsColVuc[legOrdVc[legOrdVc %in% names(obsColVuc)]]
    
    graphics::text(rep(graphics::par("usr")[2], times = length(obsColVuc)),
                   graphics::par("usr")[3] + (0.05 + 1:length(obsColVuc) * 0.05) * diff(graphics::par("usr")[3:4]),
                   adj = 1,
                   col = obsColVuc,
                   font = 2,
                   labels = names(obsColVuc),
                   pos = 2)
  }
  
}


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
