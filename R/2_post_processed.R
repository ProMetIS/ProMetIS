## Proteomics ----

.format_imputation <- function(eset) {
  
  prot_pda.df <- Biobase::pData(eset)
  prot_fda.df <- Biobase::fData(eset)
  
  ## checking that the sample names are ordered by increasing ID
  prot_samp.vi <- as.integer(substr(Biobase::sampleNames(eset), 2, 4))
  stopifnot(identical(prot_samp.vi, sort(prot_samp.vi)))
  
  ## getting imputation info
  value_origin.df <- prot_fda.df[, grep("OriginOfValueabundance",
                                        colnames(prot_fda.df), value = TRUE)]
  colnames(value_origin.df) <- gsub("_run90methode30K",
                                    "",
                                    gsub("_mgf", "",
                                         gsub("supp_OriginOfValueabundance_", "",
                                              colnames(value_origin.df))))
  
  ## re-ordering imputation info to match sample names
  if (Biobase::experimentData(eset)@title == "proteomics_liver") {
    file_to_sample.vc <- prot_pda.df[, "supp_sample name"]
    names(file_to_sample.vc) <- gsub("abundance_", "",
                                     gsub(".mgf", "",
                                          prot_pda.df[, "Sample.name"], fixed = TRUE))
    colsel.vl <- colnames(value_origin.df) %in% names(file_to_sample.vc)
    value_origin.df <- value_origin.df[, colsel.vl]
    colnames(value_origin.df) <- file_to_sample.vc[colnames(value_origin.df)]
  }
  
  value_samp.vi <- as.integer(colnames(value_origin.df), 1, 3)
  value_origin.df <- value_origin.df[, order(value_samp.vi)]
  
  stopifnot(identical(colnames(value_origin.df), as.character(prot_samp.vi)))
  colnames(value_origin.df) <- Biobase::sampleNames(eset)
  
  imputed.mi <- apply(value_origin.df, 2, DAPAR_is.MV)
  mode(imputed.mi) <- "integer"
  
  stopifnot(!any(is.na(c(imputed.mi))))
  
  colnames(imputed.mi) <- paste0("imputed_", colnames(imputed.mi))
  
  Biobase::fData(eset) <- cbind.data.frame(prot_fda.df,
                                           imputed.mi)
  
  eset
  
}


## Metabolomics ----

.filter_pcgroup <- function(input.eset,
                            cor_method.c = "pearson",
                            cor_threshold.n = 0.9) {
  
  # Authors: Antoine Gravot antoine.gravot@univ-rennes1.fr (Protocole conception)
  # and Misharl Monsoor misharl.monsoor@sb-roscoff.fr (for galaxy wrapper and R script).
  # https://github.com/workflow4metabolomics/correlation_analysis/blob/master/correlation_analysis.r
  # Adapted to ExpressionSet objects by Etienne Thevenot [2020-02-22]
  
  stopifnot("pcgroup" %in% Biobase::fvarLabels(input.eset))
  
  eset <- input.eset[stats::complete.cases(Biobase::fData(input.eset)[, "pcgroup"]), ]
  
  exprs.mn <- Biobase::exprs(eset)
  fdata.df <- Biobase::fData(eset)
  
  fdata.df[, "mean_intensity"] <- rowMeans(exprs.mn, na.rm = TRUE)
  
  order.vi <- order(fdata.df[, "pcgroup"], -fdata.df[, "mean_intensity"])
  
  exprs.mn <- exprs.mn[order.vi, ]
  fdata.df <- fdata.df[order.vi, ]
  
  exprs_cor.mn <- cor(t(exprs.mn), method = cor_method.c)
  
  pcgroups.vi <- unique(fdata.df[, "pcgroup"])
  
  redundancy.vl <- logical()
  
  for (pcgroup.i in pcgroups.vi) {
    
    pcgroup.vl <- fdata.df[, "pcgroup"] == pcgroup.i
    
    if (sum(pcgroup.vl) == 1) {
      
      redundancy.vl <- c(redundancy.vl, FALSE)
      
    } else {
      
      group_cor.mn <- exprs_cor.mn[pcgroup.vl, pcgroup.vl, drop = FALSE]
      group_cor.mn[!lower.tri(group_cor.mn)] <- NA_real_
      group_cor.vl <- apply(group_cor.mn, 1,
                            function(row.vn) {
                              row.vn <- row.vn[!is.na(row.vn)]
                              if (length(row.vn)) {
                                return(any(row.vn > cor_threshold.n))
                              } else {
                                return(FALSE)
                              }
                            })
      
      redundancy.vl <- c(redundancy.vl, group_cor.vl)
      
    }
    
  }
  
  # exprs.mn <- exprs.mn[!redundancy.vl, ]
  fdata.df <- fdata.df[!redundancy.vl, ]
  
  input.eset[rownames(fdata.df), ]
  
}



.format_metabonames <- function(metabo.mset, mice_id.df, mice_id.vc, mice_num.vc) {
  
  for (set.c in ProMetIS::metabo_sets.vc()) {
    
    eset <- metabo.mset[[set.c]]
    
    # sample names formatting and ordering
    
    pdata.df <- Biobase::pData(eset)
    
    if (grepl("(hyper|hilic)", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 11, 13)
    } else if (grepl("acqui", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    } else
      stop("Unknown metabolomics dataset name.")
    
    stopifnot(identical(sort(samp_num.vc), sort(mice_num.vc)))
    
    if (grepl("acqui", set.c)) {
      eset <- eset[, order(as.numeric(samp_num.vc))]
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    }
    
    stopifnot(identical(samp_num.vc, mice_num.vc))
    
    pdata.df <- cbind.data.frame(mice_id.df,
                                 initial_name = Biobase::sampleNames(eset),
                                 pdata.df)
    
    Biobase::sampleNames(eset) <- mice_id.vc
    Biobase::pData(eset) <- pdata.df
    
    # ordering features metadata (supplementary columns will be moved at the end
    # of the data frame with the ordermeta_mset function)
    
    fdata.df <- Biobase::fData(eset)
    
    fdata.df <- data.frame(chromato = rep(unlist(strsplit(set.c, split = "_"))[3],
                                          nrow(fdata.df)),
                           fdata.df,
                           stringsAsFactors = FALSE)
    
    supp.vi <- which(!(colnames(fdata.df) %in% ProMetIS:::.fvarLabels_first.ls()[[set.c]]))
    colnames(fdata.df)[supp.vi] <- paste0("supp_",
                                          colnames(fdata.df)[supp.vi])
    
    Biobase::fData(eset) <- fdata.df
    
    stopifnot(methods::validObject(eset))
    
    Biobase::featureNames(eset) <- fdata.df[, "namecustom"]
    
    metabo.mset <- MultiDataSet::add_eset(metabo.mset,
                                          eset,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  metabo.mset
  
}