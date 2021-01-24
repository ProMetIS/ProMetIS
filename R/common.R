## All ----

#### subsetting (MultiDataSet) ####

#' @rdname subsetting
#' @export
setMethod("subsetting", signature(x = "MultiDataSet"),
          function(x,
                   set.c = NULL,
                   genes.vc = "all",
                   sex.vc = "all",
                   tissues.vc = "all",
                   common_samples.l = FALSE,
                   na_thresh.n = 0.2,
                   var_thresh.n = .Machine$double.eps,
                   imputed_thresh.n = 0.2) {
            
            if (length(genes.vc) == 1 && genes.vc == "all")
              genes.vc <- c("WT", ProMetIS::genes.vc())
            stopifnot(all(genes.vc %in% c("WT", ProMetIS::genes.vc())))
            
            if (length(sex.vc) == 1 && sex.vc == "all")
              sex.vc <- ProMetIS::sex.vc()
            stopifnot(all(sex.vc %in% ProMetIS::sex.vc()))
            
            stopifnot(length(tissues.vc) == 1)
            
            if (tissues.vc != "all") {
              
              stopifnot(tissues.vc %in% ProMetIS::tissues.vc())
              
              x <- x[, names(x)[names(x) %in% c("preclinical",
                                                grep(tissues.vc,
                                                     ProMetIS::sets.vc(),
                                                     value = TRUE))]]
              
            }
            
            if (common_samples.l) {
              x <- MultiDataSet::commonSamples(x)
              ref.eset <- x[[names(x)[1]]]
            } else if ("preclinical" %in% names(x)) {
              ref.eset <- x[["preclinical"]]
            } else {
              sample_names.ls <- Biobase::sampleNames(x)
              all_samples.vc <- Reduce("union", sample_names.ls)
              ref_eset.i <- which.max(unlist(lapply(Biobase::sampleNames(x),
                                                    function(sample_names.vc)
                                                      sum(sample_names.vc %in% all_samples.vc))))[1]
              ref.eset <- x[[names(x)[[ref_eset.i]]]]
            }
            
            samples.vc <- Biobase::sampleNames(ref.eset)
            samples_sel.vl <- Biobase::pData(ref.eset)[, "gene"] %in% genes.vc &
              Biobase::pData(ref.eset)[, "sex"] %in% sex.vc
            
            sub.mset <- x[samples.vc[samples_sel.vl], ]
            
            sub_sets.vc <- names(sub.mset)
            
            for (set.c in sub_sets.vc) {
              
              eset <- sub.mset[[set.c]]
              
              eset <- subsetting(x = eset,
                                 set.c = set.c,
                                 genes.vc = genes.vc,
                                 sex.vc = sex.vc,
                                 tissues.vc = NULL,
                                 common_samples.l = NULL,
                                 na_thresh.n = na_thresh.n,
                                 var_thresh.n = var_thresh.n,
                                 imputed_thresh.n = imputed_thresh.n)
              
              sub.mset <- MultiDataSet::add_eset(sub.mset,
                                                 eset,
                                                 dataset.type = set.c,
                                                 GRanges = NA,
                                                 overwrite = TRUE,
                                                 warnings = FALSE)
              
            }
            
            sub.mset <- sub.mset[, sub_sets.vc]
            
            invisible(sub.mset)
            
          })


#### subsetting (ExpressionSet) ####

#' @rdname subsetting
#' @export
setMethod("subsetting", signature(x = "ExpressionSet"),
          function(x,
                   set.c,
                   genes.vc = "all",
                   sex.vc = "all",
                   tissues.vc = NULL,
                   common_samples.l = NULL,
                   na_thresh.n = 0.2,
                   var_thresh.n = .Machine$double.eps,
                   imputed_thresh.n = 0.2) {
            
            if (length(genes.vc) == 1 && genes.vc == "all")
              genes.vc <- ProMetIS::wtgenes.vc()
            stopifnot(all(genes.vc %in% ProMetIS::wtgenes.vc()))
            
            if (length(sex.vc) == 1 && sex.vc == "all")
              sex.vc <- ProMetIS::sex.vc()
            stopifnot(all(sex.vc %in% ProMetIS::sex.vc()))
            
            samples.vc <- Biobase::sampleNames(x)
            
            samples_sel.vl <- Biobase::pData(x)[, "gene"] %in% genes.vc &
              Biobase::pData(x)[, "sex"] %in% sex.vc
            
            x <- x[, samples.vc[samples_sel.vl]]
            
            filter.vi <- c(nas_zerovar = NA_integer_,
                           overimputed = NA_integer_)
            
            # NAs <= 20% and variance >= 1e-5
            na_zerovar_sel.vl <- ProMetIS:::.filter_na_zerovar(t(Biobase::exprs(x)),
                                                               na_thresh.n = na_thresh.n,
                                                               var_thresh.n = var_thresh.n)
            
            filter.vi["nas_zerovar"] <- sum(!na_zerovar_sel.vl)
            
            # proteomics: observations >= 80% in at least one condition
            if (grepl("proteomics", set.c)) {
              overimputed_sel.vl <- ProMetIS:::.filter_overimputed(eset = x,
                                                                   set.c = set.c,
                                                                   genes.vc = genes.vc,
                                                                   sex.vc = sex.vc,
                                                                   imputed_thresh.n = imputed_thresh.n)
            } else
              overimputed_sel.vl <- rep(TRUE, dim(x)["Features"])
            
            filter.vi["overimputed"] <- sum(!overimputed_sel.vl)
            
            # intersection of both conditions
            feat_sel.vl <- na_zerovar_sel.vl & overimputed_sel.vl
            
            stopifnot(length(feat_sel.vl) == dim(x)["Features"] &&
                        !any(is.na(feat_sel.vl)))
            
            x <- x[feat_sel.vl, ]
            
            # if (sum(!feat_sel.vl))
              message("Nb of discard. feat. in '", set.c, "': ",
                      paste(paste0(names(filter.vi), ": ", filter.vi), collapse = ", "))
            
            invisible(x)
            
          })


.filter_na_zerovar <- function(input.mn,
                               na_thresh.n = 0.2,
                               var_thresh.n = .Machine$double.eps) {
  
  # removing variables with > 20% NA (including 100% NA in females)
  feat_na.vn <- apply(input.mn, 2, function(feat.vn)
    sum(is.na(feat.vn)) / length(feat.vn))
  feat_notna.vl <- feat_na.vn <= na_thresh.n
  # sum(feat_notna.vl)
  
  # removing variables with variance < 1e-5
  feat_var.vn <- apply(input.mn, 2, function(feat.vn)
    var(feat.vn, na.rm = TRUE))
  feat_notzerovar.vl <- !is.na(feat_var.vn) &
    (feat_var.vn >= var_thresh.n)
  # sum(feat_notzerovar.vl)
  
  feat_sel.vl <- feat_notna.vl & feat_notzerovar.vl
  
  stopifnot(length(feat_sel.vl) == ncol(input.mn) &&
              !any(is.na(feat_sel.vl)))
  
  feat_sel.vl
  
}

.imputation_info <- function(eset,
                             set.c) {
  
  prot_pda.df <- Biobase::pData(eset)
  prot_fda.df <- Biobase::fData(eset)
  
  load(system.file("extdata/2_post_processed/metadata_supp.rdata", package = "ProMetIS"))
  
  supp_pda.df <- metadata_supp.ls[[set.c]][["pdata"]]
  supp_fda.df <- metadata_supp.ls[[set.c]][["fdata"]]
  
  stopifnot(all(rownames(prot_pda.df) %in% rownames(supp_pda.df)))
  stopifnot(all(rownames(prot_fda.df) %in% rownames(supp_fda.df)))
  
  prot_pda.df <- cbind.data.frame(prot_pda.df, supp_pda.df[rownames(prot_pda.df), , drop = FALSE])
  prot_fda.df <- cbind.data.frame(prot_fda.df, supp_fda.df[rownames(prot_fda.df), , drop = FALSE])
  
  
  ## checking that the sample names are ordered by increasing ID
  prot_samp.vi <- as.integer(substr(rownames(prot_pda.df), 2, 4))
  stopifnot(identical(prot_samp.vi, sort(prot_samp.vi)))
  
  ## getting imputation info
  value_origin.vl <- vapply(colnames(prot_fda.df), function(colname.c) {
    colname_split.vc <- unlist(strsplit(colname.c, split = "_"))
    grepl("^OriginOfValueabundance", colname.c) &
      colname_split.vc[length(colname_split.vc)] %in% rownames(prot_pda.df)
  }, FUN.VALUE = logical(1))
  value_origin.df <- prot_fda.df[, value_origin.vl]
  colnames(value_origin.df) <- gsub("_run90methode30K",
                                    "",
                                    gsub("_mgf", "",
                                         gsub("OriginOfValueabundance_", "",
                                              colnames(value_origin.df))))
  
  ## re-ordering imputation info to match sample names
  value_origin_samp.vc <- vapply(colnames(value_origin.df), function(colname.c) {
    colname_split.vc <- unlist(strsplit(colname.c, split = "_"))
    colname_split.vc[length(colname_split.vc)]}, FUN.VALUE = character(1))
  temp <- value_origin_samp.vc
  value_origin_samp.vc <- names(value_origin_samp.vc)
  names(value_origin_samp.vc) <- temp
  value_origin.df <- value_origin.df[, value_origin_samp.vc[rownames(prot_pda.df)]]
  
  stopifnot(identical(names(value_origin_samp.vc[rownames(prot_pda.df)]), Biobase::sampleNames(eset)))
  colnames(value_origin.df) <- Biobase::sampleNames(eset)
  
  imputed.mi <- apply(value_origin.df, 2, DAPAR_is.MV)
  mode(imputed.mi) <- "integer"
  
  stopifnot(!any(is.na(c(imputed.mi))))
  
  return(imputed.mi)
  
}

.filter_overimputed <- function(eset,
                                set.c,
                                genes.vc,
                                sex.vc,
                                imputed_thresh.n) {

  imputed.mi <- .imputation_info(eset = eset, set.c = set.c)
  
  stopifnot(identical(Biobase::sampleNames(eset), colnames(imputed.mi)))
  
  if (length(genes.vc) >= 2) {
    # general case: no restriction about sex
    # or e.g. LAT vs WT on males (or females) only
    
    factor.fc <- factor(Biobase::pData(eset)[, "gene"],
                        levels = ProMetIS::wtgenes.vc()[ProMetIS::wtgenes.vc() %in% genes.vc])
    
  } else if (length(genes.vc) == 1 && length(sex.vc) == 2) {
    # males vs females on LAT (or WT) only
    
    factor.fc <- factor(Biobase::pData(eset)[, "sex"],
                        levels = ProMetIS::sex.vc())
    
  } else
    stop("The corresponding imputation metric for this combination of genotype(s) and gender(s) is not currently available.",
         call. = FALSE)
  
  imputed_factor.mi <- t(apply(imputed.mi, 1, function(var.vn) {
    tapply(var.vn, factor.fc, sum)
  }))
  imputed_prop.mn <- sweep(imputed_factor.mi, 2, table(factor.fc), "/")
  
  imputed.ml <- imputed_prop.mn >= imputed_thresh.n
  imputed.vn <- rowSums(imputed.ml, na.rm = TRUE)
  
  feat_sel.vl <- imputed.vn <= 1
  
  stopifnot(length(feat_sel.vl) == dim(eset)["Features"] &&
              !any(is.na(feat_sel.vl)))
  
  feat_sel.vl
  
}


metadata_select <- function(mset,
                            step.c) { # e.g. step.c = "2_post_processed"
  
  if(is.null(names(mset))) { # preclinical expression set
    preclinical.eset <- mset
    mset <- MultiDataSet::createMultiDataSet()
    mset <- MultiDataSet::add_eset(mset,
                                   preclinical.eset,
                                   dataset.type = "preclinical",
                                   GRanges = NA,
                                   overwrite = TRUE,
                                   warnings = FALSE)
  }
  
  mset_names.vc <- names(mset)
  
  metadata_supp_file.c <- paste0("../inst/extdata/", step.c, "/metadata_supp.rdata")
  
  if(file.exists(metadata_supp_file.c)) {
    load(metadata_supp_file.c)
  } else {
    metadata_supp.ls <- vector(mode = "list", length = length(ProMetIS::sets.vc()))
    names(metadata_supp.ls) <- ProMetIS::sets.vc()
  }
  
  for (set.c in mset_names.vc) {
    
    eset <- mset[[set.c]]
    
    # sample metadata
    
    pdata.df <- Biobase::pData(eset)
    
    samplemeta.vc <- .sample_metadata_select(set.c)
    
    samplemeta.vc <- samplemeta.vc[samplemeta.vc %in% colnames(pdata.df)]
    
    samplemeta_supp.vc <- setdiff(colnames(pdata.df), samplemeta.vc)
    
    if(length(samplemeta_supp.vc)) {
      pdata_supp.df <- pdata.df[, samplemeta_supp.vc, drop = FALSE]
    } else
      pdata_supp.df <- data.frame()
    
    Biobase::pData(eset) <- pdata.df[, samplemeta.vc]
    
    # variable metadata
    
    fdata.df <- Biobase::fData(eset)
    
    variablemeta.vc <- .variable_metadata_select(fdata.df = fdata.df, set.c = set.c)
    
    variablemeta.vc <- variablemeta.vc[variablemeta.vc %in% colnames(fdata.df)]
    
    variablemeta_supp.vc <- setdiff(colnames(fdata.df), variablemeta.vc)
    
    if(length(variablemeta_supp.vc)) {
      fdata_supp.df <- fdata.df[, variablemeta_supp.vc, drop = FALSE]
    } else
      fdata_supp.df <- data.frame()
    
    Biobase::fData(eset) <- fdata.df[, variablemeta.vc]
    
    stopifnot(methods::validObject(eset))
    
    mset <- MultiDataSet::add_eset(mset,
                                   eset,
                                   dataset.type = set.c,
                                   GRanges = NA,
                                   overwrite = TRUE,
                                   warnings = FALSE)
    
    metadata_supp.ls[[set.c]] <- list(pdata = pdata_supp.df,
                                      fdata = fdata_supp.df)
 
    
  }
  
  mset <- mset[, mset_names.vc]
  
  save(metadata_supp.ls, file = metadata_supp_file.c)
  
  message("Supplementary metadata written in:\n", metadata_supp_file.c)
  
  return(invisible(mset))
  
}

.sample_metadata_select <- function(set.c) {
  
  first.vc <- c("gene",
                "mouse_nb",
                "sex")
  
  # add.vc <- ""
  # if (set.c == "preclinical")
  #   add.vc <- c("mouse_id",
  #               "genotype",
  #               "project")
  
  # metabolomics
  
  # if (grepl("metabolomics", set.c))
    # add.vc <- "initial_name"
  
  # proteomics
  
  # if (grepl("proteomics", set.c))
  #   add.vc <- "Sample.name"
  
  # first.vc <- c(first.vc, add.vc)
  
  return(first.vc)
  
}

.variable_metadata_select <- function(fdata.df, set.c) {

  ## post-processing
  
  if (set.c == "preclinical")
    varmeta.vc <- c("measurement",
                     "category")
  
  if (grepl("metabolomics", set.c)) {
    
    varmeta.vc <- c("chromato",
                     "MT",
                     "mz",
                     "rt",
                     "isotopes",
                     "adduct",
                     "pcgroup",
                     "redund_group",
                     "redund_iso_add_frag",
                     "name")
    
    if (grepl("(hyper|hilic)", set.c)) {
      
      varmeta.vc <- c(varmeta.vc,
                       c("formula",
                         "monoisotopic_mass",
                         "kegg_id",
                         "kegg_pathway_family",
                         "kegg_pathways",
                         "kegg_subpathways",
                         "chebi_id",
                         "hmdb_id",
                         "pubchem_id",
                         "inchikey",
                         "inchi"))
      
    } else if (grepl("acqui", set.c))
      varmeta.vc <- c(varmeta.vc,
                       "chebi_id",
                       "annot_level",
                       "annot_confidence")
  }
  
  if (grepl("proteomics", set.c))
    varmeta.vc <- c("accession",
                     "description",
                     "uniprot_id")

  
  ## hypothesis testing
  
  limma_col.vc <- grep("limma", colnames(fdata.df), value = TRUE)
  if (length(limma_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    limma_col.vc)
    
  
  # VIP
  
  vip_col.vc <- grep("OPLSDA_VIP-pred", colnames(fdata.df), value = TRUE)
  if (length(vip_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    vip_col.vc)
    
    # Feature selection
    
  biosign_col.vc <- grep("biosign_", colnames(fdata.df), value = TRUE)
  if (length(biosign_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    biosign_col.vc)
    
    
    # mixOmics
    
    mixomics_col.vc <- grep("mixomics_", colnames(fdata.df), value = TRUE)
    if (length(mixomics_col.vc))
      varmeta.vc <- c(varmeta.vc,
                      mixomics_col.vc)
    
 
  
  return(varmeta.vc)
  
}


#' @export
abbrev_mset <- function(mset, full_to_short.l = TRUE) {
  
  mset_names.vc <- names(mset)
  
  abbrev.mset <- MultiDataSet::createMultiDataSet()
  
  for (set.c in mset_names.vc) {
    
    eset <- mset[[set.c]]
    
    abbrev.c <- ProMetIS::sets_abbrev.vc(full_to_short.l = full_to_short.l)[set.c]
    
    exp_data <- Biobase::experimentData(eset)
    
    exp_data@title <- abbrev.c
    
    Biobase::experimentData(eset) <- exp_data
    
    stopifnot(methods::validObject(eset))
    
    abbrev.mset <- MultiDataSet::add_eset(abbrev.mset,
                                          eset,
                                          dataset.type = abbrev.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  abbrev.mset <- abbrev.mset[, ProMetIS::sets_abbrev.vc(full_to_short.l = full_to_short.l)[mset_names.vc]]
  
  invisible(abbrev.mset)
  
}

.sets_abbrev.vc <- function(full_to_short.l = TRUE) {
  
  abbrev.vc <- gsub("metabolomics", "met",
                    gsub("proteomics", "prot",
                         gsub("liver", "liv",
                              gsub("plasma", "plas",
                                   gsub("hyper", "h",
                                        gsub("acqui", "a",
                                             gsub("hilic", "hil",
                                                  gsub("-pos", "-p",
                                                       gsub("-neg", "-n",
                                                            ProMetIS::sets.vc(), fixed = TRUE),
                                                       fixed = TRUE))))))))
  full_to_short.vc <- abbrev.vc
  names(full_to_short.vc) <- ProMetIS::sets.vc()
  
  short_to_full.vc <- ProMetIS::sets.vc()
  names(short_to_full.vc) <- abbrev.vc
  
  if (full_to_short.l) {
    return(full_to_short.vc)
  } else {
    return(short_to_full.vc)
  }
  
}

## Proteomics ----

# Functions from the DAPAR R/Bioconductor package
DAPAR_is.OfType <- function(data, type) 
{
  return(type == data)
}
DAPAR_is.MV <- function(data)
{
  POV = DAPAR_is.OfType(data, "POV")
  MEC = DAPAR_is.OfType(data, "MEC")
  isNA = is.na(data)
  df <- POV | MEC | isNA
  return(df)
}

#' Computing the 'imputation' sample nb and percent (i.e. 'imputed by Prostar')
#' for each variable (for the variables with too high imputation metric in both
#' genotype conditions to be removed in the statistical analysis )
# #' @export
# imputation_metrics <- function(eset) {
#   
#   prot_pda.df <- Biobase::pData(eset)
#   prot_fda.df <- Biobase::fData(eset)
#   prot_samp.vc <- Biobase::sampleNames(eset)
#   
#   ## checking that the sample names are ordered by increasing ID
#   prot_samp.vi <- as.integer(substr(prot_samp.vc, 2, 4))
#   stopifnot(identical(prot_samp.vi, sort(prot_samp.vi)))
#   
#   ## getting imputation info
#   value_origin.df <- prot_fda.df[, grep("OriginOfValueabundance",
#                                         colnames(prot_fda.df), value = TRUE)]
#   colnames(value_origin.df) <- gsub("_run90methode30K",
#                                     "",
#                                     gsub("_mgf", "",
#                                          gsub("supp_OriginOfValueabundance_", "",
#                                               colnames(value_origin.df))))
#   
#   ## re-ordering imputation info to match sample names
#   if (Biobase::experimentData(eset)@title == "proteomics_liver") {
#     file_to_sample.vc <- prot_pda.df[, "supp_sample name"]
#     names(file_to_sample.vc) <- gsub("abundance_", "",
#                                      gsub(".mgf", "",
#                                           prot_pda.df[, "supp_Sample.name"], fixed = TRUE))
#     colsel.vl <- colnames(value_origin.df) %in% names(file_to_sample.vc)
#     value_origin.df <- value_origin.df[, colsel.vl]
#     colnames(value_origin.df) <- file_to_sample.vc[colnames(value_origin.df)]
#   }
#   
#   value_samp.vi <- as.integer(colnames(value_origin.df), 1, 3)
#   value_origin.df <- value_origin.df[, order(value_samp.vi)]
#   
#   ## getting genotype factor
#   gene.fc <- factor(prot_pda.df[, "gene"],
#                     levels = c("WT", ProMetIS::genes.vc()))
#   
#   ## getting sex factor
#   sex.fc <- factor(prot_pda.df[, "sex"],
#                    levels = ProMetIS::sex.vc())
#   
#   ## WT, LAT and MX2 imputation metric
#   imputed.mn <- t(apply(value_origin.df, 1, function(value.vc) {
#     tapply(value.vc, gene.fc, function(x) sum(DAPAR_is.MV(x)))
#   }))
#   imputed.mn <- round(sweep(imputed.mn, 2, table(gene.fc), "/"), 2)
#   colnames(imputed.mn) <- paste0("imputed_mfWLX_", colnames(imputed.mn))
#   
#   prot_fda.df <- cbind.data.frame(prot_fda.df,
#                                   imputed.mn)
#   
#   if (Biobase::experimentData(eset)@title == "proteomics_liver") {
#     
#     ## LAT and WT imputation metric by sex
#     ### M
#     mWL_sel.vl <- gene.fc %in% c("WT", "LAT") & sex.fc == "M"
#     mWL_gene.fc <- factor(gene.fc[mWL_sel.vl])
#     mWL_value.df <- value_origin.df[, mWL_sel.vl]
#     
#     mWL_imputed.mn <- t(apply(mWL_value.df, 1, function(value.vc) {
#       tapply(value.vc, mWL_gene.fc, function(x) sum(DAPAR_is.MV(x)))
#     }))
#     mWL_imputed.mn <- round(sweep(mWL_imputed.mn, 2,
#                                   table(mWL_gene.fc), "/"), 2)
#     colnames(mWL_imputed.mn) <- paste0("imputed_mWL_",
#                                        colnames(mWL_imputed.mn))
#     
#     ### F
#     fWL_sel.vl <- gene.fc %in% c("WT", "LAT") & sex.fc == "F"
#     fWL_gene.fc <- factor(gene.fc[fWL_sel.vl])
#     fWL_value.df <- value_origin.df[, fWL_sel.vl]
#     
#     fWL_imputed.mn <- t(apply(fWL_value.df, 1, function(value.vc) {
#       tapply(value.vc, fWL_gene.fc, function(x) sum(DAPAR_is.MV(x)))
#     }))
#     fWL_imputed.mn <- round(sweep(fWL_imputed.mn, 2,
#                                   table(fWL_gene.fc), "/"), 2)
#     colnames(fWL_imputed.mn) <- paste0("imputed_fWL_",
#                                        colnames(fWL_imputed.mn))
#     
#     ## M and F imputation by (LAT/WT)
#     ### WT
#     mfW_sel.vl <- gene.fc == "WT"
#     mfW_sex.fc <- factor(sex.fc[mfW_sel.vl])
#     mfW_value.df <- value_origin.df[, mfW_sel.vl]
#     
#     mfW_imputed.mn <- t(apply(mfW_value.df, 1, function(value.vc) {
#       tapply(value.vc, mfW_sex.fc, function(x) sum(DAPAR_is.MV(x)))
#     }))
#     mfW_imputed.mn <- round(sweep(mfW_imputed.mn, 2,
#                                   table(mfW_sex.fc), "/"), 2)
#     colnames(mfW_imputed.mn) <- paste0("imputed_mfW_",
#                                        colnames(mfW_imputed.mn)) 
#     
#     ### LAT
#     mfL_sel.vl <- gene.fc == "LAT"
#     mfL_sex.fc <- factor(sex.fc[mfL_sel.vl])
#     mfL_value.df <- value_origin.df[, mfL_sel.vl]
#     
#     mfL_imputed.mn <- t(apply(mfL_value.df, 1, function(value.vc) {
#       tapply(value.vc, mfL_sex.fc, function(x) sum(DAPAR_is.MV(x)))
#     }))
#     mfL_imputed.mn <- round(sweep(mfL_imputed.mn, 2,
#                                   table(mfL_sex.fc), "/"), 2)
#     colnames(mfL_imputed.mn) <- paste0("imputed_mfL_",
#                                        colnames(mfL_imputed.mn))
#     
#     
#     prot_fda.df <- cbind.data.frame(prot_fda.df,
#                                     mWL_imputed.mn,
#                                     fWL_imputed.mn,
#                                     mfW_imputed.mn,
#                                     mfL_imputed.mn)
#     
#   }
#   
#   Biobase::fData(eset) <- prot_fda.df
#   
#   eset
#   
# }