# ## Proteomics ----
# 
# .format_imputation <- function(eset) {
#   
#   prot_pda.df <- Biobase::pData(eset)
#   prot_fda.df <- Biobase::fData(eset)
#   
#   ## checking that the sample names are ordered by increasing ID
#   prot_samp.vi <- as.integer(substr(Biobase::sampleNames(eset), 2, 4))
#   stopifnot(identical(prot_samp.vi, sort(prot_samp.vi)))
#   
#   ## getting imputation info
#   value_origin.df <- prot_fda.df[, grep("OriginOfValueabundance",
#                                         colnames(prot_fda.df), value = TRUE)]
#   colnames(value_origin.df) <- gsub("_run90methode30K",
#                                     "",
#                                     gsub("_mgf", "",
#                                          gsub("OriginOfValueabundance_", "",
#                                               colnames(value_origin.df))))
#   
#   ## re-ordering imputation info to match sample names
#   if (Biobase::experimentData(eset)@title == "proteomics_liver") {
#     file_to_sample.vc <- prot_pda.df[, "sample name"]
#     names(file_to_sample.vc) <- gsub("abundance_", "",
#                                      gsub(".mgf", "",
#                                           prot_pda.df[, "Sample.name"], fixed = TRUE))
#     colsel.vl <- colnames(value_origin.df) %in% names(file_to_sample.vc)
#     value_origin.df <- value_origin.df[, colsel.vl]
#     colnames(value_origin.df) <- file_to_sample.vc[colnames(value_origin.df)]
#   }
#   
#   value_samp.vi <- as.integer(colnames(value_origin.df), 1, 3)
#   value_origin.df <- value_origin.df[, order(value_samp.vi)]
#   
#   stopifnot(identical(colnames(value_origin.df), as.character(prot_samp.vi)))
#   colnames(value_origin.df) <- Biobase::sampleNames(eset)
#   
#   imputed.mi <- apply(value_origin.df, 2, DAPAR_is.MV)
#   mode(imputed.mi) <- "integer"
#   
#   stopifnot(!any(is.na(c(imputed.mi))))
#   
#   return(imputed.mi)
#   
# }


## Metabolomics ----

# Postprocessing of metabolomics datasets
#  Imputation
#  RT filtering
#  Blank filtering
#  Pool dilution correlation
#  Signal drift correction
#  poolCV <= 0.3
#  poolCV_over_sampleCV <= 1
#  Chemical redundancy (Monnerie et al., 1999)
metabo_postprocessing <- function(metabo.mset,
                                  drift_correct.c = c("none", "pool", "sample"),
                                  .discard_pools.l = TRUE) {
  
  for (set.c in names(metabo.mset)) {
    
    print(set.c)
    
    eset <- metabo.mset[[set.c]]
    
    # Imputation (MTH Paris)
    
    if (grepl("c18hyper", set.c)) {
      
      exprs.mn <- Biobase::exprs(eset)
      
      exprs.mn[is.na(exprs.mn)] <- 1
      
      Biobase::exprs(eset) <- exprs.mn
      
    }
    
    # RT filtering (MTH Clermont)
    
    ## Converting RT in seconds
    rt.vn <- Biobase::fData(eset)[, "rt"]
    if (max(rt.vn) < 60)
      Biobase::fData(eset)[, "rt"] <- rt.vn * 60
    
    if (grepl("c18acqui", set.c))
      eset <- eset[which(Biobase::fData(eset)[, "rt"] >= 0.4 * 60 &
                           Biobase::fData(eset)[, "rt"] <= 22 * 60), ]
    
    # Blank filtering
    
    eset <- phenomis::inspecting(eset,
                                 figure.c = ifelse(.discard_pools.l, "none", "interactive"),
                                 report.c = "none")
    
    eset <- eset[which(Biobase::fData(eset)[, "blankMean_over_sampleMean"] <= 0.33), ]
    
    # Pool dilution (MTH Paris)
    
    if (grepl("(hyper|hilic)", set.c)) {
      
      eset <- eset[which(Biobase::fData(eset)[, "poolDil_cor"] >= 0.7), ]
      
      ## Setting pool1 to pool for subsequent use in the 'pool CV < 30%' filter
      
      Biobase::pData(eset)[, "sampleType"] <- gsub("pool1", "pool",
                                                   Biobase::pData(eset)[, "sampleType"])
      
    }
    
    # Discarding blanks and pool dilutions
    
    eset <- eset[, which(Biobase::pData(eset)[, "sampleType"] %in% c("sample", "pool"))]
    
    # Signal drift correction
    
    ## ordering eset according to injection order
    
    eset <- eset[, order(Biobase::pData(eset)[, "injectionOrder"])]
    
    if (grepl("acqui", set.c)) { ## resetting injection order min to 1 (MTH Clermont)
      pdata.df <- Biobase::pData(eset)
      pdata.df[, "injectionOrder"] <- pdata.df[, "injectionOrder"] - min(pdata.df[, "injectionOrder"]) + 1
      Biobase::pData(eset) <- pdata.df
    }
    
    ## keeping only the last 'pool' before the first 'sample' (to limit extrapolation)
    
    pdata.df <- Biobase::pData(eset)
    first_sample.i <- which(pdata.df[, "sampleType"] == "sample")[1]
    last_pool_before_first_sample.i <- first_sample.i - 1 
    stopifnot(pdata.df[last_pool_before_first_sample.i, "sampleType"] == "pool")
    
    eset <- eset[, seq(last_pool_before_first_sample.i, dim(eset)["Samples"], by = 1)]
    
    ## signal drift correction
    
    if (drift_correct.c != "none")
      eset <- phenomis::correcting(eset,
                                   reference.c = drift_correct.c,
                                   title.c = gsub("metabolomics_", "", set.c),
                                   span.n = 2,
                                   figure.c = ifelse(.discard_pools.l, "none", "interactive"))
    
    # NAs and variances
    
    eset <- phenomis::filtering(eset, max_na_prop.n = 0.2)
    
    # pool CV <= 0.3
    
    eset <- phenomis::inspecting(eset, report.c = "none",
                                 figure.c = ifelse(.discard_pools.l, "none", "interactive"))
    
    eset <- eset[which(Biobase::fData(eset)[, "pool_CV"] <= 0.3), ]
    
    # poolCV_over_sampleCV <= 1
    
    eset <- eset[which(Biobase::fData(eset)[, "poolCV_over_sampleCV"] <= 1), ]
    
    # Discarding pools
    
    if (.discard_pools.l)
      eset <- eset[, which(Biobase::pData(eset)[, "sampleType"] != "pool")]
    
    # Reducing chemical redundancy (Monnerie et al., 2019)
    
    eset <- phenomis::reducing(eset)
    
    eset <- eset[Biobase::fData(eset)[, "redund_is"] < 1, ]
    
    
    # Updating the eset object
    
    stopifnot(methods::validObject(eset))
    
    metabo.mset <- MultiDataSet::add_eset(metabo.mset,
                                          eset,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  return(metabo.mset)
  
}


.format_metabonames <- function(metabo.mset,
                                mice.ls) {
  # mice_id.df, mice_id.vc, mice_num.vc) {
  
  for (set.c in ProMetIS::metabo_sets.vc()) {
    
    eset <- metabo.mset[[set.c]]
    
    # sample names formatting and ordering by mouse number
    
    if (grepl("(hyper|hilic)", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 11, 13)
    } else if (grepl("acqui", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    } else
      stop("Unknown metabolomics dataset name.")
    
    stopifnot(identical(sort(samp_num.vc), sort(mice.ls[["num.vc"]])))
    
    if (grepl("acqui", set.c)) {
      eset <- eset[, order(as.numeric(samp_num.vc))]
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    }
    
    stopifnot(identical(sort(samp_num.vc), sort(mice.ls[["num.vc"]])))
    
    samp_ord.vi <- order(samp_num.vc)
    eset <- eset[, samp_ord.vi]
    samp_num.vc <- samp_num.vc[samp_ord.vi]
    
    stopifnot(identical(samp_num.vc, mice.ls[["num.vc"]]))
    
    Biobase::pData(eset) <- cbind.data.frame(mice.ls[["id.df"]],
                                             initial_name = Biobase::sampleNames(eset),
                                             Biobase::pData(eset))
    
    Biobase::sampleNames(eset) <- mice.ls[["id.vc"]]
    
    # variable metadata: adding the name of the chromatographic column

    fdata.df <- Biobase::fData(eset)

    fdata.df[, "chromato"] <- rep(unlist(strsplit(set.c, split = "_"))[3],
                                          nrow(fdata.df))

    Biobase::fData(eset) <- fdata.df

    stopifnot(methods::validObject(eset))
    
    metabo.mset <- MultiDataSet::add_eset(metabo.mset,
                                          eset,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  metabo.mset
  
}