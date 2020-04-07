testthat::context("Final aggregated datasets")


aggregated.ls <- sapply(ProMetIS::genes.vc(),
                        function(gene.c) {
                          phenomis::reading(ProMetIS::aggregated_dir.c(gene.c))
                        })


testthat::test_that("dimensions", {
  
  aggregated_dim.ls <- lapply(aggregated.ls,
                              function(gene.mset) {
                                t(sapply(names(gene.mset),
                                         function(set.c) dim(gene.mset[[set.c]])))
                              })

  # gene.c <- "MX2"
  # paste(paste(rownames(aggregated_dim.ls[[gene.c]]),
  #             apply(aggregated_dim.ls[[gene.c]], 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_dim.ls <- list(LAT = c("clinics|181|28",
                              "metabolomics_liver_c18hyper_pos|8938|28",
                              "metabolomics_liver_hilic_neg|3630|28",
                              "metabolomics_plasma_c18acqui_neg|1943|28",
                              "metabolomics_plasma_c18acqui_pos|3276|28",
                              "metabolomics_plasma_c18hyper_pos|9133|28",
                              "metabolomics_plasma_hilic_neg|3402|28",
                              "proteomics_liver|2098|28",
                              "proteomics_plasma|419|24"),
                      MX2 = c("clinics|212|29",
                              "metabolomics_liver_c18hyper_pos|8938|29",
                              "metabolomics_liver_hilic_neg|3630|29",
                              "metabolomics_plasma_c18acqui_neg|1943|29",
                              "metabolomics_plasma_c18acqui_pos|3276|29",
                              "metabolomics_plasma_c18hyper_pos|9133|29",
                              "metabolomics_plasma_hilic_neg|3402|29",
                              "proteomics_liver|2090|29",
                              "proteomics_plasma|422|25"))
  test_dim.ls <- lapply(test_dim.ls,
                        function(test_dim.vc) {
                          test_dim.mn <- sapply(test_dim.vc,
                                                function(x)
                                                  unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
                          mode(test_dim.mn) <- "integer"
                          dimnames(test_dim.mn) <- list(c("Features", "Samples"),
                                                        sapply(test_dim.vc,
                                                               function(x)
                                                                 unlist(strsplit(x, "|", fixed = TRUE))[1],
                                                               USE.NAMES = FALSE))
                          t(test_dim.mn)
                        })
  testthat::expect_identical(aggregated_dim.ls,
                             test_dim.ls)
  
})

testthat::test_that("significant", {
  
  aggregated_signif.ls <- lapply(names(aggregated.ls),
                                 function(gene.c) {
                                   gene.mset <- aggregated.ls[[gene.c]]
                                   t(sapply(names(gene.mset),
                                            function(set.c) {
                                              set_fda.df <- Biobase::fData(gene.mset[[set.c]])
                                              if (gene.c == "LAT" && set.c == "proteomics_liver") {
                                                return(apply(set_fda.df[, c("limmaM_WT.LAT_signif",
                                                                            "limmaF_WT.LAT_signif")],
                                                             2,
                                                             function(y) sum(y, na.rm = TRUE)))
                                              } else {
                                                return(c(sum(set_fda.df[, paste0("limma2ways_WT.", gene.c, "_signif")],
                                                             na.rm = TRUE),
                                                         NA))
                                              }
                                            }))
                                 })
  names(aggregated_signif.ls) <- names(aggregated.ls)
  
  # gene.c <- "MX2"
  # paste(paste(rownames(aggregated_signif.ls[[gene.c]]),
  #             apply(aggregated_signif.ls[[gene.c]], 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_signif.ls <- list(LAT = c("clinics|0|NA",
                                 "metabolomics_liver_c18hyper_pos|2412|NA",
                                 "metabolomics_liver_hilic_neg|1200|NA",
                                 "metabolomics_plasma_c18acqui_neg|0|NA",
                                 "metabolomics_plasma_c18acqui_pos|11|NA",
                                 "metabolomics_plasma_c18hyper_pos|1|NA",
                                 "metabolomics_plasma_hilic_neg|408|NA",
                                 "proteomics_liver|1|258",
                                 "proteomics_plasma|7|NA"),
                         MX2 = c("clinics|1|NA",
                                 "metabolomics_liver_c18hyper_pos|40|NA",
                                 "metabolomics_liver_hilic_neg|12|NA",
                                 "metabolomics_plasma_c18acqui_neg|1|NA",
                                 "metabolomics_plasma_c18acqui_pos|26|NA",
                                 "metabolomics_plasma_c18hyper_pos|102|NA",
                                 "metabolomics_plasma_hilic_neg|79|NA",
                                 "proteomics_liver|263|NA",
                                 "proteomics_plasma|19|NA"))
  test_signif.ls <- lapply(test_signif.ls,
                           function(test_signif.vc) {
                             test_signif.mn <- sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2:3])
                             suppressWarnings(mode(test_signif.mn) <- "integer")
                             colnames(test_signif.mn) <- sapply(test_signif.vc,
                                                                function(x) unlist(strsplit(x, "|", fixed = TRUE))[1],
                                                                USE.NAMES = FALSE)
                             t(test_signif.mn)
                           })

  testthat::expect_identical(aggregated_signif.ls,
                             test_signif.ls)
  
})

testthat::test_that("proteomics_liver", {
  
  proteomics_liver.eset <- aggregated.ls[["LAT"]][["proteomics_liver"]]
  
  proteoliv_WL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                            set.c = "proteomics_liver",
                                            genes.vc = c("WT", "LAT"))
  testthat::expect_equal(Biobase::dims(proteoliv_WL.eset )["Features", ],
                         2098)
  
  proteoliv_mWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = c("WT", "LAT"),
                                             sex.vc = "M")
  testthat::expect_equal(Biobase::dims(proteoliv_mWL.eset )["Features", ],
                         2084)
  
  proteoliv_fWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = c("WT", "LAT"),
                                             sex.vc = "F")
  testthat::expect_equal(Biobase::dims(proteoliv_fWL.eset )["Features", ],
                         2094)
  
  proteoliv_mfW.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = "WT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfW.eset )["Features", ],
                         2078)
  
  proteoliv_mfL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = "LAT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfL.eset )["Features", ],
                         2093)
  
})













testthat::context("Final aggregated data sets")

## getting the dataset fullnames
res_set.vc <- list.dirs(ProMetIS::statistics_intraomics_dir.c(),
                        full.names = FALSE,
                        recursive = FALSE)
res_set.vc <- res_set.vc[!(res_set.vc %in% c("synthesis", "tests-figures"))]

res_set.vc <- res_set.vc[!(res_set.vc %in% c("liver_metabo_annot_mth-saclay",
                                             "liver_metabo_hilic-neg_mth-saclay",
                                             "liver_metabo_c18-neg_mth-theix",
                                             "plasma_metabo_annot_mth-saclay",
                                             "plasma_metabo_hilic-neg_mth-saclay",
                                             "plasma_metabo_c18-neg_mth-theix"))]
res_set.vc <- gsub("c18-pos_", "", res_set.vc)
res_set.vc <- as.character(res_set.vc)
res_set.vc[res_set.vc == "liver_proteo__"] <- paste0("liver_proteo_", c("F", "M"), "_profi-strasbourg")


## getting the variableMetadata datasets

res_set.ls <- vector(mode = "list", length = length(res_set.vc))
names(res_set.ls) <- res_set.vc
geneset.ls <- vector(mode = "list", length = length(ProMetIS::genes.vc()))
names(geneset.ls) <- ProMetIS::genes.vc()
for (i in seq_along(geneset.ls)) {
  geneset.ls[[i]] <- res_set.ls
}
geneset.ls[["LAT"]][["liver_proteo_profi-strasbourg"]] <- NULL
geneset.ls[["MX2"]][["liver_proteo_F_profi-strasbourg"]] <- NULL
geneset.ls[["MX2"]][["liver_proteo_M_profi-strasbourg"]] <- NULL

for (gene.c in names(geneset.ls)) {
  for (set.c in names(geneset.ls[[gene.c]])) {
    file.c <- file.path(ProMetIS::cs1_pheno_result_dir.c(),
                        switch(unlist(strsplit(set.c, "_"))[2],
                               metabo = {"metabolomics"},
                               proteo = {"proteomics"},
                               phenomin = {"phenomin"}),
                        paste0(gene.c, "_", set.c, ".tsv"))
    geneset.df <- read.table(file.c,
                             comment.char = "",
                             header = TRUE,
                             quote = "",
                             row.names = 1,
                             sep = "\t",
                             stringsAsFactors = FALSE)
    geneset.ls[[gene.c]][[set.c]] <- geneset.df
  }
}

testthat::test_that("dimensions", {
  
  test_dim.an <- result_dim.an <- array(NA_integer_,
                                        dim = c(length(res_set.vc), 2, length(ProMetIS::gene.vc())),
                                        dimnames = list(res_set.vc, c("row", "col"), ProMetIS::gene.vc()))
  
  for (gene.c in ProMetIS::gene.vc()) {
    for (set.c in res_set.vc) {
      if (!is.null(geneset.ls[[gene.c]][[set.c]]))
        result_dim.an[set.c, , gene.c] <- dim(geneset.ls[[gene.c]][[set.c]])
    }
  }
  
  # test_dim.ls <- lapply(ProMetIS::gene.vc(),
  #                       function(gene.c) result_dim.an[, , gene.c])
  # lapply(test_dim.ls,
  #        function(test_dim_gene.mn)
  #        paste(paste(rownames(test_dim_gene.mn),
  #             apply(test_dim_gene.mn, 1, function(y) paste(y, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '"))
  test_dim.ls <- list(LAT = c('liver_metabo_mth-saclay|12569|102',
                              'liver_metabo_mth-theix|2579|97',
                              'liver_proteo_F_profi-strasbourg|2220|101',
                              'liver_proteo_M_profi-strasbourg|2183|96',
                              'liver_proteo_profi-strasbourg|NA|NA',
                              'phenotype_phenomin|208|39',
                              'plasma_metabo_mth-saclay|12551|96',
                              'plasma_metabo_mth-theix|5196|87',
                              'plasma_proteo_profi-toulouse|419|259'),
                      MX2 = c('liver_metabo_mth-saclay|12569|102',
                              'liver_metabo_mth-theix|2579|97',
                              'liver_proteo_F_profi-strasbourg|NA|NA',
                              'liver_proteo_M_profi-strasbourg|NA|NA',
                              'liver_proteo_profi-strasbourg|2187|279',
                              'phenotype_phenomin|208|42',
                              'plasma_metabo_mth-saclay|12551|101',
                              'plasma_metabo_mth-theix|5196|82',
                              'plasma_proteo_profi-toulouse|446|251'))
  for (gene.c in ProMetIS::gene.vc()) {
    test_dim.mn <- sapply(test_dim.ls[[gene.c]], function(x)
      unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
    suppressWarnings(mode(test_dim.mn) <- "integer")
    dimnames(test_dim.mn) <- list(c("row", "col"),
                                  sapply(test_dim.ls[[gene.c]],
                                         function(x) unlist(strsplit(x, "|", fixed = TRUE))[1],
                                         USE.NAMES = FALSE))
    
    test_dim.an[, , gene.c] <- t(test_dim.mn)
  }
  testthat::expect_identical(result_dim.an,
                             test_dim.an)
})

testthat::test_that("significant", {
  
  test_signif.mn <- result_signif.mn <- matrix(NA_integer_,
                                               nrow = length(res_set.vc),
                                               ncol = length(ProMetIS::gene.vc()),
                                               dimnames = list(res_set.vc, ProMetIS::gene.vc()))
  
  for (gene.c in ProMetIS::gene.vc()) {
    for (set.c in res_set.vc) {
      if (!is.null(geneset.ls[[gene.c]][[set.c]]))
        result_signif.mn[set.c, gene.c] <- sum(geneset.ls[[gene.c]][[set.c]][, "significant"])
    }
  }
  
  # paste(paste(rownames(result_signif.mn),
  #             apply(result_signif.mn, 1, function(y) paste(y, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_signif.vc <- c('liver_metabo_mth-saclay|3614|52',
                      'liver_metabo_mth-theix|0|0',
                      'liver_proteo_F_profi-strasbourg|281|NA',
                      'liver_proteo_M_profi-strasbourg|0|NA',
                      'liver_proteo_profi-strasbourg|NA|0',
                      'phenotype_phenomin|0|1',
                      'plasma_metabo_mth-saclay|408|187',
                      'plasma_metabo_mth-theix|11|27',
                      'plasma_proteo_profi-toulouse|7|24')
  test_signif.mn <- sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
  suppressWarnings(mode(test_signif.mn) <- "integer")
  dimnames(test_signif.mn) <- list(ProMetIS::gene.vc(),
                                   sapply(test_signif.vc,
                                          function(x) unlist(strsplit(x, "|", fixed = TRUE))[1],
                                          USE.NAMES = FALSE))
  test_signif.mn <- t(test_signif.mn)
  
  testthat::expect_identical(result_signif.mn,
                             test_signif.mn)
})   

testthat::test_that("plas_prot_toul", {
  
  protoul_fda.df <- geneset.ls[["LAT"]][["plasma_proteo_profi-toulouse"]]
  
  testthat::expect_identical(dim(protoul_fda.df),
                             as.integer(c(419, 259)))
  
  testthat::expect_identical(sum(protoul_fda.df[, "gene_WT.KO_signif"]),
                             as.integer(7))
  
  testthat::expect_identical(sum(protoul_fda.df[, "sex_M.F_signif"]),
                             as.integer(121))
  
  testthat::expect_equal(protoul_fda.df["P01746_Ig heavy chain V region .", "gene_WT.KO_fold"],
                         2.09714622,
                         tolerance = 1e-6)
  testthat::expect_equal(protoul_fda.df["P01746_Ig heavy chain V region .", "gene_WT.KO_BH"],
                         2.251392e-02,
                         tolerance = 1e-6)
  
})
