testthat::context("Datasets after intraomics statistics")

stat_intra.mset <- phenomis::reading(ProMetIS::statistics_intraomics_dir.c())

testthat::test_that("dimensions", {
  
  stat_intra_dim.mn <- t(sapply(names(stat_intra.mset),
                                function(set.c) dim(stat_intra.mset[[set.c]])))
  
  # paste(paste(rownames(stat_intra_dim.mn),
  #             apply(stat_intra_dim.mn, 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_dim.vc <- c("clinics|213|42",
                   "metabolomics_liver_c18hyper_pos|8938|42",
                   "metabolomics_liver_hilic_neg|3630|42",
                   "metabolomics_plasma_c18acqui_neg|1943|42",
                   "metabolomics_plasma_c18acqui_pos|3276|42",
                   "metabolomics_plasma_c18hyper_pos|9133|42",
                   "metabolomics_plasma_hilic_neg|3402|42",
                   "proteomics_liver|2187|42",
                   "proteomics_plasma|446|36")
  test_dim.mn <- sapply(test_dim.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
  mode(test_dim.mn) <- "integer"
  dimnames(test_dim.mn) <- list(c("Features", "Samples"),
                                sapply(test_dim.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[1], USE.NAMES = FALSE))
  test_dim.mn <- t(test_dim.mn)
  testthat::expect_identical(stat_intra_dim.mn,
                             test_dim.mn)
  
})

testthat::test_that("significant", {
  
  stat_intra_signif.mn <- t(sapply(names(stat_intra.mset),
                                   function(set.c) {
                                     set_fda.df <- Biobase::fData(stat_intra.mset[[set.c]])
                                     if (set.c == "proteomics_liver") {
                                       return(apply(set_fda.df[, c("limmaM_WT.LAT_signif",
                                                                   "limmaF_WT.LAT_signif",
                                                                   "limma2ways_WT.MX2_signif")],
                                                    2,
                                                    function(y) sum(y, na.rm = TRUE)))
                                     } else {
                                       return(c(sum(set_fda.df[, c("limma2ways_WT.LAT_signif")], na.rm = TRUE),
                                                NA,
                                                sum(set_fda.df[, c("limma2ways_WT.MX2_signif")], na.rm = TRUE)))
                                     }
                                   }))
  colnames(stat_intra_signif.mn) <- c("LAT(_M)", "LAT_F", "MX2")
  
  # paste(paste(rownames(stat_intra_signif.mn),
  #             apply(stat_intra_signif.mn, 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_signif.vc <- c("clinics|0|NA|1",
                      "metabolomics_liver_c18hyper_pos|2412|NA|40",
                      "metabolomics_liver_hilic_neg|1200|NA|12",
                      "metabolomics_plasma_c18acqui_neg|0|NA|1",
                      "metabolomics_plasma_c18acqui_pos|11|NA|26",
                      "metabolomics_plasma_c18hyper_pos|1|NA|102",
                      "metabolomics_plasma_hilic_neg|408|NA|79",
                      "proteomics_liver|1|258|263",
                      "proteomics_plasma|7|NA|19")
  test_signif.mn <- sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2:4])
  suppressWarnings(mode(test_signif.mn) <- "integer")
  dimnames(test_signif.mn) <- list(c("LAT(_M)", "LAT_F", "MX2"),
                                   sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[1], USE.NAMES = FALSE))
  test_signif.mn <- t(test_signif.mn)
  testthat::expect_identical(stat_intra_signif.mn,
                             test_signif.mn)
  
})

testthat::test_that("proteomics_liver", {
  
  proteomics_liver.eset <- stat_intra.mset[["proteomics_liver"]]
  
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
                         2120)
  
  proteoliv_fWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = c("WT", "LAT"),
                                           sex.vc = "F")
  testthat::expect_equal(Biobase::dims(proteoliv_fWL.eset )["Features", ],
                         2139)
  
  proteoliv_mfW.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = "WT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfW.eset )["Features", ],
                         2126)
  
  proteoliv_mfL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = "LAT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfL.eset )["Features", ],
                         2145)
  
})
