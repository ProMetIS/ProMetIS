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
                              "metabolomics_liver_c18hyper_pos|6209|28",
                              "metabolomics_liver_hilic_neg|2566|28",
                              "metabolomics_plasma_c18acqui_neg|1612|28",
                              "metabolomics_plasma_c18acqui_pos|6195|28",
                              "metabolomics_plasma_c18hyper_pos|5520|28",
                              "metabolomics_plasma_hilic_neg|2521|28",
                              "proteomics_liver|2098|28",
                              "proteomics_plasma|419|24"),
                      MX2 = c("clinics|212|29",
                              "metabolomics_liver_c18hyper_pos|6209|29",
                              "metabolomics_liver_hilic_neg|2566|29",
                              "metabolomics_plasma_c18acqui_neg|1612|29",
                              "metabolomics_plasma_c18acqui_pos|6195|29",
                              "metabolomics_plasma_c18hyper_pos|5520|29",
                              "metabolomics_plasma_hilic_neg|2521|29",
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
                                 "metabolomics_liver_c18hyper_pos|2133|NA",
                                 "metabolomics_liver_hilic_neg|760|NA",
                                 "metabolomics_plasma_c18acqui_neg|3|NA",
                                 "metabolomics_plasma_c18acqui_pos|8|NA",
                                 "metabolomics_plasma_c18hyper_pos|3|NA",
                                 "metabolomics_plasma_hilic_neg|648|NA",
                                 "proteomics_liver|1|258",
                                 "proteomics_plasma|7|NA"),
                         MX2 = c("clinics|1|NA",
                                 "metabolomics_liver_c18hyper_pos|91|NA",
                                 "metabolomics_liver_hilic_neg|24|NA",
                                 "metabolomics_plasma_c18acqui_neg|1|NA",
                                 "metabolomics_plasma_c18acqui_pos|115|NA",
                                 "metabolomics_plasma_c18hyper_pos|103|NA",
                                 "metabolomics_plasma_hilic_neg|54|NA",
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