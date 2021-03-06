context("etox")

test_that("get_etoxid returns correct results", {
  skip_on_cran()

  # test general
  comps <- c("Triclosan", "Glyphosate")
  o1 <- get_etoxid(comps, match = "best")
  o2 <- get_etoxid(comps, match = "all")
  o3 <- get_etoxid("Triclosan", match = "first")
  o4 <- get_etoxid("Triclosan", match = "na")
  do2 <- get_etoxid("Thiamethoxam")

  expect_is(o1, "data.frame")
  expect_is(o2, "data.frame")
  expect_is(o3, "data.frame")
  expect_is(o4, "data.frame")
  expect_is(do2, "data.frame")

  expect_equal(o1$etoxid, c("20179", "9051"))
  expect_equal(o2$etoxid, c("89236", "20179", "9051"))
  expect_equal(o3$distance, "first")
  expect_equal(do2$distance, 0)

  # only synonyms found
  expect_warning(get_etoxid("Tetracyclin"))

  # tests for the article
  data("jagst")
  ids <- get_etoxid(head(unique(jagst$substance),6), match = "best")

  expect_is(ids, "data.frame")
  expect_equal(ids$etoxid, c("8932","8494",NA,"8397","7240","7331"))
  expect_equal(ids$match,c(
    "2,4-Dimethylphenol ( 8932 )",
    "4-Chlor-2-methylphenol ( 8494 )",
    NA,
    "Atrazin ( 8397 )",
    "Benzol ( 7240 )",
    "Desethylatrazin ( 7331 )"
  ))
  expect_equal(ids$distance, c(0, 0, NA, 0, 0, 0))
  expect_equal(ids$query, c(
    "2,4-Dimethylphenol",
    "4-Chlor-2-methylphenol",
    "4-para-nonylphenol",
    "Atrazin",
    "Benzol",
    "Desethylatrazin"
  ))

})

# test_that("etox_basic returns correct results", {
#   skip_on_cran()
#
#   ids <- c("20179", "9051", "xxxxx", NA)
#   o1 <- etox_basic(ids)
#
#   expect_is(o1, 'list')
#   expect_equal(length(o1), 4)
#   expect_equal(o1[['20179']]$cas, "3380-34-5")
#   expect_equal(length(o1[['20179']]), 5)
#   expect_is(o1[['20179']]$synonyms, 'data.frame')
#   expect_true(is.na(o1[[3]]))
#   expect_true(is.na(o1[[4]]))
#})
#
#
# test_that("etox_targets returns correct results", {
#   skip_on_cran()
#
#   ids <- c("20179", "9051", "xxxxx", NA)
#   o1 <- etox_targets(ids)
#
#   expect_is(o1, 'list')
#   expect_equal(length(o1), 4)
#   expect_equal(o1[['20179']]$res$Substance[1], "Triclosan")
#   expect_equal(ncol(o1[['20179']]$res), 33)
#   expect_is(o1[['20179']]$res, 'data.frame')
#   expect_true(is.na(o1[[3]]))
#   expect_true(is.na(o1[[4]]))
# })

# test_that("etox_tests returns correct results", {
#   skip_on_cran()
#
#   ids <- c("20179", "9051", "xxxxx", NA)
#   o1 <- etox_tests(ids)
#
#   expect_is(o1, 'list')
#   expect_equal(length(o1), 4)
#   expect_equal(o1[['20179']]$res$Substance[1], "Triclosan")
#   expect_equal(ncol(o1[['20179']]$res), 41)
#   expect_is(o1[['20179']]$res, 'data.frame')
#   expect_true(is.na(o1[[3]]))
#   expect_true(is.na(o1[[4]]))
# })
#
#
# test_that("etox integration tests", {
#   skip_on_cran()
#
#   comps <- c('Triclosan', 'Glyphosate', 'xxxx')
#   ids_b <- get_etoxid(comps, match = 'best')
#   ids_a <- get_etoxid(comps, match = 'all')
#
#   # etox_*() can handle only vector inputs (so using match = 'all' does not work)
#   expect_error(etox_basic(ids_a))
#   expect_error(etox_targets(ids_a))
#   expect_error(etox_tests(ids_a))
#
#
#   int1 <- etox_basic(ids_b$etoxid)
#   int2 <- etox_targets(ids_b$etoxid)
#   int3 <- etox_tests(ids_b$etoxid)
#
#   expect_is(int1, 'list')
#   expect_equal(length(int1), 3)
#   expect_equal(int1[['20179']]$cas, "3380-34-5")
#   expect_equal(length(int1[['20179']]), 5)
#   expect_is(int1[['20179']]$synonyms, 'data.frame')
#   expect_true(is.na(int1[[3]]))
#
#   expect_is(int2, 'list')
#   expect_equal(length(int2), 3)
#   expect_equal(int2[['20179']]$res$Substance[1], "Triclosan")
#   expect_equal(ncol(int2[['20179']]$res), 33)
#   expect_is(int2[['20179']]$res, 'data.frame')
#   expect_true(is.na(int2[[3]]))
#
#   expect_is(int3, 'list')
#   expect_equal(length(int3), 3)
#   expect_equal(int3[['20179']]$res$Substance[1], "Triclosan")
#   expect_equal(ncol(int3[['20179']]$res), 41)
#   expect_is(int3[['20179']]$res, 'data.frame')
#   expect_true(is.na(int3[[3]]))
# })