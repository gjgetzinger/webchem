context("utils")


test_that("is.cas() returns correct results", {
  skip_on_cran()

  expect_true(is.cas('64-17-5'))
  expect_false(is.cas('64175'))
  expect_false(is.cas('4-17-5'))
  expect_false(is.cas('64-177-6'))
  expect_false(is.cas('64-17-55'))
  expect_false(is.cas('64-17-6'))
  expect_error(is.cas(c('64-17-5', '64-17-5')))
})

test_that("as.cas() handles properly formatted CAS",{
  skip_on_cran()

  expect_identical(as.cas("64-17-5"), "64-17-5")
  expect_silent(as.cas("64-17-5"))
})

test_that("is.inchikey() returns correct results", {
  skip_on_cran()

  expect_true(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-N'))
  expect_false(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA'))
  expect_false(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-5'))
  expect_false(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-n'))
  expect_false(is.inchikey('BQJCRHHNABKAKU/KBQPJGBKSA/N'))
  expect_false(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKXA-N'))
  expect_false(is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSB-N'))
  expect_error(is.inchikey(c('BQJCRHHNABKAKU-KBQPJGBKSA-N', 'BQJCRHHNABKAKU-KBQPJGBKSA-N')))

  g <- is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-N', type = 'chemspider')
  b <- is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA', type = 'chemspider')

  expect_true(g)
  expect_false(b)

  expect_error(is.inchikey(c('BQJCRHHNABKAKU-KBQPJGBKSA', 'BQJCRHHNABKAKU-KBQPJGBKSA-N'), type = 'chemspider'))
  expect_error(is.inchikey_cs(c('BQJCRHHNABKAKU-KBQPJGBKSA', 'BQJCRHHNABKAKU-KBQPJGBKSA-N')))
  expect_error(is.inchikey_format(c('BQJCRHHNABKAKU-KBQPJGBKSA', 'BQJCRHHNABKAKU-KBQPJGBKSA-N')))
  expect_equal(length(g), 1)
})


test_that("is.smiles() returns correct results", {
  skip_on_cran()

  expect_true(is.smiles('Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl'))
  expect_false(is.smiles('Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)ClWWX'))
  expect_error(is.smiles(c('Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl', 'Clc1ccc(cc1)C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl')))
})

test_that("extr_num() returns correct results", {
  skip_on_cran()

  expect_equal(extr_num("Melting Pt : -44.6 deg C"), -44.6)
  expect_equal(extr_num("Melting Pt : 44.6 deg C"), 44.6)
  expect_equal(extr_num("Melting Pt : 446 deg C"), 446)
})

test_that("as.cas() returns correct reults", {
  skip_on_cran()
  expect_equal(as.cas(58082), "58-08-2")
  expect_equal(as.cas(123456789), NA)
  expect_identical(as.cas(c(761659, 123456789, "hexenol")), c("761-65-9", NA, NA))
})

