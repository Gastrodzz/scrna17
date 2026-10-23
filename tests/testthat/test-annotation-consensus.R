# Tests the N-method annotation consensus helper from modules/annotation.R
ok <- tryCatch({
  suppressWarnings(suppressMessages(
    source(file.path("..", "..", "modules", "annotation.R"))))
  exists(".build_annotation_consensus")
}, error = function(e) FALSE)

test_that("3-method majority consensus resolves agreement, majority, and abstention", {
  skip_if_not(ok, "annotation.R could not be sourced in this environment")
  sctype  <- c("T cell","T cell","B cell","B cell","Unknown","Unknown")
  panglao <- c("T cell","T cell","B cell","B cell","NK cell","NK cell")
  singler <- c("T cell","T cell","Plasma","Plasma","NK cell","NK cell")
  r <- .build_annotation_consensus(
    list(ScType = sctype, PanglaoDB = panglao, clustifyr = singler),
    list(ScType = c(.9,.9,.8,.8,0,0),
         PanglaoDB = c(.7,.7,.7,.7,.6,.6),
         clustifyr = c(.8,.8,.95,.95,.7,.7)))
  expect_equal(r$consensus, c("T cell","T cell","B cell","B cell","NK cell","NK cell"))
  expect_equal(r$source[1], "all")                       # all 3 agree
  expect_equal(r$source[3], "ScType+PanglaoDB")          # 2 vs 1 majority
  expect_equal(r$source[5], "PanglaoDB+clustifyr")         # ScType abstained (Unknown)
  expect_equal(r$disagreement, c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE))
})

test_that("all-Unknown yields Unknown with source none", {
  skip_if_not(ok, "annotation.R could not be sourced in this environment")
  r <- .build_annotation_consensus(
    list(A = c("Unknown","Unknown"), B = c("Unknown","Unknown")),
    list(A = c(0,0), B = c(0,0)))
  expect_equal(r$consensus, c("Unknown","Unknown"))
  expect_equal(r$source, c("none","none"))
})

test_that("two-way tie is broken by the higher-scoring method", {
  skip_if_not(ok, "annotation.R could not be sourced in this environment")
  r <- .build_annotation_consensus(
    list(A = "B cell", B = "T cell"),
    list(A = 0.6, B = 0.9))
  expect_equal(r$consensus, "T cell")   # B had the higher score
})
