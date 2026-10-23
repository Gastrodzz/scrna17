source(file.path("..", "..", "modules", "ai_models.R"))

test_that("fallback list is returned when no api key", {
  ch <- list_claude_models(api_key = NULL, include_thinking = FALSE, force = TRUE)
  expect_true(is.character(ch))
  expect_true(length(ch) >= 3)
  expect_true(all(grepl("^claude-", unname(ch))))
  expect_true(!is.null(names(ch)))          # display names present
})

test_that("thinking variants are appended for opus/sonnet only", {
  ch <- list_claude_models(api_key = NULL, include_thinking = TRUE, force = TRUE)
  think <- unname(ch)[grepl("-thinking$", unname(ch))]
  expect_true(length(think) >= 1)
  expect_true(all(grepl("opus|sonnet", think)))   # no haiku-thinking
  expect_false(any(grepl("haiku.*-thinking", think)))
})

test_that("get_default_model prefers a non-thinking sonnet", {
  m <- get_default_model(api_key = NULL)
  expect_true(is.character(m) && length(m) == 1)
  expect_true(grepl("sonnet", m))
  expect_false(grepl("-thinking$", m))
})

test_that("cache is reused within ttl and bypassed with force", {
  a <- list_claude_models(api_key = NULL, force = TRUE)
  b <- list_claude_models(api_key = NULL, force = FALSE)
  expect_identical(a, b)
})
