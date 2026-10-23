Sys.setenv(SCRNA_MASTER_KEY = "unit-test-master-key")   # hermetic: no keyfile
source(file.path("..", "..", "modules", "auth", "auth_crypto.R"))

test_that("encrypt/decrypt round-trips", {
  s <- "sk-ant-api03-abcDEF123"
  enc <- encrypt_secret(s)
  expect_true(startsWith(enc, "enc:v1:"))
  expect_false(identical(enc, s))
  expect_equal(decrypt_secret(enc), s)
})

test_that("empty and NULL return empty string", {
  expect_equal(encrypt_secret(""), "")
  expect_equal(encrypt_secret(NULL), "")
  expect_equal(decrypt_secret(""), "")
  expect_equal(decrypt_secret(NULL), "")
})

test_that("legacy plaintext passes through decrypt unchanged", {
  expect_equal(decrypt_secret("sk-ant-plaintext-legacy"), "sk-ant-plaintext-legacy")
})

test_that("is_encrypted detects prefix", {
  expect_true(is_encrypted(encrypt_secret("x")))
  expect_false(is_encrypted("sk-ant-plain"))
  expect_false(is_encrypted(""))
})

test_that("double-encrypt is idempotent (does not wrap twice)", {
  enc <- encrypt_secret("sk-ant-abc")
  expect_equal(encrypt_secret(enc), enc)
})
