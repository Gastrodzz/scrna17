# modules/auth/auth_crypto.R
# Symmetric encryption for secrets at rest (API keys). AES-256-CBC via openssl.
# Master key: SCRNA_MASTER_KEY env var (any string, hashed to 32 bytes) or an
# auto-generated data/.keyfile (32 random bytes, chmod 0600). Ciphertext format:
#   "enc:v1:" + base64( iv[16] || ciphertext )
# The prefix lets decrypt_secret() pass through legacy plaintext untouched.

if (!requireNamespace("openssl", quietly = TRUE))
  stop("auth_crypto.R requires the 'openssl' package.")

.crypto_key_cache <- new.env(parent = emptyenv())

get_master_key <- function() {
  if (!is.null(.crypto_key_cache$key)) return(.crypto_key_cache$key)
  env <- Sys.getenv("SCRNA_MASTER_KEY", "")
  key <- if (nzchar(env)) {
    as.raw(openssl::sha256(charToRaw(env)))         # 32 bytes
  } else {
    kf <- file.path(getwd(), "data", ".keyfile")
    if (!file.exists(kf)) {
      if (!dir.exists(dirname(kf))) dir.create(dirname(kf), recursive = TRUE)
      k <- openssl::rand_bytes(32)
      writeBin(k, kf)
      tryCatch(Sys.chmod(kf, mode = "0600"), error = function(e) NULL)
      k
    } else {
      readBin(kf, what = "raw", n = 32)
    }
  }
  .crypto_key_cache$key <- key
  key
}

is_encrypted <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && startsWith(x, "enc:v1:")
}

encrypt_secret <- function(plaintext) {
  if (is.null(plaintext) || length(plaintext) == 0) return("")
  plaintext <- as.character(plaintext)
  if (!nzchar(plaintext)) return("")
  if (is_encrypted(plaintext)) return(plaintext)   # already encrypted
  iv <- openssl::rand_bytes(16)
  ct <- openssl::aes_cbc_encrypt(charToRaw(plaintext), key = get_master_key(), iv = iv)
  paste0("enc:v1:", openssl::base64_encode(c(iv, as.raw(ct))))
}

decrypt_secret <- function(ciphertext) {
  if (is.null(ciphertext) || length(ciphertext) == 0) return("")
  ciphertext <- as.character(ciphertext)
  if (!nzchar(ciphertext)) return("")
  if (!is_encrypted(ciphertext)) return(ciphertext) # legacy plaintext
  raw <- openssl::base64_decode(sub("^enc:v1:", "", ciphertext))
  iv  <- raw[1:16]; ct <- raw[-(1:16)]
  tryCatch(rawToChar(openssl::aes_cbc_decrypt(ct, key = get_master_key(), iv = iv)),
           error = function(e) "")
}
