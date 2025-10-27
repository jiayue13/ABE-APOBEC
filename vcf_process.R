#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr / stringr / readr / tidyr
})

# ---------------- Parameter parsing (no additional dependencies required) ----------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_val <- function(key, default = NULL) {
    hit_eq <- grep(paste0("^--", key, "="), args, value = TRUE)
    if (length(hit_eq) > 0) return(sub(paste0("^--", key, "="), "", hit_eq[1]))
    hit <- which(args == paste0("--", key))
    if (length(hit) > 0 && length(args) >= hit + 1) return(args[hit + 1])
    default
  }
  # s Supports both positional parameters and --s
  s_val <- get_val("s", NA)
  if (is.na(s_val)) {
    if (length(args) >= 1 && !startsWith(args[1], "--")) s_val <- args[1]
  }
  list(
    s        = suppressWarnings(as.integer(ifelse(is.na(s_val), 1L, s_val))),
    input  = get_val("input", "."),
    output   = get_val("output",  "."),
    samples  = get_val("samples", NULL)  # Optional
  )
}

# ---------------- read VCF header (#CHROM) ----------------
read_vcf_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) {
      header <- strsplit(sub("^#", "", line), "\t", fixed = FALSE)[[1]]
      break
    }
  }
  header
}

# ---------------- Sample name cleaning ----------------
sanitize_sample <- function(x) {
  b <- basename(x)
  # Remove the path suffix
  b <- sub("\\.(bam|cram)$", "", b, ignore.case = TRUE)
  b <- sub("\\.(FWD|REV)$", "", b, ignore.case = TRUE)
  # Clear redundant symbols
  b <- gsub("[^[:alnum:]_.-]+", "_", b)
  b <- gsub("_+", "_", b)
  b
}

# ---------------- Main process ----------------
main <- function() {
  a <- parse_args()
  if (!(a$s %in% c(1L, 2L))) stop("Parameter s can only be 1 or 2 (1=FWD/neg, 2=REV/pos).")

  st     <- if (a$s == 1L) "FWD" else "REV"
  strand <- if (a$s == 1L) "neg" else "pos"

  input <- a$input
  output  <- a$output
  if (!dir.exists(output)) dir.create(output, recursive = TRUE, showWarnings = FALSE)

  vcf_path <- file.path(input, sprintf("ProQ-rABE_%s.vcf", st))
  if (!file.exists(vcf_path) && file.exists(paste0(vcf_path, ".gz"))) vcf_path <- paste0(vcf_path, ".gz")
  if (!file.exists(vcf_path)) stop("Input not found:", vcf_path, " or its .gz version")

  message("[INFO] read files:", vcf_path)

  header <- read_vcf_header(vcf_path)
  if (is.null(header)) stop("No #CHROM header row found in VCF.")
  sample_cols <- header[10:length(header)]
  if (length(sample_cols) == 0) stop("Sample column not detected in VCF.")

  # Read data (skip # lines), read all characters
  df <- readr::read_tsv(
    vcf_path,
    comment   = "#",
    col_names = header,
    col_types = readr::cols(.default = readr::col_character()),
    progress  = FALSE
  )

  # Add new strand; retain only the first base of ALT; remove core deletions
  df <- df %>%
    mutate(
      strand = strand,
      ALT    = substr(ALT, 1, 1)
    ) %>%
    drop_na(CHROM, POS, REF, ALT)

  # RNAREF / RNAALT (minus strand complementation)
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  df <- df %>%
    mutate(
      RNAREF = REF,
      RNAALT = ALT
    )
  if (strand == "neg") {
    rnaref_map <- unname(comp[df$RNAREF])
    rnaalt_map <- unname(comp[df$RNAALT])
    df$RNAREF  <- ifelse(is.na(rnaref_map), df$RNAREF, rnaref_map)
    df$RNAALT  <- ifelse(is.na(rnaalt_map), df$RNAALT, rnaalt_map)
  }

  # Sample name (prefers --samples; otherwise inferred from VCF header and cleaned)
  if (!is.null(a$samples) && nchar(a$samples) > 0) {
    user_samples <- strsplit(a$samples, ",", fixed = TRUE)[[1]] %>% trimws()
    if (length(user_samples) != length(sample_cols)) {
      warning("If the number of --samples does not match the VCF sample column, it will be inferred from the VCF header.")
      sList <- sanitize_sample(sample_cols)
    } else {
      sList <- user_samples
    }
  } else {
    sList <- sanitize_sample(sample_cols)
  }

  message("[INFO] Sample columns (VCF header): ", paste(sample_cols, collapse = ", "))
  message("[INFO] Sample name (column prefix): ", paste(sList, collapse = ", "))

# ------- Fixed-position parsing: DP = second paragraph, AD = sixth paragraph; AD takes the second comma (alt) -------
# Vector-safe: Use str_split_fixed to get the sixth paragraph, then use strsplit() to extract the "second comma paragraph"
  for (k in seq_along(sample_cols)) {
    col <- sample_cols[k]
    sp  <- sList[k]

    parts  <- stringr::str_split_fixed(df[[col]], ":", 6)
    DP_chr <- parts[, 2]
    AD_raw <- parts[, 6]

    AD_alt_chr <- vapply(
      strsplit(AD_raw, ",", fixed = TRUE),
      FUN = function(x) if (length(x) >= 2) x[2] else NA_character_,
      FUN.VALUE = character(1)
    )

    df[[paste0(sp, ".DP")]]    <- suppressWarnings(as.numeric(DP_chr))
    df[[paste0(sp, ".AD")]]    <- suppressWarnings(as.numeric(AD_alt_chr))
    # ratio
    num <- df[[paste0(sp, ".AD")]]
    den <- df[[paste0(sp, ".DP")]]
    ratio <- num / den
    ratio[is.infinite(ratio)] <- NA_real_
    df[[paste0(sp, ".ratio")]] <- ratio
  }

  # Average coverage
  dp_cols   <- paste0(sList, ".DP")
  df$avgDP  <- rowMeans(as.matrix(df[dp_cols]), na.rm = TRUE)

  # Output
  out_path <- file.path(output, sprintf("mpileup_fixstrand_%s.vcf", strand))
  message("[INFO] output: ", out_path)

  ok <- tryCatch({
    write.table(df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    TRUE
  }, error = function(e) {
    message("[WARN] Write failure directly: ", conditionMessage(e))
    FALSE
  })

  if (!ok) {
    tmp <- tempfile(fileext = ".vcf")
    write.table(df, file = tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    message("[HINT] Temporary file written：", tmp)
    message("[HINT] Please copy it manually to the target directory, for example：")
    message("       cp ", shQuote(tmp), " ", shQuote(out_path))
  } else {
    message("[OK] Finish：", out_path, "  （total ", nrow(df), " rows）")
  }
}

main()

