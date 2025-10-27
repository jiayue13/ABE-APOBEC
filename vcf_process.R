#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr / stringr / readr / tidyr
})

# ---------------- 参数解析（无需额外依赖） ----------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_val <- function(key, default = NULL) {
    hit_eq <- grep(paste0("^--", key, "="), args, value = TRUE)
    if (length(hit_eq) > 0) return(sub(paste0("^--", key, "="), "", hit_eq[1]))
    hit <- which(args == paste0("--", key))
    if (length(hit) > 0 && length(args) >= hit + 1) return(args[hit + 1])
    default
  }
  # s 既支持位置参数也支持 --s
  s_val <- get_val("s", NA)
  if (is.na(s_val)) {
    if (length(args) >= 1 && !startsWith(args[1], "--")) s_val <- args[1]
  }
  list(
    s        = suppressWarnings(as.integer(ifelse(is.na(s_val), 1L, s_val))),
    workdir  = get_val("workdir", "."),
    resdir   = get_val("resdir",  "."),
    samples  = get_val("samples", NULL)  # 逗号分隔，可选
  )
}

# ---------------- 读取 VCF 表头 (#CHROM) ----------------
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

# ---------------- 样本名清洗 ----------------
sanitize_sample <- function(x) {
  b <- basename(x)
  b <- sub("\\.(bam|cram)$", "", b, ignore.case = TRUE)
  b <- gsub("[^[:alnum:]_.-]+", "_", b)
  b
}

# ---------------- 主流程 ----------------
main <- function() {
  a <- parse_args()
  if (!(a$s %in% c(1L, 2L))) stop("参数 s 只能是 1 或 2（1=FWD/neg, 2=REV/pos）。")

  st     <- if (a$s == 1L) "FWD" else "REV"
  strand <- if (a$s == 1L) "neg" else "pos"

  workdir <- a$workdir
  resdir  <- a$resdir
  if (!dir.exists(resdir)) dir.create(resdir, recursive = TRUE, showWarnings = FALSE)

  vcf_path <- file.path(workdir, sprintf("ProQ-rABE_%s.vcf", st))
  if (!file.exists(vcf_path) && file.exists(paste0(vcf_path, ".gz"))) vcf_path <- paste0(vcf_path, ".gz")
  if (!file.exists(vcf_path)) stop("未找到输入：", vcf_path, " 或其 .gz 版本")

  message("[INFO] 读取文件: ", vcf_path)

  header <- read_vcf_header(vcf_path)
  if (is.null(header)) stop("未在 VCF 中找到 #CHROM 表头行。")
  sample_cols <- header[10:length(header)]
  if (length(sample_cols) == 0) stop("VCF 中未检测到样本列。")

  # 读取数据（跳过 # 行），全部按字符读入
  df <- readr::read_tsv(
    vcf_path,
    comment   = "#",
    col_names = header,
    col_types = readr::cols(.default = readr::col_character()),
    progress  = FALSE
  )

  # 新增 strand；ALT 仅保留首个碱基；去除核心缺失
  df <- df %>%
    mutate(
      strand = strand,
      ALT    = substr(ALT, 1, 1)
    ) %>%
    drop_na(CHROM, POS, REF, ALT)

  # RNAREF / RNAALT（负链互补）—— 使用“向量安全”的方式
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  df <- df %>%
    mutate(
      RNAREF = REF,
      RNAALT = ALT
    )
  if (strand == "neg") {
    # 先映射，返回可能带 NA 的向量；再用 ifelse 替换 NA 为原值（向量安全）
    rnaref_map <- unname(comp[df$RNAREF])
    rnaalt_map <- unname(comp[df$RNAALT])
    df$RNAREF  <- ifelse(is.na(rnaref_map), df$RNAREF, rnaref_map)
    df$RNAALT  <- ifelse(is.na(rnaalt_map), df$RNAALT, rnaalt_map)
  }

  # 样本名（优先 --samples；否则从 VCF 表头推断并清洗）
  if (!is.null(a$samples) && nchar(a$samples) > 0) {
    user_samples <- strsplit(a$samples, ",", fixed = TRUE)[[1]] %>% trimws()
    if (length(user_samples) != length(sample_cols)) {
      warning("--samples 数量与 VCF 样本列不一致，将使用 VCF 表头推断。")
      sList <- sanitize_sample(sample_cols)
    } else {
      sList <- user_samples
    }
  } else {
    sList <- sanitize_sample(sample_cols)
  }

  message("[INFO] 样本列（VCF头）：", paste(sample_cols, collapse = ", "))
  message("[INFO] 样本名（列前缀）：", paste(sList, collapse = ", "))

  # ------- 固定位置解析：DP=第二段, AD=第六段；AD 取逗号第二项（alt） -------
  # 向量安全：使用 str_split_fixed 取到第 6 段，再用 strsplit() 提取“第二个逗号段”
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

  # 平均覆盖
  dp_cols   <- paste0(sList, ".DP")
  df$avgDP  <- rowMeans(as.matrix(df[dp_cols]), na.rm = TRUE)

  # 输出
  out_path <- file.path(resdir, sprintf("mpileup_fixstrand_%s.vcf", strand))
  message("[INFO] 输出路径: ", out_path)

  ok <- tryCatch({
    write.table(df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    TRUE
  }, error = function(e) {
    message("[WARN] 直接写出失败：", conditionMessage(e))
    FALSE
  })

  if (!ok) {
    tmp <- tempfile(fileext = ".vcf")
    write.table(df, file = tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    message("[HINT] 已写入临时文件：", tmp)
    message("[HINT] 请手动复制到目标目录，例如：")
    message("       cp ", shQuote(tmp), " ", shQuote(out_path))
  } else {
    message("[OK] 完成：", out_path, "  （共 ", nrow(df), " 行）")
  }
}

main()
