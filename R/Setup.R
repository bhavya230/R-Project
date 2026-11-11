# 00_setup.R -------------------------------------------------------------
# Shared by every script. Run from the RStudio Project root.
# First time only:
#install.packages(c("tidyverse","readxl","lavaan","semTools","psych","naniar"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lavaan)
  library(semTools)
  library(psych)
  library(naniar)
})

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Construct -> item map (edit here ONLY; every script reads it).
# If the questionnaire shows items were mis-assigned, change them here.
scales <- list(
  PE   = paste0("PE",   1:20),  # performance expectancy
  CU   = paste0("CU",   1:4),   # challenges using AI
  ATU  = paste0("ATU",  1:5),   # attitudes toward using AI
  AUP  = paste0("AUP",  1:5),   # perception of using AI
  MIUA = paste0("MIUA", 1:5)    # motivation / intention to keep using AI
)
items <- unlist(scales, use.names = FALSE)

# Build lavaan measurement syntax from a named list of item vectors
make_measurement <- function(sc) {
  paste(imap_chr(sc, ~ paste0(.y, " =~ ", paste(.x, collapse = " + "))),
        collapse = "\n")
}

# One-row table of the main fit indices (works for WLSMV fits)
fit_tab <- function(fit, label) {
  x <- unclass(fitMeasures(fit, c(
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.scaled", "tli.scaled",
    "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
    "srmr"
  )))
  bind_cols(tibble(model = label), as_tibble(as.list(round(x, 3))))
}