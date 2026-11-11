# 03_structural_model.R --------------------------------------------------
# Phase 2: corrected structural model (constructs -> attitude -> motivation).
# Demographics are COVARIATES, not outcomes (the reversed direction in the
# original paper).
#
# If Phase 1 changed the construct structure, update `scales` in 00_setup.R and
# the structural syntax below.
#
# Hypotheses (H8 and H9 each cover three parameters):
#  H1  PE   -> ATU    (+)    H5  PE  -> MIUA, direct  (+)
#  H2  CU   -> ATU    (-)    H6  CU  -> MIUA, direct  (-)
#  H3  AUP  -> ATU    (?)    H7  AUP -> MIUA, direct  (?)
#  H4  ATU  -> MIUA   (+)
#  H8  indirect effects of PE, CU, AUP on MIUA through ATU
#  H9  demographics (female, private, S1) -> MIUA
# The sign for AUP is "?" until the questionnaire shows how its items are
# worded / reverse-coded.
source("R/Setup.R")

dat  <- readRDS("data/clean_indonesia.rds")
meas <- make_measurement(scales)

# 1. Models -----------------------------------------------------------------
struct_full <- '
  ATU  ~ a1*PE + a2*CU + a3*AUP + d1*female + d2*private + d3*s1
  MIUA ~ c1*PE + c2*CU + c3*AUP + b*ATU + g1*female + g2*private + g3*s1

  ind_PE  := a1*b
  ind_CU  := a2*b
  ind_AUP := a3*b
  tot_PE  := c1 + a1*b
  tot_CU  := c2 + a2*b
  tot_AUP := c3 + a3*b
'

# Full-mediation alternative: PE, CU, AUP act on MIUA only through ATU
struct_med <- '
  ATU  ~ a1*PE + a2*CU + a3*AUP + d1*female + d2*private + d3*s1
  MIUA ~ b*ATU + g1*female + g2*private + g3*s1

  ind_PE  := a1*b
  ind_CU  := a2*b
  ind_AUP := a3*b
'

model_full <- paste(meas, struct_full, sep = "\n")
model_med  <- paste(meas, struct_med,  sep = "\n")

# 2. Primary estimation: ordinal items, WLSMV -------------------------------
fit_full <- sem(model_full, data = dat, ordered = items, estimator = "WLSMV")
fit_med  <- sem(model_med,  data = dat, ordered = items, estimator = "WLSMV")

fit_tbl <- bind_rows(fit_tab(fit_full, "partial mediation (full model)"),
                     fit_tab(fit_med,  "full mediation"))
print(fit_tbl)

# Scaled chi-square difference test (are the direct paths needed?)
print(lavTestLRT(fit_full, fit_med))

# Variance explained
print(lavInspect(fit_full, "rsquare"))

# Collinearity check: latent correlations above ~.85 are a warning sign
lat <- lavInspect(fit_full, "cor.lv")
print(round(lat, 3))
if (any(abs(lat[lower.tri(lat)]) > 0.85)) warning("Latent correlation > .85")

# 3. Hypothesis table --------------------------------------------------------
hyp <- tribble(
  ~label,    ~hypothesis,                 ~expected,
  "a1",      "H1: PE -> ATU",             "+",
  "a2",      "H2: CU -> ATU",             "-",
  "a3",      "H3: AUP -> ATU",            "?",
  "b",       "H4: ATU -> MIUA",           "+",
  "c1",      "H5: PE -> MIUA (direct)",   "+",
  "c2",      "H6: CU -> MIUA (direct)",   "-",
  "c3",      "H7: AUP -> MIUA (direct)",  "?",
  "ind_PE",  "H8a: PE -> ATU -> MIUA",    "+",
  "ind_CU",  "H8b: CU -> ATU -> MIUA",    "-",
  "ind_AUP", "H8c: AUP -> ATU -> MIUA",   "?",
  "g1",      "H9a: female -> MIUA",       "?",
  "g2",      "H9b: private -> MIUA",      "?",
  "g3",      "H9c: S1 -> MIUA",           "?"
)

pe <- parameterEstimates(fit_full, standardized = TRUE, ci = TRUE) %>%
  as_tibble()
std <- standardizedSolution(fit_full, type = "std.all") %>%
  as_tibble() %>%
  select(lhs, op, rhs, std_all = est.std)

hyp_tbl <- pe %>%
  filter(label %in% hyp$label) %>%
  select(label, lhs, op, rhs, est, se, z, pvalue, ci.lower, ci.upper) %>%
  left_join(hyp, by = "label") %>%
  left_join(std, by = c("lhs", "op", "rhs")) %>%
  mutate(supported = case_when(
    expected == "+" ~ est > 0 & pvalue < .05,
    expected == "-" ~ est < 0 & pvalue < .05,
    TRUE            ~ pvalue < .05      # non-directional: just significant
  )) %>%
  arrange(match(label, hyp$label)) %>%
  select(hypothesis, expected, est, se, z, pvalue, ci.lower, ci.upper,
         std_all, supported)

print(hyp_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 20)

# Total effects (direct + indirect) for the write-up
print(pe %>% filter(label %in% c("tot_PE", "tot_CU", "tot_AUP")) %>%
        select(label, est, se, pvalue, ci.lower, ci.upper))

# 4. Robustness A: continuous treatment with MLR + FIML ----------------------
fit_mlr <- sem(model_full, data = dat, estimator = "MLR", missing = "fiml")

path_compare <- standardizedSolution(fit_full, type = "std.all") %>%
  as_tibble() %>%
  filter(op == "~") %>%
  select(lhs, rhs, wlsmv = est.std) %>%
  left_join(standardizedSolution(fit_mlr, type = "std.all") %>%
              as_tibble() %>%
              filter(op == "~") %>% select(lhs, rhs, mlr = est.std),
            by = c("lhs", "rhs")) %>%
  mutate(diff = round(wlsmv - mlr, 3))
print(path_compare)
# The signs and rough sizes should agree. If they do not, report the
# disagreement instead of hiding it.

# 5. Robustness B: bootstrapped indirect effects (slow: run overnight) -------
# Test first with bootstrap = 200, then raise to 2000.
run_boot <- FALSE
if (run_boot) {
  set.seed(2026)
  fit_boot <- sem(model_full, data = dat, estimator = "ML", missing = "fiml",
                  se = "bootstrap", bootstrap = 2000,
                  parallel = if (.Platform$OS.type == "windows") "snow" else "multicore",
                  ncpus = max(1, parallel::detectCores() - 1))
  boot_tbl <- parameterEstimates(fit_boot, boot.ci.type = "perc", ci = TRUE) %>%
    filter(op == ":=") %>%
    select(label, est, se, ci.lower, ci.upper)
  print(boot_tbl)
  write_csv(boot_tbl, "output/indirect_effects_bootstrap.csv")
}

# 6. Save outputs for Phase 3 and the Shiny app ------------------------------
write_csv(hyp_tbl,   "output/hypothesis_table.csv")
write_csv(fit_tbl,   "output/structural_fit.csv")
write_csv(path_compare, "output/wlsmv_vs_mlr_paths.csv")
saveRDS(fit_full,  "output/fit_structural_full.rds")
saveRDS(model_full, "output/model_full_syntax.rds")   # reused for multi-group