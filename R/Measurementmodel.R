# 02_measurement_model.R -------------------------------------------------
# Phase 1b: validate the five-construct structure.
# Order: reliability -> EFA on half A -> CFA on half B -> final CFA on all data
#        -> AVE / CR -> discriminant validity -> modification indices.
source("R/Setup.R")

dat <- readRDS("data/clean_indonesia.rds")
cfa_5 <- make_measurement(scales)   # the original five-construct model

# 1. Reliability of the ORIGINAL scales ------------------------------------
alpha_tbl <- imap_dfr(scales, ~ {
  a <- psych::alpha(dat[.x], check.keys = FALSE)   # no auto reverse-keying
  tibble(factor = .y, k = length(.x),
         alpha = a$total$raw_alpha, alpha_std = a$total$std.alpha)
})
print(alpha_tbl)
write_csv(alpha_tbl, "output/alpha_original_scales.csv")

# 2. Split-half: explore on A, confirm on B ---------------------------------
set.seed(2026)
idx    <- sample(seq_len(nrow(dat)), size = floor(nrow(dat) / 2))
half_a <- dat[idx, ]
half_b <- dat[-idx, ]

# 3. EFA on half A (polychoric correlations because items are 4-point) -------
pc <- psych::polychoric(half_a[items])$rho

fa.parallel(pc, n.obs = nrow(half_a), fa = "fa")   # how many factors?

efa_list <- map(4:7, ~ fa(pc, nfactors = .x, n.obs = nrow(half_a),
                          rotate = "oblimin", fm = "minres"))
names(efa_list) <- paste0("k", 4:7)

# Compare solutions: want TLI >= .90, RMSEA <= .08, interpretable factors
efa_fit <- imap_dfr(efa_list, ~ tibble(solution = .y, TLI = .x$TLI,
                                       RMSEA = .x$RMSEA[1], BIC = .x$BIC))
print(efa_fit)
write_csv(efa_fit, "output/efa_solution_comparison.csv")

# READ THE LOADINGS: do the CU and AUP items split into support / risk /
# dependency factors? Compare the wording of every item with the questionnaire.
print(efa_list$k5$loadings, cutoff = 0.30)
print(efa_list$k6$loadings, cutoff = 0.30)

# 4. CFA on half B: does the ORIGINAL structure hold on unseen data? ---------
fit_b <- cfa(cfa_5, data = half_b, ordered = items, estimator = "WLSMV")
print(fit_tab(fit_b, "original 5-factor, half B"))

# 5. Final CFA on the full sample (this is the model you report) ------------
fit_cfa <- cfa(cfa_5, data = dat, ordered = items, estimator = "WLSMV")
fit_summary <- fit_tab(fit_cfa, "original 5-factor, full sample")
print(fit_summary)
# Rough targets: CFI/TLI >= .90 (ideally .95), RMSEA <= .08, SRMR <= .08.
# A 20-item PE factor often fits poorly; see the optional parcel block below.

# 6. Loadings, AVE, composite reliability -----------------------------------
loadings_tbl <- standardizedSolution(fit_cfa) %>%
  as_tibble() %>%
  filter(op == "=~") %>%
  transmute(factor = lhs, item = rhs, loading = est.std, se, pvalue,
            weak = est.std < 0.50)
print(loadings_tbl %>% mutate(loading = round(loading, 3)), n = 39)

ave_cr <- loadings_tbl %>%
  group_by(factor) %>%
  summarise(k   = n(),
            AVE = mean(loading^2),
            CR  = sum(loading)^2 / (sum(loading)^2 + sum(1 - loading^2)),
            .groups = "drop") %>%
  left_join(alpha_tbl %>% select(factor, alpha), by = "factor") %>%
  mutate(across(c(AVE, CR, alpha), ~ round(.x, 3)))
print(ave_cr)          # targets: AVE > .50, CR > .70, alpha > .70

# Ordinal-scale omega from semTools (needs a recent semTools)
tryCatch(print(semTools::compRelSEM(fit_cfa)),
         error = function(e) message("compRelSEM unavailable: update semTools"))

# 7. Discriminant validity ---------------------------------------------------
# Fornell-Larcker: sqrt(AVE) on the diagonal must exceed the correlations
lat <- lavInspect(fit_cfa, "cor.lv")
fl  <- lat
diag(fl) <- sqrt(ave_cr$AVE[match(rownames(lat), ave_cr$factor)])
print(round(fl, 3))

# HTMT: want < .85 (strict) or < .90 (lenient). Not "between -1 and 1".
htmt_res <- tryCatch(semTools::htmt(cfa_5, dat, missing = "listwise"),
                     error = function(e) { message("htmt failed: ", e$message); NULL })
if (!is.null(htmt_res)) print(round(htmt_res, 3))

# 8. Modification indices (inspect, do not free blindly) ---------------------
mi <- modindices(fit_cfa, sort. = TRUE, maximum.number = 30) %>%
  as_tibble() %>%
  filter(op == "~~", mi > 10)
print(mi)
# Residual covariances between items with near-identical wording can be
# justified; cross-loadings and dropping items need a content-based reason.

# 9. ALTERNATIVE structure (template: fill in after reading the questionnaire)
# Only split a construct if item CONTENT supports it (e.g. connectivity and
# guidance items = institutional support; privacy and dependency worries = risk).
# scales_alt <- list(
#   PE   = scales$PE,
#   SUP  = c("CU1", "CU2"),            # <- EDIT: support-type items
#   BARR = c("CU3", "CU4"),            # <- EDIT: remaining challenge items
#   ATU  = scales$ATU,
#   RISK = c("AUP1", "AUP3"),          # <- EDIT: privacy / dependency worry items
#   AUP  = c("AUP2", "AUP4", "AUP5"),  # <- EDIT: remaining perception items
#   MIUA = scales$MIUA
# )
# fit_alt <- cfa(make_measurement(scales_alt), data = dat,
#                ordered = items, estimator = "WLSMV")
# fit_summary <- bind_rows(fit_summary, fit_tab(fit_alt, "alternative structure"))
# Models with different factor structures are usually NOT nested, so compare
# CFI / RMSEA / SRMR and interpretability, not a chi-square difference test.
# Run the EFA / half-B CFA logic on the alternative as well: do not choose a
# structure only because it fits the full sample best.

# 10. OPTIONAL: PE parcels if the 20-item factor fits badly ------------------
# Only valid if PE is unidimensional. Parcels are item means, not ordinal items.
# dat <- dat %>% mutate(
#   PEp1 = rowMeans(across(PE1:PE4),   na.rm = TRUE),
#   PEp2 = rowMeans(across(PE5:PE8),   na.rm = TRUE),
#   PEp3 = rowMeans(across(PE9:PE12),  na.rm = TRUE),
#   PEp4 = rowMeans(across(PE13:PE16), na.rm = TRUE),
#   PEp5 = rowMeans(across(PE17:PE20), na.rm = TRUE))
# then refit with PE =~ PEp1 + ... + PEp5 and drop PE1:PE20 from `ordered`.

# 11. Save outputs for later phases and the Shiny app ------------------------
write_csv(fit_summary, "output/cfa_fit.csv")
write_csv(loadings_tbl, "output/cfa_loadings.csv")
write_csv(ave_cr, "output/reliability_validity.csv")
saveRDS(fit_cfa, "output/fit_cfa.rds")