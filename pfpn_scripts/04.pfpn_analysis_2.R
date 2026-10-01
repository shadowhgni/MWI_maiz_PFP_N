# =============================================================================
# PFP-N in Malawi, analysis 2: response curves, random forests, space and time
# -----------------------------------------------------------------------------
# Input : pfpn_compiled/pfpn_unified.rds (pfpn_compile.R)
# Output: pfpn_analysis_2/{tables,figures,objects}/ and, if render_report is TRUE,
#         pfpn_analysis_2/pfpn_report.html + pfpn_report.docx (from pfpn_report.Rmd)
#
#  1. Wave coverage funnel: rows of every wave through every filtering step,
#     so a wave that disappears shows where and why
#  2. Data quality by wave and N window (10-300 kg N/ha by default)
#  3. Distributions (linear scales)
#  4. PFP-N response to N rate, Eq 02 and an exponential alternative:
#       hyperbolic   PFP = A_min + A_max / (N + N_0)
#       exponential  PFP = A_min + A_max * exp(-N / N_0)
#     fitted pooled, by wave, by region and for trials; compared by AIC and
#     cross-validated error; parameters with EA-cluster bootstrap CIs
#  5. Environmental covariates: user GeoTIFFs, elevation, WorldClim, SoilGrids,
#     CHIRPS seasonal rainfall per wave (all optional, cached)
#  6. Random forests (caret::train, method "ranger") for log PFP-N with three
#     nested predictor sets (survey / + environment / + spatial signature) and
#     three cross-validation schemes (random k-fold, spatial blocks, leave one
#     wave out); variable importance, partial dependence, prediction maps
#  7. PFP-N by wave and region (all waves shown)
#  8. Results object and report
#
# Provenance tags: [verified] checked against a source; [recall] from training
# knowledge, check before relying on it.
# Packages: tidyverse; via pkg::fun caret, ranger, iml (SHAP), multcompView
# (letters), sf, terra, geodata (optional), rmarkdown + knitr (report).
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (caret, ranger, iml, emmeans, multcomp, MatchIt, sf, terra, scales, geodata, rmarkdown, knitr)

source("00.pfpn_utils.R", local = TRUE)   # paths, update-aware downloads, Malawi check, compact letter display
options(warn = 1)   # print warnings where they occur

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------

in_rds  <- pfpn_path("pfpn_compiled", "pfpn_unified.rds")
out_dir <- pfpn_path("pfpn_analysis_2")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
obj_dir <- file.path(out_dir, "objects")
walk(c(fig_dir, tab_dir, obj_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

focus_crop <- "maize"
N_min      <- 10     # kg N/ha: below this PFP-N is unstable
N_max      <- 300    # kg N/ha: above this rates are implausible for Malawian smallholders
mad_k      <- 4      # robust outlier threshold (median/MAD on log scale, per wave)
k_folds    <- 10     # random and curve-fit cross-validation folds
n_blocks   <- 10     # spatial blocks (k-means on EA coordinates) for spatial CV
n_boot     <- 200    # EA-cluster bootstrap replicates for curve parameters
n_trees    <- 500    # trees per random forest
shap_nsim    <- 20   # Monte Carlo samples per feature for SHAP values of ALL plots
shap_check_n <- 30   # plots on which the SHAP values are cross-checked with iml::Shapley
match_km     <- 30   # trial zone: survey plots within this distance of a trial site ...
match_years  <- 2    # ... and within this many years of the trial
psm_caliper  <- 0.2  # propensity-score caliper (in SD of the logit score)
N_ref      <- 50     # kg N/ha used for prediction maps
set.seed(20260930)

# Environmental covariates
covariate_dir       <- pfpn_path("covariates")         # any GeoTIFF here becomes a predictor (name = file stem)
download_covariates <- TRUE                            # elevation, WorldClim, SoilGrids via geodata
download_chirps     <- TRUE                            # seasonal rainfall per wave
cov_cache           <- pfpn_path("pfpn_compiled", "raw", "covariates")  # downloads cached and reused
# Malawi bounding box mwi_bbox comes from pfpn_utils.R
rain_months         <- c(11, 12, 1, 2, 3, 4)           # rainy season Nov-Apr

# Reference rainy season of each wave (season END year).
# IHPS 2013 asks about the 2012/2013 rainy season [verified: AG_MOD_G_13 labels];
# the others are the season before or during fieldwork [recall: check the
# questionnaires; households interviewed late may report the following season].
wave_season <- tribble(
  ~source,               ~season_end,
  "LSMS_1003_2010-2011", 2010,
  "LSMS_2248_2013",      2013,
  "LSMS_2936_2016-2017", 2016,
  "LSMS_3818_2019-2020", 2019,
  "LSMS_8507_2024-2025", 2024
)

render_report <- TRUE
report_rmd    <- "pfpn_report.Rmd"   # next to the scripts
analysis1_rds <- pfpn_path("pfpn_analysis", "objects", "analysis1.rds")   # used by the report if present

theme_set(theme_minimal(base_size = 11))
save_table <- \(x, name) write_csv(x, file.path(tab_dir, paste0(name, ".csv")), na = "")
save_fig   <- \(p, name, w = 8, h = 5) save_fig_png(p, name, w, h)
save_fig_png <- \(p, name, w, h) {
  save_plot_data(p, name, tab_dir)   # the data behind every figure -> tables/fig_<name>.csv
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 200, bg = "white")
}
dir.create(cov_cache, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data and wave coverage funnel
# -----------------------------------------------------------------------------

wave_labels <- c("LSMS_1003_2010-2011" = "IHS3 2010/11", "LSMS_2248_2013" = "IHPS 2013",
                 "LSMS_2936_2016-2017" = "IHS4 2016/17", "LSMS_3818_2019-2020" = "IHS5 2019/20",
                 "LSMS_8507_2024-2025" = "IHS6 2024/25")

# Districts by region [recall: administrative regions of Malawi]
region_of <- c(
  chitipa = "Northern", karonga = "Northern", likoma = "Northern", mzimba = "Northern",
  nkhatabay = "Northern", rumphi = "Northern",
  dedza = "Central", dowa = "Central", kasungu = "Central", lilongwe = "Central", mchinji = "Central",
  nkhotakota = "Central", ntcheu = "Central", ntchisi = "Central", salima = "Central",
  balaka = "Southern", blantyre = "Southern", chikwawa = "Southern", chiradzulu = "Southern",
  machinga = "Southern", mangochi = "Southern", mulanje = "Southern", mwanza = "Southern",
  neno = "Southern", nsanje = "Southern", phalombe = "Southern", thyolo = "Southern", zomba = "Southern"
)

robust_z <- function(x) {
  s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - median(x, na.rm = TRUE)) / s
}

if (!file.exists(in_rds)) stop("Run pfpn_compile.R first: ", in_rds, " not found", call. = FALSE)
d0 <- readRDS(in_rds) |>
  # coordinates outside Malawi (e.g. 0, 0) are errors: discarded
  mutate(latitude  = if_else(in_malawi(longitude, latitude), latitude, NA_real_),
         longitude = if_else(is.na(latitude), NA_real_, longitude))

prep <- function(x) {
  x |>
    mutate(
      wave = factor(unname(wave_labels[source]), levels = wave_labels),
      cluster = paste(source, coalesce(as.character(ea_id), hhid)),
      dkey = district_key(adm2),
      region = factor(coalesce(unname(region_of[dkey]),
                               str_to_title(str_remove(str_to_lower(adm1), " region"))),
                      levels = c("Northern", "Central", "Southern")),
      N_window = between(N_fertilizer, N_min, N_max),
      log_pfp = if_else(pfp_n > 0, log(pfp_n), NA_real_)
    )
}

lsms_all <- d0 |> filter(program == "LSMS_MWI") |> prep()

lsms_maize <- lsms_all |>
  filter(crop == focus_crop) |>
  mutate(z = robust_z(log_pfp), .by = source) |>
  mutate(outlier = (abs(z) > mad_k) %in% TRUE)

funnel_steps <- list(
  "1 compiled (all crops)"        = \(x) x,
  "2 focus crop"                  = \(x) filter(x, crop == focus_crop),
  "3 yield available"             = \(x) filter(x, crop == focus_crop, !is.na(yield)),
  "4 fertilized (N > 0)"          = \(x) filter(x, crop == focus_crop, !is.na(yield), N_fertilizer > 0),
  "5 N within window"             = \(x) filter(x, crop == focus_crop, !is.na(yield), N_window),
  "6 not an outlier"              = \(x) filter(x, crop == focus_crop, !is.na(yield), N_window, !outlier),
  "7 georeferenced"               = \(x) filter(x, crop == focus_crop, !is.na(yield), N_window, !outlier,
                                                !is.na(latitude))
)
funnel_input <- lsms_all |>
  left_join(select(lsms_maize, source, hhid, plot_id, yield_source, outlier),
            by = c("source", "hhid", "plot_id", "yield_source")) |>
  mutate(outlier = coalesce(outlier, FALSE))
funnel <- imap(funnel_steps, \(f, nm) f(funnel_input) |> count(wave, .drop = FALSE) |> mutate(step = nm)) |>
  list_rbind() |>
  complete(wave = factor(wave_labels, levels = wave_labels), step, fill = list(n = 0)) |>
  filter(!is.na(wave)) |>
  mutate(step = factor(step, levels = names(funnel_steps)))
save_table(pivot_wider(funnel, names_from = step, values_from = n), "1_wave_coverage_funnel")

missing_waves <- funnel |> filter(step == names(funnel_steps)[1], n == 0) |> pull(wave) |> as.character()
dropped_waves <- funnel |>
  mutate(first_zero = n == 0) |>
  filter(first_zero) |>
  slice_min(as.integer(step), by = wave) |>
  filter(!wave %in% missing_waves)
if (length(missing_waves) > 0) {
  message("Waves absent from ", in_rds, ": ", paste(missing_waves, collapse = ", "),
          ". Check the zip paths in pfpn_compile.R (LSMS zip not found) and TODO rows in source_map.csv.")
}

p <- ggplot(funnel, aes(step, fct_rev(wave), fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = scales::comma(n)), size = 3) +
  scale_fill_gradient(low = "#fde0dd", high = "#3182bd", name = "rows") +
  labs(x = NULL, y = NULL, title = "Rows of each LSMS wave through each filtering step",
       subtitle = "A 0 shows where a wave is lost") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
save_fig(p, "1_wave_coverage_funnel", w = 10, h = 4)

# Analysis sets
surv <- lsms_maize |>
  filter(!is.na(pfp_n), N_window, !outlier) |>
  left_join(wave_season, by = "source") |>
  mutate(season_end = coalesce(season_end, suppressWarnings(as.integer(str_sub(date, 1, 4)))))

trials <- d0 |>
  filter(program == "carob", crop == focus_crop, !is.na(pfp_n), !is_survey %in% TRUE) |>
  prep() |>
  mutate(trial_type = if_else(on_farm %in% TRUE, "On-farm trials",
                              if_else(on_farm %in% FALSE, "On-station trials", "Trials (setting unknown)")),
         season_end = suppressWarnings(as.integer(str_sub(date, 1, 4)))) |>
  filter(N_window)

message(nrow(surv), " survey rows and ", nrow(trials), " trial rows for the PFP-N analysis")
if (nrow(surv) < 100) stop("Too few survey rows after filtering; see ", file.path(tab_dir, "1_wave_coverage_funnel.csv"))

# -----------------------------------------------------------------------------
# 2. Data quality by wave, including the N window
# -----------------------------------------------------------------------------

has_flag <- \(x, f) str_detect(coalesce(x, ""), fixed(f))
pct_among <- \(cond, among) if (sum(among %in% TRUE) == 0) NA_real_ else 100 * mean(cond[among %in% TRUE] %in% TRUE)
quality <- lsms_maize |>
  filter(!is.na(yield)) |>
  mutate(fert = N_fertilizer > 0) |>
  summarise(
    plots                  = n(),
    fertilized             = sum(fert %in% TRUE),
    pct_area_reported      = 100 * mean(area_source %in% "farmer_report"),
    pct_fert_enumerator_kg = pct_among(has_flag(fert_kg_source, "enumerator"), fert),
    pct_fert_kg_mismatch   = pct_among(fert_kg_mismatch, fert),
    pct_N_unknown          = 100 * mean(is.na(N_fertilizer) & fertilizer_used %in% TRUE),
    pct_below_N_min        = pct_among(N_fertilizer < N_min, fert),
    pct_above_N_max        = pct_among(N_fertilizer > N_max, fert),
    pct_outlier            = 100 * mean(outlier),
    pct_second_harvest_blk = 100 * mean(harvest_block %in% "alt"),
    pct_district_centroid  = 100 * mean(geo_level %in% "district_centroid"),
    .by = wave
  ) |>
  complete(wave = factor(wave_labels, levels = wave_labels)) |>
  arrange(wave) |>
  mutate(across(starts_with("pct"), \(x) round(x, 1)))
save_table(quality, "2_quality_by_wave")

# -----------------------------------------------------------------------------
# 3. Distributions (linear scales; axis cut at the 99th percentile, count shown)
# -----------------------------------------------------------------------------

dist_all <- bind_rows(
  surv |> transmute(group = as.character(wave), class = "Survey (LSMS)", pfp_n, yield, N_fertilizer),
  trials |> transmute(group = trial_type, class = "Trials", pfp_n, yield, N_fertilizer)
) |>
  mutate(group = factor(group, levels = c(wave_labels, sort(unique(trials$trial_type)))))

# Letters: one-way ANOVA + emmeans + multcomp::cld (Bonferroni; "a" = highest
# mean) drawn on the figure; pairwise Wilcoxon letters (Bonferroni) kept as a
# table for comparison. Both treat plots as independent, so with thousands of
# plots small differences come out significant.
violin_box <- function(df, var, lab, name = paste0("3_", var)) {
  cap <- quantile(df[[var]], 0.99, na.rm = TRUE)
  above <- sum(df[[var]] > cap, na.rm = TRUE)
  lt_anova <- anova_cld_letters(df[[var]], df$group)
  lt_wilcox <- cld_letters(df[[var]], df$group)
  save_table(lt_anova, paste0(name, "_letters_anova_emmeans"))
  save_table(attr(lt_anova, "anova"), paste0(name, "_anova_F_test"))
  save_table(lt_wilcox, paste0(name, "_letters_wilcoxon"))
  lt <- lt_anova |> mutate(group = factor(group, levels = levels(df$group)))
  ggplot(df, aes(group, .data[[var]], fill = class)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_point(data = lt, aes(group, emmean), inherit.aes = FALSE, shape = 23, size = 2, fill = "black") +
    geom_text(data = lt, aes(group, cap * 0.97, label = letter), inherit.aes = FALSE,
              fontface = "bold", size = 4) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(0, cap)) +
    labs(x = NULL, y = lab, fill = NULL,
         caption = paste0("Letters: one-way ANOVA + emmeans, Bonferroni (alpha 0.05); 'a' = highest mean; shared letter = no difference.\n",
                          "Diamonds: means. Wilcoxon letters: tables/", name, "_letters_wilcoxon.csv.\n",
                          "Axis cut at the 99th percentile (", round(cap), "); ", above, " values above not shown.")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
save_fig(violin_box(dist_all, "pfp_n", "PFP-N (kg grain / kg N)") +
           ggtitle(paste0("PFP-N by wave and trial type, sole ", focus_crop, ", N ", N_min, "-", N_max, " kg/ha")),
         "3_pfp_violin_box", w = 10)
save_fig(violin_box(dist_all, "yield", "Yield (kg/ha)") + ggtitle("Yield"), "3_yield_violin_box", w = 10)
save_fig(violin_box(dist_all, "N_fertilizer", "N rate (kg/ha)") + ggtitle("N rate"), "3_N_violin_box", w = 10)

# -----------------------------------------------------------------------------
# 4. PFP-N response to N: hyperbolic (Eq 02) vs exponential
# -----------------------------------------------------------------------------
# Fitted on the log scale (multiplicative errors: PFP-N is right-skewed and its
# spread shrinks as N rises). Note on interpretation of Eq 02: since
# yield = PFP-N x N, the hyperbolic form implies yield = A_min*N + A_max*N/(N + N_0):
# a linear term plus a saturating (Michaelis-Menten) term. A_min is the return
# that persists at high N, A_max/N_0 the PFP-N gain at very low N, and N_0 the N
# rate at which the saturating part reaches half its maximum.

curve_fns <- list(
  hyperbolic  = \(N, A_min, A_max, N_0) A_min + A_max / (N + N_0),
  exponential = \(N, A_min, A_max, N_0) A_min + A_max * exp(-N / N_0)
)
curve_starts <- list(
  hyperbolic  = expand_grid(A_min = c(1, 10, 25), A_max = c(300, 1500, 5000), N_0 = c(2, 20, 80)),
  exponential = expand_grid(A_min = c(5, 15, 30), A_max = c(20, 60, 150), N_0 = c(15, 50, 150))
)

fit_curve <- function(df, form, start = NULL) {
  f <- curve_fns[[form]]
  starts <- if (is.null(start)) curve_starts[[form]] else as_tibble(as.list(start))
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    fit <- tryCatch(
      nls(log(pfp_n) ~ log(f(N_fertilizer, A_min, A_max, N_0)), data = df,
          start = as.list(starts[i, ]), algorithm = "port",
          lower = c(A_min = 0, A_max = 1e-6, N_0 = 1e-3),
          control = nls.control(maxiter = 300, warnOnly = FALSE)),
      error = \(e) NULL)
    if (!is.null(fit) && (is.null(best) || deviance(fit) < deviance(best))) best <- fit
  }
  best
}

predict_curve <- \(fit, form, N) { p <- unname(coef(fit)[c("A_min", "A_max", "N_0")]); curve_fns[[form]](N, p[1], p[2], p[3]) }

# Marginal yield response implied by each form (kg grain per extra kg N), from
# yield = PFP-N x N:
#   hyperbolic   dY/dN = A_min + A_max * N_0 / (N + N_0)^2
#   exponential  dY/dN = A_min + A_max * exp(-N / N_0) * (1 - N / N_0)
# This is the quantity most Malawi studies report (e.g. Burke et al. 2020), so
# it makes our curves comparable with theirs. Cross-sectional curves are not
# causal: plots with more N differ in other ways too.
mp_fns <- list(
  hyperbolic  = \(N, A_min, A_max, N_0) A_min + A_max * N_0 / (N + N_0)^2,
  exponential = \(N, A_min, A_max, N_0) A_min + A_max * exp(-N / N_0) * (1 - N / N_0)
)
marginal_curve <- \(fit, form, N) { p <- unname(coef(fit)[c("A_min", "A_max", "N_0")]); mp_fns[[form]](N, p[1], p[2], p[3]) }
derived_at <- \(fit, form) c(
  PFP_25 = predict_curve(fit, form, 25), PFP_50 = predict_curve(fit, form, 50),
  PFP_100 = predict_curve(fit, form, 100), PFP_150 = predict_curve(fit, form, 150),
  MP_25 = marginal_curve(fit, form, 25), MP_50 = marginal_curve(fit, form, 50),
  MP_100 = marginal_curve(fit, form, 100), MP_150 = marginal_curve(fit, form, 150))

# EA-grouped k-fold CV error (PFP-N scale) for one form on one data set
cv_curve <- function(df, form, k = k_folds) {
  eas <- unique(df$cluster)
  fold <- setNames(sample(rep_len(seq_len(k), length(eas))), eas)[df$cluster]
  errs <- map(seq_len(k), \(i) {
    fit <- fit_curve(df[fold != i, ], form)
    if (is.null(fit)) return(NULL)
    test <- df[fold == i, ]
    tibble(obs = test$pfp_n, pred = predict_curve(fit, form, test$N_fertilizer))
  }) |> list_rbind()
  if (nrow(errs) == 0) return(c(cv_rmse = NA, cv_mae = NA))
  c(cv_rmse = sqrt(mean((errs$obs - errs$pred)^2)), cv_mae = mean(abs(errs$obs - errs$pred)))
}

# EA-cluster bootstrap of parameters and of PFP-N at selected N rates
boot_curve <- function(df, form, fit, B) {
  idx <- split(seq_len(nrow(df)), df$cluster)
  map(seq_len(B), \(b) {
    s <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
    bf <- fit_curve(df[s, ], form, start = coef(fit))
    if (is.null(bf)) return(NULL)
    p <- coef(bf)
    bind_cols(tibble(A_min = unname(p["A_min"]), A_max = unname(p["A_max"]), N_0 = unname(p["N_0"])),
              as_tibble(as.list(derived_at(bf, form))))
  }) |> list_rbind()
}

curve_groups <- bind_rows(
  surv |> mutate(group_type = "Pooled", group = "All survey plots"),
  surv |> mutate(group_type = "Wave", group = as.character(wave)),
  surv |> filter(!is.na(region)) |> mutate(group_type = "Region", group = as.character(region)),
  trials |> mutate(group_type = "Trials", group = trial_type, cluster = paste(dataset_id, latitude, longitude))
)

curve_results <- curve_groups |>
  group_by(group_type, group) |>
  group_map(\(df, key) {
    if (nrow(df) < 30 || n_distinct(df$N_fertilizer) < 4) return(NULL)
    map(names(curve_fns), \(form) {
      fit <- fit_curve(df, form)
      if (is.null(fit)) return(NULL)
      cv <- cv_curve(df, form)
      bt <- boot_curve(df, form, fit, if (key$group_type == "Pooled") n_boot else ceiling(n_boot / 2))
      pred <- df$pfp_n - predict_curve(fit, form, df$N_fertilizer)
      list(key = mutate(key, form = form, n = nrow(df), clusters = n_distinct(df$cluster)),
           fit = fit, cv = cv, boot = bt, rmse = sqrt(mean(pred^2)))
    }) |> compact()
  }) |>
  flatten()

curve_comp <- map(curve_results, \(r) mutate(r$key, AIC = AIC(r$fit), rmse = r$rmse,
                                             cv_rmse = r$cv["cv_rmse"], cv_mae = r$cv["cv_mae"])) |>
  list_rbind() |>
  mutate(delta_AIC = AIC - min(AIC), .by = c(group_type, group)) |>
  mutate(preferred_AIC = delta_AIC == 0,
         preferred_CV = cv_rmse == min(cv_rmse, na.rm = TRUE), .by = c(group_type, group))
save_table(curve_comp, "4_curve_form_comparison")

ci <- \(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
curve_params <- map(curve_results, \(r) {
  est <- c(coef(r$fit), derived_at(r$fit, r$key$form))
  tibble(parameter = names(est), estimate = unname(est)) |>
    mutate(ci_low  = map_dbl(parameter, \(p) if (nrow(r$boot)) ci(r$boot[[p]])[1] else NA_real_),
           ci_high = map_dbl(parameter, \(p) if (nrow(r$boot)) ci(r$boot[[p]])[2] else NA_real_),
           boot_ok = nrow(r$boot)) |>
    bind_cols(select(r$key, group_type, group, form, n))
}) |> list_rbind() |>
  relocate(group_type, group, form, n)
save_table(curve_params, "4_curve_parameters")

# Figure: binned medians and both fitted curves per group
bins <- curve_groups |>
  mutate(N_bin = cut(N_fertilizer, seq(N_min, N_max + 10, by = 10), include.lowest = TRUE)) |>
  summarise(N_mid = median(N_fertilizer), pfp_median = median(pfp_n), n = n(),
            .by = c(group_type, group, N_bin)) |>
  filter(n >= 5)
curves <- map(curve_results, \(r) {
  N <- seq(N_min, N_max, length.out = 150)
  mutate(r$key, data = list(tibble(N_fertilizer = N, pfp = predict_curve(r$fit, r$key$form, N))))
}) |> list_rbind() |> unnest(data)

group_order <- c("All survey plots", wave_labels, "Northern", "Central", "Southern", sort(unique(trials$trial_type)))
bins   <- bins   |> mutate(group = factor(group, levels = unique(c(group_order, group))))
curves <- curves |> mutate(group = factor(group, levels = unique(c(group_order, group))))

for (gt in unique(curves$group_type)) {
  p <- ggplot() +
    geom_point(data = filter(bins, group_type == gt), aes(N_mid, pfp_median, size = n), colour = "grey45", alpha = 0.7) +
    geom_line(data = filter(curves, group_type == gt), aes(N_fertilizer, pfp, colour = form), linewidth = 0.9) +
    facet_wrap(~group) +
    scale_colour_manual(values = c(hyperbolic = "#d95f02", exponential = "#1b9e77")) +
    labs(x = "N rate (kg/ha)", y = "PFP-N (kg/kg)", colour = "Form", size = "plots in\n10-kg bin",
         title = paste("PFP-N response to N rate:", str_to_lower(gt)),
         subtitle = "Points: median PFP-N per 10 kg N bin; lines: fitted Eq 02 (hyperbolic) and exponential forms")
  save_fig(p, paste0("4_curves_", str_to_lower(gt)), w = if (gt == "Wave") 11 else 9, h = if (gt == "Wave") 6 else 4.5)
}

# -----------------------------------------------------------------------------
# 5. Environmental covariates (optional, cached)
# -----------------------------------------------------------------------------

usable_tif <- \(f) file.exists(f) && file.size(f) > 0 && !is.null(tryCatch(terra::rast(f), error = \(e) NULL))
mwi_ext <- terra::ext(mwi_bbox[["xmin"]], mwi_bbox[["xmax"]], mwi_bbox[["ymin"]], mwi_bbox[["ymax"]])
crop_mwi <- \(r) terra::crop(r, mwi_ext)

cached_layer <- function(name, loader) {
  f <- file.path(cov_cache, paste0(name, ".tif"))
  if (usable_tif(f)) return(terra::rast(f))
  r <- tryCatch(loader(), error = \(e) { message("Covariate ", name, " not available: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  r <- crop_mwi(r[[1]])
  names(r) <- name
  terra::writeRaster(r, f, overwrite = TRUE)
  terra::rast(f)
}

static_layers <- list()
user_tifs <- list.files(covariate_dir, "\\.tif$", full.names = TRUE)
for (f in user_tifs) {
  nm <- make.names(tools::file_path_sans_ext(basename(f)))
  static_layers[[nm]] <- crop_mwi(terra::rast(f)[[1]])
}
if (download_covariates && requireNamespace("geodata", quietly = TRUE)) {
  gpath <- file.path(cov_cache, "geodata")
  # geodata function names and variables [recall: check ?geodata if a layer is missing]
  static_layers$elevation <- cached_layer("elevation", \() geodata::elevation_30s(country = "MWI", path = gpath))
  bio <- tryCatch(geodata::worldclim_country(country = "MWI", var = "bio", path = gpath), error = \(e) NULL)
  if (!is.null(bio)) {
    static_layers$temp_mean <- cached_layer("temp_mean", \() bio[[grep("bio_1$", names(bio))]])
    static_layers$rain_clim <- cached_layer("rain_clim", \() bio[[grep("bio_12$", names(bio))]])
  }
  for (v in c("phh2o", "soc", "nitrogen", "clay")) {
    static_layers[[paste0("soil_", v)]] <- cached_layer(paste0("soil_", v),
      \() geodata::soil_world(var = v, depth = 5, path = gpath))
  }
}
static_layers <- compact(static_layers)

# CHIRPS monthly rainfall for Africa [recall: URL pattern of CHIRPS-2.0 africa_monthly]
chirps_url <- \(y, m) sprintf("https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs/chirps-v2.0.%d.%02d.tif.gz", y, m)
gunzip_file <- function(gz, out) {
  con <- gzfile(gz, "rb"); on.exit(close(con))
  o <- file(out, "wb"); on.exit(close(o), add = TRUE)
  repeat { b <- readBin(con, "raw", 1e7); if (length(b) == 0) break; writeBin(b, o) }
}
# The downloaded .tif.gz is kept (update-aware fetch); the Malawi crop is rebuilt
# whenever the download is newer than it.
chirps_month <- function(y, m) {
  gz  <- file.path(cov_cache, "chirps_raw", sprintf("chirps-v2.0.%d.%02d.tif.gz", y, m))
  out <- file.path(cov_cache, sprintf("chirps_%d_%02d.tif", y, m))
  ok <- tryCatch({ fetch(chirps_url(y, m), gz); TRUE },
                 error = \(e) { message("CHIRPS ", y, "-", m, " not available: ", conditionMessage(e)); FALSE })
  if (!ok) return(if (usable_tif(out)) terra::rast(out) else NULL)
  if (outdated(out, gz) || !usable_tif(out)) {
    tif <- tempfile(fileext = ".tif")
    gunzip_file(gz, tif)
    r <- crop_mwi(terra::rast(tif))
    r[r < 0] <- NA
    names(r) <- sprintf("chirps_%d_%02d", y, m)
    terra::writeRaster(r, out, overwrite = TRUE)
  }
  terra::rast(out)
}
season_rain <- function(end_year) {
  ym <- tibble(m = rain_months, y = if_else(rain_months >= 9, end_year - 1L, end_year))
  layers <- map2(ym$y, ym$m, chirps_month)
  if (any(map_lgl(layers, is.null))) return(NULL)
  terra::app(terra::rast(layers), sum)
}

seasons <- sort(unique(c(surv$season_end, trials$season_end)))
rain_layers <- if (download_chirps) {
  set_names(map(seasons, \(y) season_rain(as.integer(y))), seasons) |> compact()
} else list()

extract_at <- function(r, lon, lat) terra::extract(r, cbind(lon, lat))[, 1]

add_covariates <- function(df) {
  for (nm in names(static_layers)) df[[nm]] <- extract_at(static_layers[[nm]], df$longitude, df$latitude)
  if (length(rain_layers) > 0) {
    df$rain_season <- NA_real_
    for (y in names(rain_layers)) {
      i <- which(df$season_end == as.integer(y) & !is.na(df$latitude))
      if (length(i)) df$rain_season[i] <- extract_at(rain_layers[[y]], df$longitude[i], df$latitude[i])
    }
  }
  df
}

env_vars <- c(names(static_layers), if (length(rain_layers) > 0) "rain_season")
message("Environmental covariates available: ", if (length(env_vars)) paste(env_vars, collapse = ", ") else "none")

# -----------------------------------------------------------------------------
# 6. Random forests with caret (ranger), three CV schemes
# -----------------------------------------------------------------------------

# Spatial signature: coordinates plus oblique geographic coordinates (the
# coordinates rotated by several angles), which let trees split along
# directions other than north-south and east-west.
add_spatial_signature <- function(df, angles = c(30, 60, 120, 150)) {
  x <- df$longitude - mean(mwi_bbox[c("xmin", "xmax")])
  y <- df$latitude  - mean(mwi_bbox[c("ymin", "ymax")])
  for (a in angles) df[[paste0("ogc_", a)]] <- sqrt(x^2 + y^2) * cos(a * pi / 180 - atan2(y, x))
  df
}

rf_data <- surv |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  add_covariates() |>
  add_spatial_signature() |>
  mutate(
    N_rate        = N_fertilizer,
    year          = season_end,
    plot_area_ha  = plot_area / 1e4,
    share_planted = coalesce(intercrop_fraction, 1),
    manure        = as.integer(OM_used %in% TRUE),
    urea          = as.integer(str_detect(str_to_lower(coalesce(fertilizer_type, "")), "urea")),
    compound      = as.integer(str_detect(str_to_lower(coalesce(fertilizer_type, "")), "23:21|npk|compound|dap")),
    centroid      = as.integer(geo_level %in% "district_centroid")
  )

survey_vars  <- c("N_rate", "year", "plot_area_ha", "share_planted", "manure", "urea", "compound")
spatial_vars <- c("longitude", "latitude", paste0("ogc_", c(30, 60, 120, 150)))
env_ok <- env_vars[map_lgl(env_vars, \(v) mean(!is.na(rf_data[[v]])) > 0.8)]

predictor_sets <- list(
  "Survey"                        = survey_vars,
  "Survey + environment"          = c(survey_vars, env_ok),
  "Survey + environment + space"  = c(survey_vars, env_ok, spatial_vars)
)
if (length(env_ok) == 0) predictor_sets[["Survey + environment"]] <- NULL

# complete cases on all predictors used by the largest set (median-impute rare gaps)
all_vars <- unique(unlist(predictor_sets))
rf_data <- rf_data |>
  mutate(across(all_of(all_vars), \(x) if_else(is.na(x), median(x, na.rm = TRUE), as.numeric(x)))) |>
  filter(!is.na(log_pfp))
message("Random forest data: ", nrow(rf_data), " plots, ", n_distinct(rf_data$cluster), " EAs")

# CV folds (training-row indices, as caret expects)
group_folds <- function(groups, k) {
  ug <- unique(groups)
  f <- setNames(sample(rep_len(seq_len(k), length(ug))), ug)[groups]
  set_names(map(seq_len(k), \(i) which(f != i)), paste0("Fold", seq_len(k)))
}
ea_xy <- rf_data |> summarise(longitude = mean(longitude), latitude = mean(latitude), .by = cluster)
blocks <- kmeans(ea_xy[, c("longitude", "latitude")], centers = min(n_blocks, nrow(ea_xy) - 1), nstart = 25)$cluster
rf_data$block <- setNames(blocks, ea_xy$cluster)[rf_data$cluster]

cv_schemes <- list(
  "Random k-fold"        = caret::createFolds(rf_data$log_pfp, k = k_folds, returnTrain = TRUE),
  "Spatial blocks"       = group_folds(rf_data$block, n_distinct(rf_data$block)),
  "Leave one wave out"   = group_folds(as.character(rf_data$wave), n_distinct(rf_data$wave))
)
if (n_distinct(rf_data$wave) < 2) cv_schemes[["Leave one wave out"]] <- NULL

p <- ggplot(rf_data |> distinct(cluster, .keep_all = TRUE), aes(longitude, latitude, colour = factor(block))) +
  geom_point(size = 1) + coord_equal() +
  labs(colour = "block", title = "Spatial CV blocks (k-means on EA coordinates)", x = NULL, y = NULL)
save_fig(p, "6_spatial_cv_blocks", w = 5, h = 8)

train_rf <- function(vars, index, grid) {
  caret::train(
    x = as.data.frame(rf_data[, vars]), y = rf_data$log_pfp,
    method = "ranger", metric = "RMSE", tuneGrid = grid,
    trControl = caret::trainControl(method = "cv", index = index, savePredictions = "final"),
    num.trees = n_trees, importance = "permutation")
}

rf_runs <- imap(predictor_sets, \(vars, set_name) {
  message("Random forest: ", set_name, " (", length(vars), " predictors)")
  p <- length(vars)
  grid <- expand.grid(mtry = unique(pmax(1, round(c(1/3, 1/2, 2/3) * p))),
                      splitrule = "variance", min.node.size = c(5, 25))
  # tune on spatial blocks (the honest scheme), then evaluate every scheme with the chosen settings
  tuned <- train_rf(vars, cv_schemes[["Spatial blocks"]], grid)
  best <- tuned$bestTune
  evals <- imap(cv_schemes, \(idx, scheme) {
    fit <- if (scheme == "Spatial blocks") tuned else train_rf(vars, idx, best)
    pr <- fit$pred
    tibble(predictors = set_name, cv_scheme = scheme,
           rmse_log = sqrt(mean((pr$obs - pr$pred)^2)),
           r2_log = cor(pr$obs, pr$pred)^2,
           rmse_pfp = sqrt(mean((exp(pr$obs) - exp(pr$pred))^2)),
           mae_pfp = mean(abs(exp(pr$obs) - exp(pr$pred))),
           r2_pfp = cor(exp(pr$obs), exp(pr$pred))^2)
  }) |> list_rbind()
  list(set = set_name, vars = vars, model = tuned, best = best, evals = evals)
})

rf_metrics <- map(rf_runs, "evals") |> list_rbind() |>
  mutate(predictors = factor(predictors, levels = names(predictor_sets)),
         cv_scheme = factor(cv_scheme, levels = names(cv_schemes)))
save_table(rf_metrics, "6_rf_cv_metrics")
save_table(imap(rf_runs, \(r, nm) mutate(r$best, predictors = nm)) |> list_rbind(), "6_rf_tuning")

p <- ggplot(rf_metrics, aes(cv_scheme, r2_log, fill = predictors)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", r2_log)), position = position_dodge(0.8), vjust = -0.3, size = 3) +
  labs(x = NULL, y = expression(R^2~"(held-out, log PFP-N)"), fill = NULL,
       title = "Random forest skill by predictor set and cross-validation scheme",
       subtitle = "Random k-fold over-states skill when plots in the same place share information")
save_fig(p, "6_rf_cv_metrics")

final <- rf_runs[[length(rf_runs)]]
imp <- caret::varImp(final$model, scale = FALSE)$importance |>
  rownames_to_column("variable") |>
  rename(importance = Overall) |>
  mutate(group = case_when(variable %in% survey_vars ~ "survey", variable %in% spatial_vars ~ "space",
                           .default = "environment")) |>
  arrange(desc(importance))
save_table(imp, "6_rf_importance")
p <- ggplot(imp, aes(importance, fct_reorder(variable, importance), fill = group)) +
  geom_col() +
  labs(x = "Permutation importance (increase in MSE, log PFP-N)", y = NULL, fill = NULL,
       title = paste("Variable importance:", final$set))
save_fig(p, "6_rf_importance", h = 6)

# SHAP values for ALL plots (Monte Carlo Shapley, mc_shapley() in pfpn_utils.R).
# Each predictor's contribution to each plot's prediction on the log scale, so
# exp(SHAP) is its multiplicative effect on PFP-N. The SHAP-based partial
# dependence of a predictor is exp(mean prediction + its SHAP value), plotted
# against the predictor for every plot. A random subset is cross-checked
# against iml::Shapley.
X_all <- as.data.frame(rf_data[, final$vars])
pred_log <- \(d) predict(final$model, newdata = d)
message("SHAP values for all ", nrow(X_all), " plots (", shap_nsim, " Monte Carlo samples per predictor)")
phi <- mc_shapley(pred_log, X_all, nsim = shap_nsim)
base_log <- mean(pred_log(X_all))

shap <- as_tibble(phi) |>
  mutate(row = row_number()) |>
  pivot_longer(-row, names_to = "variable", values_to = "phi") |>
  mutate(value = map2_dbl(row, variable, \(i, v) X_all[i, v]))
save_table(mutate(shap, cluster = rf_data$cluster[row], wave = as.character(rf_data$wave[row])), "6_shap_values_all_plots")

shap_imp <- shap |>
  summarise(mean_abs_shap = mean(abs(phi)), .by = variable) |>
  arrange(desc(mean_abs_shap)) |>
  mutate(group = case_when(variable %in% survey_vars ~ "survey", variable %in% spatial_vars ~ "space",
                           .default = "environment"))
save_table(shap_imp, "6_shap_importance")

# cross-check with iml on a random subset
shap_check <- NULL
if (requireNamespace("iml", quietly = TRUE) && shap_check_n > 0) {
  predictor <- iml::Predictor$new(final$model, data = X_all, y = rf_data$log_pfp,
                                  predict.function = \(model, newdata) predict(model, newdata = newdata))
  rows <- sort(sample.int(nrow(X_all), min(shap_check_n, nrow(X_all))))
  iml_phi <- map(rows, \(i) {
    sh <- iml::Shapley$new(predictor, x.interest = X_all[i, , drop = FALSE], sample.size = 50)
    tibble(row = i, variable = sh$results$feature, phi_iml = sh$results$phi)
  }) |> list_rbind()
  shap_check <- iml_phi |> left_join(shap, by = c("row", "variable"))
  save_table(shap_check, "6_shap_check_vs_iml")
  shap_check <- tibble(plots = length(rows), correlation = cor(shap_check$phi, shap_check$phi_iml))
  message("SHAP cross-check with iml::Shapley on ", length(rows), " plots: r = ", round(shap_check$correlation, 3))
}

p <- shap |>
  mutate(scaled = (value - min(value)) / (max(value) - min(value) + 1e-12), .by = variable) |>
  mutate(variable = factor(variable, levels = rev(shap_imp$variable))) |>
  ggplot(aes(phi, variable, colour = scaled)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_jitter(height = 0.2, width = 0, size = 0.4, alpha = 0.4) +
  scale_colour_viridis_c(name = "predictor\nvalue", breaks = c(0, 1), labels = c("low", "high")) +
  labs(x = "SHAP value (contribution to log PFP-N)", y = NULL, title = "SHAP summary",
       subtitle = paste(nrow(X_all), "plots; exp(SHAP) = multiplicative effect on PFP-N"))
save_fig(p, "6_shap_summary", h = 6)

# SHAP-based partial dependence (all plots)
pd_vars <- unique(c("N_rate", "year", head(setdiff(shap_imp$variable, c("N_rate", "year", spatial_vars)), 4)))
shap_pd <- shap |>
  filter(variable %in% pd_vars) |>
  mutate(pfp = exp(base_log + phi), variable = factor(variable, levels = pd_vars))
save_table(shap_pd, "6_shap_partial_dependence")

# classic partial dependence kept as a table for comparison
pd <- map(pd_vars, \(v) {
  vals <- if (v == "year") sort(unique(rf_data$year)) else
    unique(quantile(rf_data[[v]], seq(0.02, 0.98, length.out = 25), na.rm = TRUE))
  X <- X_all
  tibble(variable = v, value = vals,
         pfp = map_dbl(vals, \(x) { X[[v]] <- x; exp(mean(predict(final$model, newdata = X))) }))
}) |> list_rbind()
save_table(pd, "6_rf_partial_dependence_classic")

pooled_curve <- curves |> filter(group_type == "Pooled")
p_pdN <- ggplot() +
  geom_point(data = filter(shap_pd, variable == "N_rate"), aes(value, pfp), alpha = 0.15, size = 0.5, colour = "grey30") +
  geom_smooth(data = filter(shap_pd, variable == "N_rate"), aes(value, pfp, linetype = "Random forest (SHAP-based)"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, colour = "black") +
  geom_line(data = pooled_curve, aes(N_fertilizer, pfp, colour = form), linewidth = 0.8) +
  scale_colour_manual(values = c(hyperbolic = "#d95f02", exponential = "#1b9e77")) +
  labs(x = "N rate (kg/ha)", y = "PFP-N (kg/kg)", colour = "Eq fit", linetype = NULL,
       title = "PFP-N vs N rate: SHAP-based partial dependence and fitted curves",
       subtitle = "Points: every plot, exp(mean prediction + SHAP of N rate)")
save_fig(p_pdN, "6_rf_pd_N_vs_curves")

pd_other <- filter(shap_pd, variable != "N_rate")
pd_smooth <- pd_other |> filter(n_distinct(value) >= 10, .by = variable)   # smooth continuous predictors only
pd_means  <- pd_other |> filter(n_distinct(value) < 10, .by = variable) |>
  summarise(pfp = mean(pfp), .by = c(variable, value))                    # mean per value otherwise
p <- ggplot(pd_other, aes(value, pfp)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(data = pd_smooth, method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = FALSE, colour = "firebrick") +
  geom_point(data = pd_means, colour = "firebrick", size = 2.5) +
  geom_hline(yintercept = exp(base_log), linetype = 2, colour = "grey50") +
  facet_wrap(~variable, scales = "free_x") +
  labs(x = NULL, y = "PFP-N (kg/kg), SHAP-based", title = "SHAP-based partial dependence of PFP-N",
       subtitle = "Every plot: exp(mean prediction + SHAP of the predictor); dashed: mean prediction")
save_fig(p, "6_rf_partial_dependence", w = 10, h = 6)

# Spatial signature: summed SHAP of coordinates and oblique coordinates, mean per EA
sp_space <- shap |>
  filter(variable %in% spatial_vars) |>
  summarise(phi_space = sum(phi), .by = row) |>
  mutate(cluster = rf_data$cluster[row], longitude = X_all$longitude[row], latitude = X_all$latitude[row]) |>
  summarise(phi_space = mean(phi_space), longitude = mean(longitude), latitude = mean(latitude), plots = n(),
            .by = cluster)
if (nrow(sp_space) > 0) {
  p <- ggplot(sp_space, aes(longitude, latitude, colour = exp(phi_space))) +
    geom_point(size = 1.6) + coord_equal() +
    scale_colour_gradient2(midpoint = 1, low = "#b2182b", mid = "grey85", high = "#2166ac", name = "x PFP-N") +
    labs(x = NULL, y = NULL, title = "Spatial signature (SHAP)",
         subtitle = "Mean per EA, all plots: where location raises (blue) or lowers (red) PFP-N")
  save_fig(p, "6_shap_spatial_signature", w = 5, h = 8)
}

# Prediction maps: PFP-N at N_ref kg/ha for each wave (needs covariates as rasters)
map_ok <- length(env_ok) == length(env_vars) || length(env_ok) == 0
if (map_ok) {
  grid <- expand_grid(longitude = seq(min(rf_data$longitude), max(rf_data$longitude), by = 0.05),
                      latitude  = seq(min(rf_data$latitude),  max(rf_data$latitude),  by = 0.05))
  near <- lengths(sf::st_is_within_distance(
    sf::st_as_sf(grid, coords = c("longitude", "latitude"), crs = 4326),
    sf::st_as_sf(distinct(rf_data, longitude, latitude), coords = c("longitude", "latitude"), crs = 4326),
    dist = 25000)) > 0
  grid <- grid[near, ]
  typical <- rf_data |> summarise(across(c(plot_area_ha, share_planted), median),
                                  urea = 1L, compound = 1L, manure = 0L)
  pred_maps <- map(sort(unique(rf_data$year)), \(yr) {
    g <- grid |> mutate(year = yr, season_end = yr, N_rate = N_ref) |> bind_cols(typical) |>
      add_covariates() |> add_spatial_signature()
    g <- g |> filter(if_all(all_of(final$vars), \(x) !is.na(x)))
    g$pfp <- exp(predict(final$model, newdata = as.data.frame(g[, final$vars])))
    g |> select(longitude, latitude, year, pfp)
  }) |> list_rbind() |>
    left_join(distinct(select(rf_data, year, wave)), by = "year")
  save_table(pred_maps, "6_rf_prediction_grid")
  p <- ggplot(pred_maps, aes(longitude, latitude, fill = pfp)) +
    geom_raster() + coord_equal() + facet_wrap(~wave, nrow = 1) +
    scale_fill_viridis_c(name = "PFP-N") +
    labs(x = NULL, y = NULL, title = paste0("Predicted PFP-N at ", N_ref, " kg N/ha (random forest, ", final$set, ")"),
         subtitle = "Typical plot; within 25 km of surveyed EAs") +
    theme(axis.text = element_blank(), panel.grid = element_blank())
  save_fig(p, "6_rf_prediction_maps", w = 12, h = 6)
} else {
  message("Prediction maps skipped: some covariates are not available as rasters")
}

# -----------------------------------------------------------------------------
# 7. PFP-N by wave and region (all waves shown)
# -----------------------------------------------------------------------------

cluster_boot_ci <- function(x, cl, B = n_boot) {
  if (length(x) < 10) return(c(NA_real_, NA_real_))
  idx <- split(seq_along(x), cl)
  reps <- replicate(B, median(x[unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)]))
  unname(quantile(reps, c(0.025, 0.975)))
}

wave_tbl <- surv |>
  summarise(plots = n(), EAs = n_distinct(cluster), median_pfp_n = median(pfp_n),
            ci = list(cluster_boot_ci(pfp_n, cluster)), median_yield = median(yield),
            median_N = median(N_fertilizer), .by = wave) |>
  unnest_wider(ci, names_sep = "_") |>
  rename(ci_low = ci_1, ci_high = ci_2) |>
  complete(wave = factor(wave_labels, levels = wave_labels), fill = list(plots = 0L)) |>
  arrange(wave)
save_table(wave_tbl, "7_pfp_by_wave")

region_tbl <- surv |>
  filter(!is.na(region)) |>
  summarise(plots = n(), median_pfp_n = median(pfp_n), ci = list(cluster_boot_ci(pfp_n, cluster)),
            .by = c(wave, region)) |>
  unnest_wider(ci, names_sep = "_") |>
  rename(ci_low = ci_1, ci_high = ci_2) |>
  complete(wave = factor(wave_labels, levels = wave_labels), region) |>
  arrange(region, wave)
save_table(region_tbl, "7_pfp_by_region_wave")

p <- ggplot(wave_tbl, aes(wave, median_pfp_n, ymin = ci_low, ymax = ci_high)) +
  geom_pointrange(na.rm = TRUE) +
  geom_text(aes(y = 0, label = paste0("n=", plots)), vjust = 0, size = 3, colour = "grey40") +
  scale_x_discrete(drop = FALSE) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "Median PFP-N (kg/kg)", title = "PFP-N by survey wave",
       subtitle = paste0("Sole ", focus_crop, ", N ", N_min, "-", N_max, " kg/ha; 95% EA-cluster bootstrap CI; empty wave = no data"))
save_fig(p, "7_pfp_by_wave")

p <- ggplot(filter(region_tbl, !is.na(region)), aes(wave, median_pfp_n, ymin = ci_low, ymax = ci_high,
                                                    colour = region, group = region)) +
  geom_line(position = position_dodge(0.3), na.rm = TRUE) +
  geom_pointrange(position = position_dodge(0.3), na.rm = TRUE) +
  scale_x_discrete(drop = FALSE) + expand_limits(y = 0) +
  labs(x = NULL, y = "Median PFP-N (kg/kg)", colour = NULL, title = "PFP-N by region and wave")
save_fig(p, "7_pfp_by_region_wave")

# -----------------------------------------------------------------------------
# 8. Survey plots near trials: PFP-N patterns and agronomic use efficiency
# -----------------------------------------------------------------------------
# Trial zone: survey plots within match_km of a trial site AND within
# match_years of the trial's season (IHS6 plots are at district centroids and
# are left out of the zones).
#
# AUE-N = (yield with N - yield without N) / N applied   (kg grain per kg N)
#   Trials : each fertilized plot against the mean yield of the N = 0 plots of
#            the same trial site and season (other treatments, e.g. variety, are
#            pooled, so this is approximate).
#   Surveys: no plot is observed with and without N, so each fertilized plot is
#            matched to an unfertilized plot with a similar propensity score
#            (MatchIt: logistic propensity model, nearest neighbour, with
#            replacement, caliper psm_caliper SD). Within the trial zones the
#            match is exact on the zone (same trial site); across Malawi it is
#            exact on the wave. Matching balances OBSERVED characteristics only:
#            if farmers who fertilize also have better soils or management, the
#            survey AUE-N is biased upwards.

# District boundaries for the map (optional: geodata)
bnd <- if (requireNamespace("geodata", quietly = TRUE)) {
  tryCatch(sf::st_as_sf(geodata::gadm("MWI", level = 1, path = file.path(cov_cache, "geodata"))),
           error = \(e) NULL)
} else NULL

trials_all <- d0 |>
  filter(program == "carob", crop == focus_crop, !is.na(yield), yield > 0, !is.na(N_fertilizer),
         !is_survey %in% TRUE) |>
  prep() |>
  mutate(latitude  = if_else(in_malawi(longitude, latitude), latitude, NA_real_),
         longitude = if_else(is.na(latitude), NA_real_, longitude),
         trial_type = if_else(on_farm %in% TRUE, "On-farm trials",
                              if_else(on_farm %in% FALSE, "On-station trials", "Trials (setting unknown)")),
         season_end = suppressWarnings(as.integer(str_sub(date, 1, 4))),
         site = paste(dataset_id, round(latitude, 3), round(longitude, 3), season_end))

# Survey plots with a yield, either unfertilized (controls) or within the N window
surv_ae <- lsms_maize |>
  filter(!is.na(yield), yield > 0, !is.na(N_fertilizer), N_fertilizer == 0 | N_window) |>
  mutate(z_y = robust_z(log(yield)), .by = source) |>
  filter(abs(z_y) <= mad_k) |>
  left_join(wave_season, by = "source") |>
  mutate(season_end = coalesce(season_end, suppressWarnings(as.integer(str_sub(date, 1, 4)))),
         treat = as.integer(N_fertilizer > 0))

sites <- trials_all |>
  filter(!is.na(latitude), !is.na(season_end)) |>
  summarise(trial_plots = n(), trial_plots_N0 = sum(N_fertilizer == 0),
            .by = c(site, dataset_id, trial_type, latitude, longitude, season_end))
geo_sv <- surv_ae |>
  filter(!is.na(latitude), !is.na(longitude), !geo_level %in% "district_centroid") |>
  mutate(sv_row = row_number())

zone_ok <- nrow(sites) > 0 && nrow(geo_sv) > 0
if (zone_ok) {
  sites_sf <- sf::st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  sv_sf    <- sf::st_as_sf(geo_sv, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  hits <- sf::st_is_within_distance(sites_sf, sv_sf, dist = match_km * 1000)
  pairs <- tibble(site_i = rep(seq_len(nrow(sites)), lengths(hits)), sv_row = unlist(hits)) |>
    mutate(site = sites$site[site_i], site_year = sites$season_end[site_i],
           sv_year = geo_sv$season_end[sv_row]) |>
    filter(abs(site_year - sv_year) <= match_years)
  pairs$dist_km <- if (nrow(pairs)) as.numeric(sf::st_distance(sites_sf[pairs$site_i, ], sv_sf[pairs$sv_row, ],
                                                                by_element = TRUE)) / 1000 else numeric()
  nearest <- pairs |> slice_min(dist_km, by = sv_row, with_ties = FALSE)
  geo_sv <- geo_sv |>
    left_join(select(nearest, sv_row, zone_site = site, dist_km), by = "sv_row") |>
    left_join(select(sites, zone_site = site, zone_trial_type = trial_type), by = "zone_site") |>
    mutate(in_zone = !is.na(zone_site))
  zone_ok <- any(geo_sv$in_zone)
}

if (zone_ok) {
  zone_tbl <- geo_sv |>
    filter(in_zone) |>
    summarise(survey_plots = n(), fertilized = sum(treat == 1), unfertilized = sum(treat == 0),
              EAs = n_distinct(cluster), median_dist_km = median(dist_km),
              waves = paste(sort(unique(as.character(wave))), collapse = "; "), .by = zone_site) |>
    right_join(sites, by = c("zone_site" = "site")) |>
    mutate(across(c(survey_plots, fertilized, unfertilized, EAs), \(x) coalesce(x, 0L))) |>
    select(site = zone_site, dataset_id, trial_type, season_end, latitude, longitude, trial_plots, trial_plots_N0,
           survey_plots, fertilized, unfertilized, EAs, median_dist_km, waves) |>
    arrange(dataset_id, season_end)
  save_table(zone_tbl, "8_trial_zones")

  # Map of trial sites, their 30 km zones and the survey EAs inside / outside
  buffers <- sf::st_buffer(sites_sf, dist = match_km * 1000)
  ea_pts <- geo_sv |> summarise(longitude = mean(longitude), latitude = mean(latitude),
                                in_zone = any(in_zone), .by = cluster)
  p <- ggplot() +
    { if (!is.null(bnd)) geom_sf(data = bnd, fill = "grey95", colour = "grey70", linewidth = 0.2) } +
    geom_sf(data = buffers, aes(fill = trial_type), alpha = 0.25, colour = NA) +
    geom_point(data = ea_pts, aes(longitude, latitude, colour = in_zone), size = 0.7) +
    geom_point(data = sites, aes(longitude, latitude, shape = trial_type), size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey55"), labels = c(`TRUE` = "inside a zone", `FALSE` = "outside"),
                        name = "Survey EAs") +
    labs(fill = paste0(match_km, " km zone"), shape = "Trial site", x = NULL, y = NULL,
         title = paste0("Trial sites and ", match_km, " km zones"),
         subtitle = paste0("Survey plots count as inside when also within ", match_years, " years of the trial")) +
    coord_sf()
  save_fig(p, "8_trial_zones_map", w = 6, h = 9)

  # PFP-N: inside vs outside the zones vs trials (same figure type as section 3)
  zone_pfp <- bind_rows(
    geo_sv |> filter(treat == 1, !is.na(pfp_n), !outlier) |>
      transmute(group = if_else(in_zone, paste0("Survey: within ", match_km, " km of trials"), "Survey: outside trial zones"),
                class = "Survey (LSMS)", pfp_n, yield, N_fertilizer, wave = as.character(wave)),
    trials |> transmute(group = trial_type, class = "Trials", pfp_n, yield, N_fertilizer, wave = NA_character_)
  ) |>
    mutate(group = factor(group, levels = c("Survey: outside trial zones", paste0("Survey: within ", match_km, " km of trials"),
                                            sort(unique(trials$trial_type)))))
  save_fig(violin_box(zone_pfp, "pfp_n", "PFP-N (kg grain / kg N)", name = "8_pfp_zones") +
             ggtitle(paste0("PFP-N near trials vs the rest of Malawi, sole ", focus_crop)), "8_pfp_zones", w = 9)

  # Do the wave patterns near trials follow the national pattern?
  zone_wave <- geo_sv |>
    filter(treat == 1, !is.na(pfp_n), !outlier) |>
    summarise(plots = n(), median_pfp_n = median(pfp_n), .by = c(wave, in_zone)) |>
    mutate(where = if_else(in_zone, "within trial zones", "outside trial zones")) |>
    select(-in_zone) |>
    bind_rows(wave_tbl |> filter(plots > 0) |> transmute(wave, plots, median_pfp_n, where = "all Malawi (section 7)")) |>
    complete(wave = factor(wave_labels, levels = wave_labels), where)
  save_table(zone_wave, "8_pfp_by_wave_zones")
  zone_test <- geo_sv |>
    filter(treat == 1, !is.na(pfp_n), !outlier) |>
    summarise(n_inside = sum(in_zone), n_outside = sum(!in_zone),
              p_wilcoxon = if (sum(in_zone) >= 3 && sum(!in_zone) >= 3)
                suppressWarnings(wilcox.test(pfp_n[in_zone], pfp_n[!in_zone])$p.value) else NA_real_,
              .by = wave)
  save_table(zone_test, "8_pfp_zone_vs_outside_tests")
  p <- ggplot(filter(zone_wave, !is.na(median_pfp_n)), aes(wave, median_pfp_n, colour = where, group = where)) +
    geom_line() + geom_point(aes(size = plots)) +
    scale_x_discrete(drop = FALSE) + expand_limits(y = 0) +
    labs(x = NULL, y = "Median PFP-N (kg/kg)", colour = NULL, size = "plots",
         title = "PFP-N by wave: near trials vs the rest of Malawi")
  save_fig(p, "8_pfp_by_wave_zones")
}

# AUE-N in trials
ae_trial <- trials_all |>
  mutate(y0 = mean(yield[N_fertilizer == 0]), n0 = sum(N_fertilizer == 0), .by = site) |>
  filter(N_fertilizer > 0, n0 > 0, between(N_fertilizer, N_min, N_max)) |>
  mutate(aue = (yield - y0) / N_fertilizer,
         n_class = cut(N_fertilizer, c(0, 50, 100, Inf), labels = c("1-50", "51-100", ">100")))
save_table(ae_trial, "8_aue_trial_plots")

boot_ci <- function(x, fun, B = 500) {
  if (length(x) < 5) return(c(NA_real_, NA_real_))
  unname(quantile(replicate(B, fun(sample(x, replace = TRUE))), c(0.025, 0.975), na.rm = TRUE))
}
boot_ratio_ci <- function(dy, n, B = 500) {
  if (length(dy) < 5) return(c(NA_real_, NA_real_))
  unname(quantile(replicate(B, { i <- sample.int(length(dy), replace = TRUE); sum(dy[i]) / sum(n[i]) }),
                  c(0.025, 0.975), na.rm = TRUE))
}

# AUE-N in surveys by propensity-score matching
psm_aue <- function(df, label, exact_var) {
  covs <- intersect(c("log_area", "manure", "latitude", "longitude", env_ok), names(df))
  df <- df |>
    mutate(log_area = log(plot_area), manure = as.integer(OM_used %in% TRUE)) |>
    filter(if_all(all_of(covs), \(x) !is.na(x)), !is.na(.data[[exact_var]]))
  # strata (zone or wave) need both fertilized and unfertilized plots
  ok <- df |> summarise(t = sum(treat == 1), c = sum(treat == 0), .by = all_of(exact_var)) |> filter(t > 0, c > 0)
  df <- df |> semi_join(ok, by = exact_var) |> mutate(stratum = factor(.data[[exact_var]]))
  info <- tibble(analysis = label, treated = sum(df$treat == 1), controls = sum(df$treat == 0),
                 strata = n_distinct(df$stratum))
  if (min(info$treated, info$controls) < 20) {
    return(list(info = mutate(info, possible = FALSE,
                              reason = "fewer than 20 fertilized or unfertilized plots with all covariates"),
                pairs = NULL, balance = NULL))
  }
  f <- reformulate(covs, "treat")
  m <- MatchIt::matchit(f, data = as.data.frame(df), method = "nearest", distance = "glm",
                        replace = TRUE, caliper = psm_caliper, std.caliper = TRUE,
                        exact = if (nlevels(df$stratum) > 1) ~ stratum else NULL)
  md <- as_tibble(MatchIt::get_matches(m))
  pairs <- md |>
    summarise(y1 = yield[treat == 1][1], y0 = yield[treat == 0][1], N = N_fertilizer[treat == 1][1],
              stratum = first(as.character(stratum)), .by = subclass) |>
    filter(!is.na(y1), !is.na(y0)) |>
    mutate(aue = (y1 - y0) / N, analysis = label)
  sm <- summary(m)
  balance <- tibble(variable = rownames(sm$sum.all), smd_before = sm$sum.all[, "Std. Mean Diff."],
                    smd_after = sm$sum.matched[rownames(sm$sum.all), "Std. Mean Diff."], analysis = label)
  list(info = mutate(info, possible = TRUE, reason = NA_character_, pairs = nrow(pairs),
                     treated_unmatched = info$treated - n_distinct(md$subclass),
                     max_abs_smd_after = max(abs(balance$smd_after), na.rm = TRUE)),
       pairs = pairs, balance = balance)
}

psm_data <- surv_ae |> filter(!is.na(latitude), !is.na(longitude)) |> add_covariates()
psm_runs <- list()
if (zone_ok) {
  zone_data <- psm_data |>
    inner_join(select(filter(geo_sv, in_zone), source, hhid, plot_id, yield_source, zone_site),
               by = c("source", "hhid", "plot_id", "yield_source"))
  psm_runs$zones <- psm_aue(zone_data, paste0("Survey PSM within ", match_km, " km trial zones"), "zone_site")
}
psm_runs$malawi <- psm_aue(psm_data |> mutate(wave_chr = as.character(wave)), "Survey PSM, all of Malawi", "wave_chr")

psm_info    <- map(psm_runs, "info") |> list_rbind()
psm_pairs   <- map(psm_runs, "pairs") |> compact() |> list_rbind()
psm_balance <- map(psm_runs, "balance") |> compact() |> list_rbind()
save_table(psm_info, "8_psm_feasibility")
if (nrow(psm_pairs)) save_table(psm_pairs, "8_psm_matched_pairs")
if (nrow(psm_balance)) save_table(psm_balance, "8_psm_balance")

summ_aue <- function(aue, dy, N, label, n_units) {
  ci_med <- boot_ci(aue, median); ci_rat <- boot_ratio_ci(dy, N)
  tibble(estimate_of = label, n = n_units, median_aue = median(aue), median_ci_low = ci_med[1], median_ci_high = ci_med[2],
         ratio_aue = sum(dy) / sum(N), ratio_ci_low = ci_rat[1], ratio_ci_high = ci_rat[2])
}
aue_tbl <- bind_rows(
  ae_trial |> group_by(trial_type) |>
    group_map(\(x, k) summ_aue(x$aue, x$yield - x$y0, x$N_fertilizer, k$trial_type, nrow(x))) |> list_rbind(),
  psm_pairs |> group_by(analysis) |>
    group_map(\(x, k) summ_aue(x$aue, x$y1 - x$y0, x$N, k$analysis, nrow(x))) |> list_rbind()
)
mp <- curve_params |> filter(group_type == "Pooled", parameter %in% c("MP_50", "MP_100"))
if (nrow(mp)) {
  aue_tbl <- bind_rows(aue_tbl, mp |> transmute(estimate_of = paste0("Curve marginal response, ", form, ", ",
                                                                     str_remove(parameter, "MP_"), " kg N/ha"),
                                                n = n, median_aue = estimate, median_ci_low = ci_low, median_ci_high = ci_high))
}
save_table(aue_tbl, "8_aue_estimates")

if (nrow(aue_tbl)) {
  p <- aue_tbl |>
    mutate(estimate_of = fct_rev(fct_inorder(estimate_of))) |>
    ggplot(aes(median_aue, estimate_of, xmin = median_ci_low, xmax = median_ci_high)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange() +
    labs(x = "AUE-N (kg grain per kg N): median, 95% bootstrap CI", y = NULL,
         title = "Agronomic use efficiency of N: trials vs matched survey plots",
         subtitle = "Survey AUE-N from propensity-score matched fertilized/unfertilized plots; curves: marginal response")
  save_fig(p, "8_aue_estimates", w = 9, h = 4.5)
}
if (nrow(ae_trial)) {
  p <- ggplot(ae_trial, aes(n_class, aue, fill = trial_type)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.size = 0.6) +
    labs(x = "N rate (kg/ha)", y = "AUE-N (kg grain / kg N)", fill = NULL,
         title = "AUE-N in trials by N rate", subtitle = "Each fertilized plot vs the N = 0 plots of its trial site and season")
  save_fig(p, "8_aue_trials_by_N")
}
if (nrow(psm_balance)) {
  p <- psm_balance |>
    pivot_longer(c(smd_before, smd_after), names_to = "when", values_to = "smd") |>
    mutate(when = factor(if_else(when == "smd_before", "before matching", "after matching"),
                         levels = c("before matching", "after matching"))) |>
    ggplot(aes(abs(smd), variable, colour = when)) +
    geom_vline(xintercept = 0.1, linetype = 2) + geom_point(size = 2) +
    facet_wrap(~analysis) +
    labs(x = "|standardised mean difference|", y = NULL, colour = NULL,
         title = "Covariate balance, fertilized vs unfertilized plots", subtitle = "Below 0.1 = well balanced")
  save_fig(p, "8_psm_balance", w = 10)
}

# -----------------------------------------------------------------------------
# 9. Results object and report
# -----------------------------------------------------------------------------

results <- list(
  created = Sys.time(), in_rds = in_rds, focus_crop = focus_crop,
  settings = list(N_min = N_min, N_max = N_max, mad_k = mad_k, k_folds = k_folds, n_blocks = n_blocks,
                  n_boot = n_boot, n_trees = n_trees, N_ref = N_ref, psm_caliper = psm_caliper, match_km = match_km, match_years = match_years, shap_nsim = shap_nsim),
  n_survey = nrow(surv), n_trials = nrow(trials),
  funnel = funnel, missing_waves = missing_waves, dropped_waves = dropped_waves,
  quality = quality,
  curve_comp = curve_comp, curve_params = curve_params,
  env_vars = env_vars, env_used = env_ok, predictor_sets = predictor_sets,
  rf_n = nrow(rf_data), rf_eas = n_distinct(rf_data$cluster), rf_blocks = n_distinct(rf_data$block),
  rf_metrics = rf_metrics, rf_importance = imp, rf_final_set = final$set, rf_pd = pd,
  shap_importance = shap_imp, shap_n = nrow(X_all), shap_nsim = shap_nsim, shap_check = shap_check,
  match_km = match_km, match_years = match_years, zone_ok = zone_ok,
  zone_tbl = if (zone_ok) zone_tbl else NULL, zone_wave = if (zone_ok) zone_wave else NULL,
  zone_test = if (zone_ok) zone_test else NULL, aue_tbl = aue_tbl, psm_info = psm_info,
  wave_tbl = wave_tbl, region_tbl = region_tbl,
  fig_dir = normalizePath(fig_dir), tab_dir = normalizePath(tab_dir)
)
saveRDS(results, file.path(obj_dir, "results.rds"))
saveRDS(final$model, file.path(obj_dir, "rf_final_model.rds"))

if (render_report) {
  if (!file.exists(report_rmd)) {
    message("Report template not found: ", report_rmd)
  } else if (!rmarkdown::pandoc_available()) {
    message("Pandoc not available; open ", report_rmd, " in RStudio and knit it")
  } else {
    for (fmt in c("html_document", "word_document")) {
      rmarkdown::render(report_rmd, output_format = fmt, output_dir = out_dir,
                        params = list(results = normalizePath(file.path(obj_dir, "results.rds")),
                                      analysis1 = if (file.exists(analysis1_rds)) normalizePath(analysis1_rds) else ""),
                        envir = new.env(), quiet = TRUE)
    }
    message("Report: ", file.path(out_dir, "pfpn_report.html"), " and .docx")
  }
}
message("Done. Outputs in ", out_dir)

# ==============================================================================
# END
# ==============================================================================