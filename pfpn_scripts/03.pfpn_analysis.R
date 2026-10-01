# =============================================================================
# How does PFP-N vary across space and time in Malawi?
# -----------------------------------------------------------------------------
# Input : pfpn_compiled/pfpn_unified.rds (from pfpn_compile.R)
# Output: pfpn_analysis/{tables,figures}/ and pfpn_analysis/summary.md
#
# Part 1  Data robustness (by quality and by source)
#   1.1 overview, quality flags and robust outliers per source
#   1.2 plot area: reported vs GPS (and why PFP-N is immune to area error)
#   1.3 fertilizer quantities: enumerator kg vs quantity x unit
#   1.4 harvest: second question block (IHS3, IHPS 2013); crop cut vs farmer report (LCAS)
#   1.5 trials as benchmarks for survey data
#       - distributions by source type (survey, on-farm trial, on-station trial)
#       - attainable-yield frontier from trials (quantile regression) and the share
#         of survey plots above it at the same N rate
#       - survey EAs matched in space and time to trial sites, compared by N class
# Part 2  PFP-N across space and time (LSMS sole maize, outliers removed)
#   2.1 national medians per wave with EA-cluster bootstrap CIs
#   2.2 district x wave medians and maps
#   2.3 robust GAMs (scaled-t) with a spatial smooth, wave effects and EA random
#       effects, for log PFP-N, log yield and log N (decomposition), and log PFP-N
#       at a common N rate
#   2.4 variance components: district / EA / household / plot
#   2.5 sensitivity of wave effects to data-quality choices
#
# Trials are NOT ground truth: they are researcher-managed, often better
# managed, on other sites and years. They are used as a plausibility benchmark
# (attainable yield at a given N) and for like-for-like comparisons where sites
# and survey EAs are close in space and time.
#
# Packages: tidyverse; via pkg::fun mgcv, glmmTMB, quantreg, sf, splines,
# patchwork, and (optional, for district maps) geodata + terra.
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (mgcv, glmmTMB, quantreg, sf, splines, patchwork, geodata)

source("00.pfpn_utils.R", local = TRUE)   # paths, Malawi coordinate check, compact letter display

in_rds  <- pfpn_path("pfpn_compiled", "pfpn_unified.rds")
out_dir <- pfpn_path("pfpn_analysis")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
walk(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# Settings
focus_crop  <- "maize"   # PFP-N is crop specific; maize receives most fertilizer in Malawi
low_N       <- 10        # kg N/ha; PFP-N from smaller rates is unstable
mad_k       <- 4         # robust z (median/MAD, log scale) beyond which a row is an outlier
match_km    <- 30        # survey EA within this distance of a trial site ...
match_years <- 2         # ... and within this many years
n_boot      <- 500       # EA-cluster bootstrap replicates
min_rows    <- 30        # minimum rows for a section / model to run
boundary_cache <- pfpn_path("pfpn_compiled", "raw")   # where geodata stores GADM files
set.seed(20260930)

theme_set(theme_minimal(base_size = 11))
save_table <- \(x, name) write_csv(x, file.path(tab_dir, paste0(name, ".csv")), na = "")
save_fig   <- \(p, name, w = 8, h = 5) {
  save_plot_data(p, name, tab_dir)   # the data behind every figure -> tables/fig_<name>.csv
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
findings <- character()
note <- \(...) findings <<- c(findings, paste0(...))
pct  <- \(x) sprintf("%.0f%%", 100 * mean(x, na.rm = TRUE))
skip <- \(what, why) message("Skipped ", what, ": ", why)

# -----------------------------------------------------------------------------
# 0. Data
# -----------------------------------------------------------------------------

wave_labels <- c("LSMS_1003_2010-2011" = "IHS3 2010/11", "LSMS_2248_2013" = "IHPS 2013",
                 "LSMS_2936_2016-2017" = "IHS4 2016/17", "LSMS_3818_2019-2020" = "IHS5 2019/20",
                 "LSMS_8507_2024-2025" = "IHS6 2024/25")

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

d <- d0 |>
  filter(crop == focus_crop) |>
  mutate(
    source_type = case_when(
      program == "LSMS_MWI"                          ~ "Survey: LSMS farmer report",
      program == "LCAS" & yield_source == "crop_cut" ~ "Survey: LCAS crop cut",
      program == "LCAS"                              ~ "Survey: LCAS farmer report",
      program == "carob" & is_survey %in% TRUE       ~ "Survey: Carob",
      program == "carob" & on_farm %in% TRUE         ~ "Trial: on-farm",
      program == "carob" & on_farm %in% FALSE        ~ "Trial: on-station",
      .default                                       = "Trial: setting unknown"),
    source_class = if_else(str_starts(source_type, "Trial"), "trial", "survey"),
    year     = suppressWarnings(as.integer(str_sub(as.character(date), 1, 4))),
    wave     = factor(unname(wave_labels[source]), levels = wave_labels),
    cluster  = paste(source, coalesce(as.character(ea_id), hhid)),   # EA (LSMS) or household
    hh       = paste(source, hhid),
    dkey     = district_key(adm2),
    log_y    = if_else(yield > 0, log(yield), NA_real_),
    log_n    = if_else(N_fertilizer > 0, log(N_fertilizer), NA_real_),
    log_pfp  = if_else(pfp_n > 0, log(pfp_n), NA_real_),
    n_class  = cut(N_fertilizer, c(-Inf, 0, 50, 100, Inf), labels = c("0", "1-50", "51-100", ">100"))
  ) |>
  mutate(across(c(log_y, log_n, log_pfp), robust_z, .names = "z_{.col}"), .by = c(source_type, source)) |>
  mutate(outlier = (abs(z_log_y) > mad_k) %in% TRUE | (abs(z_log_n) > mad_k) %in% TRUE |
                   (abs(z_log_pfp) > mad_k) %in% TRUE)

if (nrow(d) == 0) stop("No rows for crop '", focus_crop, "' in ", in_rds)
message(nrow(d), " ", focus_crop, " rows: ",
        paste(names(table(d$source_type)), table(d$source_type), sep = " = ", collapse = "; "))

if ("survey_weight" %in% names(d)) {
  note("Survey weights are present but not used by this script; results describe the sample.")
} else {
  note("No survey weights in the compiled data: results describe the sampled plots, not the ",
       "population of Malawian maize plots.")
}

# =============================================================================
# PART 1. DATA ROBUSTNESS
# =============================================================================

# 1.1 Overview, quality flags, outliers ----------------------------------------

overview <- d |>
  summarise(
    rows          = n(),
    plots         = n_distinct(paste(source, hhid, plot_id)),
    with_pfp_n    = sum(!is.na(pfp_n)),
    zero_N        = sum(N_fertilizer %in% 0),
    N_unknown     = sum(is.na(N_fertilizer)),
    outliers      = sum(outlier),
    out_of_range  = sum(!is.na(qc_out_of_range)),
    median_yield  = median(yield, na.rm = TRUE),
    median_N_pos  = median(N_fertilizer[N_fertilizer > 0], na.rm = TRUE),
    median_pfp_n  = median(pfp_n, na.rm = TRUE),
    .by = c(source_class, source_type, source)
  ) |>
  arrange(source_class, source_type, source)
save_table(overview, "1_1_overview_by_source")

flag_shares <- d |>
  select(source_type, source, quality_flags) |>
  mutate(row = row_number()) |>
  separate_longer_delim(quality_flags, "; ") |>
  filter(!is.na(quality_flags)) |>
  count(source_type, source, flag = quality_flags) |>
  left_join(count(d, source_type, source, name = "rows"), by = c("source_type", "source")) |>
  mutate(pct_rows = round(100 * n / rows, 1)) |>
  arrange(source, desc(pct_rows))
save_table(flag_shares, "1_1_quality_flags_by_source")

note("Robust outliers (|median/MAD z| > ", mad_k, " on log yield, log N or log PFP-N, per source): ",
     sum(d$outlier), " of ", nrow(d), " rows (", pct(d$outlier), ").")

# 1.2 Plot area ----------------------------------------------------------------
# For sole crops, yield and N are both per planted hectare, so
# PFP-N = harvest kg / N kg and the area cancels: area errors bias yield and N
# rates but not PFP-N. Checked below on rows with both areas.

area_pairs <- d |>
  filter(source_class == "survey", plot_area_gps_m2 > 0, plot_area_reported_m2 > 0) |>
  distinct(source, hhid, plot_id, .keep_all = TRUE) |>
  mutate(
    log_ratio   = log(plot_area_reported_m2 / plot_area_gps_m2),
    gps_decile  = ntile(plot_area_gps_m2, 10),
    # yield and PFP-N recomputed with the reported instead of the GPS area
    yield_rep   = yield * plot_area_gps_m2 / plot_area_reported_m2,
    n_rep       = N_fertilizer * plot_area_gps_m2 / plot_area_reported_m2,
    pfp_rep     = if_else(n_rep > 0, yield_rep / n_rep, NA_real_)
  )

if (nrow(area_pairs) >= min_rows) {
  area_tbl <- area_pairs |>
    summarise(n = n(),
              median_gps_m2 = median(plot_area_gps_m2),
              median_reported_over_gps = exp(median(log_ratio)),
              q25 = exp(quantile(log_ratio, 0.25)), q75 = exp(quantile(log_ratio, 0.75)),
              .by = gps_decile) |>
    arrange(gps_decile)
  save_table(area_tbl, "1_2_reported_vs_gps_area_by_size_decile")

  area_fit <- quantreg::rq(log_ratio ~ log(plot_area_gps_m2), tau = 0.5, data = area_pairs)
  slope <- coef(area_fit)[2]
  pfp_change <- with(area_pairs, max(abs(pfp_rep / pfp_n - 1), na.rm = TRUE))

  p <- ggplot(area_pairs, aes(plot_area_gps_m2, exp(log_ratio))) +
    geom_bin_2d(bins = 60) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_quantile(quantiles = c(0.25, 0.5, 0.75), formula = y ~ x, colour = "firebrick") +
    scale_x_log10() + scale_y_log10() + scale_fill_viridis_c(trans = "log10") +
    labs(x = "GPS area (m2)", y = "Reported / GPS area",
         title = "Farmer-reported vs GPS plot area",
         subtitle = "Median regression (red); 1 = no error")
  save_fig(p, "1_2_reported_vs_gps_area")

  note("Plot area (", nrow(area_pairs), " plots with GPS and reported area): median reported/GPS = ",
       round(exp(median(area_pairs$log_ratio)), 2), "; median-regression slope on log GPS area = ",
       round(slope, 2), if (slope < 0) " (small plots are over-reported, large plots under-reported)." else ".",
       " Yield and N rates inherit this error; PFP-N does not (largest change when swapping areas: ",
       signif(100 * pfp_change, 2), "%).")
} else skip("1.2 area comparison", "too few plots with both GPS and reported area")

# 1.3 Fertilizer quantities ----------------------------------------------------

fert_q <- d |>
  filter(source_class == "survey", fertilizer_used %in% TRUE) |>
  summarise(fertilized_rows  = n(),
            pct_enumerator_kg = round(100 * mean(str_detect(coalesce(fert_kg_source, ""), "enumerator")), 1),
            pct_kg_mismatch   = round(100 * mean(fert_kg_mismatch %in% TRUE), 1),
            pct_N_unknown     = round(100 * mean(is.na(N_fertilizer)), 1),
            pct_N_below_low   = round(100 * mean(N_fertilizer > 0 & N_fertilizer < low_N, na.rm = TRUE), 1),
            .by = source)
save_table(fert_q, "1_3_fertilizer_quantity_quality")

pfp_by_kg_source <- d |>
  filter(source_class == "survey", !is.na(pfp_n)) |>
  mutate(kg_basis = case_when(str_detect(coalesce(fert_kg_source, ""), ";") ~ "mixed",
                              .default = coalesce(fert_kg_source, "none")),
         kg_mismatch = fert_kg_mismatch %in% TRUE) |>
  summarise(n = n(), median_pfp_n = median(pfp_n),
            q25 = quantile(pfp_n, 0.25), q75 = quantile(pfp_n, 0.75),
            .by = c(source, kg_basis, kg_mismatch))
save_table(pfp_by_kg_source, "1_3_pfp_by_fertilizer_kg_basis")

# 1.4 Harvest measurement ------------------------------------------------------

alt_sources <- d |> filter(harvest_block %in% "alt") |> distinct(source) |> pull()
if (length(alt_sources) > 0) {
  hb <- d |> filter(source %in% alt_sources, !is.na(harvest_block), !is.na(log_y))
  hb_tbl <- hb |>
    summarise(n = n(), median_yield = median(yield), median_pfp_n = median(pfp_n, na.rm = TRUE),
              .by = c(source, harvest_block))
  hb_test <- hb |>
    summarise(p_wilcoxon_log_yield = if (n_distinct(harvest_block) == 2)
                suppressWarnings(wilcox.test(log_y ~ harvest_block)$p.value) else NA_real_,
              .by = source)
  save_table(left_join(hb_tbl, hb_test, by = "source"), "1_4_harvest_block_comparison")
}

cc_pairs <- d |>
  filter(program == "LCAS", !is.na(yield)) |>
  select(source, hhid, plot_id, yield_source, yield) |>
  pivot_wider(names_from = yield_source, values_from = yield, values_fn = first) |>
  filter(if_all(any_of(c("crop_cut", "farmer_report")), \(x) !is.na(x)))
if (all(c("crop_cut", "farmer_report") %in% names(cc_pairs)) && nrow(cc_pairs) >= min_rows) {
  cc_pairs <- cc_pairs |> mutate(log_ratio = log(farmer_report / crop_cut))
  save_table(summarise(cc_pairs, n = n(), median_farmer_over_cropcut = exp(median(log_ratio)),
                       q25 = exp(quantile(log_ratio, 0.25)), q75 = exp(quantile(log_ratio, 0.75)),
                       .by = source), "1_4_farmer_report_vs_crop_cut")
  note("LCAS farmer-reported vs crop-cut yield (", nrow(cc_pairs), " plots): median ratio ",
       round(exp(median(cc_pairs$log_ratio)), 2), ".")
} else skip("1.4 crop cut vs farmer report", "no LCAS plots with both yield measures")

# 1.5 Trials as benchmarks -----------------------------------------------------

dist_tbl <- d |>
  pivot_longer(c(yield, N_fertilizer, pfp_n), names_to = "variable", values_to = "value") |>
  filter(!is.na(value), !(variable == "N_fertilizer" & value == 0)) |>
  summarise(n = n(), q10 = quantile(value, 0.10), q25 = quantile(value, 0.25),
            median = median(value), q75 = quantile(value, 0.75), q90 = quantile(value, 0.90),
            .by = c(variable, source_type)) |>
  arrange(variable, source_type)
save_table(dist_tbl, "1_5_distributions_by_source_type")

vd <- d |> filter(!is.na(pfp_n))
cap <- quantile(vd$pfp_n, 0.99)
letters_st <- cld_letters(vd$pfp_n, vd$source_type) |> rename(source_type = group)
save_table(letters_st, "1_5_pfp_by_source_type_letters")
p <- ggplot(vd, aes(source_type, pfp_n)) +
  geom_violin(fill = "grey90", colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_text(data = letters_st, aes(source_type, cap * 0.97, label = letter), inherit.aes = FALSE,
            fontface = "bold", size = 4) +
  coord_cartesian(ylim = c(0, cap)) +
  labs(x = NULL, y = "PFP-N (kg grain / kg N)", title = paste("PFP-N by source type,", focus_crop),
       caption = paste0("Letters: pairwise Wilcoxon tests, Bonferroni-adjusted (alpha 0.05); shared letter = no difference.\n",
                        "Axis cut at the 99th percentile (", round(cap), "); ", sum(vd$pfp_n > cap), " values above not shown.")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
save_fig(p, "1_5_pfp_by_source_type", w = 9)

tr <- d |> filter(source_class == "trial", !is.na(yield), !is.na(N_fertilizer))
sv <- d |> filter(source_class == "survey", !is.na(yield), !is.na(N_fertilizer))

# Attainable-yield frontier from trials: yield quantiles as a smooth function of N
if (nrow(tr) >= min_rows && n_distinct(tr$N_fertilizer) >= 3) {
  n_max <- as.numeric(quantile(tr$N_fertilizer, 0.99))
  df_sp <- min(3, n_distinct(tr$N_fertilizer) - 1)
  taus  <- c(0.5, 0.9, 0.95)
  tr_fit <- filter(tr, N_fertilizer <= n_max)
  frontier <- quantreg::rq(yield ~ splines::bs(N_fertilizer, df = df_sp, Boundary.knots = c(0, n_max)),
                           tau = taus, data = tr_fit)

  sv_in <- sv |> filter(between(N_fertilizer, 0, n_max))
  pred  <- predict(frontier, newdata = sv_in)
  sv_in <- sv_in |> mutate(above_trial_q50 = yield > pred[, 1], above_trial_q90 = yield > pred[, 2],
                           above_trial_q95 = yield > pred[, 3])
  frontier_tbl <- sv_in |>
    summarise(n = n(),
              pct_above_trial_median = round(100 * mean(above_trial_q50), 1),
              pct_above_trial_q90    = round(100 * mean(above_trial_q90), 1),
              pct_above_trial_q95    = round(100 * mean(above_trial_q95), 1),
              .by = c(source_type, source, area_source))
  save_table(frontier_tbl, "1_5_survey_vs_trial_frontier")

  # which survey quality flags go with implausibly high yields?
  flag_frontier <- sv_in |>
    mutate(row = row_number()) |>
    separate_longer_delim(quality_flags, "; ") |>
    mutate(quality_flags = coalesce(quality_flags, "(no flag)")) |>
    summarise(n = n(), pct_above_trial_q95 = round(100 * mean(above_trial_q95), 1), .by = quality_flags) |>
    arrange(desc(pct_above_trial_q95))
  save_table(flag_frontier, "1_5_frontier_exceedance_by_quality_flag")

  grid_n <- tibble(N_fertilizer = seq(0, n_max, length.out = 100))
  curves <- as_tibble(predict(frontier, newdata = grid_n), .name_repair = \(x) paste0("q", taus * 100)) |>
    bind_cols(grid_n) |>
    pivot_longer(-N_fertilizer, names_to = "quantile", values_to = "yield")
  p <- ggplot() +
    geom_bin_2d(data = sv_in, aes(N_fertilizer, yield), bins = 50) +
    geom_point(data = tr_fit, aes(N_fertilizer, yield), shape = 21, colour = "black", fill = "white",
               size = 1.2, alpha = 0.7) +
    geom_line(data = curves, aes(N_fertilizer, yield, linetype = quantile), colour = "firebrick", linewidth = 0.8) +
    scale_fill_viridis_c(trans = "log10", name = "survey plots") +
    scale_y_log10() +
    labs(x = "N applied (kg/ha)", y = "Yield (kg/ha, log scale)", linetype = "trial quantile",
         title = "Survey plots against the trial yield frontier",
         subtitle = "Open circles: trial plots; red: trial yield quantiles given N; survey outliers included")
  save_fig(p, "1_5_survey_vs_trial_frontier")

  note("Trial frontier (", nrow(tr_fit), " trial rows, N 0-", round(n_max), " kg/ha): ",
       pct(sv_in$above_trial_q95), " of survey plots exceed the trials' 95th-percentile yield at the ",
       "same N (", pct(sv_in$above_trial_q50), " exceed the trial median). Shares well above 5% point ",
       "to over-reported harvests, under-reported N or unit-conversion errors, unless trials in ",
       "Malawi are not representative of farm conditions.")
} else skip("1.5 trial frontier", "fewer than 30 trial rows or fewer than 3 N levels")

# Survey EAs close in space and time to trial sites
sites <- tr |>
  filter(!is.na(latitude), !is.na(longitude), !is.na(year)) |>
  distinct(dataset_id, latitude, longitude, year, source_type) |>
  mutate(site = row_number())
sv_geo <- sv |>
  filter(!is.na(latitude), !is.na(longitude), !is.na(year), !geo_level %in% "district_centroid") |>
  mutate(sv_row = row_number())

if (nrow(sites) > 0 && nrow(sv_geo) > 0) {
  hits <- sf::st_is_within_distance(
    sf::st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326),
    sf::st_as_sf(sv_geo, coords = c("longitude", "latitude"), crs = 4326),
    dist = match_km * 1000)
  pairs <- tibble(site = rep(sites$site, lengths(hits)), sv_row = unlist(hits)) |>
    left_join(select(sites, site, dataset_id, site_year = year, trial_type = source_type), by = "site") |>
    left_join(select(sv_geo, sv_row, sv_year = year, n_class, yield, pfp_n), by = "sv_row") |>
    filter(abs(site_year - sv_year) <= match_years)

  if (nrow(pairs) > 0) {
    trial_site <- tr |>
      inner_join(select(sites, site, dataset_id, latitude, longitude, year),
                 by = c("dataset_id", "latitude", "longitude", "year")) |>
      summarise(trial_n = n(), trial_median_yield = median(yield),
                trial_median_pfp = median(pfp_n, na.rm = TRUE), .by = c(site, n_class))
    survey_site <- pairs |>
      summarise(survey_n = n(), survey_median_yield = median(yield),
                survey_median_pfp = median(pfp_n, na.rm = TRUE), .by = c(site, trial_type, n_class))
    matched <- inner_join(survey_site, trial_site, by = c("site", "n_class")) |>
      mutate(yield_ratio = survey_median_yield / trial_median_yield,
             pfp_ratio   = survey_median_pfp / trial_median_pfp)
    save_table(matched, "1_5_matched_sites_by_N_class")
    matched_sum <- matched |>
      summarise(sites = n_distinct(site), survey_plots = sum(survey_n),
                median_yield_ratio = median(yield_ratio, na.rm = TRUE),
                median_pfp_ratio = median(pfp_ratio, na.rm = TRUE), .by = c(trial_type, n_class)) |>
      arrange(trial_type, n_class)
    save_table(matched_sum, "1_5_matched_summary")
    note("Matched comparison (survey EAs within ", match_km, " km and ", match_years, " years of ",
         n_distinct(matched$site), " trial sites): median survey/trial yield ratio by N class = ",
         paste0(str_remove(matched_sum$trial_type, "Trial: "), " trials, N ", matched_sum$n_class, ": ",
                round(matched_sum$median_yield_ratio, 2), collapse = "; "),
         ". Ratios below 1 are expected (yield gap); ratios above 1 are a warning sign.")
  } else skip("1.5 matched comparison", "no survey EA within the distance and year window")
} else skip("1.5 matched comparison", "no georeferenced, dated trial sites or survey EAs")

# =============================================================================
# PART 2. PFP-N ACROSS SPACE AND TIME (LSMS, sole maize, outliers removed)
# =============================================================================

a <- d |>
  filter(program == "LSMS_MWI", !is.na(log_pfp), !outlier, !is.na(wave)) |>
  mutate(wave_f = droplevels(wave), cluster_f = factor(cluster), hh_f = factor(hh),
         district_f = factor(dkey))
if (nrow(a) < min_rows) stop("Too few LSMS rows with PFP-N for Part 2 (", nrow(a), ")")
ref_wave <- levels(a$wave_f)[1]

# 2.1 National medians per wave with EA-cluster bootstrap ---------------------

cluster_boot_ci <- function(x, cl, B) {
  idx <- split(seq_along(x), cl)
  reps <- replicate(B, {
    s <- sample.int(length(idx), length(idx), replace = TRUE)
    median(x[unlist(idx[s], use.names = FALSE)])
  })
  unname(quantile(reps, c(0.025, 0.975)))
}

wave_tbl <- a |>
  summarise(plots = n(), EAs = n_distinct(cluster),
            median_pfp_n = median(pfp_n),
            ci = list(cluster_boot_ci(pfp_n, cluster, n_boot)),
            median_yield = median(yield), median_N = median(N_fertilizer),
            .by = wave_f) |>
  unnest_wider(ci, names_sep = "_") |>
  rename(pfp_ci_low = ci_1, pfp_ci_high = ci_2) |>
  arrange(wave_f)
save_table(wave_tbl, "2_1_pfp_by_wave")

p1 <- ggplot(wave_tbl, aes(wave_f, median_pfp_n, ymin = pfp_ci_low, ymax = pfp_ci_high)) +
  geom_pointrange() +
  labs(x = NULL, y = "Median PFP-N (kg/kg)", title = "PFP-N by survey wave",
       subtitle = paste0("Sole ", focus_crop, ", fertilized plots; 95% EA-cluster bootstrap CI"))
p2 <- wave_tbl |>
  select(wave_f, `Yield (kg/ha)` = median_yield, `N (kg/ha)` = median_N) |>
  pivot_longer(-wave_f) |>
  ggplot(aes(wave_f, value, group = name)) + geom_line() + geom_point() +
  facet_wrap(~name, scales = "free_y") + labs(x = NULL, y = "Median")
save_fig(patchwork::wrap_plots(p1, p2, ncol = 1), "2_1_pfp_by_wave", h = 7)

note("National median PFP-N by wave: ",
     paste0(wave_tbl$wave_f, " ", round(wave_tbl$median_pfp_n, 1), " [", round(wave_tbl$pfp_ci_low, 1),
            "-", round(wave_tbl$pfp_ci_high, 1), "]", collapse = "; "), " kg/kg.")

# 2.2 District x wave ----------------------------------------------------------

district_tbl <- a |>
  summarise(plots = n(), EAs = n_distinct(cluster), median_pfp_n = median(pfp_n),
            q25 = quantile(pfp_n, 0.25), q75 = quantile(pfp_n, 0.75),
            district = first(adm2), .by = c(district_f, wave_f)) |>
  mutate(reliable = plots >= 20 & EAs >= 3) |>
  arrange(district_f, wave_f)
save_table(district_tbl, "2_2_pfp_by_district_wave")

get_boundaries <- function() {
  if (!requireNamespace("geodata", quietly = TRUE)) return(NULL)
  v <- tryCatch(geodata::gadm("MWI", level = 1, path = boundary_cache), error = \(e) NULL)
  if (is.null(v)) return(NULL)
  sf::st_as_sf(v) |> mutate(district_f = district_key(NAME_1))
}
bnd <- get_boundaries()

if (!is.null(bnd)) {
  map_df <- bnd |>
    select(district_f) |>
    left_join(district_tbl |> mutate(district_f = as.character(district_f)), by = "district_f",
              relationship = "one-to-many") |>
    filter(!is.na(wave_f)) |>
    mutate(median_pfp_n = if_else(reliable, median_pfp_n, NA_real_))
  p <- ggplot() +
    geom_sf(data = bnd, fill = "grey92", colour = "white", linewidth = 0.2) +
    geom_sf(data = map_df, aes(fill = median_pfp_n), colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(trans = "log10", na.value = "grey80", name = "median\nPFP-N") +
    facet_wrap(~wave_f, nrow = 1) +
    labs(title = "Median PFP-N by district and wave",
         subtitle = "Grey: fewer than 20 plots or 3 EAs") +
    theme(axis.text = element_blank(), panel.grid = element_blank())
  save_fig(p, "2_2_pfp_district_maps", w = 12, h = 6)
  unmatched <- setdiff(unique(as.character(district_tbl$district_f)), bnd$district_f)
  if (length(unmatched) > 0) message("Districts not matched to GADM names: ", paste(unmatched, collapse = ", "))
} else skip("2.2 district maps", "geodata not installed or GADM download failed (tables still written)")

# 2.3 Robust GAMs: space, time, and decomposition ------------------------------
# log PFP-N = log yield - log N, so wave and spatial effects on PFP-N split
# exactly into a yield part and an N-rate part. Scaled-t errors down-weight
# remaining extreme values. The EA random effect absorbs within-EA clustering.

g <- a |> filter(!is.na(latitude), !is.na(longitude))
k_sp <- max(10, min(100, floor(n_distinct(paste(g$latitude, g$longitude)) / 3)))

fit_gam <- function(formula, data) {
  tryCatch(mgcv::bam(formula, family = mgcv::scat(), data = data, method = "fREML", discrete = TRUE),
           error = \(e) { message("GAM failed: ", conditionMessage(e)); NULL })
}

if (nrow(g) >= min_rows && nlevels(droplevels(g$wave_f)) >= 1) {
  g <- g |> mutate(wave_f = droplevels(wave_f), cluster_f = droplevels(cluster_f))
  multi_wave <- nlevels(g$wave_f) > 1
  rhs <- paste(if (multi_wave) "wave_f +" else "",
               "s(longitude, latitude, k = k_sp) + s(cluster_f, bs = 're')")
  m_pfp   <- fit_gam(as.formula(paste("log_pfp ~", rhs)), g)
  m_y     <- fit_gam(as.formula(paste("log_y ~", rhs)), g)
  m_n     <- fit_gam(as.formula(paste("log_n ~", rhs)), g)
  m_pfp_N <- fit_gam(as.formula(paste("log_pfp ~", rhs, "+ s(log_n, k = 5)")), g)   # at a common N rate
  m_space <- fit_gam(as.formula(paste("log_pfp ~ s(longitude, latitude, k = k_sp) + s(cluster_f, bs = 're')")), g)
  m_time  <- if (multi_wave) fit_gam(log_pfp ~ wave_f + s(cluster_f, bs = "re"), g) else NULL

  wave_effects <- function(m, label) {
    if (is.null(m) || !multi_wave) return(NULL)
    pt <- summary(m)$p.table
    pt <- pt[str_starts(rownames(pt), "wave_f"), , drop = FALSE]
    tibble(model = label, wave = str_remove(rownames(pt), "^wave_f"),
           pct_vs_ref = 100 * (exp(pt[, 1]) - 1),
           ci_low = 100 * (exp(pt[, 1] - 1.96 * pt[, 2]) - 1),
           ci_high = 100 * (exp(pt[, 1] + 1.96 * pt[, 2]) - 1))
  }
  decomp <- bind_rows(wave_effects(m_pfp, "log PFP-N"), wave_effects(m_y, "log yield"),
                      wave_effects(m_n, "log N rate"), wave_effects(m_pfp_N, "log PFP-N at common N")) |>
    mutate(reference_wave = ref_wave)
  save_table(decomp, "2_3_wave_effects_decomposition")

  if (nrow(decomp) > 0) {
    p <- ggplot(decomp, aes(pct_vs_ref, wave, xmin = ci_low, xmax = ci_high, colour = model)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_pointrange(position = position_dodge(width = 0.6)) +
      labs(x = paste("% difference vs", ref_wave, "(robust GAM, spatial smooth + EA effect)"), y = NULL,
           colour = NULL, title = "What drives PFP-N change over time: yield or N rate?")
    save_fig(p, "2_3_wave_effects_decomposition")
  }

  dev_tbl <- tibble(
    model = c("space + time", "space only", "time only", "space + time + N rate"),
    fit = list(m_pfp, m_space, m_time, m_pfp_N)) |>
    mutate(deviance_explained = map_dbl(fit, \(m) if (is.null(m)) NA_real_ else summary(m)$dev.expl),
           AIC = map_dbl(fit, \(m) if (is.null(m)) NA_real_ else AIC(m))) |>
    select(-fit)
  save_table(dev_tbl, "2_3_space_vs_time_deviance")

  # Spatial pattern: multiplicative deviation from the national level
  if (!is.null(m_pfp)) {
    pts_sf <- sf::st_as_sf(distinct(g, longitude, latitude), coords = c("longitude", "latitude"), crs = 4326)
    grid <- expand_grid(longitude = seq(min(g$longitude), max(g$longitude), by = 0.05),
                        latitude  = seq(min(g$latitude),  max(g$latitude),  by = 0.05))
    grid_sf <- sf::st_as_sf(grid, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
    near <- lengths(sf::st_is_within_distance(grid_sf, pts_sf, dist = 25000)) > 0   # no extrapolation
    if (!is.null(bnd)) near <- near & lengths(sf::st_intersects(grid_sf, sf::st_union(bnd))) > 0
    grid <- grid[near, ] |>
      mutate(wave_f = factor(ref_wave, levels = levels(g$wave_f)),
             cluster_f = factor(levels(g$cluster_f)[1], levels = levels(g$cluster_f)))
    tm <- predict(m_pfp, newdata = grid, type = "terms", terms = "s(longitude,latitude)")
    grid$spatial_factor <- exp(tm[, 1])
    p <- ggplot(grid, aes(longitude, latitude, fill = spatial_factor)) +
      geom_raster() +
      { if (!is.null(bnd)) geom_sf(data = bnd, inherit.aes = FALSE, fill = NA, colour = "grey30", linewidth = 0.2) } +
      scale_fill_gradient2(trans = "log", midpoint = 0, low = "#b2182b", high = "#2166ac",
                           breaks = c(0.5, 0.75, 1, 1.33, 2), name = "x national\nlevel") +
      coord_sf() +
      labs(title = "Spatial pattern of PFP-N (robust GAM smooth)",
           subtitle = "Multiplicative deviation from the national level, all waves pooled",
           x = NULL, y = NULL)
    save_fig(p, "2_3_pfp_spatial_smooth", w = 5, h = 8)
    note("Spatial pattern: the GAM smooth ranges from x", round(min(grid$spatial_factor), 2),
         " to x", round(max(grid$spatial_factor), 2), " the national PFP-N level across Malawi.")
  }

  # Wave-specific spatial patterns (does the geography change over time?)
  if (multi_wave) {
    k_by <- max(10, floor(k_sp / 2))
    m_by <- fit_gam(log_pfp ~ wave_f + s(longitude, latitude, by = wave_f, k = k_by) +
                      s(cluster_f, bs = "re"), g)
    if (!is.null(m_by) && !is.null(m_pfp)) {
      aic_by <- AIC(m_by) - AIC(m_pfp)
      save_table(tibble(model = c("common spatial pattern", "wave-specific spatial patterns"),
                        AIC = c(AIC(m_pfp), AIC(m_by))), "2_3_spatial_pattern_stability")
      note("Spatial pattern stability: wave-specific smooths change AIC by ", round(aic_by, 1),
           if (aic_by < -10) " (the geography of PFP-N changes between waves)." else
             " (no strong evidence that the geography changes between waves).")
    }
  }

  gam_summaries <- list(pfp = m_pfp, yield = m_y, N = m_n, pfp_at_common_N = m_pfp_N) |>
    compact() |>
    map(\(m) capture.output(summary(m)))
  writeLines(unlist(imap(gam_summaries, \(x, nm) c(paste("=====", nm), x, ""))),
             file.path(tab_dir, "2_3_gam_summaries.txt"))

  if (!is.null(m_pfp) && multi_wave) {
    dd <- decomp |> filter(model %in% c("log PFP-N", "log yield", "log N rate"))
    note("Wave effects vs ", ref_wave, " (robust GAM): ",
         paste0(dd$model, " / ", dd$wave, ": ", sprintf("%+.0f%%", dd$pct_vs_ref), collapse = "; "),
         ". Because log PFP-N = log yield - log N, these show whether PFP-N changes come from yields ",
         "or from N rates.")
  }
} else skip("2.3 GAMs", "too few georeferenced LSMS rows")

# 2.4 Variance components ------------------------------------------------------

vc_data <- a |> filter(!is.na(district_f)) |> droplevels()
if (nrow(vc_data) >= min_rows) {
  vc_formula <- if (nlevels(vc_data$wave_f) > 1)
    log_pfp ~ wave_f + (1 | district_f) + (1 | cluster_f) + (1 | hh_f) else
    log_pfp ~ (1 | district_f) + (1 | cluster_f) + (1 | hh_f)
  vc_fit <- tryCatch(glmmTMB::glmmTMB(vc_formula, data = vc_data, family = gaussian()),
                     error = \(e) { message("Variance components failed: ", conditionMessage(e)); NULL })
  if (!is.null(vc_fit)) {
    vc <- glmmTMB::VarCorr(vc_fit)$cond
    fixed_var <- var(predict(vc_fit, re.form = NA))
    vc_tbl <- tibble(
      level = c("wave (fixed)", "district", "EA within district", "household within EA", "plot (residual)"),
      variance = c(fixed_var, as.numeric(vc$district_f), as.numeric(vc$cluster_f),
                   as.numeric(vc$hh_f), sigma(vc_fit)^2)) |>
      mutate(share_pct = round(100 * variance / sum(variance), 1))
    save_table(vc_tbl, "2_4_variance_components")
    note("Variance of log PFP-N: ", paste0(vc_tbl$level, " ", vc_tbl$share_pct, "%", collapse = ", "),
         ". Most variation between plots within the same place points to farm management and ",
         "measurement error rather than geography or year.")
  }
}

# 2.5 Sensitivity of wave effects to data-quality choices ----------------------

base <- d |>
  filter(program == "LSMS_MWI", !is.na(log_pfp), !is.na(wave)) |>
  mutate(wave_f = droplevels(wave), cluster_f = factor(cluster))

scenarios <- list(
  "Main: outliers removed"               = \(x) filter(x, !outlier),
  "Outliers kept"                        = \(x) x,
  "N rate >= low_N"                      = \(x) filter(x, !outlier, N_fertilizer >= low_N),
  "Enumerator kg only"                   = \(x) filter(x, !outlier, fert_kg_source %in% "enumerator_kg"),
  "No fertilizer kg mismatch"            = \(x) filter(x, !outlier, !fert_kg_mismatch %in% TRUE),
  "Main harvest block only"              = \(x) filter(x, !outlier, !harvest_block %in% "alt"),
  "Stand code observed (not inferred)"   = \(x) filter(x, !outlier, crop_stand %in% "sole"),
  "No quality flags at all"              = \(x) filter(x, !outlier, is.na(quality_flags))
)
names(scenarios) <- str_replace(names(scenarios), "low_N", paste(low_N, "kg/ha"))

sens <- imap(scenarios, \(f, nm) {
  x <- f(base) |> droplevels()
  if (nrow(x) < min_rows || nlevels(x$wave_f) < 2 || !ref_wave %in% levels(x$wave_f)) {
    return(tibble(scenario = nm, rows = nrow(x), wave = NA_character_, pct_vs_ref = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_))
  }
  x$wave_f <- relevel(x$wave_f, ref = ref_wave)
  m <- tryCatch(glmmTMB::glmmTMB(log_pfp ~ wave_f + (1 | cluster_f), data = x, family = gaussian()),
                error = \(e) NULL)
  if (is.null(m)) return(tibble(scenario = nm, rows = nrow(x)))
  ct <- summary(m)$coefficients$cond
  ct <- ct[str_starts(rownames(ct), "wave_f"), , drop = FALSE]
  tibble(scenario = nm, rows = nrow(x), wave = str_remove(rownames(ct), "^wave_f"),
         pct_vs_ref = 100 * (exp(ct[, 1]) - 1),
         ci_low = 100 * (exp(ct[, 1] - 1.96 * ct[, 2]) - 1),
         ci_high = 100 * (exp(ct[, 1] + 1.96 * ct[, 2]) - 1))
}) |> list_rbind() |> mutate(reference_wave = ref_wave)
save_table(sens, "2_5_sensitivity_wave_effects")

sens_plot <- filter(sens, !is.na(wave))
if (nrow(sens_plot) > 0) {
  p <- ggplot(sens_plot, aes(pct_vs_ref, fct_rev(factor(scenario, levels = names(scenarios))),
                             xmin = ci_low, xmax = ci_high)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(size = 0.3) +
    facet_wrap(~wave, nrow = 1) +
    labs(x = paste("% difference in PFP-N vs", ref_wave), y = NULL,
         title = "Are wave differences robust to data-quality choices?",
         subtitle = "Linear mixed model on log PFP-N with EA random effect; empty rows: wave absent in that subset")
  save_fig(p, "2_5_sensitivity_wave_effects", w = 11, h = 5)

  spread <- sens_plot |> summarise(range = max(pct_vs_ref) - min(pct_vs_ref), .by = wave)
  note("Sensitivity: across ", n_distinct(sens_plot$scenario), " data-quality scenarios, wave effects vary by ",
       paste0(spread$wave, " ", round(spread$range), " points", collapse = "; "),
       ". Conclusions that hold in every scenario are robust to the known data problems.")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

dir.create(file.path(out_dir, "objects"), showWarnings = FALSE)
saveRDS(list(created = Sys.time(), findings = findings, n_rows = nrow(d),
             fig_dir = normalizePath(fig_dir), tab_dir = normalizePath(tab_dir)),
        file.path(out_dir, "objects", "analysis1.rds"))

writeLines(c(
  paste0("# PFP-N of sole ", focus_crop, " in Malawi: data robustness and variation in space and time"),
  "",
  paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " from ", in_rds, "."),
  "Tables in `tables/`, figures in `figures/`. Key findings computed in this run:",
  "",
  paste0("- ", findings),
  "",
  "## How to read the robustness part",
  "",
  "- PFP-N (harvest kg / N kg on the same planted area) is immune to plot-area error, which is the",
  "  best documented weakness of survey yields; it remains exposed to errors in harvest quantity,",
  "  harvest unit conversion, fertilizer quantity and fertilizer N content.",
  "- Trials are a benchmark, not ground truth: they show what yields are attainable at a given N",
  "  rate under researcher management. Survey plots far above the trial frontier are more likely",
  "  measurement errors than exceptional farmers.",
  "",
  "## How to read the space-time part",
  "",
  "- Wave effects and spatial patterns are from robust GAMs (scaled-t errors) with an EA random",
  "  effect; the yield and N-rate models show which of the two drives PFP-N differences.",
  "- 'PFP-N at common N' removes the mechanical effect of N rate (PFP-N falls as N rises).",
  "- IHS6 plots are placed at district centroids (no EA coordinates released), so they add little",
  "  to the spatial smooth; IHS3 areas are farmer-reported, which affects yield and N rate but not PFP-N."
), file.path(out_dir, "summary.md"))
message("Done: see ", file.path(out_dir, "summary.md"))

# ==============================================================================
# END
# ==============================================================================