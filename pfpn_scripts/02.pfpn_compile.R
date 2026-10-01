# =============================================================================
# PFP-N compile: LSMS Malawi (.dta zips) + LCAS + Carob -> one terminag dataset
# -----------------------------------------------------------------------------
# Unit of analysis : plot-crop x yield source (farmer_report / crop_cut / ...)
# Output columns   : terminag names where a term exists (yield kg/ha,
#                    N_fertilizer kg/ha, plot_area m2, ...), plus a few extras
#                    (program, source, yield_source, area_source, geo_level,
#                    fert_kg_source, organic covariates, pfp_n, qc_*).
#
# Scope: mineral N only (organic fertilizer kept as covariates OM_used,
# organic_qty, organic_unit); IHS2 dropped; coordinates = EA coordinates where
# the wave provides them, district centroids (flagged) for IHS6.
#
# Workflow
#   1. Run pfpn_variable_lookup.R: it checks the per-wave variable spec against
#      each DDI and writes pfpn_config/source_map.csv. Review spec_check.csv and
#      replace any TODO rows.
#   2. Set the zip paths below and run this script.
#   3. Each run syncs the conversion tables in pfpn_config/ with the labels
#      found in the data. Rows with an empty factor show up in
#      pfpn_compiled/missing_conversions.csv. Fill them in and re-run.
#
# Provenance tags: [verified] checked against a source; [recall] from training
# knowledge, confirm before relying on it.
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (haven, readxl, terra, geodata)

source("00.pfpn_utils.R", local = TRUE)   # paths, update-aware downloads, zip extraction, Malawi check

lookup_dir <- pfpn_path("pfpn_lookup")     # outputs of pfpn_variable_lookup.R
cfg_dir    <- pfpn_path("pfpn_config")     # mapping + conversion tables you curate
out_dir    <- pfpn_path("pfpn_compiled")
walk(c(cfg_dir, out_dir, file.path(out_dir, "raw")), dir.create,
     recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. Configuration: sources
#    `source` must match the source column of source_map.csv.
#    date = year of the reference rainy season harvest; adjust if you prefer
#    the survey year.
# -----------------------------------------------------------------------------

# Where to get the LSMS data
# -----------------------------------------------------------------------------
# The Malawi Integrated Household Surveys (IHS) are free public-use data from
# the World Bank Microdata Library. Downloading needs a (free) account and
# accepting the terms of use; then, on each study page, open "Get Microdata"
# and download the Stata version. Keep the zip as downloaded: this script
# reads the .dta files directly from it.
#   IHS3 2010-11  https://microdata.worldbank.org/catalog/1003
#   IHS4 2016-17  https://microdata.worldbank.org/catalog/2936
#   IHS5 2019-20  https://microdata.worldbank.org/catalog/3818
#   IHS6 2024-25  https://microdata.worldbank.org/catalog/8507
#   IHPS 2013     https://microdata.worldbank.org/catalog/2248
#                 (Integrated Household Panel Survey 2010-2013, short-term panel)
# The variable documentation of each wave is at <study page>/data-dictionary
# (e.g. https://microdata.worldbank.org/catalog/8507/data-dictionary).
# IHS2 2004-05 (catalog 2307) is not used: no plot-level harvest.
# IHPS 2010-2013: households of 204 of IHS3's 768 EAs, re-interviewed in 2013
# (rainy season 2012/13). Only its 2013 round (*_13 files) is used; its 2010
# round (*_10 files) repeats IHS3 households already read from catalog 1003.
# -----------------------------------------------------------------------------

# Provide the path to the ZIP file for each LSMS wave (to be modified depending on local)
lsms_mwi_2010_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2010/MWI_2010_IHS-III_v01_M_STATA8.zip'
lsms_mwi_2013_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2013/MWI_2010-2013_IHPS_v01_M_Stata.zip'
lsms_mwi_2016_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2016/MWI_2016_IHS-IV_v04_M_STATA14.zip'
lsms_mwi_2019_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2019/MWI_2019_IHS-V_v06_M_Stata.zip'
lsms_mwi_2024_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2024/MWI_2024-2025_IHS-VI_v01_M_STATA14.zip'

lsms_waves <- tribble(
  ~source,                ~zip,               ~date,  ~season, ~country,
  "LSMS_1003_2010-2011",  lsms_mwi_2010_zip,  "2010", "wet",   "Malawi",
  "LSMS_2936_2016-2017",  lsms_mwi_2016_zip,  "2016", "wet",   "Malawi",
  "LSMS_3818_2019-2020",  lsms_mwi_2019_zip,  "2019", "wet",   "Malawi",
  "LSMS_8507_2024-2025",  lsms_mwi_2024_zip,  "2024", "wet",   "Malawi",
  "LSMS_2248_2013",       lsms_mwi_2013_zip,  "2013", "wet",   "Malawi"   # IHPS 2013 round
)

# Waves whose zip is not found are skipped, and listed so the skip is visible.
# Relative paths are resolved from the working directory.
missing_zip <- lsms_waves |> filter(!file.exists(zip))
if (nrow(missing_zip) > 0) {
  message("LSMS zip not found, wave skipped (working directory: ", getwd(), "):\n",
          paste0("  ", missing_zip$source, " -> ", missing_zip$zip, collapse = "\n"))
}
lsms_waves <- lsms_waves |> filter(file.exists(zip))

# Waves without EA coordinates: fill latitude/longitude with the centroid of
# the household's district (adm2), flagged geo_level = "district_centroid".
centroid_sources <- c("LSMS_8507_2024-2025")
centroid_country <- "MWI"   # GADM code; level 1 polygons = districts [recall, checked below]

# Crop stand. PFP-N is computed for sole crops only; intercrops are dropped.
#   LSMS rows without a stand code are classed from the number of crops on the plot.
#   unknown_stand: rows from sources with no stand information at all (LCAS as
#   mapped, Carob datasets without intercrop columns):
#     "include_flagged" keeps them, crop_stand = "unknown"; "exclude" drops them.
unknown_stand <- "include_flagged"

# Flag (not drop) PFP-N computed from very small N rates, where the ratio explodes
low_N_threshold <- 10   # kg N/ha

# One row per LCAS dataset (csv or xlsx), already renamed to standard LCAS
# names with rename_lcas.R. All use the mapping rows with source == "LCAS".
lcas_datasets <- tribble(
  ~source,        ~path,                        ~date,  ~season, ~country
  # "LCAS_eth_2023", "data/lcas/lcas_eth_2023.csv", "2023", "wet",  "Ethiopia"
)

# Carob collections to include [verified: URL pattern from caramba source]
carob_groups    <- c("agronomy", "survey")
carob_countries <- 'Malawi'         # e.g. c("Malawi", "Zambia"); NULL keeps all
carob_url <- function(group) sprintf("https://geodata.ucdavis.edu/carob/carob_%s_latest-cc.zip", group)

# This is just a selection of CAROB datasets
# carob_mwi_list <- c('doi_10.25502_20180814_0923_HJ',
#                     'doi_10.25502_20180814_1554_HJ',
#                     'doi_10.34725_DVN_25746', 
#                     'doi_10.5061_dryad.fg15tg2', 
#                     'doi_10.7910_DVN_UNLRGC',
#                     'doi_10.7910_DVN_GXUNAZ',
#                     'doi_10.7910_DVN_QLJUY7',
#                     'doi_10.7910_DVN_UJIPSW',
#                     'doi_10.7910_DVN_BUPNF4',
#                     'hdl_11529_10868',
#                     'doi_10.25502_egvm-8g12_d',
#                     'doi_10.7910_DVN_O71ROU',
#                     'doi_10.25502_syv0-ye87_d',
#                     'doi_10.25502_HNKM-Y645_D',
#                     'doi_10.71682_10549300',
#                     'doi_10.71682_10549314',
#                     'doi_10.7910_DVN_1T4Q3F',
#                     'doi_10.7910_DVN_NZW56Q',
#                     'doi_10.7910_DVN_1A6WMD',
#                     'hdl_11529_10549178',
#                     'doi_10.7910_DVN_ILVSXJ'
#                     )

# Helpers (ensure, fetch, unzip_once, is_zip, usable, in_malawi, district_key) come from pfpn_utils.R

# -----------------------------------------------------------------------------
# 1. Roles: what each mapped source variable means
# -----------------------------------------------------------------------------

roles <- tribble(
  ~role,               ~level,     ~kind,
  "hhid",              "key",      "key",
  "plot_id",           "key",      "key",
  "crop",              "plotcrop", "cat",   # joins on the label, not the code
  "season",            "plotcrop", "cat",
  "harvest_qty",       "plotcrop", "num",
  "harvest_qty_kg",    "plotcrop", "num",   # harvest already in kg (after `multiply`), no unit variable
  "harvest_qty_alt",   "plotcrop", "num",   # second harvest block (IHS3 ag_g09); used where main is empty
  "harvest_unit_alt",  "plotcrop", "cat",
  "harvest_condition_alt", "plotcrop", "cat",
  "harvest_unit",      "plotcrop", "cat",
  "harvest_condition", "plotcrop", "cat",
  "conv_region",       "plotcrop", "cat",
  "yield_reported",    "plotcrop", "num",   # already per ha; use `multiply` to reach kg/ha
  "cropcut_yield",     "plotcrop", "num",   # already per ha; use `multiply` to reach kg/ha
  "yield_moisture",    "plotcrop", "num",
  "crop_cut",          "plotcrop", "cat",
  "intercrop_fraction","plotcrop", "cat",   # numeric or labelled; handled below
  "intercropped",      "plotcrop", "cat",
  "plot_area_gps",     "plot",     "num",   # `multiply` REQUIRED (unit -> m2)
  "plot_area_reported","plot",     "num",
  "plot_area_reported_m2", "plot",  "num",  # reported area already in one unit (use `multiply` -> m2)
  "area_unit",         "plot",     "cat",
  "fertilizer_used",   "plot",     "cat",
  "irrigated",         "plot",     "cat",
  "fertilizer_type",   "fert",     "cat",
  "fertilizer_amount", "fert",     "num",
  "fertilizer_unit",   "fert",     "cat",
  "fertilizer_kg",     "fert",     "num",   # enumerator-computed total kg; preferred over amount x unit
  "OM_used",           "plot",     "cat",   # organic fertilizer: covariate only
  "organic_qty",       "plot",     "num",
  "organic_unit",      "plot",     "cat",
  "whole_plot_planted","plotcrop", "cat",   # yes -> intercrop_fraction = 1
  "intercrop_percent", "plotcrop", "num",   # % of plot under the crop (IHS6 ag_g11_2)
  "ea_id",             "hh",       "cat",
  "adm1",              "hh",       "cat",
  "adm2",              "hh",       "cat",
  "adm3",              "hh",       "cat",
  "latitude",          "hh",       "num",
  "longitude",         "hh",       "num"
)

# -----------------------------------------------------------------------------
# 2. Source map (written by pfpn_variable_lookup.R)
# -----------------------------------------------------------------------------
# Columns: source, file, var, role, app, multiply, fert_type_fixed, note
#   app             fertilizer application / product index (1, 2, ...) for fert roles
#   multiply        numeric factor applied to num roles (e.g. 1000 for t/ha -> kg/ha,
#                   4046.856 for acres -> m2 on plot_area_gps)
#   fert_type_fixed fertilizer name when a column holds one product only
#   file            LSMS: data file name (extension ignored); LCAS: ignored

map_file <- file.path(cfg_dir, "source_map.csv")
if (!file.exists(map_file)) stop("Run pfpn_variable_lookup.R first: it writes ", map_file)

src_map <- read_csv(map_file, col_types = cols(.default = "c"), show_col_types = FALSE) |>
  filter(!is.na(role))

todo <- src_map |> filter(is.na(var) | str_detect(var, "^(TODO|\\?)"))
if (nrow(todo) > 0) {
  message(nrow(todo), " source_map rows still TODO (skipped):")
  print(select(todo, source, file, role, app, any_of("note")), n = Inf)
}

src_map <- src_map |>
  filter(!is.na(var), !str_detect(var, "^(TODO|\\?)")) |>
  mutate(app = as.integer(app), multiply = as.numeric(multiply)) |>
  left_join(roles, by = "role")

bad_roles <- src_map |> filter(is.na(level)) |> distinct(role)
if (nrow(bad_roles) > 0) stop("Unknown roles in source_map.csv: ", paste(bad_roles$role, collapse = ", "))

no_mult <- src_map |> filter(role == "plot_area_gps", is.na(multiply))
if (nrow(no_mult) > 0) stop("Set `multiply` (unit -> m2) for plot_area_gps rows: ",
                            paste(no_mult$source, no_mult$var, collapse = "; "))

# Column name each mapped variable gets after reading
src_map <- src_map |>
  mutate(
    new_name = case_when(
      level == "fert" ~ paste0(role, "__app", coalesce(app, 1L)),
      kind == "key"   ~ paste0(role, "__k", row_number()),
      .default = role
    ),
    .by = c(source, file)
  )

dups <- src_map |> count(source, file, new_name) |> filter(n > 1)
if (nrow(dups) > 0) stop("A non-key role is mapped twice within one file: ",
                         paste(dups$source, dups$file, dups$new_name, collapse = "; "))

# -----------------------------------------------------------------------------
# 3. Readers
# -----------------------------------------------------------------------------

# Labelled (Stata) columns -> label text for categorical roles, numbers otherwise
decode <- function(x, kind) {
  if (kind == "num") return(suppressWarnings(as.numeric(haven::zap_labels(x))))
  if (kind == "key") return(as.character(haven::zap_labels(x)))
  if (haven::is.labelled(x)) as.character(haven::as_factor(x, levels = "default"))
  else as.character(x)
}

apply_map <- function(df, m) {
  lower <- setNames(names(df), str_to_lower(names(df)))
  m <- m |> mutate(actual = unname(lower[str_to_lower(var)]))
  missing <- m |> filter(is.na(actual))
  if (nrow(missing) > 0) warning("Variables not found: ",
                                 paste(missing$source, missing$file, missing$var, collapse = "; "))
  m <- m |> filter(!is.na(actual))
  out <- pmap(m, \(actual, new_name, kind, multiply, ...) {
    v <- decode(df[[actual]], kind)
    if (kind == "num" && !is.na(multiply)) v <- v * multiply
    v
  })
  names(out) <- m$new_name
  as_tibble(out)
}

read_lsms_file <- function(zip, file, m) {
  entries <- utils::unzip(zip, list = TRUE)$Name
  # DDI file names may carry .NSDstat (IHS2/IHS3) or no extension; the zip holds .dta
  file_key <- str_to_lower(tools::file_path_sans_ext(file))
  hit <- entries[str_to_lower(tools::file_path_sans_ext(basename(entries))) == file_key &
                 str_detect(str_to_lower(entries), "\\.dta$")]
  if (length(hit) == 0) {
    warning(sprintf("%s: %s.dta not found in zip", basename(zip), file))
    return(NULL)
  }
  # Some zips hold the same module twice (IHS3: Full_Sample/ and Panel/).
  # Use the full cross-section; the panel households are covered by IHPS.
  if (length(hit) > 1) {
    pick <- c(hit[str_detect(hit, "(?i)full_sample")], hit[!str_detect(hit, "(?i)panel")], hit)[1]
    message(basename(zip), ": ", length(hit), " copies of ", file, ".dta; using ", pick)
    hit <- pick
  }
  # extracted once to pfpn_compiled/raw/lsms/<zip name>/ and reused afterwards;
  # re-extracted only if the copy on disk cannot be read
  exdir <- file.path(out_dir, "raw", "lsms", tools::file_path_sans_ext(basename(zip)))
  path  <- unzip_once(zip, exdir, files = hit[1])
  nm <- tryCatch(names(haven::read_dta(path, n_max = 0)), error = \(e) NULL)
  if (is.null(nm)) {
    message("Extracted copy unreadable, extracting again: ", basename(path))
    utils::unzip(zip, files = hit[1], exdir = exdir, overwrite = TRUE)
    nm <- names(haven::read_dta(path, n_max = 0))
  }
  keep <- nm[str_to_lower(nm) %in% str_to_lower(m$var)]
  haven::read_dta(path, col_select = all_of(keep)) |> apply_map(m)
}

read_lcas_file <- function(path, m) {
  df <- if (str_detect(str_to_lower(path), "\\.xlsx?$")) readxl::read_excel(path)
        else read_csv(path, show_col_types = FALSE, guess_max = 1e5)
  apply_map(df, m)
}

# Combine repeated key columns (e.g. gardenid + plotid) into one key
unite_keys <- function(df) {
  for (k in c("hhid", "plot_id", "crop")) {
    cols <- names(df)[str_detect(names(df), paste0("^", k, "__k"))]
    if (length(cols) > 0) df <- df |> unite(!!k, all_of(cols), sep = "_", na.rm = TRUE)
  }
  df
}

# Split one source's tables into hh / plot / plotcrop / fert levels
split_levels <- function(tables, m, src_meta) {
  lvl_cols <- function(lvl) m |> filter(level == lvl) |> pull(role) |> unique()
  key_by <- list(hh = "hhid", plot = c("hhid", "plot_id"),
                 plotcrop = c("hhid", "plot_id", "crop"), fert = c("hhid", "plot_id"))

  pieces <- map(names(key_by), \(lvl) {
    wanted <- lvl_cols(lvl)
    parts <- tables |>
      map(\(t) {
        cols <- if (lvl == "fert") names(t)[str_detect(names(t), "__app\\d+$")]
                else intersect(wanted, names(t))
        keys <- intersect(key_by[[lvl]], names(t))
        if (length(cols) == 0 || !all(key_by[[lvl]] %in% names(t))) return(NULL)
        t |> select(all_of(c(keys, cols))) |> distinct()
      }) |>
      compact()
    if (length(parts) == 0) return(NULL)
    reduce(parts, \(a, b) full_join(a, b, by = key_by[[lvl]], relationship = "many-to-many"))
  })
  names(pieces) <- names(key_by)
  pieces
}

# Wide fert columns (role__appN) -> long, one row per application / product
fert_long <- function(fert, m) {
  if (is.null(fert)) return(NULL)
  fixed <- m |> filter(!is.na(fert_type_fixed)) |> transmute(app = coalesce(app, 1L), fert_type_fixed)
  fert |>
    pivot_longer(matches("__app\\d+$"), names_to = c(".value", "app"),
                 names_pattern = "(.*)__app(\\d+)$") |>
    mutate(app = as.integer(app)) |>
    left_join(distinct(fixed), by = "app") |>
    ensure(c("fertilizer_type", "fertilizer_amount", "fertilizer_kg")) |>
    mutate(
      from_fixed      = is.na(fertilizer_type) & !is.na(fert_type_fixed),
      fertilizer_type = coalesce(fertilizer_type, fert_type_fixed),
      has_amount      = !is.na(fertilizer_amount) | !is.na(fertilizer_kg)
    ) |>
    # a product-specific column only counts when it holds an amount
    filter(if_else(from_fixed, has_amount, has_amount | !is.na(fertilizer_type))) |>
    select(-has_amount) |>
    select(-fert_type_fixed, -from_fixed)
}

read_source <- function(source, reader_arg, reader, src_meta) {
  m <- src_map |> filter(source == !!src_meta$map_source)
  if (nrow(m) == 0) { warning("No mapping rows for ", src_meta$map_source); return(NULL) }
  tables <- m |>
    group_split(file) |>
    map(\(mf) reader(reader_arg, unique(mf$file), mf)) |>
    compact() |>
    map(unite_keys) |>
    map(\(t) {                       # LCAS: one plot per respondent
      if (!"plot_id" %in% names(t)) t$plot_id <- "plot1"
      if (!"hhid" %in% names(t)) t$hhid <- as.character(seq_len(nrow(t)))
      t
    })
  lv <- split_levels(tables, m, src_meta)
  lv$fert <- fert_long(lv$fert, m)
  lv <- map(lv, \(x) if (is.null(x)) NULL else mutate(x, source = source))
  lv$meta <- src_meta
  lv
}

lsms <- lsms_waves |>
  pmap(\(source, zip, date, season, country) {
    read_source(source, zip, \(z, f, mf) read_lsms_file(z, f, mf),
                tibble(source, date, season, country, program = "LSMS_MWI", map_source = source))
  })

lcas <- lcas_datasets |>
  pmap(\(source, path, date, season, country) {
    read_source(source, path, \(p, f, mf) read_lcas_file(p, mf),
                tibble(source, date, season, country, program = "LCAS", map_source = "LCAS"))
  })

survey <- c(lsms, lcas) |> compact()
if (length(survey) == 0) {
  stop("No LSMS or LCAS source could be read. Check the zip paths in `lsms_waves` ",
       "(working directory: ", getwd(), ") and the paths in `lcas_datasets`.", call. = FALSE)
}
stack  <- function(lvl) survey |> map(lvl) |> compact() |> list_rbind()

hh       <- stack("hh")
plot     <- stack("plot")
plotcrop <- stack("plotcrop")

# Second harvest block (IHS3): fill quantity, unit and condition together where
# the main block is empty, so a unit never gets paired with the other block's quantity
if (nrow(plotcrop) > 0) {
  plotcrop <- plotcrop |>
    ensure(c("harvest_qty", "harvest_unit", "harvest_condition",
             "harvest_qty_alt", "harvest_unit_alt", "harvest_condition_alt")) |>
    mutate(
      harvest_block     = case_when(!is.na(harvest_qty) ~ "main",
                                    !is.na(harvest_qty_alt) ~ "alt", .default = NA_character_),
      use_alt           = harvest_block %in% "alt",
      harvest_qty       = if_else(use_alt, as.numeric(harvest_qty_alt), as.numeric(harvest_qty)),
      harvest_unit      = if_else(use_alt, harvest_unit_alt, harvest_unit),
      harvest_condition = if_else(use_alt, harvest_condition_alt, harvest_condition)
    ) |>
    select(-use_alt, -harvest_qty_alt, -harvest_unit_alt, -harvest_condition_alt)
  if (all(is.na(plotcrop$harvest_condition))) plotcrop$harvest_condition <- NULL
  if (all(is.na(plotcrop$harvest_unit)))      plotcrop$harvest_unit <- NULL
}
fert     <- stack("fert")
meta     <- stack("meta") |> select(-any_of("map_source"))

# -----------------------------------------------------------------------------
# 4. Conversion tables: sync with labels observed in the data
# -----------------------------------------------------------------------------

# ci = TRUE: compare label keys case- and space-insensitively ("UREA" == "urea")
ci_key <- function(x) str_to_lower(str_squish(x))

sync_table <- function(name, observed, value_col, prefill = NULL, ci = FALSE) {
  path <- file.path(cfg_dir, name)
  keys <- setdiff(names(observed), c(value_col, "basis"))
  if (ci) observed <- observed[!duplicated(ci_key(do.call(paste, observed[keys]))), ]
  observed <- distinct(observed)
  if (!value_col %in% names(observed)) observed[[value_col]] <- NA_real_
  if (!"basis" %in% names(observed)) observed$basis <- NA_character_
  if (!is.null(prefill)) observed <- prefill(observed)
  if (file.exists(path)) {
    old <- read_csv(path, col_types = cols(.default = "c"), show_col_types = FALSE) |>
      mutate(across(all_of(value_col), as.numeric))
    new <- if (ci) {
      observed |>
        filter(!ci_key(do.call(paste, observed[keys])) %in% ci_key(do.call(paste, old[keys])))
    } else anti_join(observed, old, by = keys)
    new <- new |> mutate(across(all_of(keys), as.character))
    if (nrow(new) > 0) message(name, ": ", nrow(new), " new rows appended")
    tbl <- bind_rows(old, new)
  } else {
    message(name, ": created with ", nrow(observed), " rows")
    tbl <- observed
  }
  write_csv(tbl, path, na = "")
  tbl
}

# Unit label -> kg (weight parsed from the label; crop condition NOT applied)
prefill_kg <- function(df, label_col) {
  l <- str_to_lower(df[[label_col]])
  num <- suppressWarnings(as.numeric(str_match(l, "(\\d+(?:\\.\\d+)?)\\s*(?:kg|kilo)")[, 2]))
  guess <- case_when(
    !is.na(num)                                ~ num,
    str_detect(l, "^\\s*(kgs?|kilo ?grams?)\\b") ~ 1,
    str_detect(l, "^\\s*grams?\\b")             ~ 0.001,
    str_detect(l, "^\\s*(tonnes?|tons?)\\b")      ~ 1000,
    .default = NA_real_
  )
  df |> mutate(kg_per_unit = coalesce(kg_per_unit, guess),
               basis = if_else(is.na(basis) & !is.na(guess), "parsed from label", basis))
}

# Area unit label -> m2. Unit definitions are exact; local units must be filled in.
prefill_m2 <- function(df) {
  l <- str_to_lower(df$area_unit)
  guess <- case_when(
    str_detect(l, "acre")                   ~ 4046.8564224,
    str_detect(l, "hectare|\\bha\\b")        ~ 10000,
    str_detect(l, "square met|\\bm2\\b|sq\\.? ?m") ~ 1,
    .default = NA_real_
  )
  df |> mutate(m2_per_unit = coalesce(m2_per_unit, guess),
               basis = if_else(is.na(basis) & !is.na(guess), "unit definition", basis))
}

# Crop-share label -> fraction (simple fractions, halves/quarters, percentages)
prefill_frac <- function(df) {
  l <- str_to_lower(df$intercrop_fraction)
  fr <- str_match(l, "(\\d)\\s*/\\s*(\\d)")
  pct <- suppressWarnings(as.numeric(str_match(l, "(\\d+(?:\\.\\d+)?)\\s*%")[, 2]))
  guess <- case_when(
    str_detect(l, "less than|more than|between") ~ NA_real_,   # ranges: decide yourself
    !is.na(fr[, 2]) ~ as.numeric(fr[, 2]) / as.numeric(fr[, 3]),
    !is.na(pct)     ~ pct / 100,
    str_detect(l, "almost all|whole|entire|all") ~ 1,
    str_detect(l, "half")    ~ 0.5,
    str_detect(l, "quarter") ~ 0.25,
    .default = suppressWarnings(as.numeric(l))
  )
  df |> mutate(fraction = coalesce(fraction, if_else(guess > 1, guess / 100, guess)),
               basis = if_else(is.na(basis) & !is.na(fraction), "parsed from label", basis))
}

cols_present <- function(df, cols) intersect(cols, names(df))

harvest_keys <- cols_present(plotcrop, c("crop", "harvest_unit", "harvest_condition", "conv_region"))
harvest_conv <- if ("harvest_qty" %in% names(plotcrop)) {
  plotcrop |>
    filter(!is.na(harvest_qty)) |>
    distinct(source, across(all_of(harvest_keys))) |>
    sync_table("harvest_unit_kg.csv", observed = _, "kg_per_unit",
               prefill = \(d) if ("harvest_unit" %in% names(d)) prefill_kg(d, "harvest_unit") else d)
} else NULL

fert_unit_conv <- if (!is.null(fert) && "fertilizer_unit" %in% names(fert)) {
  fert |> distinct(source, fertilizer_unit) |> filter(!is.na(fertilizer_unit)) |>
    sync_table("fert_unit_kg.csv", observed = _, "kg_per_unit",
               prefill = \(d) prefill_kg(d, "fertilizer_unit"))
} else NULL

area_conv <- if ("area_unit" %in% names(plot)) {
  plot |> distinct(source, area_unit) |> filter(!is.na(area_unit)) |>
    sync_table("area_unit_m2.csv", observed = _, "m2_per_unit", prefill = prefill_m2)
} else NULL

frac_conv <- if ("intercrop_fraction" %in% names(plotcrop)) {
  plotcrop |> distinct(source, intercrop_fraction) |> filter(!is.na(intercrop_fraction)) |>
    sync_table("intercrop_fraction.csv", observed = _, "fraction", prefill = prefill_frac)
} else NULL

# Fertilizer N and crop names start from the lookup script's stubs
seed_from_lookup <- function(cfg_name, lookup_name) {
  dest <- file.path(cfg_dir, cfg_name)
  src  <- file.path(lookup_dir, lookup_name)
  if (!file.exists(dest) && file.exists(src)) invisible(file.copy(src, dest))
}
seed_from_lookup("fertilizer_n.csv", "fertilizer_n_content_stub.csv")
seed_from_lookup("crop_map.csv", "crop_map_stub.csv")

fert_n <- if (!is.null(fert) && "fertilizer_type" %in% names(fert)) {
  fert |> distinct(value_label = fertilizer_type) |> filter(!is.na(value_label)) |>
    sync_table("fertilizer_n.csv", observed = _, "N_pct", ci = TRUE)
} else NULL

crop_map <- if ("crop" %in% names(plotcrop)) {
  plotcrop |> distinct(value_label = crop) |> filter(!is.na(value_label)) |>
    sync_table("crop_map.csv", observed = _, "fresh_moisture", ci = TRUE) |>
    ensure(c("terminag_crop", "yield_part"))   # fill terminag_crop for unmatched labels
} else NULL

# -----------------------------------------------------------------------------
# 5. Compute: area, harvest, yields, fertilizer N
# -----------------------------------------------------------------------------

to_logical <- function(x) {
  l <- str_to_lower(str_squish(as.character(x)))
  case_when(
    str_detect(l, "^(yes|y|true|t|1)$") ~ TRUE,
    str_detect(l, "^(no|n|false|f|0|2)$") ~ FALSE,   # LSMS/ODK often code 2 = no [recall]
    .default = NA
  )
}

# Plot area (m2): GPS first, farmer-reported second
plot_std <- plot |>
  ensure(c("plot_area_gps", "plot_area_reported", "plot_area_reported_m2", "area_unit",
           "fertilizer_used", "irrigated",
           "OM_used", "organic_qty", "organic_unit")) |>
  mutate(plot_area_gps = as.numeric(plot_area_gps),
         plot_area_reported = as.numeric(plot_area_reported)) |>
  left_join(if (is.null(area_conv)) tibble(source = character(), area_unit = character(), m2_per_unit = numeric())
            else select(area_conv, source, area_unit, m2_per_unit),
            by = c("source", "area_unit"), na_matches = "never") |>
  mutate(
    area_rep_m2 = coalesce(plot_area_reported * m2_per_unit, as.numeric(plot_area_reported_m2)),
    plot_area   = coalesce(if_else(plot_area_gps > 0, plot_area_gps, NA_real_), area_rep_m2),
    area_source = case_when(!is.na(plot_area_gps) & plot_area_gps > 0 ~ "gps",
                            !is.na(area_rep_m2) ~ "farmer_report",
                            .default = NA_character_),
    fertilizer_used = to_logical(fertilizer_used),
    irrigated       = to_logical(irrigated),
    # organic fertilizer: covariates only, not part of N_fertilizer
    OM_used         = coalesce(to_logical(OM_used), if_else(as.numeric(organic_qty) > 0, TRUE, NA)),
    organic_qty     = as.numeric(organic_qty)
  ) |>
  mutate(plot_area_gps_m2 = if_else(plot_area_gps > 0, plot_area_gps, NA_real_),
         plot_area_reported_m2 = area_rep_m2) |>
  select(source, hhid, plot_id, plot_area, area_source, plot_area_gps_m2, plot_area_reported_m2,
         fertilizer_used, irrigated, OM_used, organic_qty, organic_unit)

# Fertilizer: kg product -> kg N, summed per plot
fert_plot <- if (is.null(fert) || nrow(fert) == 0) {
  tibble(source = character(), hhid = character(), plot_id = character())
} else {
  fert |>
    ensure(c("fertilizer_unit", "fertilizer_type", "fertilizer_kg")) |>
    mutate(fertilizer_amount = as.numeric(fertilizer_amount),
           fertilizer_kg     = as.numeric(fertilizer_kg)) |>
    left_join(if (is.null(fert_unit_conv)) tibble(source = character(), fertilizer_unit = character(), kg_per_unit = numeric())
              else select(fert_unit_conv, source, fertilizer_unit, kg_per_unit),
              by = c("source", "fertilizer_unit"), na_matches = "never") |>
    mutate(
      kg_from_unit = fertilizer_amount * if_else(is.na(fertilizer_unit), 1, kg_per_unit),
      kg_enum      = if_else(fertilizer_kg > 0, fertilizer_kg, NA_real_),
      kg_product   = coalesce(kg_enum, kg_from_unit),
      kg_source    = case_when(!is.na(kg_enum) ~ "enumerator_kg",
                               !is.na(kg_from_unit) ~ "qty_x_unit",
                               .default = NA_character_),
      # >10 % disagreement between the two when both exist
      kg_mismatch  = !is.na(kg_enum) & !is.na(kg_from_unit) & kg_from_unit > 0 &
                     abs(kg_enum - kg_from_unit) / kg_from_unit > 0.1
    ) |>
    mutate(fert_key = ci_key(fertilizer_type)) |>
    left_join(if (is.null(fert_n)) tibble(fert_key = character(), N_pct = numeric(), terminag_fertilizer = character())
              else fert_n |> transmute(fert_key = ci_key(value_label), N_pct, across(any_of("terminag_fertilizer"))) |>
                     distinct(fert_key, .keep_all = TRUE),
              by = "fert_key", na_matches = "never") |>
    mutate(kg_N = kg_product * N_pct / 100) |>
    summarise(
      fert_kg         = sum(kg_product, na.rm = TRUE),
      N_kg            = sum(kg_N, na.rm = TRUE),
      n_complete      = all(!is.na(kg_N) | coalesce(kg_product, 0) == 0),
      fert_kg_source  = paste(sort(unique(na.omit(kg_source))), collapse = "; "),
      fert_kg_mismatch = any(kg_mismatch),
      fertilizer_type = paste(sort(unique(na.omit(coalesce(terminag_fertilizer, fertilizer_type)))),
                              collapse = "; "),
      .by = c(source, hhid, plot_id)
    )
}

# Crop names -> terminag
crop_std <- if (is.null(crop_map)) NULL else
  crop_map |> transmute(crop_key = ci_key(value_label), terminag_crop, across(any_of("yield_part"))) |>
  distinct(crop_key, .keep_all = TRUE)

# Plot-crop: harvest -> kg; area share; yields per source
pc <- plotcrop |>
  ensure(c("crop", "harvest_qty", "harvest_qty_kg", "harvest_unit", "harvest_condition", "conv_region",
           "yield_reported", "cropcut_yield", "yield_moisture", "crop_cut",
           "intercrop_fraction", "intercrop_percent", "whole_plot_planted", "intercropped", "season")) |>
  mutate(across(c(harvest_qty, harvest_qty_kg, yield_reported, cropcut_yield, yield_moisture,
                  intercrop_percent), as.numeric))

if (!is.null(harvest_conv)) {
  pc <- pc |> left_join(select(harvest_conv, source, all_of(harvest_keys), kg_per_unit),
                        by = c("source", harvest_keys), na_matches = "never")
} else pc$kg_per_unit <- NA_real_

frac_num <- suppressWarnings(as.numeric(pc$intercrop_fraction))
if (!is.null(frac_conv)) {
  pc <- pc |> left_join(select(frac_conv, source, intercrop_fraction, fraction),
                        by = c("source", "intercrop_fraction"), na_matches = "never")
} else pc$fraction <- NA_real_

pc <- pc |>
  mutate(
    # priority: % of plot (numeric) > "whole plot planted" = yes > share label table
    intercrop_fraction = coalesce(
      if_else(between(intercrop_percent, 0, 100), intercrop_percent / 100, NA_real_),
      if_else(to_logical(whole_plot_planted) %in% TRUE, 1, NA_real_),
      if_else(frac_num > 1, frac_num / 100, frac_num),
      fraction
    ),
    intercropped = case_when(
      str_detect(str_to_lower(intercropped), "inter|mix") ~ TRUE,
      str_detect(str_to_lower(intercropped), "pure|mono|sole") ~ FALSE,
      .default = to_logical(intercropped)
    ),
    harvest_kg = coalesce(harvest_qty * kg_per_unit, harvest_qty_kg)
  ) |>
  # Crop stand: stand code first; for LSMS plots without one, the number of
  # crops listed on the plot in module G decides
  mutate(n_crops_plot = n_distinct(crop), .by = c(source, hhid, plot_id)) |>
  mutate(
    crop_stand = case_when(
      intercropped %in% TRUE  ~ "intercrop",
      intercropped %in% FALSE ~ "sole",
      str_starts(source, "LSMS") & n_crops_plot == 1 ~ "sole (inferred: only crop on plot)",
      str_starts(source, "LSMS") & n_crops_plot > 1  ~ "intercrop (inferred: several crops on plot)",
      .default = "unknown"),
    intercropped = if_else(crop_stand == "unknown", NA, str_starts(crop_stand, "intercrop"))
  ) |>
  left_join(plot_std, by = c("source", "hhid", "plot_id")) |>
  left_join(fert_plot, by = c("source", "hhid", "plot_id")) |>
  mutate(
    # Sole crops: yield AND N are expressed per planted hectare (plot area x share
    # planted), so a partly planted plot does not dilute the N rate
    crop_area_ha   = plot_area * coalesce(intercrop_fraction, 1) / 1e4,
    yield_farmer   = coalesce(harvest_kg / crop_area_ha, yield_reported),
    fertilizer_amount = if_else(fertilizer_used %in% FALSE & is.na(fert_kg), 0, fert_kg / crop_area_ha),
    N_fertilizer   = case_when(
      fertilizer_used %in% FALSE & (is.na(fert_kg) | fert_kg == 0) ~ 0,
      n_complete %in% TRUE ~ N_kg / crop_area_ha,
      .default = NA_real_
    ),
    fertilizer_used = coalesce(fertilizer_used, fert_kg > 0)
  )

if (!is.null(crop_std)) {
  pc <- pc |> mutate(crop_key = ci_key(crop)) |>
    left_join(crop_std, by = "crop_key", na_matches = "never") |>
    mutate(crop_label = crop, crop = coalesce(terminag_crop, str_to_lower(crop))) |>
    select(-terminag_crop, -crop_key)
}

hh_std <- if (nrow(hh) == 0) tibble(source = character(), hhid = character()) else
  hh |> ensure(c("ea_id", "adm1", "adm2", "adm3", "latitude", "longitude")) |>
  mutate(across(c(latitude, longitude), as.numeric)) |>
  summarise(across(c(ea_id, adm1, adm2, adm3, latitude, longitude), \(x) first(na.omit(x))),
            .by = c(source, hhid)) |>
  # coordinates outside Malawi (e.g. 0, 0) are data errors: discarded and flagged
  mutate(coords_invalid = !is.na(latitude) & !in_malawi(longitude, latitude),
         latitude  = if_else(coords_invalid, NA_real_, latitude),
         longitude = if_else(coords_invalid, NA_real_, longitude)) |>
  mutate(geo_level = case_when(
    !is.na(latitude) & str_starts(source, "LSMS") ~ "EA (modified coordinates)",
    !is.na(latitude)                               ~ "plot GPS",
    .default = NA_character_),
    geo_uncertainty = NA_real_)
if (any(hh_std$coords_invalid %in% TRUE)) {
  message("Coordinates outside Malawi discarded for ", sum(hh_std$coords_invalid), " households: ",
          paste(names(table(hh_std$source[hh_std$coords_invalid])), table(hh_std$source[hh_std$coords_invalid]),
                sep = " = ", collapse = "; "))
}

# District centroids for waves without EA coordinates. Built once from GADM and
# saved to pfpn_config/district_centroids.csv (edit or replace it freely).
# geo_uncertainty = radius (m) of a circle with the district's area.
district_centroids <- function() {
  path <- file.path(cfg_dir, "district_centroids.csv")
  if (!file.exists(path)) {
    if (!requireNamespace("geodata", quietly = TRUE)) {
      message("geodata not installed: no district centroids (install it, or provide ", path, ")")
      return(NULL)
    }
    v <- tryCatch(geodata::gadm(centroid_country, level = 1, path = file.path(out_dir, "raw")),
                  error = \(e) NULL)
    if (is.null(v)) { message("GADM download failed: no district centroids this run"); return(NULL) }
    cen <- terra::crds(terra::centroids(v, inside = TRUE))
    tibble(adm_name = v$NAME_1, longitude = cen[, 1], latitude = cen[, 2],
           geo_uncertainty = sqrt(terra::expanse(v, unit = "m") / pi),
           basis = "GADM level 1, point inside polygon (terra::centroids)") |>
      write_csv(path)
  }
  read_csv(path, show_col_types = FALSE)
}

cents <- if (any(hh_std$source %in% centroid_sources)) district_centroids() else NULL
if (!is.null(cents)) {
  cents <- cents |> mutate(dkey = district_key(adm_name))
  need  <- hh_std$source %in% centroid_sources & is.na(hh_std$latitude)
  matched <- tibble(adm2 = hh_std$adm2[need]) |>
    mutate(dkey = district_key(adm2)) |>
    left_join(select(cents, dkey, adm_name, c_lat = latitude, c_lon = longitude, c_unc = geo_uncertainty),
              by = "dkey", na_matches = "never", relationship = "many-to-one")
  hh_std$latitude[need]        <- matched$c_lat
  hh_std$longitude[need]       <- matched$c_lon
  hh_std$geo_uncertainty[need] <- matched$c_unc
  hh_std$geo_level[need]       <- if_else(is.na(matched$c_lat), NA_character_, "district_centroid")
  distinct(matched, adm2, adm_name) |>
    write_csv(file.path(out_dir, "district_centroid_matches.csv"), na = "")
  unmatched <- matched |> filter(is.na(adm_name)) |> distinct(adm2)
  if (nrow(unmatched) > 0) {
    message("Districts without a centroid match (edit district_centroids.csv): ",
            paste(unmatched$adm2, collapse = ", "))
  }
}

survey_unified <- pc |>
  left_join(hh_std, by = c("source", "hhid")) |>
  left_join(meta |> rename(season_default = season), by = "source") |>
  mutate(season = coalesce(season, season_default)) |>
  pivot_longer(c(yield_farmer, cropcut_yield), names_to = "yield_source", values_to = "yield") |>
  filter(!is.na(yield)) |>
  mutate(
    yield_source   = recode(yield_source, yield_farmer = "farmer_report", cropcut_yield = "crop_cut"),
    crop_cut       = yield_source == "crop_cut",
    yield_moisture = if_else(yield_source == "crop_cut", yield_moisture, NA_real_),
    is_survey = TRUE, on_farm = TRUE, dataset_id = source
  )

# -----------------------------------------------------------------------------
# 6. Carob (already terminag)
# -----------------------------------------------------------------------------

keep_cols <- c("dataset_id", "country", "adm1", "adm2", "adm3", "latitude", "longitude",
               "date", "planting_date", "season", "crop", "yield_part", "intercropped",
               "intercrop_fraction", "plot_area", "irrigated", "fertilizer_used",
               "fertilizer_type", "fertilizer_amount", "N_fertilizer", "crop_cut",
               "yield", "yield_moisture", "is_survey", "on_farm", "trial_id", "plot_id",
               "hhid", "treatment", "geo_uncertainty", "intercrops")

read_carob <- function(group) {
  zf <- tryCatch(fetch(carob_url(group), file.path(out_dir, "raw", basename(carob_url(group))),
                       check = is_zip),
                 error = \(e) { warning("carob ", group, ": download failed"); NULL })
  if (is.null(zf)) return(NULL)
  ff <- unzip_once(zf, file.path(out_dir, "raw", paste0("carob_", group)))
  ff <- ff[str_detect(ff, "\\.csv$") & !str_detect(str_to_lower(basename(ff)), "meta|long")]
  if (length(ff) == 0) { warning("carob ", group, ": no data csv in zip"); return(NULL) }
  map(ff, \(f) read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)) |>
    list_rbind() |>
    select(any_of(keep_cols)) |>
    mutate(program = "carob", source = paste0("carob_", group))
}

carob <- map(carob_groups, read_carob) |> compact() |> list_rbind()

if (nrow(carob) > 0) {
  carob <- carob |>
    ensure(c("geo_uncertainty", "yield", "N_fertilizer", "fertilizer_amount", "plot_area", "intercrop_fraction",
             "yield_moisture", "latitude", "longitude", "crop_cut", "is_survey", "on_farm",
             "intercropped", "fertilizer_used", "irrigated", "plot_id", "hhid", "trial_id",
             "intercrops")) |>
    mutate(
      across(c(yield, N_fertilizer, fertilizer_amount, plot_area, intercrop_fraction,
               yield_moisture, latitude, longitude), \(x) suppressWarnings(as.numeric(x))),
      across(c(crop_cut, is_survey, on_farm, intercropped, fertilizer_used, irrigated), to_logical),
      yield_source = case_when(crop_cut %in% TRUE ~ "crop_cut",
                               is_survey %in% TRUE ~ "farmer_report",
                               .default = "experiment"),
      plot_id = coalesce(plot_id, trial_id),
      geo_uncertainty = suppressWarnings(as.numeric(geo_uncertainty)),
      coords_invalid = !is.na(latitude) & !in_malawi(longitude, latitude),
      latitude  = if_else(coords_invalid, NA_real_, latitude),
      longitude = if_else(coords_invalid, NA_real_, longitude),
      geo_level = if_else(!is.na(latitude), "as reported (carob)", NA_character_),
      crop_stand = case_when(
        intercropped %in% TRUE | (!is.na(intercrops) & !str_to_lower(intercrops) %in% c("", "none", "no")) ~ "intercrop",
        intercropped %in% FALSE ~ "sole",
        .default = "unknown")
    ) |>
    filter(!is.na(yield))
  if (!is.null(carob_countries)) carob <- carob |> filter(country %in% carob_countries)
}

# -----------------------------------------------------------------------------
# 7. Unify, PFP-N, QC against terminag ranges
# -----------------------------------------------------------------------------

final_cols <- c("program", "source", "dataset_id", "country", "adm1", "adm2", "adm3",
                "latitude", "longitude", "geo_level", "geo_uncertainty", "ea_id",
                "date", "season", "hhid", "plot_id", "crop",
                "crop_label", "yield_part", "crop_stand", "intercropped", "intercrop_fraction", "plot_area",
                "area_source", "plot_area_gps_m2", "plot_area_reported_m2", "irrigated", "fertilizer_used", "fertilizer_type",
                "fertilizer_amount", "N_fertilizer", "n_complete", "fert_kg_source",
                "fert_kg_mismatch", "OM_used", "organic_qty", "organic_unit", "harvest_block", "yield_source",
                "crop_cut", "yield", "yield_moisture", "pfp_n", "is_survey", "on_farm",
                "treatment", "quality_flags")

all_rows <- bind_rows(
  survey_unified |> mutate(across(any_of(c("date", "plot_id", "hhid")), as.character)),
  carob
) |>
  ensure(c("crop_stand", "area_source", "intercrop_fraction", "fert_kg_mismatch", "fert_kg_source",
           "n_complete", "harvest_block", "geo_level", "yield_source", "coords_invalid"))

stand_keep <- c("sole", "sole (inferred: only crop on plot)",
                if (unknown_stand == "include_flagged") "unknown")

stand_summary <- all_rows |>
  count(program, source, crop_stand) |>
  mutate(kept = crop_stand %in% stand_keep)
write_csv(stand_summary, file.path(out_dir, "crop_stand_summary.csv"), na = "")

flag <- \(cond, label) if_else(cond %in% TRUE, label, NA_character_)

unified <- all_rows |>
  filter(crop_stand %in% stand_keep) |>
  mutate(
    pfp_n = if_else(N_fertilizer > 0, yield / N_fertilizer, NA_real_),
    quality_flags = pmap_chr(
      list(
        flag(area_source == "farmer_report", "area_farmer_reported"),
        flag(intercrop_fraction < 1, "partly_planted"),
        flag(str_detect(crop_stand, "inferred"), "stand_inferred"),
        flag(crop_stand == "unknown", "stand_unknown"),
        flag(fert_kg_mismatch, "fert_kg_mismatch"),
        flag(str_detect(fert_kg_source, "qty_x_unit"), "fert_kg_from_unit_conversion"),
        flag(fertilizer_used %in% TRUE & is.na(N_fertilizer), "N_unknown_product"),
        flag(N_fertilizer > 0 & N_fertilizer < low_N_threshold, paste0("N_below_", low_N_threshold)),
        flag(harvest_block == "alt", "harvest_second_block"),
        flag(geo_level == "district_centroid", "location_district_centroid, coords_outside_malawi_discarded"),
        flag(coords_invalid, "coords_outside_malawi_discarded")
      ),
      \(...) { f <- na.omit(c(...)); if (length(f) == 0) NA_character_ else paste(f, collapse = "; ") }
    )
  ) |>
  select(any_of(final_cols))

# terminag valid ranges (reuse the lookup script's download if present)
tg_zip <- file.path(lookup_dir, "raw", "terminag-main.zip")
if (!usable(tg_zip, is_zip)) tg_zip <- fetch("https://github.com/controvoc/terminag/archive/refs/heads/main.zip",
                                             file.path(out_dir, "raw", "terminag-main.zip"), check = is_zip)
tg_dir <- file.path(out_dir, "raw", "terminag")
invisible(unzip_once(tg_zip, tg_dir))
tg_all <- list.files(file.path(tg_dir, "terminag-main", "variables"), "\\.csv$", full.names = TRUE) |>
  map(\(f) read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)) |>
  list_rbind() |>
  distinct(name, .keep_all = TRUE)
tg_vars <- tg_all |>
  mutate(valid_min = suppressWarnings(as.numeric(valid_min)),
         valid_max = suppressWarnings(as.numeric(valid_max))) |>
  filter(name %in% names(unified), !is.na(valid_min) | !is.na(valid_max))

# terminag's plot_area range (1-350 m2 in the Sep 2026 version [verified]) is meant
# for trial plots; farmers' fields exceed it, so that check is skipped for surveys.
qc_skip_survey <- c("plot_area")

qc <- map(seq_len(nrow(tg_vars)), \(i) {
  v <- tg_vars$name[i]; x <- unified[[v]]
  out <- (!is.na(tg_vars$valid_min[i]) & x < tg_vars$valid_min[i]) |
         (!is.na(tg_vars$valid_max[i]) & x > tg_vars$valid_max[i])
  if (v %in% qc_skip_survey) out <- out & !(unified$is_survey %in% TRUE)
  if_else(out %in% TRUE, v, NA_character_)
})
unified$qc_out_of_range <- if (length(qc) == 0) NA_character_ else
  pmap_chr(qc, \(...) { f <- na.omit(c(...)); if (length(f) == 0) NA_character_ else paste(f, collapse = "; ") })

# -----------------------------------------------------------------------------
# 8. Data dictionary and labelled outputs
#    .rds : every column carries attributes label, description, unit, derivation
#           (label is the attribute read by RStudio's viewer, haven and labelled);
#           the full dictionary and dataset metadata are attributes of the table.
#    .csv : pfpn_unified.csv (data) + pfpn_unified_dictionary.csv (dictionary).
# -----------------------------------------------------------------------------

# Per-column documentation. `unit` = NA falls back to terminag's unit.
col_doc <- tribble(
  ~variable,          ~label,                                   ~unit,          ~derivation,
  "program",          "Data programme",                         NA,             "LSMS_MWI, LCAS or carob",
  "source",           "Source dataset / survey wave",           NA,             "LSMS: LSMS_<catalog id>_<years>; LCAS: dataset name; Carob: carob_<collection>",
  "dataset_id",       "Dataset identifier",                     NA,             "Carob dataset_id; for surveys equal to source",
  "country",          "Country",                                NA,             "LSMS/LCAS: from configuration; Carob: as published",
  "adm1",             "Administrative level 1",                 NA,             "LSMS: region (hh_mod_a_filt)",
  "adm2",             "Administrative level 2",                 NA,             "LSMS: district (hh_mod_a_filt; IHS3 hh_a01)",
  "adm3",             "Administrative level 3",                 NA,             "Not mapped for surveys; Carob as published",
  "latitude",         "Latitude",                               "degrees",      "LSMS: EA coordinates (modified/displaced by the data producer); IHS6: district centroid; LCAS: plot GPS; Carob: as published. See geo_level",
  "longitude",        "Longitude",                              "degrees",      "As latitude",
  "geo_level",        "Spatial level of the coordinates",       NA,             "EA (modified coordinates), district_centroid, plot GPS, or as reported (carob)",
  "geo_uncertainty",  "Uncertainty of the coordinates",         "m",            "District centroids: radius of a circle with the district's area; otherwise NA or as published",
  "ea_id",            "Enumeration area identifier",            NA,             "LSMS hh_mod_a_filt ea_id; useful for clustering",
  "date",             "Year of the reference season",           NA,             "From configuration (lsms_waves / lcas_datasets); Carob as published",
  "season",           "Season",                                 NA,             "wet (LSMS rainy-season modules) unless mapped otherwise",
  "hhid",             "Household identifier",                   NA,             "LSMS case_id (IHS3-6) or y2_hhid (IHPS 2013); LCAS respondent id (row number if not mapped); Carob as published",
  "plot_id",          "Plot identifier",                        NA,             "LSMS gardenid_plotid (IHS4-6), visit_plotid (IHS3) or plot ID ag_c00/ag_d00/ag_g00 (IHPS 2013); Carob plot_id or trial_id",
  "crop",             "Crop (terminag name)",                   NA,             "Source crop label mapped to terminag with crop_map.csv; unmatched labels kept in lower case",
  "crop_label",       "Crop as labelled in the source",         NA,             "Original crop label before mapping",
  "yield_part",       "Harvested part the yield refers to",     NA,             "From crop_map.csv (terminag values_crop) or Carob",
  "crop_stand",       "Crop stand",                             NA,             "sole; sole (inferred: only crop on plot); unknown. Intercrops are dropped",
  "intercropped",     "Intercropped",                           NA,             "FALSE for kept rows; NA where the stand is unknown",
  "intercrop_fraction","Share of the plot planted with the crop","fraction 0-1", "Priority: % planted (IHS6) > whole plot planted = yes > share category (intercrop_fraction.csv). Used for the planted area",
  "plot_area",        "Plot area",                              "m2",           "GPS area if > 0, otherwise farmer-reported area converted with area_unit_m2.csv; see area_source",
  "area_source",      "Source of the plot area",                NA,             "gps or farmer_report",
  "plot_area_gps_m2", "GPS-measured plot area",                 "m2",           "LSMS GPS area converted to m2 (NA if not measured); kept to assess reported-area bias",
  "plot_area_reported_m2", "Farmer-reported plot area",         "m2",           "Reported area converted with area_unit_m2.csv; kept to assess reported-area bias",
  "irrigated",        "Irrigated",                              NA,             "Not mapped for surveys; Carob as published",
  "fertilizer_used",  "Inorganic fertilizer used",              NA,             "LSMS ag_d38; if missing, TRUE when fertilizer kg > 0",
  "fertilizer_type",  "Inorganic fertilizer product(s)",        NA,             "Products applied on the plot (terminag names where matched), separated by '; '",
  "fertilizer_amount","Inorganic fertilizer applied",           "kg/ha",        "Sum of product kg over applications / planted area (plot_area x intercrop_fraction)",
  "N_fertilizer",     "Mineral N applied",                      "kg/ha",        "Sum over applications of product kg x N% (fertilizer_n.csv) / planted area; 0 if no fertilizer; NA if any product has unknown N%",
  "n_complete",       "N content known for all products",       NA,             "FALSE when a product used on the plot has no N% in fertilizer_n.csv",
  "fert_kg_source",   "Basis of fertilizer kg",                 NA,             "enumerator_kg (enumerator total, preferred) and/or qty_x_unit (quantity x unit weight)",
  "fert_kg_mismatch", "Enumerator kg differs from qty x unit",  NA,             "TRUE when both exist and differ by > 10%",
  "OM_used",          "Organic fertilizer used",                NA,             "LSMS ag_d36 (or organic quantity > 0). Covariate only; not included in N_fertilizer",
  "organic_qty",      "Organic fertilizer quantity",            "see organic_unit", "LSMS ag_d37a; not converted to kg",
  "organic_unit",     "Unit of organic fertilizer quantity",    NA,             "LSMS ag_d37b label",
  "harvest_block",    "Harvest question block used",            NA,             "main, or alt (IHS3 / IHPS 2013 ag_g09 block, used where ag_g13 is empty)",
  "yield_source",     "Yield measurement method",               NA,             "farmer_report, crop_cut or experiment (Carob trials). Do not pool without accounting for method",
  "crop_cut",         "Yield from crop cut",                    NA,             "TRUE when yield_source is crop_cut",
  "yield",            "Yield",                                  "kg/ha",        "Surveys: harvest kg (harvest_unit_kg.csv) / planted area, or reported/crop-cut yield; Carob as published. Moisture basis in yield_moisture where known",
  "yield_moisture",   "Moisture content of the yield",          "%",            "Crop cuts where recorded; unknown for farmer reports",
  "pfp_n",            "Partial factor productivity of N",       "kg yield / kg N", "yield / N_fertilizer, sole crops only; NA when N_fertilizer is 0 or unknown",
  "is_survey",        "Survey data",                            NA,             "TRUE for LSMS and LCAS; Carob as published",
  "on_farm",          "On-farm",                                NA,             "TRUE for LSMS and LCAS; Carob as published",
  "treatment",        "Experimental treatment",                 NA,             "Carob trials only",
  "quality_flags",    "Data-quality flags",                     NA,             "Flags separated by '; ' (see data_quality_statement.md): area_farmer_reported, partly_planted, stand_inferred, stand_unknown, fert_kg_mismatch, fert_kg_from_unit_conversion, N_unknown_product, N_below_<threshold>, harvest_second_block, location_district_centroid, coords_outside_malawi_discarded",
  "qc_out_of_range",  "Values outside terminag valid range",    NA,             "Names of variables outside terminag valid_min/valid_max (plot_area not checked for surveys)"
)

describe_values <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  if (is.numeric(x)) return(paste(signif(min(x), 4), "to", signif(max(x), 4)))
  u <- sort(unique(as.character(x)))
  # " | " separates categories; some values themselves contain "; "
  if (length(u) > 15) paste(length(u), "distinct values") else paste(u, collapse = " | ")
}

dictionary <- tibble(variable = names(unified)) |>
  left_join(col_doc, by = "variable") |>
  left_join(tg_all |> transmute(variable = name, terminag_unit = unit,
                                terminag_description = description),
            by = "variable") |>
  mutate(
    terminag_term = variable %in% tg_all$name,
    label         = coalesce(label, variable),
    unit          = coalesce(unit, na_if(terminag_unit, "")),
    type          = map_chr(variable, \(v) switch(class(unified[[v]])[1],
                                                  numeric = "numeric", integer = "integer",
                                                  logical = "logical", "character")),
    values        = map_chr(variable, \(v) describe_values(unified[[v]])),
    n_non_missing = map_int(variable, \(v) sum(!is.na(unified[[v]]))),
    pct_missing   = round(100 * (1 - n_non_missing / nrow(unified)), 1)
  ) |>
  select(variable, label, type, unit, terminag_term, terminag_description, derivation,
         values, n_non_missing, pct_missing)

undocumented <- dictionary |> filter(is.na(derivation)) |> pull(variable)
if (length(undocumented) > 0) warning("Columns without documentation in col_doc: ",
                                      paste(undocumented, collapse = ", "))

# Attach labels to the columns (rds)
unified_labelled <- unified
for (i in seq_len(nrow(dictionary))) {
  v <- dictionary$variable[i]
  attr(unified_labelled[[v]], "label")       <- dictionary$label[i]
  attr(unified_labelled[[v]], "description") <- coalesce(dictionary$terminag_description[i], dictionary$label[i])
  attr(unified_labelled[[v]], "unit")        <- dictionary$unit[i]
  attr(unified_labelled[[v]], "derivation")  <- dictionary$derivation[i]
}
attr(unified_labelled, "title")      <- "PFP-N: yield and mineral N for sole crops, Malawi LSMS + LCAS + Carob (terminag names)"
attr(unified_labelled, "created")    <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
attr(unified_labelled, "sources")    <- sort(unique(unified$source))
attr(unified_labelled, "dictionary") <- dictionary
attr(unified_labelled, "notes")      <- c(
  "One row per plot-crop x yield_source; sole crops only; mineral N only (organic = covariates).",
  "Conversion tables used: pfpn_config/*.csv at the time of this run.",
  "Data quality: see data_quality_statement.md and the quality_flags column."
)

write_csv(unified, file.path(out_dir, "pfpn_unified.csv"), na = "")
write_csv(dictionary, file.path(out_dir, "pfpn_unified_dictionary.csv"), na = "")
saveRDS(unified_labelled, file.path(out_dir, "pfpn_unified.rds"))
message("Wrote pfpn_unified.csv + pfpn_unified_dictionary.csv, and pfpn_unified.rds (labelled)")

missing_rows <- function(tbl, table, value_col, key_cols) {
  if (is.null(tbl)) return(NULL)
  key_cols <- intersect(key_cols, names(tbl))
  tbl |>
    filter(is.na(.data[[value_col]])) |>
    transmute(table, source = if ("source" %in% names(tbl)) source else NA_character_,
              key = do.call(paste, c(across(all_of(key_cols)), sep = " | ")), fill_column = value_col)
}
missing_conv <- bind_rows(
  missing_rows(harvest_conv,   "harvest_unit_kg.csv",    "kg_per_unit", harvest_keys),
  missing_rows(fert_unit_conv, "fert_unit_kg.csv",       "kg_per_unit", "fertilizer_unit"),
  missing_rows(area_conv,      "area_unit_m2.csv",       "m2_per_unit", "area_unit"),
  missing_rows(frac_conv,      "intercrop_fraction.csv", "fraction",    "intercrop_fraction"),
  missing_rows(fert_n,         "fertilizer_n.csv",       "N_pct",       "value_label")
)
write_csv(missing_conv, file.path(out_dir, "missing_conversions.csv"), na = "")

summary_tbl <- unified |>
  summarise(
    rows          = n(),
    plots         = n_distinct(paste(hhid, plot_id)),
    with_yield    = sum(!is.na(yield)),
    with_N        = sum(!is.na(N_fertilizer)),
    zero_N        = sum(N_fertilizer %in% 0),
    with_pfp_n    = sum(!is.na(pfp_n)),
    median_pfp_n  = median(pfp_n, na.rm = TRUE),
    qc_flagged    = sum(!is.na(qc_out_of_range)),
    .by = c(program, source, yield_source)
  )
write_csv(summary_tbl, file.path(out_dir, "summary_by_source.csv"), na = "")

print(summary_tbl, n = Inf)

# -----------------------------------------------------------------------------
# 9. Data quality statement (pfpn_compiled/data_quality_statement.md)
#    Fixed text on error sources + indicators computed from this run.
# -----------------------------------------------------------------------------

pct <- \(x) if (length(x) == 0 || all(is.na(x))) "-" else sprintf("%.0f%%", 100 * mean(x, na.rm = TRUE))
has <- \(x, pattern) str_detect(coalesce(x, ""), pattern)

dq <- unified |>
  summarise(
    rows              = as.character(n()),
    `area GPS`        = pct(area_source == "gps"),
    `area reported`   = pct(area_source == "farmer_report"),
    `partly planted`  = pct(intercrop_fraction < 1),
    `stand inferred/unknown` = pct(has(crop_stand, "inferred|unknown")),
    `fert kg enumerator` = pct(if_else(fertilizer_used %in% TRUE, has(fert_kg_source, "enumerator"), NA)),
    `fert kg mismatch >10%` = pct(if_else(fertilizer_used %in% TRUE, fert_kg_mismatch %in% TRUE, NA)),
    `N unknown (fertilised)` = pct(if_else(fertilizer_used %in% TRUE, is.na(N_fertilizer), NA)),
    `N < threshold`   = pct(if_else(N_fertilizer > 0, N_fertilizer < low_N_threshold, NA)),
    `2nd harvest block` = pct(harvest_block == "alt"),
    `district centroid` = pct(geo_level == "district_centroid"),
    `PFP-N median [IQR]` = if (all(is.na(pfp_n))) "-" else
      sprintf("%.1f [%.1f-%.1f]", median(pfp_n, na.rm = TRUE),
              quantile(pfp_n, 0.25, na.rm = TRUE), quantile(pfp_n, 0.75, na.rm = TRUE)),
    .by = c(source, yield_source)
  )

md_table <- function(d) {
  d <- mutate(d, across(everything(), as.character))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

dropped <- stand_summary |> filter(!kept) |> summarise(n = sum(n), .by = source)

writeLines(c(
  "# Data quality statement: PFP-N (sole crops)",
  "",
  paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), ". PFP-N = yield (kg/ha) / mineral N applied ",
         "(kg N/ha), both per planted hectare, computed for sole crops only; plots without mineral N ",
         "have no PFP-N. Organic inputs are recorded as covariates and not counted as N."),
  "",
  "## N applied (denominator)",
  "",
  "- **Quantities are farmer recall** of product applied, often in local units (bags, pails, plates).",
  "  The enumerator's total in kg is used where given; otherwise quantity x unit, with unit",
  "  weights read from labels or entered in `fert_unit_kg.csv`. Where both exist and differ",
  "  by more than 10 %, the row is flagged `fert_kg_mismatch`.",
  "- **N content is nominal**: grade on the label (e.g. 23:21:0+4S) or terminag's typical value",
  "  for the product. Actual products are not analysed, and adulterated, mixed or mislabelled",
  "  fertilizer is not detected. Generic answers ('NPK', 'Other') have no N content, so",
  "  those plots have unknown N (`N_unknown_product`) rather than an understated N.",
  "- **Plot attribution**: N is recorded per plot and application; farmers who spread one",
  "  purchase over several plots may misallocate it between plots.",
  "- **Mineral N only**: where manure or compost was applied (`OM_used`), PFP-N credits all",
  "  yield to mineral N and so overstates its productivity.",
  "",
  "## Yield (numerator)",
  "",
  "- **Farmer-reported harvest** in non-standard units is converted with `harvest_unit_kg.csv`.",
  "  Label-derived weights (e.g. '50 KG BAG' = 50 kg) ignore crop condition (shelled vs",
  "  unshelled, fresh vs dry); the official IHS conversion factors should replace them.",
  "  Moisture content of farmer-reported harvests is unknown.",
  "- **Area**: GPS area where available, otherwise farmer-reported area (`area_farmer_reported`;",
  "  all of IHS3). Survey methods work generally finds self-reported areas and harvests",
  "  differ systematically from GPS-measured areas and crop cuts [recall: LSMS",
  "  methodological studies; not verified here], so yields from reported areas are less reliable.",
  "- **Partly planted plots** (`partly_planted`) depend on a reported share of the plot",
  "  (categorical in most waves), which adds error to both yield and N per hectare.",
  "- **Crop cuts** (LCAS) measure a small sampled area and scale up; they avoid recall and",
  "  unit problems but are sensitive to where the quadrats fall and to moisture adjustment.",
  "",
  "## PFP-N (ratio)",
  "",
  "- Errors in yield and N combine in the ratio, and small N rates make it unstable:",
  paste0("  PFP-N from N below ", low_N_threshold, " kg N/ha is flagged (`N_below_", low_N_threshold,
         "`). Report medians/IQR or trim before comparing sources."),
  "- Farmer-report, crop-cut and experimental (Carob) yields are kept as separate rows",
  "  (`yield_source`) and should not be pooled without accounting for their different biases.",
  "- Sole crops only: intercropped plots were dropped; LSMS plots without a stand code were",
  "  classed from the number of crops on the plot (`stand_inferred`). Sources without any stand",
  paste0("  information were ", if (unknown_stand == "include_flagged") "kept and flagged `stand_unknown`." else "dropped."),
  "- Location: LSMS EA coordinates are deliberately displaced by the data producer; IHS6 rows",
  "  use district centroids (`location_district_centroid`) and suit only coarse spatial analysis.",
  "",
  "## Indicators from this run",
  "",
  md_table(dq),
  "",
  "Rows dropped as intercrop (or unknown stand, if excluded):",
  "",
  if (nrow(dropped) == 0) "none" else md_table(dropped)
), file.path(out_dir, "data_quality_statement.md"))
message("Data quality statement: ", file.path(out_dir, "data_quality_statement.md"))
if (nrow(missing_conv) > 0) {
  message("\n", nrow(missing_conv), " conversion rows still empty -> ",
          file.path(out_dir, "missing_conversions.csv"), ". Fill them in ", cfg_dir, " and re-run.")
}

# ==============================================================================
# END
# ==============================================================================
