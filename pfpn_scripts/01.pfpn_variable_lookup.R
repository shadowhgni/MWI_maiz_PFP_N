# =============================================================================
# PFP-N variable specification: LSMS Malawi (IHS3-IHS6) + LCAS -> terminag
# -----------------------------------------------------------------------------
# Replaces the keyword-tagging version. Variables are now listed explicitly per
# wave (from the World Bank data dictionaries) and checked against each wave's
# DDI codebook before any data are read:
#   - does the variable exist in that file?
#   - does its label match what the role expects (e.g. "QUANTITY")?
#   - how many valid observations does the DDI report?
# Outputs (pfpn_lookup/):
#   spec_check.csv             every spec row with its check result
#   data_dictionary_<src>.csv  full variable list of the files used, per wave
#   value_labels.csv           codes of the categorical variables used
#   fertilizer_n_content_stub.csv, crop_map_stub.csv   seeds for the compile step
#   pfpn_config/source_map.csv the mapping read by pfpn_compile.R (written only
#                              if absent; otherwise source_map_generated.csv)
#
# Scope decisions: plot-crop level; mineral N only (organic = covariate);
# IHS2 (2004-05) dropped (no plot-level harvest); EA coordinates where the
# wave provides them, district centroids for IHS6 (handled in pfpn_compile.R).
#
# Status column in the spec:
#   verified  name and label seen in the online data dictionary
#             (microdata.worldbank.org/catalog/<id>/data-dictionary) or the DDI
#   expected  same name as a verified wave; confirmed or rejected by the check
#
# Online dictionaries checked: IHS6 ag_mod_d, ag_mod_g, hh_mod_a_filt; IHS5 ag_mod_d;
# IHS4 AG_MOD_D, AG_MOD_G; IHS3 AG_MOD_C, AG_MOD_G; IHPS 2013 AG_MOD_C_13,
# AG_MOD_D_13, AG_MOD_G_13. Other rows: names seen in the DDI-based lookup of the
# earlier script version, or expected (checked against the DDI when this runs).
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (xml2, readxl, terra)

source("00.pfpn_utils.R", local = TRUE)   # paths, update-aware downloads, file checks

out_dir <- pfpn_path("pfpn_lookup")
cfg_dir <- pfpn_path("pfpn_config")
raw_dir <- file.path(out_dir, "raw")
walk(c(raw_dir, cfg_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

normalise <- function(x) x |> str_to_lower() |> str_replace_all("[-_/]", " ") |> str_squish()
file_stem <- function(x) str_to_lower(tools::file_path_sans_ext(basename(x)))

ACRE_M2 <- 4046.8564224

# -----------------------------------------------------------------------------
# 1. LSMS specification, one block per wave
#    check: regex (case-insensitive) the DDI label must match
# -----------------------------------------------------------------------------

sr <- function(file, var, role, check = NA, app = NA, multiply = NA,
               status = "verified", note = NA) {
  tibble(file, var, role, check, app = as.integer(app), multiply = as.numeric(multiply),
         status, note)
}

keys_gp <- function(f) bind_rows(sr(f, "case_id", "hhid", "household|\\bhh\\b"),
                                 sr(f, "gardenid", "plot_id", "garden"),
                                 sr(f, "plotid", "plot_id", "plot"))

# IHS4, IHS5, IHS6 share names. `st` gives the status of blocks not seen online.
spec_ihs456 <- function(C, D, G, A, geo = NULL, st = list()) {
  s <- \(x) st[[x]] %||% "verified"
  bind_rows(
    keys_gp(C), keys_gp(D), keys_gp(G),
    sr(C, "ag_c04a", "plot_area_reported", "self reported area"),
    sr(C, "ag_c04b", "area_unit", "unit", status = s("c04b")),
    sr(C, "ag_c04c", "plot_area_gps", "gps", multiply = ACRE_M2,
       note = "assumed acres: label states no unit"),
    sr(D, "ag_d38",  "fertilizer_used",   "inorganic fertili"),
    sr(D, "ag_d39a", "fertilizer_type",   "type of inorganic", 1),
    sr(D, "ag_d39b", "fertilizer_amount", "quantity", 1, status = s("fert")),
    sr(D, "ag_d39c", "fertilizer_unit",   "unit",     1, status = s("fert")),
    sr(D, "ag_d39d", "fertilizer_kg",     "kg",       1, status = s("fert"),
       note = "ENUMERATOR: TOTAL KGs"),
    sr(D, "ag_d39g", "fertilizer_type",   "type of inorganic", 2),
    sr(D, "ag_d39h", "fertilizer_amount", "quantity", 2, status = s("fert")),
    sr(D, "ag_d39i", "fertilizer_unit",   "unit",     2, status = s("fert")),
    sr(D, "ag_d39j", "fertilizer_kg",     "kg",       2, status = s("fert"),
       note = "ENUMERATOR: TOTAL KGs"),
    sr(D, "ag_d36",  "OM_used",      "organic fertili", status = s("fert")),
    sr(D, "ag_d37a", "organic_qty",  "organic fertili"),
    sr(D, "ag_d37b", "organic_unit", "unit", status = s("fert")),
    sr(G, "crop_code", "crop", "crop"),
    sr(G, "ag_g01",  "intercropped",       "crop stand"),
    sr(G, "ag_g02",  "whole_plot_planted", "entire area", status = s("share")),
    sr(G, "ag_g03",  "intercrop_fraction", "how much of the", status = s("share")),
    sr(G, "ag_g13a", "harvest_qty",       "harvest.*quant"),
    sr(G, "ag_g13b", "harvest_unit",      "harvest.*unit"),
    sr(G, "ag_g13c", "harvest_condition", "s/u"),
    sr(A, "case_id", "hhid", "household|\\bhh\\b"),
    sr(A, "ea_id",   "ea_id", "ea"),
    sr(A, "region",  "adm1", "region"),
    sr(A, "district", "adm2", "district"),
    geo
  )
}

geo_rows <- function(f, lat, lon) bind_rows(sr(f, "case_id", "hhid", "household|\\bhh\\b"),
                                            sr(f, lat, "latitude",  "latitude"),
                                            sr(f, lon, "longitude", "longitude"))

lsms_spec <- bind_rows(
  # IHS6 2024-25 [ag_mod_d, ag_mod_g, hh_mod_a_filt verified online; no EA coordinates released]
  spec_ihs456("ag_mod_c", "ag_mod_d", "ag_mod_g", "hh_mod_a_filt",
              st = list(c04b = "expected")) |>
    bind_rows(sr("ag_mod_g", "ag_g11_2", "intercrop_percent", "percentage")) |>
    mutate(source = "LSMS_8507_2024-2025", catalog_id = "8507"),
  # IHS5 2019-20 [ag_mod_d verified online; share block expected (same in IHS4 and IHS6)]
  spec_ihs456("ag_mod_c", "ag_mod_d", "ag_mod_g", "hh_mod_a_filt",
              geo = geo_rows("householdgeovariables_ihs5", "ea_lat_mod", "ea_lon_mod"),
              st = list(c04b = "expected", share = "expected")) |>
    mutate(source = "LSMS_3818_2019-2020", catalog_id = "3818"),
  # IHS4 2016-17 [AG_MOD_D, AG_MOD_G verified online]
  spec_ihs456("AG_MOD_C", "AG_MOD_D", "AG_MOD_G", "HH_MOD_A_FILT",
              geo = geo_rows("HouseholdGeovariablesIHS4", "lat_modified", "lon_modified")) |>
    mutate(source = "LSMS_2936_2016-2017", catalog_id = "2936"),
  # IHS3 2010-11 [AG_MOD_C, AG_MOD_G verified online; AG_MOD_D from DDI lookup]
  # - no GPS plot area released; ag_c04c = the reported area expressed in acres
  # - panel households carry `visit`; it is part of the plot key
  # - harvest asked in two blocks (ag_g09 after the panel/cross-section filter
  #   ag_g06-g08, ag_g13 later); mapped as main + alternative and combined
  bind_rows(
    sr("AG_MOD_C", "case_id", "hhid", "household|\\bhh\\b"), sr("AG_MOD_C", "visit", "plot_id", "visit"),
    sr("AG_MOD_C", "ag_c00", "plot_id", "plot"),
    sr("AG_MOD_D", "case_id", "hhid", "household|\\bhh\\b"), sr("AG_MOD_D", "visit", "plot_id", "visit",
       status = "expected"),
    sr("AG_MOD_D", "ag_d00", "plot_id", "plot"),
    sr("AG_MOD_G", "case_id", "hhid", "household|\\bhh\\b"), sr("AG_MOD_G", "visit", "plot_id", "visit"),
    sr("AG_MOD_G", "ag_g0b", "plot_id", "plot"),
    sr("AG_MOD_C", "ag_c04c", "plot_area_reported_m2", "area in acres", multiply = ACRE_M2,
       note = "reported area in acres; no GPS area in IHS3"),
    sr("AG_MOD_D", "ag_d38",  "fertilizer_used",   "inorganic fertili"),
    sr("AG_MOD_D", "ag_d39a", "fertilizer_type",   "type", 1),
    sr("AG_MOD_D", "ag_d39b", "fertilizer_amount", "quant", 1),
    sr("AG_MOD_D", "ag_d39c", "fertilizer_unit",   "unit", 1),
    sr("AG_MOD_D", "ag_d39d", "fertilizer_kg",     "kgs", 1),
    sr("AG_MOD_D", "ag_d39f", "fertilizer_type",   "type", 2),
    sr("AG_MOD_D", "ag_d39g", "fertilizer_amount", "quant", 2),
    sr("AG_MOD_D", "ag_d39h", "fertilizer_unit",   "unit", 2),
    sr("AG_MOD_D", "ag_d39i", "fertilizer_kg",     "kgs", 2),
    sr("AG_MOD_D", "ag_d36",  "OM_used", "organic fertili", status = "expected"),
    sr("AG_MOD_D", "ag_d37a", "organic_qty", "organic fertili"),
    sr("AG_MOD_D", "ag_d37b", "organic_unit", "unit"),
    sr("AG_MOD_G", "ag_g0d",  "crop", "crop code"),
    sr("AG_MOD_G", "ag_g01",  "intercropped", "crop stand"),
    sr("AG_MOD_G", "ag_g02",  "whole_plot_planted", "entire area"),
    sr("AG_MOD_G", "ag_g03",  "intercrop_fraction", "how much of the"),
    sr("AG_MOD_G", "ag_g13a", "harvest_qty", "harvest.*quant"),
    sr("AG_MOD_G", "ag_g13b", "harvest_unit", "harvest.*unit"),
    sr("AG_MOD_G", "ag_g13c", "harvest_condition", "s/u"),
    sr("AG_MOD_G", "ag_g09a", "harvest_qty_alt", "harvest.*quant"),
    sr("AG_MOD_G", "ag_g09b", "harvest_unit_alt", "harvest.*unit"),
    sr("AG_MOD_G", "ag_g09c", "harvest_condition_alt", "s/u"),
    sr("HH_MOD_A_FILT", "case_id", "hhid", "household|\\bhh\\b"),
    sr("HH_MOD_A_FILT", "ea_id", "ea_id", "ea"),
    sr("HH_MOD_A_FILT", "hh_a01", "adm2", "district"),
    geo_rows("HouseholdGeovariables", "lat_modified", "lon_modified")
  ) |> mutate(source = "LSMS_1003_2010-2011", catalog_id = "1003"),
  # IHPS 2013 (catalog 2248): the 2013 round of the panel of IHS3 households
  # (204 of IHS3's 768 EAs, re-interviewed Apr-Dec 2013, rainy season 2012/13).
  # Only the *_13 files are used: the *_10 files repeat IHS3 households already
  # in catalog 1003. [AG_MOD_C_13, AG_MOD_D_13, AG_MOD_G_13 verified online]
  # - household key y2_hhid; plot IDs ag_c00 / ag_d00 / ag_g00
  # - module G codes differ from IHS3: ag_g00 = plot, ag_g0b = crop code
  # - fertilizer block follows IHS3: a-e 1st application, f-j 2nd (ag_d39f = type,
  #   verified from its codes); ag_d39d / ag_d39i = total kg, inferred from their
  #   values (0-1500, mostly 25/50/100). Labels are truncated in the DDI, so only
  #   the common prefix is checked.
  # - ag_c04c = GPS area in acres (verified label)
  bind_rows(
    sr("AG_MOD_C_13", "y2_hhid", "hhid", "household|\\bhh\\b"), sr("AG_MOD_C_13", "ag_c00", "plot_id", "plot"),
    sr("AG_MOD_D_13", "y2_hhid", "hhid", "household|\\bhh\\b"), sr("AG_MOD_D_13", "ag_d00", "plot_id", "plot"),
    sr("AG_MOD_G_13", "y2_hhid", "hhid", "household|\\bhh\\b"), sr("AG_MOD_G_13", "ag_g00", "plot_id", "plot"),
    sr("AG_MOD_C_13", "ag_c04a", "plot_area_reported", "area"),
    sr("AG_MOD_C_13", "ag_c04b", "area_unit", "unit"),
    sr("AG_MOD_C_13", "ag_c04c", "plot_area_gps", "gps", multiply = ACRE_M2),
    sr("AG_MOD_D_13", "ag_d36",  "OM_used",      "organic fertili"),
    sr("AG_MOD_D_13", "ag_d37a", "organic_qty",  "organic fertili"),
    sr("AG_MOD_D_13", "ag_d37b", "organic_unit", "organic fertili"),
    sr("AG_MOD_D_13", "ag_d38",  "fertilizer_used",   "inorganic fertili"),
    sr("AG_MOD_D_13", "ag_d39a", "fertilizer_type",   "inorganic fertili", 1, status = "expected"),
    sr("AG_MOD_D_13", "ag_d39b", "fertilizer_amount", "inorganic fertili", 1, status = "expected"),
    sr("AG_MOD_D_13", "ag_d39c", "fertilizer_unit",   "inorganic fertili", 1, status = "expected"),
    sr("AG_MOD_D_13", "ag_d39d", "fertilizer_kg",     "inorganic fertili", 1,
       note = "kg: inferred from values (0-1500, mostly 25/50/100)"),
    sr("AG_MOD_D_13", "ag_d39f", "fertilizer_type",   "inorganic fertili", 2),
    sr("AG_MOD_D_13", "ag_d39g", "fertilizer_amount", "inorganic fertili", 2, status = "expected"),
    sr("AG_MOD_D_13", "ag_d39h", "fertilizer_unit",   "inorganic fertili", 2, status = "expected"),
    sr("AG_MOD_D_13", "ag_d39i", "fertilizer_kg",     "inorganic fertili", 2, status = "expected"),
    sr("AG_MOD_G_13", "ag_g0b",  "crop", "crop code"),
    sr("AG_MOD_G_13", "ag_g01",  "intercropped", "crop stand"),
    sr("AG_MOD_G_13", "ag_g02",  "whole_plot_planted", "entire area"),
    sr("AG_MOD_G_13", "ag_g03",  "intercrop_fraction", "how much of the"),
    sr("AG_MOD_G_13", "ag_g13a", "harvest_qty", "harvest"),
    sr("AG_MOD_G_13", "ag_g13b", "harvest_unit", "harvest"),
    sr("AG_MOD_G_13", "ag_g13c", "harvest_condition", "harvest"),
    sr("AG_MOD_G_13", "ag_g09a", "harvest_qty_alt", "harvest"),
    sr("AG_MOD_G_13", "ag_g09b", "harvest_unit_alt", "harvest"),
    sr("AG_MOD_G_13", "ag_g09c", "harvest_condition_alt", "harvest"),
    sr("HH_MOD_A_FILT_13", "y2_hhid", "hhid", "household|\\bhh\\b"),
    sr("HH_MOD_A_FILT_13", "ea_id", "ea_id", "ea", status = "expected"),
    sr("HH_MOD_A_FILT_13", "region", "adm1", "region", status = "expected"),
    sr("HH_MOD_A_FILT_13", "district", "adm2", "district", status = "expected"),
    sr("HouseholdGeovariables_IHPS_13", "y2_hhid", "hhid", "household|\\bhh\\b"),
    sr("HouseholdGeovariables_IHPS_13", "LAT_DD_MOD", "latitude", "latitude", status = "expected",
       note = "name from recall; replace with the latitude variable if the check reports it missing"),
    sr("HouseholdGeovariables_IHPS_13", "LON_DD_MOD", "longitude", "longitude", status = "expected",
       note = "name from recall; replace with the longitude variable if the check reports it missing")
  ) |> mutate(source = "LSMS_2248_2013", catalog_id = "2248")
)

# -----------------------------------------------------------------------------
# 2. Read DDI codebooks (cached)
# -----------------------------------------------------------------------------

ddi_url <- function(id) sprintf("https://microdata.worldbank.org/metadata/export/%s/ddi", id)
node_text <- function(nodes, xpath) xml2::xml_find_first(nodes, xpath) |> xml2::xml_text() |> str_squish()

read_ddi <- function(id) {
  path <- file.path(raw_dir, sprintf("ddi_%s.xml", id))
  path <- tryCatch(fetch(ddi_url(id), path, check = is_ddi), error = \(e) {
    stop("No usable DDI codebook for catalog ", id, " (", conditionMessage(e), "). Download the ",
         "DDI/XML export from https://microdata.worldbank.org/catalog/", id,
         " and save it as ", path, call. = FALSE)
  })
  doc <- path |>
    xml2::read_xml() |> xml2::xml_ns_strip()
  files <- xml2::xml_find_all(doc, "//fileDscr")
  vars  <- xml2::xml_find_all(doc, "//dataDscr/var")
  vars_tbl <- tibble(
    file_id  = xml2::xml_attr(vars, "files") |> word(1),
    order    = seq_along(vars),
    var      = xml2::xml_attr(vars, "name"),
    label    = node_text(vars, "./labl"),
    question = node_text(vars, "./qstn/qstnLit"),
    valid_n  = suppressWarnings(as.numeric(node_text(vars, "./sumStat[@type='vald']"))),
    node     = lapply(seq_along(vars), \(i) vars[[i]])
  )
  tibble(file_id = xml2::xml_attr(files, "ID"), file = node_text(files, "./fileTxt/fileName")) |>
    right_join(vars_tbl, by = "file_id") |>
    mutate(catalog_id = id, file_key = file_stem(file))
}

extract_cats <- function(node) {
  cg <- xml2::xml_find_all(node, "./catgry")
  tibble(value = node_text(cg, "./catValu"), value_label = node_text(cg, "./labl"))
}

ddi <- unique(lsms_spec$catalog_id) |> map(read_ddi) |> list_rbind()

# -----------------------------------------------------------------------------
# 3. Check the LSMS spec against the DDI
# -----------------------------------------------------------------------------

check_spec <- function(spec, dict) {
  spec |>
    mutate(file_key = file_stem(file), var_key = str_to_lower(var)) |>
    left_join(dict |> transmute(catalog_id, file_key, var_key = str_to_lower(var),
                                ddi_var = var, ddi_label = label, valid_n),
              by = c("catalog_id", "file_key", "var_key")) |>
    mutate(result = case_when(
      is.na(ddi_var)                                           ~ "missing",
      is.na(check)                                             ~ "ok",
      str_detect(str_to_lower(coalesce(ddi_label, "")), coalesce(check, ".")) ~ "ok",
      .default = "label_mismatch"
    )) |>
    select(-file_key, -var_key)
}

lsms_check <- check_spec(lsms_spec, ddi)

# Full dictionaries of the files used, for browsing
ddi |>
  semi_join(lsms_spec |> transmute(catalog_id, file_key = file_stem(file)) |> distinct(),
            by = c("catalog_id", "file_key")) |>
  left_join(distinct(lsms_spec, catalog_id, source), by = "catalog_id") |>
  select(source, file, order, var, label, question, valid_n) |>
  group_by(source) |>
  group_walk(\(d, k) write_csv(d, file.path(out_dir, paste0("data_dictionary_", k$source, ".csv")), na = ""))

# Browse helper: show_vars("3818", "ag_mod_g", "^ag_g0")
show_vars <- function(id, file_regex, name_regex = ".") {
  ddi |>
    filter(catalog_id == id, str_detect(file_key, str_to_lower(file_regex)), str_detect(var, name_regex)) |>
    select(file, var, label, valid_n) |>
    print(n = Inf, width = 200)
}

# -----------------------------------------------------------------------------
# 4. LCAS specification, checked against the module xlsforms
# -----------------------------------------------------------------------------

lcas_base <- "https://systems-agronomy.github.io/lcas/module_documentation/"

lcas_spec <- bind_rows(
  sr("05_respondent", "TODO", "hhid", status = "todo",
     note = "respondent/instance id of your LCAS export (e.g. ODK KEY); not in the xlsform"),
  sr("05_respondent", "crop_name", "crop", "crop"),
  sr("06_land", "surveyed_plot", "plot_area_reported", "largest", note = "in local land unit llu"),
  sr("06_land", "llu", "area_unit", "land unit"),
  sr("08_fertility", "apply_minfert", "fertilizer_used", "mineral fertili"),
  sr("08_fertility", "fym_applied_qty", "organic_qty", "fym", note = "quintal per application"),
  sr("09_detailed_fertility", "amt_dap_basal",    "fertilizer_amount", "dap",    1, 1),
  sr("09_detailed_fertility", "amt_npk_basal",    "fertilizer_amount", "npk",    2, 1),
  sr("09_detailed_fertility", "amt_npks_basal",   "fertilizer_amount", "npks",   3, 1),
  sr("09_detailed_fertility", "amt_urea_basal",   "fertilizer_amount", "urea",   4, 1),
  sr("09_detailed_fertility", "amt_mop_basal",    "fertilizer_amount", "mop",    5, 1),
  sr("09_detailed_fertility", "amt_ssp_basal",    "fertilizer_amount", "ssp",    6, 1),
  sr("09_detailed_fertility", "amt_tsp_basal",    "fertilizer_amount", "tsp",    7, 1),
  sr("09_detailed_fertility", "amt_znso4_basal",  "fertilizer_amount", "zinc",   8, 1),
  sr("09_detailed_fertility", "amt_gypsum_basal", "fertilizer_amount", "gypsum", 9, 1),
  sr("09_detailed_fertility", "amt_boron_basal",  "fertilizer_amount", "boron", 10, 1),
  sr("13_harvest", "total_production_lp", "harvest_qty_kg", "quintal", multiply = 100),
  sr("04_cropcut", "cropcut_done", "crop_cut", "crop cut"),
  sr("04_cropcut", "grainYield_tonPerHa", "cropcut_yield", NA, multiply = 1000,
     note = "t/ha -> kg/ha; calculate field without label"),
  sr("18_geolocation", "latitude",  "latitude",  NA, status = "expected", note = "plot GPS, not EA"),
  sr("18_geolocation", "longitude", "longitude", NA, note = "plot GPS, not EA")
) |>
  mutate(source = "LCAS", catalog_id = NA_character_,
         fert_type_fixed = c(DAP = "DAP", NPK = "NPK", NPKS = "NPKS", UREA = "Urea", MOP = "MoP",
                             SSP = "SSP", TSP = "TSP", ZNSO4 = "ZnSO4", GYPSUM = "Gypsum",
                             BORON = "Boron")[str_to_upper(str_match(var, "^amt_(.+)_basal$")[, 2])])

lcas_xlsx_links <- function(module) {
  page <- paste0(lcas_base, module, ".html")
  xml2::read_html(page) |> xml2::xml_find_all("//a[contains(@href, '.xlsx')]") |>
    xml2::xml_attr("href") |> xml2::url_absolute(page) |> unique()
}

pick_label <- function(df) {
  col <- names(df)[str_detect(str_to_lower(names(df)), "^label")][1]
  if (is.na(col)) rep(NA_character_, nrow(df)) else df[[col]]
}

# Cached as raw/lcas_<module>__<file>.xlsx; used when the site is unreachable
read_xlsform <- function(module) {
  urls <- tryCatch(lcas_xlsx_links(module), error = \(e) character())
  paths <- if (length(urls) > 0) {
    map_chr(urls, \(u) fetch(u, file.path(raw_dir, paste0("lcas_", module, "__", basename(u))),
                             check = is_xlsx))
  } else {
    list.files(raw_dir, paste0("^lcas_", module, "__.*\\.xlsx$"), full.names = TRUE)
  }
  map(paths, \(p) {
    s <- readxl::read_excel(p, sheet = "survey", col_types = "text")
    tibble(module, form = basename(p), var = s$name, label = str_squish(pick_label(s)))
  }) |> list_rbind()
}

lcas_fields <- unique(lcas_spec$file) |> map(read_xlsform) |> list_rbind()
if (nrow(lcas_fields) == 0) {
  # no form could be read (site unreachable and nothing cached): keep the columns
  message("LCAS xlsforms not available: LCAS rows are not checked this run")
  lcas_fields <- tibble(module = character(), form = character(), var = character(), label = character())
}

lcas_check <- lcas_spec |>
  left_join(lcas_fields |> filter(!is.na(var)) |> distinct(module, var, .keep_all = TRUE) |>
              transmute(file = module, var, ddi_var = var, ddi_label = label),
            by = c("file", "var")) |>
  mutate(valid_n = NA_real_,
         result = case_when(
           status == "todo"    ~ "todo",
           nrow(lcas_fields) == 0 ~ "not_checked (xlsforms unavailable)",
           is.na(ddi_var)      ~ "missing",
           is.na(check)        ~ "ok",
           str_detect(str_to_lower(coalesce(ddi_label, "")), coalesce(check, ".")) ~ "ok",
           .default = "label_mismatch"))

# -----------------------------------------------------------------------------
# 5. Combined check report and source_map
# -----------------------------------------------------------------------------

spec_check <- bind_rows(lsms_check, lcas_check) |>
  select(source, file, var, role, app, multiply, any_of("fert_type_fixed"), status, result,
         ddi_label, valid_n, check, note)

write_csv(spec_check, file.path(out_dir, "spec_check.csv"), na = "")

message("\nSpec check (rows per source x result):")
print(count(spec_check, source, result) |> pivot_wider(names_from = result, values_from = n, values_fill = 0))
problems <- spec_check |> filter(!result %in% c("ok", "todo"))
if (nrow(problems) > 0) {
  message("\nRows needing attention:")
  print(select(problems, source, file, var, role, result, ddi_label), n = Inf, width = 200)
}

source_map <- spec_check |>
  mutate(var = if_else(result %in% c("missing", "todo"), "TODO", var),
         note = case_when(result == "missing" ~ paste("not in DDI/xlsform:", coalesce(note, "")),
                          result == "label_mismatch" ~ paste("CHECK label:", ddi_label),
                          .default = note)) |>
  transmute(source, file, var, role, app, multiply, fert_type_fixed, note)

map_path <- file.path(cfg_dir, "source_map.csv")
if (file.exists(map_path)) map_path <- file.path(out_dir, "source_map_generated.csv")
write_csv(source_map, map_path, na = "")
message("\nsource map written to ", map_path)

# -----------------------------------------------------------------------------
# 6. terminag + seeds for fertilizer N and crop names
# -----------------------------------------------------------------------------

load_terminag <- function() {
  zip   <- fetch("https://github.com/controvoc/terminag/archive/refs/heads/main.zip",
                 file.path(raw_dir, "terminag-main.zip"), check = is_zip)
  exdir <- file.path(raw_dir, "terminag")
  unzip_once(zip, exdir)
  root  <- file.path(exdir, "terminag-main")
  read_dir <- function(dir) {
    files <- list.files(file.path(root, dir), "\\.csv$", full.names = TRUE)
    names(files) <- basename(files) |> str_remove_all(paste0("^", dir, "_|\\.csv$"))
    map(files, \(f) read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE))
  }
  list(variables = read_dir("variables") |> list_rbind() |> distinct(name, .keep_all = TRUE),
       values = read_dir("values"))
}
tg <- load_terminag()

# Value labels of the categorical LSMS variables used
cat_roles <- c("crop", "fertilizer_type", "fertilizer_unit", "harvest_unit", "harvest_condition",
               "harvest_unit_alt", "harvest_condition_alt",
               "area_unit", "intercropped", "intercrop_fraction", "organic_unit", "adm2")
value_labels <- lsms_check |>
  filter(result != "missing", role %in% cat_roles) |>
  mutate(file_key = file_stem(file), var_key = str_to_lower(var)) |>
  inner_join(ddi |> transmute(catalog_id, file_key, var_key = str_to_lower(var), node),
             by = c("catalog_id", "file_key", "var_key")) |>
  mutate(cats = map(node, extract_cats)) |>
  select(source, file, var, role, cats) |>
  unnest(cats)
write_csv(value_labels, file.path(out_dir, "value_labels.csv"), na = "")

match_vocab <- function(labels, keys) {
  keys <- keys |> filter(!is.na(key), key != "") |> mutate(key_ns = str_remove_all(key, " "))
  map(labels, \(lab) {
    hit <- keys |> filter(str_detect(lab, paste0("\\b", str_escape(key), "\\b"))) |>
      slice_max(nchar(key), n = 1, with_ties = FALSE)
    if (nrow(hit) == 0) hit <- keys |>
      filter(nchar(key_ns) >= 5, str_detect(str_remove_all(lab, " "), fixed(key_ns))) |>
      slice_max(nchar(key_ns), n = 1, with_ties = FALSE)
    select(hit, -key_ns)
  })
}

fert_tbl <- tg$values$fertilizer_type |>
  mutate(across(any_of(c("N", "P", "K")), \(x) suppressWarnings(as.numeric(x))))
fert_keys <- bind_rows(transmute(fert_tbl, terminag_fertilizer = name, key = normalise(name)),
                       transmute(fert_tbl, terminag_fertilizer = name, key = normalise(description))) |>
  distinct()
grade_rx <- "(\\d{1,2}(?:\\.\\d+)?)\\s*[:\\-]\\s*(\\d{1,2}(?:\\.\\d+)?)\\s*[:\\-]\\s*(\\d{1,2}(?:\\.\\d+)?)"

fert_labels <- bind_rows(
  value_labels |> filter(role == "fertilizer_type") |> distinct(value_label),
  tibble(value_label = na.omit(unique(lcas_spec$fert_type_fixed)))
) |> filter(!is.na(value_label), value_label != "") |> distinct()

fert_n_stub <- fert_labels |>
  mutate(lab = normalise(value_label),
         grade_N = as.numeric(str_match(str_to_lower(value_label), grade_rx)[, 2]),
         hit = match_vocab(lab, fert_keys)) |>
  unnest(hit, keep_empty = TRUE) |>
  left_join(select(fert_tbl, terminag_fertilizer = name, tg_N = N),
            by = "terminag_fertilizer", na_matches = "never") |>
  mutate(N_pct = coalesce(grade_N, tg_N),
         N_pct_basis = case_when(!is.na(grade_N) ~ "parsed from N:P:K grade",
                                 !is.na(tg_N) ~ "terminag values_fertilizer_type",
                                 .default = "FILL IN")) |>
  select(value_label, terminag_fertilizer, N_pct, N_pct_basis)
# Corrections reviewed by hand for Malawi labels (kept in code so that a fresh
# run from empty folders reproduces them)
fert_n_stub <- fert_n_stub |>
  mutate(lab = str_to_lower(value_label),
         fix = case_when(
           str_detect(lab, "compound d|d[- ]?compound") ~ "D-compound",
           str_detect(lab, "^mop$|muriate")           ~ "MOP",
           str_detect(lab, "boron")                   ~ "boron",
           str_detect(lab, "^none$")                  ~ "none",
           .default = NA_character_),
         N_pct = case_when(fix == "D-compound" ~ 10, fix %in% c("MOP", "boron", "none") ~ 0, .default = N_pct),
         N_pct_basis = case_when(
           fix == "D-compound" ~ "terminag values_fertilizer_type (D-compound; label word order differs)",
           fix == "MOP" ~ "muriate of potash (KCl) contains no N",
           fix == "boron" ~ "boron fertilizers contain no N",
           fix == "none" ~ "no fertilizer",
           str_detect(lab, "^npks$") ~ "terminag generic N, P, K, S blend - CHECK the grade actually used",
           .default = N_pct_basis),
         terminag_fertilizer = coalesce(if_else(fix == "D-compound", "D-compound", NA_character_), terminag_fertilizer)) |>
  select(-lab, -fix)
write_csv(fert_n_stub, file.path(out_dir, "fertilizer_n_content_stub.csv"), na = "")

crop_tbl  <- tg$values$crop
crop_keys <- bind_rows(transmute(crop_tbl, terminag_crop = name, key = normalise(name)),
                       crop_tbl |> transmute(terminag_crop = name, altname) |>
                         separate_longer_delim(altname, ";") |>
                         transmute(terminag_crop, key = normalise(altname))) |> distinct()
crop_map_stub <- value_labels |> filter(role == "crop") |> distinct(value_label) |>
  mutate(lab = normalise(value_label), hit = match_vocab(lab, crop_keys)) |>
  unnest(hit, keep_empty = TRUE) |>
  left_join(select(crop_tbl, terminag_crop = name, yield_part = harvest, fresh_moisture, crop_group = group),
            by = "terminag_crop", na_matches = "never") |>
  select(value_label, terminag_crop, yield_part, fresh_moisture, crop_group)
# Malawi crop labels that the automatic match misses (reviewed by hand)
crop_fixes <- tribble(
  ~pattern,                   ~terminag_crop,       ~basis,
  "^rise",                    "rice",               "label typo of RICE",
  "^ground bean\\(nzama",     "bambara groundnut",  "nzama = bambara groundnut (Chichewa) [recall]",
  "^soyabean$",               "soybean",            "spelling variant",
  "^macadamia$",              "macadamia nut",      "terminag name",
  "^beans$",                  "common bean",        "ASSUMED Phaseolus; check",
  "^peas$",                   "pea",                "ASSUMED Pisum; check"
)
fix_crop <- function(lab) {
  hit <- crop_fixes |> filter(str_detect(str_to_lower(lab), pattern))
  if (nrow(hit) == 0) NA_character_ else hit$terminag_crop[1]
}
crop_map_stub <- crop_map_stub |>
  mutate(fixed = map_chr(value_label, fix_crop),
         basis = case_when(!is.na(terminag_crop) ~ "terminag match",
                           !is.na(fixed) ~ map_chr(value_label, \(l) crop_fixes$basis[str_detect(str_to_lower(l), crop_fixes$pattern)][1]),
                           .default = "FILL IN"),
         terminag_crop = coalesce(terminag_crop, fixed)) |>
  select(-fixed, -yield_part, -fresh_moisture, -crop_group) |>
  left_join(select(crop_tbl, terminag_crop = name, yield_part = harvest, fresh_moisture, crop_group = group),
            by = "terminag_crop", na_matches = "never")
write_csv(crop_map_stub, file.path(out_dir, "crop_map_stub.csv"), na = "")

# Roles whose terminag term should exist
role_terms <- c(hhid = "hhid", plot_id = "plot_id", crop = "crop", fertilizer_used = "fertilizer_used",
                fertilizer_type = "fertilizer_type", fertilizer_amount = "fertilizer_amount",
                intercropped = "intercropped", intercrop_fraction = "intercrop_fraction",
                OM_used = "OM_used", adm1 = "adm1", adm2 = "adm2", latitude = "latitude",
                longitude = "longitude", crop_cut = "crop_cut", cropcut_yield = "yield")
missing_terms <- setdiff(role_terms, tg$variables$name)
if (length(missing_terms) > 0) warning("terminag terms not found: ", paste(missing_terms, collapse = ", "))
message("\nDone. Review ", file.path(out_dir, "spec_check.csv"), " before running pfpn_compile.R")



# ==============================================================================
# END
# ==============================================================================
