# =============================================================================
# 00.pfpn_utils.R: helpers shared by all PFP-N scripts (sourced at their top)
# -----------------------------------------------------------------------------
# - pfpn_path(): every output folder sits under one root folder, set with
#   options(pfpn.root = "...") (run_all.R does this); default "." = working dir
# - fetch(): download once, then re-download only when the remote file changed
#   (ETag, else HTTP Last-Modified, else size); offline or without that
#   information, the existing file is kept
# - unzip_once(): extract zip entries once; re-extract when the zip is newer
# - in_malawi(): coordinates inside Malawi's bounding box ((0, 0) etc. fail)
# =============================================================================

# clean global environment
rm(list = ls())

# load core packages
library(tidyverse)
# Other packages used via pkg::fun (xml2, readxl, terra, sf, multcompView, emmeans, multcomp)

pfpn_path <- function(...) file.path(getOption("pfpn.root", "."), ...)

ensure <- function(df, cols) { for (c in cols) if (!c %in% names(df)) df[[c]] <- NA; df }

# Malawi bounding box (slightly padded); used to reject impossible coordinates
mwi_bbox  <- c(xmin = 32.6, xmax = 36.0, ymin = -17.2, ymax = -9.3)
in_malawi <- function(lon, lat) {
  !is.na(lon) & !is.na(lat) &
    lon >= mwi_bbox[["xmin"]] & lon <= mwi_bbox[["xmax"]] &
    lat >= mwi_bbox[["ymin"]] & lat <= mwi_bbox[["ymax"]]
}

# "Lilongwe City" -> "lilongwe", "Mzimba North" -> "mzimba", "Nkhata Bay" -> "nkhatabay"
# Keys (after district_key()) that differ from GADM district names
district_aliases <- c(
  mzuzu       = "mzimba",     # Mzuzu City lies in Mzimba District [recall]
  blanytyre   = "blantyre",   # misspelling of Blantyre in some waves
  dzalekacamp = "dowa",       # Dzaleka refugee camp lies in Dowa District [recall]
  zombanon    = "zomba"       # "Zomba Non-City"
)
district_key <- function(x) {
  k <- x |> str_to_lower() |>
    str_remove_all("\\b(city|boma|urban|rural|district|north|south)\\b") |>
    str_remove_all("[^a-z]")
  coalesce(unname(district_aliases[k]), k)
}

# --- file checks --------------------------------------------------------------
is_zip  <- \(f) nrow(utils::unzip(f, list = TRUE)) > 0
is_xml  <- \(f) inherits(xml2::read_xml(f), "xml_document")
is_xlsx <- \(f) length(readxl::excel_sheets(f)) > 0
is_tif  <- \(f) !is.null(terra::rast(f))
# A DDI codebook, not an HTML page that happens to parse as XML (error or sign-in page)
is_ddi  <- \(f) {
  doc <- xml2::xml_ns_strip(xml2::read_xml(f))
  xml2::xml_name(doc) == "codeBook" && length(xml2::xml_find_all(doc, "//dataDscr/var")) > 0
}
usable <- \(f, check = \(f) TRUE) file.exists(f) && file.size(f) > 0 &&
  isTRUE(tryCatch(check(f), error = \(e) FALSE))

# --- update-aware download ----------------------------------------------------
# pfpn.refresh option: "if_updated" (default), "never" (always reuse), "always"
# pfpn.download option: download function (url, destfile); lets tests run offline

http_date <- function(x) {
  # "Wed, 23 Sep 2026 07:30:00 GMT", parsed without depending on the locale
  m <- str_match(x, "(\\d{1,2}) ([A-Za-z]{3}) (\\d{4}) (\\d{2}):(\\d{2}):(\\d{2})")
  if (is.na(m[1, 1])) return(as.POSIXct(NA))
  mon <- match(str_to_lower(m[1, 3]), c("jan", "feb", "mar", "apr", "may", "jun",
                                        "jul", "aug", "sep", "oct", "nov", "dec"))
  ISOdatetime(as.integer(m[1, 4]), mon, as.integer(m[1, 2]), as.integer(m[1, 5]),
              as.integer(m[1, 6]), as.integer(m[1, 7]), tz = "GMT")
}

remote_info <- function(url) {
  h <- tryCatch(suppressWarnings(curlGetHeaders(url, timeout = 20)), error = \(e) NULL)
  if (is.null(h) || isTRUE(attr(h, "status") >= 400)) return(NULL)
  # only the final response counts (redirects come first, e.g. GitHub's 302)
  starts <- grep("^HTTP/", h)
  if (length(starts)) h <- h[tail(starts, 1):length(h)]
  field <- \(k) {
    l <- grep(paste0("^", k, ":"), h, ignore.case = TRUE, value = TRUE)
    if (length(l)) str_trim(sub("^[^:]+:", "", tail(l, 1))) else NA_character_
  }
  modified <- http_date(field("last-modified"))
  served   <- http_date(field("date"))
  # a Last-Modified equal to the response time means the file is generated on
  # request (no real modification date): ignore it
  if (!is.na(modified) && !is.na(served) && abs(as.numeric(served - modified, units = "secs")) < 120) modified <- as.POSIXct(NA)
  size <- suppressWarnings(as.numeric(field("content-length")))
  list(modified = modified, size = if (isTRUE(size > 0)) size else NA_real_, etag = field("etag"))
}

download_to <- function(url, dest) {
  dl <- getOption("pfpn.download", \(url, destfile) utils::download.file(url, destfile, mode = "wb", quiet = TRUE))
  tmp <- paste0(dest, ".part")
  dl(url, tmp)
  tmp
}

fetch <- function(url, dest, check = \(f) TRUE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  mode <- getOption("pfpn.refresh", "if_updated")
  have <- usable(dest, check)

  if (have && mode == "never") return(dest)
  etag_file <- paste0(dest, ".etag")
  if (have && mode == "if_updated") {
    info <- remote_info(url)
    if (is.null(info)) return(dest)                          # offline / no headers: keep
    changed <- if (!is.na(info$etag)) {
      if (!file.exists(etag_file)) { writeLines(info$etag, etag_file); FALSE }   # first check: record it
      else info$etag != readLines(etag_file, warn = FALSE)[1]
    } else if (!is.na(info$modified)) {
      info$modified > file.mtime(dest) + 60
    } else if (!is.na(info$size)) {
      info$size != file.size(dest)
    } else FALSE                                             # nothing to compare: keep
    if (!changed) return(dest)                               # up to date
    message("Remote file has changed, updating: ", basename(dest))
  }
  if (file.exists(dest) && !have) {
    kept <- paste0(dest, ".unusable_", format(Sys.time(), "%Y%m%d%H%M%S"))
    file.rename(dest, kept)
    message("Existing file not usable (kept as ", basename(kept), "); downloading ", basename(dest))
  }
  tmp <- tryCatch(download_to(url, dest), error = \(e) NULL)
  if (is.null(tmp) || !usable(tmp, check)) {
    if (!is.null(tmp)) unlink(tmp)
    if (have) { message("Update failed, keeping the existing ", basename(dest)); return(dest) }
    stop("Download failed or file not usable: ", url, call. = FALSE)
  }
  if (have) file.rename(dest, paste0(dest, ".previous"))     # one backup of the replaced version
  file.rename(tmp, dest)
  info <- remote_info(url)
  if (!is.null(info) && !is.na(info$modified)) Sys.setFileTime(dest, info$modified)
  if (!is.null(info) && !is.na(info$etag)) writeLines(info$etag, etag_file)
  dest
}

# --- zip extraction -----------------------------------------------------------
# Extract entries that are missing, or older than the zip itself (zip replaced).
# Extracted files get the extraction time, so later comparisons are meaningful.
unzip_once <- function(zip, exdir, files = NULL) {
  entries <- utils::unzip(zip, list = TRUE)$Name
  entries <- entries[!str_ends(entries, "/")]
  if (!is.null(files)) entries <- intersect(entries, files)
  targets <- file.path(exdir, entries)
  stale <- !file.exists(targets) | file.mtime(targets) < file.mtime(zip)
  if (any(stale)) {
    utils::unzip(zip, files = entries[stale], exdir = exdir, overwrite = TRUE)
    Sys.setFileTime(targets[stale], Sys.time())
  }
  targets
}

# Output of a step is out of date when any input is newer (used for derived files)
outdated <- function(output, inputs) {
  !file.exists(output) || any(file.mtime(inputs[file.exists(inputs)]) > file.mtime(output))
}

# --- compact letter display for group comparisons -----------------------------
# Pairwise Wilcoxon rank-sum tests with Bonferroni adjustment; groups sharing a
# letter do not differ at `alpha`. A symmetric p-value matrix is passed to
# multcompView so group names containing "-" are safe. Groups with < 3 values
# get no letter.
# One-way ANOVA, estimated marginal means (emmeans) and compact letter display
# (multcomp::cld) with Bonferroni adjustment. reversed = TRUE gives "a" to the
# highest mean, so letters decrease with the mean. Returns the group means and
# letters, with the ANOVA F test attached as attribute "anova".
anova_cld_letters <- function(value, group, alpha = 0.05) {
  ok <- !is.na(value) & !is.na(group)
  df <- tibble(value = value[ok], group = droplevels(as.factor(group[ok])))
  df <- df |> filter(group %in% names(which(table(group) >= 3))) |> mutate(group = droplevels(group))
  if (nlevels(df$group) < 2) return(tibble(group = character(), letter = character()))
  fit <- lm(value ~ group, data = df)
  em  <- emmeans::emmeans(fit, ~ group)
  cl  <- multcomp::cld(em, adjust = "bonferroni", Letters = letters, reversed = TRUE, alpha = alpha)
  out <- as_tibble(as.data.frame(cl)) |>
    transmute(group = as.character(group), emmean, SE, lower_CL = lower.CL, upper_CL = upper.CL,
              letter = str_trim(.group))
  a <- anova(fit)
  attr(out, "anova") <- tibble(df_group = a$Df[1], df_resid = a$Df[2], F = a$`F value`[1], p_value = a$`Pr(>F)`[1])
  out
}

# Write the data behind a ggplot (and every layer, and patchwork panels) to
# tables/fig_<name>.csv (fig_<name>_2.csv, ... when layers use different data)
save_plot_data <- function(p, name, dir) {
  plots <- if (inherits(p, "patchwork")) c(list(p), p$patches$plots) else list(p)
  dfs <- list()
  for (q in plots) {
    if (!inherits(q, "ggplot")) next
    for (x in c(list(q$data), lapply(q$layers, \(l) l$data))) {
      if (!is.data.frame(x) || nrow(x) == 0) next
      if (inherits(x, "sf")) x <- sf::st_drop_geometry(x)
      x <- as_tibble(x) |> select(where(\(col) !is.list(col)))
      if (!any(map_lgl(dfs, \(d) identical(d, x)))) dfs[[length(dfs) + 1]] <- x
    }
  }
  for (i in seq_along(dfs)) {
    suffix <- if (i == 1) "" else paste0("_", i)
    write_csv(dfs[[i]], file.path(dir, paste0("fig_", name, suffix, ".csv")), na = "")
  }
  invisible(length(dfs))
}

# Monte Carlo Shapley values for every row of X (the sampling estimator of
# Strumbelj & Kononenko 2014, as used by iml::Shapley, vectorised over rows):
# for each feature j and simulation, a random feature order and a random
# background row; features before j come from the plot, the rest from the
# background row; phi_j = mean of f(with j) - f(without j). pred(df) must
# return numeric predictions. Returns an n x p matrix.
mc_shapley <- function(pred, X, nsim = 20) {
  X <- as.data.frame(X)
  n <- nrow(X); p <- ncol(X)
  Xm  <- as.matrix(X)
  phi <- matrix(0, n, p, dimnames = list(NULL, names(X)))
  for (s in seq_len(nsim)) {
    ranks <- t(apply(matrix(runif(n * p), n, p), 1, rank))
    Wm <- Xm[sample.int(n, n, replace = TRUE), , drop = FALSE]
    for (j in seq_len(p)) {
      before <- ranks < ranks[, j]
      plus  <- ifelse(before | col(before) == j, Xm, Wm)
      minus <- ifelse(before, Xm, Wm)
      both  <- as.data.frame(rbind(plus, minus)); names(both) <- names(X)
      f <- pred(both)
      phi[, j] <- phi[, j] + (f[seq_len(n)] - f[n + seq_len(n)]) / nsim
    }
  }
  phi
}

cld_letters <- function(value, group, alpha = 0.05) {
  ok <- !is.na(value) & !is.na(group)
  value <- value[ok]
  group <- droplevels(as.factor(group[ok]))
  keep <- names(which(table(group) >= 3))
  if (length(keep) < 2) return(tibble(group = character(), letter = character()))
  sel <- group %in% keep
  value <- value[sel]; group <- droplevels(group[sel])
  pw <- suppressWarnings(pairwise.wilcox.test(value, group, p.adjust.method = "bonferroni", exact = FALSE))
  lev <- levels(group)
  pm <- matrix(1, length(lev), length(lev), dimnames = list(lev, lev))
  low <- pw$p.value
  for (i in rownames(low)) for (j in colnames(low)) if (!is.na(low[i, j])) pm[i, j] <- pm[j, i] <- low[i, j]
  lt <- multcompView::multcompLetters(pm, threshold = alpha)$Letters
  tibble(group = names(lt), letter = unname(lt))
}

# ==============================================================================
# END
# ==============================================================================