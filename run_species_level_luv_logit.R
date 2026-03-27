#!/usr/bin/env Rscript

# Species-level logistic analysis of routine trichromacy vs frugivory + lifetime UV exposure.
#
# Data integration assumptions:
# - Main dataset contains species-level records derived from PanTHERIA + EltonTraits.
# - AnAge longevity is supplied in tab-delimited format (default: anage_data.txt).
#
# Usage:
#   Rscript run_species_level_luv_logit.R [main_data.csv] [anage_data.txt]

suppressPackageStartupMessages({
  required_packages <- c("readr", "dplyr", "stringr", "tibble")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing package(s): ", paste(missing_packages, collapse = ", "),
      "\nInstall with install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))."
    )
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
})

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[[1]] else "diurnal_primate_trichromacy_lifespan_master_v0_3_anage.csv"
anage_file <- if (length(args) >= 2) args[[2]] else "anage_data.txt"

safe_numeric <- function(x) {
  readr::parse_number(as.character(x), na = c("", "NA", "NaN", "unknown", "-999", "-998"))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

get_abs_latitude <- function(df) {
  nms <- names(df)
  if ("abs_latitude_mid_dd" %in% nms) return(safe_numeric(df$abs_latitude_mid_dd))
  if ("abs_latitude_mid_pantheria" %in% nms) return(safe_numeric(df$abs_latitude_mid_pantheria))
  if ("lat_mid_dd_pantheria" %in% nms) return(abs(safe_numeric(df$lat_mid_dd_pantheria)))
  if ("lat_midrange_dd_pantheria" %in% nms) return(abs(safe_numeric(df$lat_midrange_dd_pantheria)))
  if (all(c("lat_min_dd_pantheria", "lat_max_dd_pantheria") %in% nms)) {
    lat_mid <- (safe_numeric(df$lat_min_dd_pantheria) + safe_numeric(df$lat_max_dd_pantheria)) / 2
    return(abs(lat_mid))
  }
  stop("No usable latitude field found (expected abs latitude or range midpoint columns).")
}

classify_canopy_multiplier <- function(x) {
  x <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  out <- rep(NA_real_, length(x))
  is_canopy <- stringr::str_detect(x, "canopy|arboreal|tree")
  is_non_canopy <- stringr::str_detect(x, "ground|terrestrial|open")
  out[is_canopy] <- 0.1
  out[is_non_canopy] <- 1.0
  out
}

load_anage_longevity <- function(path) {
  if (!file.exists(path)) stop("AnAge file not found: ", path)
  an <- readr::read_tsv(path, show_col_types = FALSE)
  needed <- c("Genus", "Species", "Maximum longevity (yrs)")
  miss <- setdiff(needed, names(an))
  if (length(miss) > 0) stop("AnAge is missing required columns: ", paste(miss, collapse = ", "))

  an %>%
    transmute(
      scientific_name = stringr::str_trim(stringr::str_c(Genus, Species, sep = " ")),
      anage_longevity_years = safe_numeric(`Maximum longevity (yrs)`)
    ) %>%
    filter(!is.na(scientific_name), scientific_name != "") %>%
    group_by(scientific_name) %>%
    summarise(anage_longevity_years = suppressWarnings(max(anage_longevity_years, na.rm = TRUE)), .groups = "drop") %>%
    mutate(anage_longevity_years = ifelse(is.infinite(anage_longevity_years), NA_real_, anage_longevity_years))
}

if (!file.exists(data_file)) stop("Main dataset not found: ", data_file)
main_df <- readr::read_csv(data_file, show_col_types = FALSE)

required <- c("routine_trichromacy", "diet_pct_fruit_eltontraits", "foraging_stratum_eltontraits")
missing_required <- setdiff(required, names(main_df))
if (length(missing_required) > 0) {
  stop("Main dataset missing required column(s): ", paste(missing_required, collapse = ", "))
}

if (!("scientific_name" %in% names(main_df))) {
  if (all(c("genus", "species_epithet") %in% names(main_df))) {
    main_df <- main_df %>% mutate(scientific_name = stringr::str_trim(paste(genus, species_epithet)))
  } else {
    stop("Need scientific_name OR genus + species_epithet.")
  }
}

anage_lookup <- load_anage_longevity(anage_file)

if (!("max_longevity_years_pantheria" %in% names(main_df))) main_df$max_longevity_years_pantheria <- NA_real_
if (!("included_in_core_non_nocturnal_build" %in% names(main_df))) main_df$included_in_core_non_nocturnal_build <- NA
if (!("order" %in% names(main_df))) main_df$order <- NA_character_

analysis_df <- main_df %>%
  mutate(
    order = stringr::str_to_lower(as.character(order)),
    routine_trichromacy = safe_numeric(routine_trichromacy),
    diet_pct_fruit_eltontraits = safe_numeric(diet_pct_fruit_eltontraits),
    max_longevity_years_pantheria = safe_numeric(max_longevity_years_pantheria),
    included_in_core_non_nocturnal_build = as.character(included_in_core_non_nocturnal_build)
  ) %>%
  left_join(anage_lookup, by = "scientific_name") %>%
  mutate(
    longevity_years = dplyr::coalesce(max_longevity_years_pantheria, anage_longevity_years),
    longevity_years = if_else(longevity_years > 0, longevity_years, NA_real_),
    log_longevity = log(longevity_years),
    abs_latitude = get_abs_latitude(cur_data_all()),
    uv = cos(abs_latitude * pi / 180),
    arboreal = classify_canopy_multiplier(foraging_stratum_eltontraits),
    arboreal = dplyr::coalesce(arboreal, 1.0),
    uv_arb =  uv * arboreal,
    luv_arb =  uv_arb * longevity_years,
    z_frugivory = zscore(diet_pct_fruit_eltontraits),
    z_uv_arb = zscore(uv_arb),
    z_luv_arb = zscore(luv_arb),
    z_lifetime = zscore(longevity_years)
  ) %>%
  filter(
    order == "primates" | is.na(order),
    tolower(included_in_core_non_nocturnal_build) %in% c("true", "t", "1", "yes", "y") |
      is.na(included_in_core_non_nocturnal_build)
  ) %>%
  distinct(scientific_name, .keep_all = TRUE) %>%
  filter(!is.na(routine_trichromacy), !is.na(z_frugivory), !is.na(uv_arb))

if (nrow(analysis_df) < 10) {
  stop("Too few complete species for logistic regression after filtering (n = ", nrow(analysis_df), ").")
}

model_full <- glm(routine_trichromacy ~ z_frugivory + z_luv_arb + z_lifetime, data = analysis_df, family = binomial())
model_frug <- glm(routine_trichromacy ~ z_frugivory, data = analysis_df, family = binomial())
model_luv <- glm(routine_trichromacy ~ z_luv_arb, data = analysis_df, family = binomial())
model_lifetime <- glm(routine_trichromacy ~ z_lifetime, data = analysis_df, family = binomial())
cat("\n===== DATA SUMMARY =====\n")
cat("Species analyzed:", nrow(analysis_df), "\n")
cat("Routine trichromat species:", sum(analysis_df$routine_trichromacy == 1, na.rm = TRUE), "\n")
cat("\n===== FULL BINOMIAL LOGIT =====\n")
print(summary(model_full))
cat("\n===== REDUCED MODELS =====\n")
print(summary(model_frug))
print(summary(model_luv))
print(summary(model_lifetime))
cat("\n===== AIC COMPARISON =====\n")
print(AIC(model_full, model_frug, model_uv, model_lifetime))

coef_tbl <- tibble::tibble(
  term = names(coef(model_full)),
  log_odds = as.numeric(coef(model_full)),
  odds_ratio = exp(log_odds)
)
cat("\n===== FULL MODEL COEFFICIENTS =====\n")
print(coef_tbl)

analysis_export <- analysis_df %>%
  select(
    scientific_name,
    routine_trichromacy,
    diet_pct_fruit_eltontraits,
    z_frugivory,
    longevity_years,
    log_longevity,
    abs_latitude,
    uv,
    foraging_stratum_eltontraits,
    arboreal,
    uv_arb,
    luv_arb,
    z_uv_arb
  )

readr::write_csv(analysis_export, "species_level_luv_analysis_dataset.csv")
readr::write_csv(coef_tbl, "species_level_luv_full_model_coefficients.csv")

cat("\nSaved outputs:\n")
cat(" - species_level_luv_analysis_dataset.csv\n")
cat(" - species_level_luv_full_model_coefficients.csv\n")
