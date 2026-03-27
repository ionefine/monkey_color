## run_frugivory_lifetime_light_phylo_updated.R
## Fits logistic and phylogenetic logistic models for routine trichromacy.
##
## Required input dataset columns:
##   - routine_trichromacy
##   - diet_pct_fruit_eltontraits
##   - either scientific_name OR genus + species_epithet


##   - one latitude option used by get_abs_latitude()
##
## This script now merges AnAge longevity directly from anage_data.txt
## (tab-delimited; uses Genus + Species and Maximum longevity (yrs)).
install.packages("ape")
suppressPackageStartupMessages({
  required_packages <- c("readr", "dplyr", "stringr", "ggplot2", "ape", "phylolm")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them first, e.g. install.packages(c(",
      paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
      "))."
    )
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
})

args <- commandArgs(trailingOnly = TRUE)
setwd("C:/Users/Ione Fine/OneDrive - UW/Documents/code/monkey_color")
arg_data_file  <- if (length(args) >= 1) args[[1]] else Sys.getenv("DATA_FILE", unset = "")
arg_anage_file <- if (length(args) >= 2) args[[2]] else Sys.getenv("ANAGE_FILE", unset = "")
arg_tree_file  <- if (length(args) >= 3) args[[3]] else Sys.getenv("TREE_FILE", unset = "")

data_file  <- if (nzchar(arg_data_file)) arg_data_file else "diurnal_primate_trichromacy_lifespan_master_v0_3_anage.csv"
anage_file <- if (nzchar(arg_anage_file)) arg_anage_file else "anage_data.txt"
tree_file  <- if (nzchar(arg_tree_file)) arg_tree_file else "primate_tree.tre"

resolve_tree_file <- function(path) {
  if (nzchar(path) && file.exists(path)) {
    return(path)
  }

  if (nzchar(path) && !file.exists(path)) {
    message("Tree file not found at provided/default path: ", path)
  }

  tree_pattern <- "\\.(tre|tree|nwk|nex|nexus)(\\.gz)?$"
  search_dirs <- c(".", "data", "trees", "input", "inputs")
  search_dirs <- unique(search_dirs[dir.exists(search_dirs)])

  candidates <- unlist(
    lapply(
      search_dirs,
      function(dir_path) list.files(
        path = dir_path,
        pattern = tree_pattern,
        ignore.case = TRUE,
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    use.names = FALSE
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  if (length(candidates) == 0) {
    candidates <- list.files(pattern = tree_pattern, ignore.case = TRUE, full.names = TRUE)
    candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  }

  if (length(candidates) == 1) {
    message("Using detected tree file: ", candidates[[1]])
    return(candidates[[1]])
  }

  if (length(candidates) > 1) {
    message(
      "Multiple tree files detected (",
      paste(candidates, collapse = ", "),
      "). Pass one explicitly as arg 3 or TREE_FILE."
    )
  } else {
    message("No phylogenetic tree file found; phylogenetic models will be skipped.")
    message(
      "Expected a tree file such as 'primate_tree.tre' (or .tree/.nwk/.nex/.nexus) ",
      "in the working directory or in ./data, ./trees, ./input, or ./inputs."
    )
  }

  NA_character_
}

tree_file <- resolve_tree_file(tree_file)

zscore_local <- function(x) {
  x <- as.numeric(x)
  x[x <= -998] <- NA_real_
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_numeric <- function(x) {
  readr::parse_number(as.character(x), na = c("", "NA", "NaN", "unknown"))
}

get_abs_latitude <- function(df) {
  nms <- names(df)
  if ("abs_latitude_mid_pantheria" %in% nms) return(as.numeric(df$abs_latitude_mid_pantheria))
  if ("abs_latitude_mid_dd_pantheria" %in% nms) return(as.numeric(df$abs_latitude_mid_dd_pantheria))
  if ("lat_midrange_dd_pantheria" %in% nms) return(abs(as.numeric(df$lat_midrange_dd_pantheria)))
  if (all(c("lat_min_dd_pantheria", "lat_max_dd_pantheria") %in% nms)) {
    lat_mid <- (as.numeric(df$lat_min_dd_pantheria) + as.numeric(df$lat_max_dd_pantheria)) / 2
    return(abs(lat_mid))
  }
  stop("No usable latitude variable found.")
}

recode_canopy_multiplier <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  out <- rep(NA_real_, length(x))
  is_canopy <- str_detect(x, "canopy|arboreal|tree")
  is_open   <- str_detect(x, "ground|terrestrial|open")
  out[is_canopy] <- 0.3
  out[is_open]   <- 1.0
  out
}

recode_canopy_class <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  out <- rep(NA_character_, length(x))
  is_canopy <- str_detect(x, "canopy|arboreal|tree")
  is_open   <- str_detect(x, "ground|terrestrial|open")
  out[is_canopy] <- "canopy"
  out[is_open]   <- "non_canopy"
  out
}

read_main_data <- function(path) {
  if (!file.exists(path)) {
    stop("Main dataset not found: ", path,
         "\nPass a file path as the first argument, e.g.\n",
         "Rscript run_frugivory_lifetime_light_phylo_updated.R your_data.csv")
  }
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("tsv", "txt")) {
    read_tsv(path, show_col_types = FALSE)
  } else {
    read_csv(path, show_col_types = FALSE)
  }
}

assert_required_columns <- function(df) {
  required <- c("routine_trichromacy", "diet_pct_fruit_eltontraits")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Main dataset is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

read_anage_data <- function(path) {
  if (!file.exists(path)) {
    stop("AnAge file not found: ", path)
  }
  anage_raw <- read_tsv(path, show_col_types = FALSE)

  required <- c("Genus", "Species", "Maximum longevity (yrs)")
  missing_cols <- setdiff(required, names(anage_raw))
  if (length(missing_cols) > 0) {
    stop("AnAge file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  anage_raw %>%
    transmute(
      scientific_name = str_trim(str_c(Genus, Species, sep = " ")),
      anage_max_longevity_years_from_txt = safe_numeric(`Maximum longevity (yrs)`)
    ) %>%
    filter(!is.na(scientific_name), scientific_name != "") %>%
    group_by(scientific_name) %>%
    summarise(
      anage_max_longevity_years_from_txt = suppressWarnings(max(anage_max_longevity_years_from_txt, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      anage_max_longevity_years_from_txt = ifelse(
        is.infinite(anage_max_longevity_years_from_txt),
        NA_real_,
        anage_max_longevity_years_from_txt
      )
    )
}

dat <- read_main_data(data_file)
assert_required_columns(dat)
anage_lookup <- read_anage_data(anage_file)

if (all(c("genus", "species_epithet") %in% names(dat))) {
  dat <- dat %>% mutate(scientific_name = str_trim(paste(genus, species_epithet)))
} else if (!"scientific_name" %in% names(dat)) {
  stop("Need either scientific_name or genus + species_epithet in main dataset.")
}

dat <- dat %>%
  left_join(anage_lookup, by = "scientific_name")

if (!"anage_max_longevity_years" %in% names(dat)) {
  dat$anage_max_longevity_years <- NA_real_
}

if (!"foraging_stratum_eltontraits" %in% names(dat)) {
  dat$foraging_stratum_eltontraits <- NA_character_
}
if (!"max_longevity_years_pantheria" %in% names(dat)) {
  dat$max_longevity_years_pantheria <- NA_real_
}

dat <- dat %>%
  mutate(
    routine_trichromacy = as.numeric(routine_trichromacy),
    diet_pct_fruit_eltontraits = safe_numeric(diet_pct_fruit_eltontraits),
    max_longevity_years_pantheria = safe_numeric(max_longevity_years_pantheria),
    anage_max_longevity_years = coalesce(
      safe_numeric(anage_max_longevity_years),
      anage_max_longevity_years_from_txt
    )
  )

if ("max_longevity_years_combined" %in% names(dat)) {
  dat <- dat %>%
    mutate(
      longevity_years_analysis = coalesce(
        safe_numeric(max_longevity_years_combined),
        max_longevity_years_pantheria,
        anage_max_longevity_years
      )
    )
} else {
  dat <- dat %>%
    mutate(longevity_years_analysis = coalesce(max_longevity_years_pantheria, anage_max_longevity_years))
}

dat <- dat %>%
  mutate(longevity_years_analysis = if_else(longevity_years_analysis > 0, longevity_years_analysis, NA_real_))

dat$abs_latitude <- get_abs_latitude(dat)

dat <- dat %>%
  mutate(
    uv_raw = cos(abs_latitude * pi / 180),
    log_longevity = log(longevity_years_analysis),
    canopy_multiplier = recode_canopy_multiplier(foraging_stratum_eltontraits),
    canopy_class = recode_canopy_class(foraging_stratum_eltontraits),
    canopy_multiplier = coalesce(canopy_multiplier, 1.0),
    canopy_class = coalesce(canopy_class, "unknown"),
    lifetime_light_raw = log_longevity * uv_raw * canopy_multiplier,
    z_fruit = zscore_local(diet_pct_fruit_eltontraits),
    z_lifetime_light = zscore_local(lifetime_light_raw)
  )

D <- dat %>%
  filter(
    !is.na(scientific_name),
    !is.na(routine_trichromacy),
    !is.na(z_fruit),
    !is.na(z_lifetime_light)
  ) %>%
  distinct(scientific_name, .keep_all = TRUE)

if (nrow(D) < 10) {
  stop("Too few complete-case species after filtering (n = ", nrow(D), "). Check input columns and missingness.")
}

message("Complete-case species for baseline model: ", nrow(D))
message("Species with longevity for analysis: ", sum(!is.na(dat$longevity_years_analysis)))
message("Species matched to AnAge TXT longevity: ", sum(!is.na(dat$anage_max_longevity_years_from_txt)))

glm_main <- glm(
  routine_trichromacy ~ z_fruit + z_lifetime_light,
  data = D,
  family = binomial()
)

glm_fruit_only <- glm(routine_trichromacy ~ z_fruit, data = D, family = binomial())
glm_light_only <- glm(routine_trichromacy ~ z_lifetime_light, data = D, family = binomial())

cat("\n===== STANDARD LOGISTIC MODELS =====\n")
print(summary(glm_main))
print(summary(glm_fruit_only))
print(summary(glm_light_only))

cat("\nAIC comparison:\n")
print(AIC(glm_main, glm_fruit_only, glm_light_only))

cat("\nOdds ratios (main model):\n")
print(exp(cbind(Estimate = coef(glm_main), confint.default(glm_main))))

if (is.na(tree_file) || !file.exists(tree_file)) {
  message("Skipping phylogenetic models because no usable tree file was found.")
} else {
  tree <- read.tree(tree_file)
  keep <- intersect(tree$tip.label, D$scientific_name)
  if (length(keep) < 20) {
    warning("Very few species overlap between tree and data; skipping phylogenetic models.")
  } else {
    tree2 <- drop.tip(tree, setdiff(tree$tip.label, keep))
    D_phy <- D %>%
      filter(scientific_name %in% keep) %>%
      arrange(match(scientific_name, tree2$tip.label))
    stopifnot(identical(tree2$tip.label, D_phy$scientific_name))

    message("Species in phylogenetic model: ", nrow(D_phy))

    phy_main <- phyloglm(
      routine_trichromacy ~ z_fruit + z_lifetime_light,
      data = D_phy,
      phy = tree2,
      method = "logistic_MPLE",
      btol = 50
    )

    phy_fruit_only <- phyloglm(
      routine_trichromacy ~ z_fruit,
      data = D_phy,
      phy = tree2,
      method = "logistic_MPLE",
      btol = 50
    )

    phy_light_only <- phyloglm(
      routine_trichromacy ~ z_lifetime_light,
      data = D_phy,
      phy = tree2,
      method = "logistic_MPLE",
      btol = 50
    )

    cat("\n===== PHYLOGENETIC LOGISTIC MODELS =====\n")
    print(summary(phy_main))
    print(summary(phy_fruit_only))
    print(summary(phy_light_only))

    phy_compare <- tibble(
      model = c("phy_main", "phy_fruit_only", "phy_light_only"),
      logLik = c(logLik(phy_main), logLik(phy_fruit_only), logLik(phy_light_only)),
      AIC = c(AIC(phy_main), AIC(phy_fruit_only), AIC(phy_light_only))
    )

    cat("\nPhylogenetic model comparison:\n")
    print(phy_compare)

    cat("\nApproximate odds ratios (phylogenetic main model):\n")
    print(exp(coef(phy_main)))

    pred_grid <- expand.grid(
      z_lifetime_light = seq(min(D_phy$z_lifetime_light), max(D_phy$z_lifetime_light), length.out = 200),
      z_fruit = c(-1, 0, 1)
    )

    pred_grid$pred_prob <- predict(phy_main, newdata = pred_grid, type = "response")
    pred_grid$fruit_level <- factor(
      pred_grid$z_fruit,
      levels = c(-1, 0, 1),
      labels = c("Low frugivory (-1 SD)", "Mean frugivory", "High frugivory (+1 SD)")
    )

    p <- ggplot(pred_grid, aes(x = z_lifetime_light, y = pred_prob, color = fruit_level)) +
      geom_line(linewidth = 1) +
      labs(
        x = "Standardized lifetime light exposure",
        y = "Predicted P(routine trichromacy)",
        color = "",
        title = "Phylogenetic logistic regression",
        subtitle = "Routine trichromacy predicted by frugivory and lifetime light exposure"
      ) +
      theme_minimal(base_size = 12)

    print(p)
    ggsave("phylogenetic_trichromacy_predictions_updated.png", p, width = 8, height = 5, dpi = 300)

    analysis_subset <- D_phy %>%
      select(any_of(c(
        "scientific_name",
        "routine_trichromacy",
        "longevity_years_analysis",
        "max_longevity_years_pantheria",
        "anage_max_longevity_years",
        "anage_max_longevity_years_from_txt",
        "anage_match_type",
        "anage_match_name",
        "diet_pct_fruit_eltontraits",
        "z_fruit",
        "log_longevity",
        "abs_latitude",
        "uv_raw",
        "foraging_stratum_eltontraits",
        "canopy_class",
        "canopy_multiplier",
        "lifetime_light_raw",
        "z_lifetime_light"
      )))

    write_csv(analysis_subset, "frugivory_lifetime_light_phylo_analysis_subset_updated.csv")

    cat("\nSaved:\n")
    cat("  frugivory_lifetime_light_phylo_analysis_subset_updated.csv\n")
    cat("  phylogenetic_trichromacy_predictions_updated.png\n")
  }
}
