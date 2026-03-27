#!/usr/bin/env Rscript

# Combined species-level analysis of routine trichromacy vs frugivory +
# lifetime UV exposure. This script runs:
#   1) standard (non-phylogenetic) logistic regressions (same specification
#      as run_species_level_luv_logit.R), and
#   2) phylogenetic logistic regressions with tree matching/remapping.
#
# Usage:
#   Rscript phylo_code.R [main_data.csv] [anage_data.txt] [tree_file] [remap_file]

suppressPackageStartupMessages({
  required_packages <- c("readr", "dplyr", "stringr", "tibble", "ape", "phylolm")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing package(s): ", paste(missing_packages, collapse = ", "),
      "\nInstall with install.packages(c(",
      paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
      "))."
    )
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
})

args <- commandArgs(trailingOnly = TRUE)
data_file  <- if (length(args) >= 1) args[[1]] else "diurnal_primate_trichromacy_lifespan_master_v0_3_anage.csv"
anage_file <- if (length(args) >= 2) args[[2]] else "anage_data.txt"
tree_file  <- if (length(args) >= 3) args[[3]] else "consensusTree_10kTrees_Primates_Version3.nex"
remap_file <- if (length(args) >= 4) args[[4]] else "remapping.csv"

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

clean_species_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\.", " ")
  x <- stringr::str_replace_all(x, "-", " ")
  x <- stringr::str_replace_all(x, "\\(.*?\\)", " ")
  x <- stringr::str_replace_all(x, "\\[.*?\\]", " ")
  x <- stringr::str_replace(x, "\\b(sp|species|cf|aff)\\b\\.?$", "")
  x <- stringr::str_replace(x, "\\b(voucher|sample|isolate|strain)\\b.*$", "")
  x <- stringr::str_squish(x)
  x <- stringr::str_to_lower(x)
  
  toks <- stringr::str_split(x, "\\s+", simplify = TRUE)
  first <- toks[, 1]
  second <- toks[, 2]
  
  out <- ifelse(first != "" & second != "", paste(first, second), x)
  out <- stringr::str_squish(out)
  out[out == ""] <- NA_character_
  out
}

clean_relation_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\.", " ")
  x <- stringr::str_replace_all(x, "-", " ")
  x <- stringr::str_replace_all(x, "\\(.*?\\)", " ")
  x <- stringr::str_replace_all(x, "\\[.*?\\]", " ")
  x <- stringr::str_squish(x)
  x <- stringr::str_to_lower(x)
  x[x == ""] <- NA_character_
  x
}

display_species_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  x
}

load_phylo_tree <- function(path) {
  if (!file.exists(path)) stop("Tree file not found: ", path)
  
  ext <- tolower(tools::file_ext(path))
  
  tree <- if (ext %in% c("nwk", "newick", "tre")) {
    ape::read.tree(path)
  } else if (ext %in% c("nex", "nexus")) {
    ape::read.nexus(path)
  } else {
    stop(
      "Unsupported tree file extension: ", ext,
      "\nUse one of: .nwk, .newick, .tre, .nex, .nexus"
    )
  }
  
  if (inherits(tree, "multiPhylo")) {
    message("Tree file contains multiple trees; using the first tree.")
    tree <- tree[[1]]
  }
  
  if (!inherits(tree, "phylo")) stop("Loaded tree is not a valid 'phylo' object.")
  if (is.null(tree$tip.label) || length(tree$tip.label) == 0) stop("Tree has no tip labels.")
  if (is.null(tree$edge.length)) stop("Tree has no branch lengths. phyloglm requires branch lengths.")
  
  tree
}

plot_and_save_tree <- function(tree, pdf_file, png_file, width = 12, height = 16, cex_tip = NULL) {
  n_tips <- length(tree$tip.label)
  
  if (is.null(cex_tip)) {
    cex_tip <- if (n_tips <= 40) {
      0.9
    } else if (n_tips <= 80) {
      0.7
    } else if (n_tips <= 150) {
      0.5
    } else {
      0.35
    }
  }
  
  grDevices::pdf(pdf_file, width = width, height = height)
  par(mar = c(1, 1, 3, 1))
  plot.phylo(tree, type = "phylogram", cex = cex_tip, no.margin = TRUE)
  title(main = paste0("Pruned/grafted primate tree (n = ", n_tips, " tips)"))
  grDevices::dev.off()
  
  grDevices::png(png_file, width = 1800, height = 2400, res = 220)
  par(mar = c(1, 1, 3, 1))
  plot.phylo(tree, type = "phylogram", cex = cex_tip, no.margin = TRUE)
  title(main = paste0("Pruned/grafted primate tree (n = ", n_tips, " tips)"))
  grDevices::dev.off()
}

load_remap_table <- function(path) {
  if (!file.exists(path)) {
    message("Remapping file not found at: ", path, " ; proceeding without remaps.")
    return(tibble::tibble(
      scientific_name_clean = character(),
      tip_relation = character(),
      other_name = character(),
      action = character()
    ))
  }
  
  mp <- suppressWarnings(readr::read_csv(path, locale = readr::locale(encoding = "Latin1"), show_col_types = FALSE))
  
  needed <- c("scientific_name_clean", "tip_relation", "other name")
  miss <- setdiff(needed, names(mp))
  if (length(miss) > 0) {
    stop("Remapping file is missing required columns: ", paste(miss, collapse = ", "))
  }
  
  mp %>%
    transmute(
      scientific_name_clean = clean_species_name(scientific_name_clean),
      tip_relation = clean_relation_name(tip_relation),
      other_name = clean_species_name(`other name`)
    ) %>%
    mutate(
      action = dplyr::case_when(
        !is.na(other_name) & other_name != "" ~ "merge",
        !is.na(tip_relation) & tip_relation != "" ~ "graft",
        TRUE ~ "none"
      )
    ) %>%
    filter(!is.na(scientific_name_clean), scientific_name_clean != "") %>%
    distinct(scientific_name_clean, .keep_all = TRUE)
}

apply_merge_map <- function(x_clean, remap_tbl) {
  if (nrow(remap_tbl) == 0) return(x_clean)
  merge_tbl <- remap_tbl %>% filter(action == "merge")
  if (nrow(merge_tbl) == 0) return(x_clean)
  
  idx <- match(x_clean, merge_tbl$scientific_name_clean)
  replacement <- merge_tbl$other_name[idx]
  dplyr::coalesce(replacement, x_clean)
}

find_anchor_node <- function(tree, relation_clean) {
  relation_clean <- clean_relation_name(relation_clean)
  tip_clean <- clean_species_name(tree$tip.label)
  
  exact_tip_idx <- which(tip_clean == relation_clean)
  if (length(exact_tip_idx) == 1) {
    return(list(type = "tip", tip_index = exact_tip_idx, node = NULL, relation = relation_clean))
  }
  if (length(exact_tip_idx) > 1) {
    stop("Multiple exact tip matches found for relation: ", relation_clean)
  }
  
  genus <- strsplit(relation_clean, "\\s+")[[1]][1]
  if (is.na(genus) || genus == "") {
    stop("Could not parse genus from tip_relation: ", relation_clean)
  }
  
  tip_genus <- vapply(strsplit(tip_clean, "\\s+"), `[`, character(1), 1)
  genus_tips <- which(tip_genus == genus)
  
  if (length(genus_tips) == 0) {
    stop("No tips found in tree for tip_relation genus: ", relation_clean)
  }
  
  if (length(genus_tips) == 1) {
    return(list(type = "tip", tip_index = genus_tips, node = NULL, relation = relation_clean))
  }
  
  mrca_node <- ape::getMRCA(tree, genus_tips)
  if (is.null(mrca_node) || is.na(mrca_node)) {
    stop("Could not compute MRCA for genus relation: ", relation_clean)
  }
  
  list(type = "node", tip_index = NULL, node = mrca_node, relation = relation_clean)
}

graft_tip_at_relation <- function(tree, new_tip_label, relation_clean,
                                  branch_fraction = 0.5,
                                  min_branch = 1e-6) {
  anchor <- find_anchor_node(tree, relation_clean)
  
  new_tip <- list(
    edge = matrix(c(2L, 1L), nrow = 1, ncol = 2),
    tip.label = new_tip_label,
    Nnode = 1L,
    edge.length = min_branch
  )
  class(new_tip) <- "phylo"
  
  if (anchor$type == "tip") {
    tip_idx <- anchor$tip_index
    
    parent_row <- which(tree$edge[, 2] == tip_idx)
    if (length(parent_row) != 1) {
      stop("Could not identify unique parent edge for tip: ", tree$tip.label[tip_idx])
    }
    
    old_len <- tree$edge.length[parent_row]
    if (is.na(old_len) || old_len <= 0) old_len <- min_branch * 4
    
    split_len <- max(old_len * branch_fraction, min_branch)
    
    tree2 <- ape::bind.tree(
      x = tree,
      y = new_tip,
      where = tip_idx,
      position = split_len
    )
    
    return(tree2)
  }
  
  if (anchor$type == "node") {
    tree2 <- ape::bind.tree(
      x = tree,
      y = new_tip,
      where = anchor$node,
      position = 0
    )
    
    return(tree2)
  }
  
  stop("Unsupported anchor type for grafting.")
}

apply_graft_plan <- function(tree, graft_tbl) {
  if (nrow(graft_tbl) == 0) {
    return(list(tree = tree, audit = tibble::tibble()))
  }
  
  audit_rows <- list()
  out_tree <- tree
  
  for (i in seq_len(nrow(graft_tbl))) {
    sp <- graft_tbl$scientific_name_clean[i]
    rel <- graft_tbl$tip_relation[i]
    before_n <- length(out_tree$tip.label)
    
    res <- tryCatch({
      new_tree <- graft_tip_at_relation(out_tree, new_tip_label = sp, relation_clean = rel)
      list(ok = TRUE, tree = new_tree, error = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, tree = out_tree, error = conditionMessage(e))
    })
    
    out_tree <- res$tree
    after_n <- length(out_tree$tip.label)
    
    audit_rows[[i]] <- tibble::tibble(
      scientific_name_clean = sp,
      tip_relation = rel,
      grafted = isTRUE(res$ok) && after_n == before_n + 1,
      error_message = res$error
    )
  }
  
  list(tree = out_tree, audit = dplyr::bind_rows(audit_rows))
}

deduplicate_tree_by_clean_species <- function(tree) {
  tip_tbl <- tibble::tibble(
    tree_tip_raw = tree$tip.label,
    tree_tip_clean = clean_species_name(tree$tip.label)
  )
  
  dup_tbl <- tip_tbl %>%
    group_by(tree_tip_clean) %>%
    summarise(
      n_raw_tips = dplyr::n(),
      raw_tips = paste(tree_tip_raw, collapse = "; "),
      .groups = "drop"
    ) %>%
    filter(!is.na(tree_tip_clean), n_raw_tips > 1)
  
  keep_tbl <- tip_tbl %>%
    filter(!is.na(tree_tip_clean)) %>%
    group_by(tree_tip_clean) %>%
    slice(1) %>%
    ungroup()
  
  keep_raw <- keep_tbl$tree_tip_raw
  drop_raw <- setdiff(tree$tip.label, keep_raw)
  
  tree_dedup <- if (length(drop_raw) > 0) ape::drop.tip(tree, drop_raw) else tree
  
  match_idx <- match(tree_dedup$tip.label, keep_tbl$tree_tip_raw)
  tree_dedup$tip.label <- keep_tbl$tree_tip_clean[match_idx]
  
  list(
    tree = tree_dedup,
    duplicate_table = dup_tbl,
    keep_table = keep_tbl,
    dropped_raw_tips = tibble::tibble(tree_tip_raw_dropped = drop_raw)
  )
}

resolve_binary_with_target_priority <- function(x, original_names, final_name) {
  keep <- !is.na(x)
  x <- x[keep]
  original_names <- original_names[keep]
  
  if (length(x) == 0) return(NA_real_)
  
  ux <- unique(x)
  if (length(ux) == 1) return(as.numeric(ux))
  
  target_idx <- which(original_names == final_name)
  if (length(target_idx) >= 1) {
    target_vals <- unique(x[target_idx])
    if (length(target_vals) == 1) return(as.numeric(target_vals))
  }
  
  tab <- table(x)
  max_n <- max(tab)
  winners <- as.numeric(names(tab)[tab == max_n])
  
  if (length(winners) == 1) return(winners)
  
  NA_real_
}

resolve_binary_conflict <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(FALSE)
  length(unique(x)) > 1
}

resolve_categorical_agree <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  if (length(ux) == 1) return(ux)
  NA_character_
}

resolve_categorical_conflict <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) <= 1) return(FALSE)
  length(unique(x)) > 1
}

resolve_numeric_mean <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

collapse_merged_species <- function(df) {
  has_col <- function(nm) nm %in% names(df)
  
  df %>%
    group_by(scientific_name_clean) %>%
    summarise(
      scientific_name_raw = dplyr::first(na.omit(scientific_name_raw)),
      scientific_name_display = dplyr::first(na.omit(scientific_name_display)),
      scientific_name_clean_original = dplyr::first(na.omit(scientific_name_clean_original)),
      
      merged_from_species = paste(sort(unique(scientific_name_clean_original[!is.na(scientific_name_clean_original)])), collapse = "; "),
      n_rows_collapsed = dplyr::n(),
      n_unique_input_species = dplyr::n_distinct(scientific_name_clean_original, na.rm = TRUE),
      
      routine_trichromacy = resolve_binary_with_target_priority(
        routine_trichromacy,
        scientific_name_clean_original,
        dplyr::first(scientific_name_clean)
      ),
      routine_trichromacy_conflict = resolve_binary_conflict(routine_trichromacy),
      
      diet_pct_fruit_eltontraits = resolve_numeric_mean(diet_pct_fruit_eltontraits),
      diet_pct_fruit_n = sum(!is.na(diet_pct_fruit_eltontraits)),
      diet_pct_fruit_min = suppressWarnings(min(diet_pct_fruit_eltontraits, na.rm = TRUE)),
      diet_pct_fruit_max = suppressWarnings(max(diet_pct_fruit_eltontraits, na.rm = TRUE)),
      
      max_longevity_years_pantheria = {
        pan <- max_longevity_years_pantheria
        pan <- pan[!is.na(pan)]
        if (length(pan) == 0) NA_real_ else mean(pan)
      },
      max_longevity_pantheria_n = sum(!is.na(max_longevity_years_pantheria)),
      
      included_in_core_non_nocturnal_build = resolve_categorical_agree(included_in_core_non_nocturnal_build),
      included_in_core_conflict = resolve_categorical_conflict(included_in_core_non_nocturnal_build),
      
      order = resolve_categorical_agree(order),
      order_conflict = resolve_categorical_conflict(order),
      
      foraging_stratum_eltontraits = resolve_categorical_agree(foraging_stratum_eltontraits),
      foraging_stratum_conflict = resolve_categorical_conflict(foraging_stratum_eltontraits),
      
      abs_latitude_mid_dd = if (has_col("abs_latitude_mid_dd")) resolve_numeric_mean(abs_latitude_mid_dd) else NA_real_,
      abs_latitude_mid_pantheria = if (has_col("abs_latitude_mid_pantheria")) resolve_numeric_mean(abs_latitude_mid_pantheria) else NA_real_,
      lat_mid_dd_pantheria = if (has_col("lat_mid_dd_pantheria")) resolve_numeric_mean(lat_mid_dd_pantheria) else NA_real_,
      lat_midrange_dd_pantheria = if (has_col("lat_midrange_dd_pantheria")) resolve_numeric_mean(lat_midrange_dd_pantheria) else NA_real_,
      lat_min_dd_pantheria = if (has_col("lat_min_dd_pantheria")) resolve_numeric_mean(lat_min_dd_pantheria) else NA_real_,
      lat_max_dd_pantheria = if (has_col("lat_max_dd_pantheria")) resolve_numeric_mean(lat_max_dd_pantheria) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      diet_pct_fruit_min = ifelse(is.infinite(diet_pct_fruit_min), NA_real_, diet_pct_fruit_min),
      diet_pct_fruit_max = ifelse(is.infinite(diet_pct_fruit_max), NA_real_, diet_pct_fruit_max)
    )
}

fit_phyloglm_robust <- function(formula, data, phy, model_name = "model") {
  settings_tbl <- tibble::tribble(
    ~method,          ~btol, ~log.alpha.bound,
    "logistic_MPLE",   10,    2.0,
    "logistic_MPLE",    8,    2.0,
    "logistic_MPLE",   10,    1.5,
    "logistic_IG10",   10,    2.0
  )
  
  audit_rows <- vector("list", nrow(settings_tbl))
  fit_obj <- NULL
  last_error <- NULL
  
  for (i in seq_len(nrow(settings_tbl))) {
    cfg <- settings_tbl[i, ]
    
    message(
      "Trying ", model_name,
      " with method=", cfg$method,
      ", btol=", cfg$btol,
      ", log.alpha.bound=", cfg$log.alpha.bound
    )
    
    fit_try <- tryCatch(
      {
        suppressWarnings(
          phylolm::phyloglm(
            formula = formula,
            data = data,
            phy = phy,
            method = cfg$method,
            btol = cfg$btol,
            log.alpha.bound = cfg$log.alpha.bound
          )
        )
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    
    ok <- !is.null(fit_try)
    
    audit_rows[[i]] <- tibble::tibble(
      model_name = model_name,
      attempt = i,
      method = cfg$method,
      btol = cfg$btol,
      log.alpha.bound = cfg$log.alpha.bound,
      success = ok,
      error_message = if (ok) NA_character_ else last_error
    )
    
    if (ok) {
      attr(fit_try, "fit_settings") <- cfg
      fit_obj <- fit_try
      break
    }
  }
  
  audit_tbl <- dplyr::bind_rows(audit_rows)
  
  list(
    fit = fit_obj,
    audit = audit_tbl,
    last_error = last_error
  )
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
    main_df <- main_df %>%
      mutate(scientific_name = stringr::str_trim(paste(genus, species_epithet)))
  } else {
    stop("Need scientific_name OR genus + species_epithet.")
  }
}

anage_lookup_raw <- load_anage_longevity(anage_file)

# ----- Remap/collapse data build for phylogenetic analysis -----
remap_tbl <- load_remap_table(remap_file)
readr::write_csv(remap_tbl, "remapping_cleaned.csv")

main_df <- main_df %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name_display = display_species_name(scientific_name),
    scientific_name_clean_original = clean_species_name(scientific_name),
    scientific_name_clean = apply_merge_map(scientific_name_clean_original, remap_tbl)
  )

merge_remap_audit <- main_df %>%
  select(
    scientific_name_raw,
    scientific_name_display,
    scientific_name_clean_original,
    scientific_name_clean
  ) %>%
  distinct() %>%
  mutate(was_merged = scientific_name_clean_original != scientific_name_clean)

readr::write_csv(merge_remap_audit, "merge_remap_audit.csv")

main_df_collapsed <- collapse_merged_species(main_df)
readr::write_csv(main_df_collapsed, "merged_species_collapse_audit.csv")

anage_lookup <- anage_lookup_raw %>%
  mutate(
    scientific_name_clean_original = clean_species_name(scientific_name),
    scientific_name_clean = apply_merge_map(scientific_name_clean_original, remap_tbl)
  ) %>%
  group_by(scientific_name_clean) %>%
  summarise(anage_longevity_years = resolve_numeric_mean(anage_longevity_years), .groups = "drop")

if (!("max_longevity_years_pantheria" %in% names(main_df_collapsed))) main_df_collapsed$max_longevity_years_pantheria <- NA_real_
if (!("included_in_core_non_nocturnal_build" %in% names(main_df_collapsed))) main_df_collapsed$included_in_core_non_nocturnal_build <- NA
if (!("order" %in% names(main_df_collapsed))) main_df_collapsed$order <- NA_character_

analysis_df <- main_df_collapsed %>%
  mutate(
    order = stringr::str_to_lower(as.character(order)),
    routine_trichromacy = safe_numeric(routine_trichromacy),
    diet_pct_fruit_eltontraits = safe_numeric(diet_pct_fruit_eltontraits),
    max_longevity_years_pantheria = safe_numeric(max_longevity_years_pantheria),
    included_in_core_non_nocturnal_build = as.character(included_in_core_non_nocturnal_build)
  ) %>%
  left_join(anage_lookup, by = "scientific_name_clean") %>%
  mutate(
    longevity_years = dplyr::coalesce(max_longevity_years_pantheria, anage_longevity_years),
    longevity_years = if_else(longevity_years > 0, longevity_years, NA_real_),
    log_longevity = log(longevity_years),
    abs_latitude = get_abs_latitude(cur_data_all()),
    uv = cos(abs_latitude * pi / 180),
    arboreal = classify_canopy_multiplier(foraging_stratum_eltontraits),
    arboreal = dplyr::coalesce(arboreal, 1.0),
    uv_arb = uv * arboreal,
    luv_arb = uv_arb * longevity_years,
    z_frugivory = zscore(diet_pct_fruit_eltontraits),
    z_luv_arb = zscore(luv_arb),
    z_lifetime = zscore(longevity_years)
  ) %>%
  filter(
    order == "primates" | is.na(order),
    tolower(included_in_core_non_nocturnal_build) %in% c("true", "t", "1", "yes", "y") |
      is.na(included_in_core_non_nocturnal_build)
  ) %>%
  filter(!is.na(scientific_name_clean), !is.na(routine_trichromacy), !is.na(z_frugivory), !is.na(uv_arb))

if (nrow(analysis_df) < 10) {
  stop("Too few complete species for phylogenetic analysis after filtering (n = ", nrow(analysis_df), ").")
}

phy_raw <- load_phylo_tree(tree_file)

tree_tip_audit_before <- tibble::tibble(
  tree_tip_raw = phy_raw$tip.label,
  tree_tip_display = display_species_name(phy_raw$tip.label),
  tree_tip_clean = clean_species_name(phy_raw$tip.label)
)
readr::write_csv(tree_tip_audit_before, "tree_tip_audit_before_grafting.csv")

dedup_result <- deduplicate_tree_by_clean_species(phy_raw)
phy_base <- dedup_result$tree

readr::write_csv(dedup_result$duplicate_table, "tree_clean_name_duplicates_before_pruning.csv")
readr::write_csv(dedup_result$keep_table, "tree_tip_keep_table_after_dedup.csv")
readr::write_csv(dedup_result$dropped_raw_tips, "tree_tip_dropped_during_dedup.csv")

tree_tip_audit_dedup <- tibble::tibble(
  tree_tip_raw = phy_base$tip.label,
  tree_tip_display = display_species_name(phy_base$tip.label),
  tree_tip_clean = clean_species_name(phy_base$tip.label)
)
readr::write_csv(tree_tip_audit_dedup, "tree_tip_audit_after_dedup_before_grafting.csv")

graft_tbl <- remap_tbl %>%
  filter(action == "graft") %>%
  filter(scientific_name_clean %in% analysis_df$scientific_name_clean)

graft_result <- apply_graft_plan(phy_base, graft_tbl)
phy_augmented <- graft_result$tree
graft_audit <- graft_result$audit
readr::write_csv(graft_audit, "tree_graft_audit.csv")

failed_grafts <- graft_audit %>% filter(!grafted)
if (nrow(failed_grafts) > 0) {
  message("Some grafts failed. See tree_graft_audit.csv for details.")
}

tree_tip_audit_after <- tibble::tibble(
  tree_tip_raw = phy_augmented$tip.label,
  tree_tip_display = display_species_name(phy_augmented$tip.label),
  tree_tip_clean = clean_species_name(phy_augmented$tip.label)
)
readr::write_csv(tree_tip_audit_after, "tree_tip_audit_after_grafting.csv")

if (anyDuplicated(phy_augmented$tip.label)) {
  dupes <- unique(phy_augmented$tip.label[duplicated(phy_augmented$tip.label)])
  stop(
    "Duplicate tip labels remain after grafting: ",
    paste(dupes, collapse = ", "),
    "\nThis usually means a grafted species name already exists in the deduplicated tree."
  )
}

tree_tip_table <- tibble::tibble(
  tree_tip_raw = phy_augmented$tip.label,
  tree_tip_display = display_species_name(phy_augmented$tip.label),
  tree_tip_clean = clean_species_name(phy_augmented$tip.label)
) %>%
  filter(!is.na(tree_tip_clean)) %>%
  distinct(tree_tip_clean, .keep_all = TRUE)

data_species_table <- analysis_df %>%
  select(
    scientific_name_raw,
    scientific_name_display,
    scientific_name_clean_original,
    scientific_name_clean
  ) %>%
  distinct()

match_table <- data_species_table %>%
  left_join(tree_tip_table, by = c("scientific_name_clean" = "tree_tip_clean")) %>%
  mutate(in_tree = !is.na(tree_tip_raw))

data_only_species <- match_table %>%
  filter(!in_tree) %>%
  transmute(
    scientific_name_raw,
    scientific_name_display,
    scientific_name_clean_original,
    scientific_name_clean
  )

tree_only_species <- tree_tip_table %>%
  anti_join(data_species_table, by = c("tree_tip_clean" = "scientific_name_clean")) %>%
  transmute(
    tree_tip_raw,
    tree_tip_display,
    tree_tip_clean
  )

matched_species <- match_table %>%
  filter(in_tree) %>%
  transmute(
    scientific_name_raw,
    scientific_name_display,
    scientific_name_clean_original,
    scientific_name_clean,
    tree_tip_raw,
    tree_tip_display
  )

readr::write_csv(match_table, "species_name_match_table.csv")
readr::write_csv(data_only_species, "species_in_data_not_in_tree.csv")
readr::write_csv(tree_only_species, "species_in_tree_not_in_data.csv")
readr::write_csv(matched_species, "species_matched_between_data_and_tree.csv")

cat("\n===== TREE / DATA MATCH SUMMARY =====\n")
cat("Species in filtered dataset:", nrow(data_species_table), "\n")
cat("Unique cleaned tips in tree after grafting:", nrow(tree_tip_table), "\n")
cat("Matched species:", nrow(matched_species), "\n")
cat("Species in data but not tree:", nrow(data_only_species), "\n")
cat("Species in tree but not data:", nrow(tree_only_species), "\n")
cat("Merge remaps loaded:", sum(remap_tbl$action == "merge"), "\n")
cat("Graft remaps loaded:", sum(remap_tbl$action == "graft"), "\n")
cat("Successful grafts:", sum(graft_audit$grafted, na.rm = TRUE), "\n")

if (nrow(data_only_species) > 0) {
  cat("\nFirst species in data but not tree:\n")
  print(utils::head(data_only_species, 20))
}

if (nrow(tree_only_species) > 0) {
  cat("\nFirst species in tree but not data:\n")
  print(utils::head(tree_only_species, 20))
}

if (nrow(matched_species) < 10) {
  stop("Too few matched species between dataset and tree (n = ", nrow(matched_species), ").")
}

tree_tips_to_keep_raw <- matched_species$tree_tip_raw

phy_pruned <- ape::drop.tip(
  phy_augmented,
  setdiff(phy_augmented$tip.label, tree_tips_to_keep_raw)
)

phy_pruned$tip.label <- clean_species_name(phy_pruned$tip.label)

if (anyDuplicated(phy_pruned$tip.label)) {
  dupes <- unique(phy_pruned$tip.label[duplicated(phy_pruned$tip.label)])
  stop(
    "Duplicate cleaned tip labels remain after pruning: ",
    paste(dupes, collapse = ", ")
  )
}

analysis_phylo_df <- analysis_df %>%
  filter(scientific_name_clean %in% phy_pruned$tip.label) %>%
  slice(match(phy_pruned$tip.label, scientific_name_clean))

if (!identical(analysis_phylo_df$scientific_name_clean, phy_pruned$tip.label)) {
  stop("Aligned data rows do not exactly match pruned tree tip labels.")
}

analysis_phylo_df <- analysis_phylo_df %>%
  mutate(routine_trichromacy = as.integer(routine_trichromacy))

analysis_phylo_df <- as.data.frame(analysis_phylo_df)
rownames(analysis_phylo_df) <- analysis_phylo_df$scientific_name_clean

if (!all(phy_pruned$tip.label %in% rownames(analysis_phylo_df))) {
  stop("Tree tips do not all match data rownames.")
}
if (!all(rownames(analysis_phylo_df) %in% phy_pruned$tip.label)) {
  stop("Data rownames do not all match tree tips.")
}

cat("\n===== MODEL INPUT DIAGNOSTICS =====\n")
cat("Rows in analysis_phylo_df:", nrow(analysis_phylo_df), "\n")

diag_tbl <- dplyr::as_tibble(analysis_phylo_df, rownames = NA) %>%
  mutate(
    miss_routine = is.na(routine_trichromacy),
    miss_z_frug = is.na(z_frugivory),
    miss_z_luv_arb = is.na(z_luv_arb),
    miss_z_lifetime = is.na(z_lifetime),
    complete_full_model = !miss_routine & !miss_z_frug & !miss_z_luv_arb & !miss_z_lifetime
  )

cat("Complete cases for full model:", sum(diag_tbl$complete_full_model), "\n")
cat("Missing routine_trichromacy:", sum(diag_tbl$miss_routine), "\n")
cat("Missing z_frugivory:", sum(diag_tbl$miss_z_frug), "\n")
cat("Missing z_luv_arb:", sum(diag_tbl$miss_z_luv_arb), "\n")
cat("Missing z_lifetime:", sum(diag_tbl$miss_z_lifetime), "\n")

readr::write_csv(diag_tbl, "model_input_diagnostics.csv")
readr::write_csv(
  diag_tbl %>% filter(!complete_full_model),
  "model_rows_dropped_for_missingness.csv"
)

if (sum(diag_tbl$complete_full_model) < 10) {
  cat("\nFirst rows with missing model inputs:\n")
  print(
    diag_tbl %>%
      filter(!complete_full_model) %>%
      select(
        scientific_name_clean,
        merged_from_species,
        routine_trichromacy,
        z_frugivory,
        z_luv_arb,
        z_lifetime,
        routine_trichromacy_conflict,
        foraging_stratum_conflict
      ) %>%
      head(30)
  )
  stop("Too few complete cases for phylogenetic logistic regression after missing-data filtering.")
}

# Restrict both analyses to the same complete-case species set.
analysis_phylo_df <- analysis_phylo_df[diag_tbl$complete_full_model, , drop = FALSE]
phy_pruned <- ape::keep.tip(phy_pruned, rownames(analysis_phylo_df))
analysis_phylo_df <- analysis_phylo_df[match(phy_pruned$tip.label, rownames(analysis_phylo_df)), , drop = FALSE]

if (!all(analysis_phylo_df$routine_trichromacy %in% c(0L, 1L))) {
  stop("routine_trichromacy must be coded as 0/1 for phylogenetic logistic regression.")
}

# ----- Non-phylogenetic logistic regressions on the same species as phylo -----
analysis_nonphy_df <- tibble::as_tibble(analysis_phylo_df, rownames = NA)

model_full <- glm(routine_trichromacy ~ z_frugivory + z_luv_arb + z_lifetime, data = analysis_nonphy_df, family = binomial())
model_frug <- glm(routine_trichromacy ~ z_frugivory, data = analysis_nonphy_df, family = binomial())
model_luv <- glm(routine_trichromacy ~ z_luv_arb, data = analysis_nonphy_df, family = binomial())
model_lifetime <- glm(routine_trichromacy ~ z_lifetime, data = analysis_nonphy_df, family = binomial())

cat("\n===== NON-PHYLOGENETIC DATA SUMMARY (PHYLO-MATCHED SPECIES) =====\n")
cat("Species analyzed:", nrow(analysis_nonphy_df), "\n")
cat("Routine trichromat species:", sum(analysis_nonphy_df$routine_trichromacy == 1, na.rm = TRUE), "\n")
cat("\n===== FULL BINOMIAL LOGIT =====\n")
print(summary(model_full))
cat("\n===== REDUCED MODELS =====\n")
print(summary(model_frug))
print(summary(model_luv))
print(summary(model_lifetime))
cat("\n===== AIC COMPARISON =====\n")
print(AIC(model_full, model_frug, model_luv, model_lifetime))

coef_tbl_nonphy <- tibble::tibble(
  term = names(coef(model_full)),
  log_odds = as.numeric(coef(model_full)),
  odds_ratio = exp(log_odds)
)

analysis_export_nonphy <- analysis_nonphy_df %>%
  select(
    scientific_name_clean,
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
    z_luv_arb,
    z_lifetime
  )

readr::write_csv(analysis_export_nonphy, "species_level_luv_analysis_dataset.csv")
readr::write_csv(coef_tbl_nonphy, "species_level_luv_full_model_coefficients.csv")

plot_and_save_tree(
  tree = phy_pruned,
  pdf_file = "species_level_luv_pruned_tree.pdf",
  png_file = "species_level_luv_pruned_tree.png"
)

ape::write.tree(phy_pruned, file = "species_level_luv_pruned_tree.nwk")

full_res <- fit_phyloglm_robust(
  routine_trichromacy ~ z_frugivory + z_luv_arb + z_lifetime,
  data = analysis_phylo_df,
  phy = phy_pruned,
  model_name = "full"
)

frug_res <- fit_phyloglm_robust(
  routine_trichromacy ~ z_frugivory,
  data = analysis_phylo_df,
  phy = phy_pruned,
  model_name = "frugivory_only"
)

luv_res <- fit_phyloglm_robust(
  routine_trichromacy ~ z_luv_arb,
  data = analysis_phylo_df,
  phy = phy_pruned,
  model_name = "luv_only"
)

lifetime_res <- fit_phyloglm_robust(
  routine_trichromacy ~ z_lifetime,
  data = analysis_phylo_df,
  phy = phy_pruned,
  model_name = "lifetime_only"
)

fit_audit_tbl <- dplyr::bind_rows(
  full_res$audit,
  frug_res$audit,
  luv_res$audit,
  lifetime_res$audit
)
readr::write_csv(fit_audit_tbl, "phyloglm_fit_audit.csv")

if (is.null(full_res$fit)) {
  stop("Full phylogenetic model failed under all attempted settings. See phyloglm_fit_audit.csv")
}
if (is.null(frug_res$fit)) {
  stop("Frugivory-only phylogenetic model failed under all attempted settings. See phyloglm_fit_audit.csv")
}
if (is.null(luv_res$fit)) {
  stop("LUV-only phylogenetic model failed under all attempted settings. See phyloglm_fit_audit.csv")
}
if (is.null(lifetime_res$fit)) {
  stop("Lifetime-only phylogenetic model failed under all attempted settings. See phyloglm_fit_audit.csv")
}

model_full_phy <- full_res$fit
model_frug_phy <- frug_res$fit
model_luv_phy <- luv_res$fit
model_lifetime_phy <- lifetime_res$fit

cat("\n===== DATA SUMMARY =====\n")
cat("Species analyzed after tree matching:", nrow(analysis_phylo_df), "\n")
cat("Routine trichromat species:", sum(analysis_phylo_df$routine_trichromacy == 1, na.rm = TRUE), "\n")

cat("\n===== FIT SETTINGS USED =====\n")
print(attr(model_full_phy, "fit_settings"))
print(attr(model_frug_phy, "fit_settings"))
print(attr(model_luv_phy, "fit_settings"))
print(attr(model_lifetime_phy, "fit_settings"))

cat("\n===== FULL PHYLOGENETIC LOGISTIC MODEL =====\n")
print(summary(model_full_phy))

cat("\n===== REDUCED PHYLOGENETIC MODELS =====\n")
print(summary(model_frug_phy))
print(summary(model_luv_phy))
print(summary(model_lifetime_phy))

coef_tbl_phy <- tibble::tibble(
  term = names(coef(model_full_phy)),
  log_odds = as.numeric(coef(model_full_phy)),
  odds_ratio = exp(log_odds)
)

cat("\n===== FULL PHYLOGENETIC MODEL COEFFICIENTS =====\n")
print(coef_tbl_phy)

analysis_export <- tibble::as_tibble(analysis_phylo_df) %>%
  select(
    scientific_name_raw,
    scientific_name_display,
    scientific_name_clean_original,
    scientific_name_clean,
    merged_from_species,
    n_rows_collapsed,
    n_unique_input_species,
    routine_trichromacy,
    routine_trichromacy_conflict,
    diet_pct_fruit_eltontraits,
    longevity_years,
    log_longevity,
    abs_latitude,
    uv,
    foraging_stratum_eltontraits,
    foraging_stratum_conflict,
    arboreal,
    luv_arb,
    z_frugivory,
    z_luv_arb,
    z_lifetime
  )

readr::write_csv(analysis_export, "species_level_luv_phylo_analysis_dataset.csv")
readr::write_csv(coef_tbl_phy, "species_level_luv_phylo_full_model_coefficients.csv")

cat("\nSaved outputs:\n")
cat(" - species_level_luv_analysis_dataset.csv\n")
cat(" - species_level_luv_full_model_coefficients.csv\n")
cat(" - species_level_luv_phylo_analysis_dataset.csv\n")
cat(" - species_level_luv_phylo_full_model_coefficients.csv\n")
cat(" - species_level_luv_pruned_tree.nwk\n")
cat(" - species_level_luv_pruned_tree.pdf\n")
cat(" - species_level_luv_pruned_tree.png\n")
cat(" - remapping_cleaned.csv\n")
cat(" - merge_remap_audit.csv\n")
cat(" - merged_species_collapse_audit.csv\n")
cat(" - tree_tip_audit_before_grafting.csv\n")
cat(" - tree_clean_name_duplicates_before_pruning.csv\n")
cat(" - tree_tip_keep_table_after_dedup.csv\n")
cat(" - tree_tip_dropped_during_dedup.csv\n")
cat(" - tree_tip_audit_after_dedup_before_grafting.csv\n")
cat(" - tree_graft_audit.csv\n")
cat(" - tree_tip_audit_after_grafting.csv\n")
cat(" - species_name_match_table.csv\n")
cat(" - species_in_data_not_in_tree.csv\n")
cat(" - species_in_tree_not_in_data.csv\n")
cat(" - species_matched_between_data_and_tree.csv\n")
cat(" - model_input_diagnostics.csv\n")
cat(" - model_rows_dropped_for_missingness.csv\n")
cat(" - phyloglm_fit_audit.csv\n")
