## run_frugivory_lifetime_light_phylo_updated.R
## Uses updated CSV with AnAge merged in.
##
## Outcome:
##   routine_trichromacy (0/1)
##
## Predictors:
##   z_fruit
##   z_lifetime_light
##
## Lifetime light exposure =
##   log(longevity) * UV * canopy_multiplier
##
## UV = cos(abs(latitude) * pi/180)
## canopy_multiplier = 0.3 for canopy species
## canopy_multiplier = 1.0 for non-canopy species
##
## Longevity source logic:
##   if max_longevity_years_combined exists, use it
##   else use coalesce(max_longevity_years_pantheria, anage_max_longevity_years)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ape)
  library(phylolm)
})

data_file <- "diurnal_primate_trichromacy_lifespan_master_v0_3_anage.csv"
tree_file <- "primate_tree.tre"   # change to your actual tree file

zscore_local <- function(x) {
  x <- as.numeric(x)
  x[x <= -998] <- NA_real_
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
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

dat <- read_csv(data_file, show_col_types = FALSE)

if (all(c("genus", "species_epithet") %in% names(dat))) {
  dat <- dat %>% mutate(scientific_name = str_trim(paste(genus, species_epithet)))
} else if (!("scientific_name" %in% names(dat))) {
  stop("Need either scientific_name or genus + species_epithet.")
}

dat <- dat %>%
  mutate(
    routine_trichromacy = as.numeric(routine_trichromacy),
    max_longevity_years_pantheria = as.numeric(max_longevity_years_pantheria),
    anage_max_longevity_years = as.numeric(anage_max_longevity_years)
  )

if ("max_longevity_years_combined" %in% names(dat)) {
  dat <- dat %>% mutate(longevity_years_analysis = as.numeric(max_longevity_years_combined))
} else {
  dat <- dat %>% mutate(longevity_years_analysis = coalesce(max_longevity_years_pantheria,
                                                            anage_max_longevity_years))
}

dat$abs_latitude <- get_abs_latitude(dat)

dat <- dat %>%
  mutate(
    uv_raw = cos(abs_latitude * pi / 180),
    log_longevity = log(longevity_years_analysis),
    canopy_multiplier = recode_canopy_multiplier(foraging_stratum_eltontraits),
    canopy_class = recode_canopy_class(foraging_stratum_eltontraits),
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

message("Complete-case species for baseline model: ", nrow(D))
message("Species with longevity for analysis: ", sum(!is.na(dat$longevity_years_analysis)))

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

tree <- read.tree(tree_file)
keep <- intersect(tree$tip.label, D$scientific_name)
if (length(keep) < 20) stop("Very few species overlap between tree and data. Check tip labels.")

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
  select(
    scientific_name,
    routine_trichromacy,
    longevity_years_analysis,
    max_longevity_years_pantheria,
    anage_max_longevity_years,
    anage_match_type,
    anage_match_name,
    diet_pct_fruit_eltontraits,
    z_fruit,
    log_longevity,
    abs_latitude,
    uv_raw,
    foraging_stratum_eltontraits,
    canopy_class,
    canopy_multiplier,
    lifetime_light_raw,
    z_lifetime_light
  )

write_csv(analysis_subset, "frugivory_lifetime_light_phylo_analysis_subset_updated.csv")

cat("\nSaved:\n")
cat("  frugivory_lifetime_light_phylo_analysis_subset_updated.csv\n")
cat("  phylogenetic_trichromacy_predictions_updated.png\n")
