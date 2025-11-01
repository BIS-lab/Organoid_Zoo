library(tidyverse)
library(drc)        # LL.4
library(pracma)     # trapz
library(pheatmap)
library(RColorBrewer)

# ---------------------------
# Core DSS(LL.4) function
# ---------------------------
# Following Yadav et al. (2014):
# 1) Convert viability(%) -> inhibition(%) = 100 - viability
# 2) Fit LL.4 on viability vs. dose and densely sample predictions on the log10 dose axis
# 3) Integrate only the positive part of (inhibition - thr) over the log10 dose range, then normalize by (100 - thr)
# 4) Rescale to 0–50
compute_dss_ll4 <- function(conc, viability,
                            thr = 10,        # activity threshold (%) (commonly 10)
                            grid_n = 200,    # number of points for prediction grid
                            fallback_nonparam = TRUE) {

  stopifnot(length(conc) == length(viability))
  # Keep valid pairs
  ok <- is.finite(conc) & is.finite(viability)
  conc <- conc[ok]; viability <- viability[ok]
  if (length(unique(conc)) < 3 || length(conc) < 4) return(NA_real_)

  # Inhibition (%)
  inhib <- pmin(pmax(100 - viability, 0), 100)

  # LL.4 is typically stable when modeling a decreasing viability curve: y ~ dose
  # (One could also model inhibition as increasing, but viability is more common.)
  df_fit <- tibble(
    dose = conc,
    resp = viability
  )

  # Attempt LL.4 fit
  fit <- tryCatch(
    drm(resp ~ dose, data = df_fit, fct = LL.4(names = c("b","c","d","e"))),
    error = function(e) NULL
  )

  # Fallback: non-parametric DSS approximation if LL.4 fails
  if (is.null(fit)) {
    if (!fallback_nonparam) return(NA_real_)
    x_log10 <- log10(conc)
    xr <- range(x_log10)
    if (!is.finite(diff(xr)) || diff(xr) <= 0) return(NA_real_)
    excess <- pmax((100 - viability) - thr, 0) # (inhib - thr)
    area <- trapz(x_log10[order(x_log10)], excess[order(x_log10)])
    dss01 <- area / (diff(xr) * (100 - thr))
    return(50 * max(min(dss01, 1), 0))
  }

  # Prediction grid within the observed concentration range (log-spaced; natural log)
  minC <- min(conc, na.rm = TRUE)
  maxC <- max(conc, na.rm = TRUE)
  if (!is.finite(minC) || !is.finite(maxC) || minC <= 0 || maxC <= minC) return(NA_real_)
  grid <- exp(seq(log(minC), log(maxC), length.out = grid_n))
  grid_log10 <- log10(grid)

  # Predict viability (%)
  pred_viab <- tryCatch(
    predict(fit, newdata = data.frame(dose = grid)),
    error = function(e) rep(NA_real_, length(grid))
  )

  # If prediction fails, fallback to non-parametric DSS
  if (all(!is.finite(pred_viab))) {
    if (!fallback_nonparam) return(NA_real_)
    x_log10 <- log10(conc)
    xr <- range(x_log10)
    if (!is.finite(diff(xr)) || diff(xr) <= 0) return(NA_real_)
    excess <- pmax((100 - viability) - thr, 0)
    area <- trapz(x_log10[order(x_log10)], excess[order(x_log10)])
    dss01 <- area / (diff(xr) * (100 - thr))
    return(50 * max(min(dss01, 1), 0))
  }

  # Predicted inhibition (%)
  pred_inhib <- pmin(pmax(100 - pred_viab, 0), 100)

  # DSS integration
  excess <- pmax(pred_inhib - thr, 0)
  width  <- diff(range(grid_log10))
  if (!is.finite(width) || width <= 0 || all(excess <= 0)) return(0)
  area <- trapz(grid_log10, excess)  # integrate over log10 concentration axis
  dss01 <- area / (width * (100 - thr))
  dss <- 50 * max(min(dss01, 1), 0)  # scale to 0–50
  as.numeric(dss)
}

# ===========================
# Data loading & preprocessing
# ===========================
cols_div <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)

df <- read_csv("/vol/data/drug_screen/drug_screening_simple_251016.csv", col_types = cols()) %>%
  mutate(
    `Concentration (µM)` = str_replace_all(`Concentration (µM)`, ",", "."),
    Concentration = as.numeric(`Concentration (µM)`),
    Drug = factor(Drug, levels = unique(Drug))
  )

species_cols <- setdiff(names(df), c("Drug", "Concentration (µM)", "Concentration"))

# Aggregate replicates → species-level mean
eps <- 1e-12
df <- df %>% mutate(Conc_log10 = log10(pmax(Concentration, eps)))

long_avg <- df %>%
  pivot_longer(all_of(species_cols), names_to = "Sample", values_to = "Viability") %>%
  mutate(
    # If Sample names have an underscore, the prefix is treated as Species
    Species = if_else(str_detect(Sample, "_"),
                      str_extract(Sample, "^[^_]+"),
                      Sample)
  ) %>%
  filter(!is.na(Viability), !is.na(Concentration), Concentration > 0) %>%
  group_by(Drug, Species, Concentration) %>%
  summarise(Viability = mean(Viability, na.rm = TRUE), .groups = "drop")

# ===========================
# Compute DSS(LL.4)
# ===========================
# For each Drug–Species curve, fit LL.4 and summarize into DSS
dss_df <- long_avg %>%
  arrange(Drug, Species, Concentration) %>%
  group_by(Drug, Species) %>%
  summarise(
    DSS = compute_dss_ll4(Concentration, Viability, thr = 10, grid_n = 200),
    .groups = "drop"
  )

# Pivot to matrix: rows=Drug, cols=Species
dss_mat <- dss_df %>%
  pivot_wider(names_from = Species, values_from = DSS) %>%
  arrange(Drug) %>%
  column_to_rownames("Drug") %>%
  as.matrix()

# Column order (desired species order)
species_order <- c(
  "mouse", "rat", "spiny", "hmouse","guineapig", "human",
  "cynomolgus", "marmoset", "dog", "ferret", "pig", "cattle",
  "goat", "steppe", "serotine", "blackhouse","sakhalin", "chick"
)
species_order <- intersect(species_order, colnames(dss_mat))
dss_mat <- dss_mat[, species_order, drop = FALSE]

# (Optional) Fix drug row order
drugs <- c(
  "Lorlatinib","Cabozantinib","Gilteritinib","Erdafitinib","Afatinib","Osimertinib","Lapatinib","Pazopanib","Axitinib","Linsitinib",
  "Bosutinib","Ibrutinib",
  "LY3009120", "Encorafenib","Sorafenib","Trametinib","Tipifarnib","Lonafarnib",
  "Alpelisib","Everolimus",
  "Pimozide","Ruxolitinib",
  "ICG-001","Enzalutamid",
  "Palbociclib","Abemaciclib","AMG 900","Alisertib","Venetoclax",
  "Rucaparib","Romidepsin","Panobinostat","Onalespib","Ganetespib",
  "Thiomyristoyl","Elesclomol"
)
keep <- intersect(drugs, rownames(dss_mat))
if (length(keep) > 0) dss_mat <- dss_mat[keep, , drop = FALSE]

# ===========================
# Heatmap (DSS: 0–50; higher = more sensitive)
# ===========================
bk <- seq(0, 50, length.out = length(cols_div) + 1)

options(repr.plot.width = 8, repr.plot.height = 8)
p <- pheatmap(
  dss_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = cols_div,
  breaks = bk,
  angle_col = 45,
  main = "DSS (LL.4 fit; species-averaged replicates; 0–50)",
  na_col = "grey95",
  border_color = NA
)
ggsave(
  filename = "/vol/output/2025final/drug_screening_DSS_heatmap.pdf",
  plot = p,
  width = 8,
  height = 8,
  dpi = 300
) # Fig 7b


##### Dose–response curves for selected rodents and drugs
rodents <- c("mouse","rat","spiny","hmouse","guineapig")
chosen_drugs <- c("Panobinostat","Tipifarnib","Osimertinib")

long_raw <- df %>%
  pivot_longer(all_of(species_cols), names_to = "Sample", values_to = "Viability") %>%
  mutate(
    Species = if_else(str_detect(Sample, "_"), str_extract(Sample, "^[^_]+"), Sample),
    Rep     = str_extract(Sample, "(?<=_)\\d+")
  ) %>%
  filter(Species %in% rodents,
         Drug %in% chosen_drugs,
         is.finite(Concentration), Concentration > 0) %>%
  mutate(
    Species = factor(Species, levels = rodents),
    Drug    = factor(Drug,    levels = chosen_drugs)
  )

stats_line <- long_raw %>%
  group_by(Drug, Species, Concentration) %>%
  summarise(
    n    = sum(is.finite(Viability)),
    mean = mean(Viability, na.rm = TRUE),
    sd   = sd(Viability, na.rm = TRUE),
    se   = if_else(n > 1, sd / sqrt(n), NA_real_),
    .groups = "drop"
  )

# custom x축 라벨 (2×5^n)
label_2x5n <- function(x) {
  n <- log(x / 2, 5)
  n_round <- round(n, 1)
  paste0("2×5^", n_round)
}

# plot
options(repr.plot.width = 16, repr.plot.height = 9)
p <- ggplot(stats_line, aes(x = Concentration, y = mean)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 2, color = "blue") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0, color = "blue", linewidth = 0.5, na.rm = TRUE) +
  geom_hline(yintercept = 100, linetype = 3, color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 50,  linetype = 1, color = "grey40", linewidth = 0.5) +
  facet_grid(Drug ~ Species, scales = "free_x", switch = "y") +
  scale_x_log10(
    breaks = function(x) sort(unique(stats_line$Concentration)),
    labels = label_2x5n
  ) +
  coord_cartesian(ylim = c(0, 150)) +
  labs(
    title = "Dose–response curves (mean ± SE)",
    x = "Concentration (µM)",
    y = "Relative viability (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.x = element_text(size = 24, face = "bold"),  # species 이름 (열 제목)
    strip.text.y = element_text(size = 24, face = "bold"),  # 약명 (행 제목)
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(p, filename="/vol/output/2025final/dss_rodent_AUCplot_15_251020.pdf", width = 16, height = 9) # Fig 7e