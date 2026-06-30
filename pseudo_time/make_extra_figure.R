## Nature-style compound supplementary figure (2 x 2 macro)
## Claim: heavy/light peptide abundances diverge in a coordinated, biologically
##        interpretable way along the Basal -> Secretory principal-curve pseudotime.
## Each macro cell is itself a stacked sub-composition.
## Layout: 183 mm x 205 mm.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(svglite)
  library(ragg)
})

# ---- Paths --------------------------------------------------------------
dat_dir    <- "/Users/christina/Desktop/ProteinResearch/Demo/dat"
sup_dir    <- "/Users/christina/Desktop/ProteinResearch/Demo/pseudo_time/supplementary_figure"
out_dir    <- file.path(sup_dir, "figure_output")
cache_path <- file.path(sup_dir, "extra_pipeline_cache.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================================
# Slingshot pipeline (cached: delete extra_pipeline_cache.rds to force rerun)
# =========================================================================
if (file.exists(cache_path)) {
  cache <- readRDS(cache_path)
} else {
  suppressPackageStartupMessages({
    library(QuantQC)
    library(slingshot)
  })
  meta <- read.csv(file.path(
    dat_dir, "04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv"
  ), row.names = 1)
  meta_sub <- meta |> filter(Cell_Type %in% c("Basal", "Secratory"))

  prot <- as.matrix(read.csv(file.path(
    dat_dir, "04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv"
  ), row.names = 1))
  prot_sub <- prot[, intersect(colnames(prot), meta_sub$ID)]

  min_cells <- 200
  prot_filt <- prot_sub[rowSums(!is.na(prot_sub)) >= min_cells, ]
  prot_filt <- QuantQC::Normalize_reference_vector_log(prot_filt)
  prot_imp  <- t(apply(prot_filt, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE); x
  }))

  cor_mat   <- cor(prot_imp, use = "pairwise.complete.obs")
  pca_res   <- eigen(cor_mat)
  pc_scores <- pca_res$vectors
  rownames(pc_scores) <- colnames(prot_imp)

  pc_df <- as.data.frame(pc_scores[, 1:5])
  colnames(pc_df) <- paste0("PC", 1:5)
  pc_df$ID        <- rownames(pc_df)
  pc_df$Cell_Type <- meta_sub$Cell_Type[match(pc_df$ID, meta_sub$ID)]

  pc_mat <- as.matrix(pc_df[, paste0("PC", 1:5)])
  rownames(pc_mat) <- pc_df$ID
  sce <- slingshot(pc_mat, clusterLabels = factor(pc_df$Cell_Type),
                   start.clus = "Basal", end.clus = "Secratory", df = 3)
  pc_df$pseudotime <- slingPseudotime(sce)[, 1]

  krt5_uniprot <- "Q922U2"
  agr2_uniprot <- "O88312"
  df_marker <- data.frame(
    ID   = meta_sub$ID,
    Krt5 = as.numeric(prot_sub[krt5_uniprot, meta_sub$ID]),
    Arg2 = as.numeric(prot_sub[agr2_uniprot, meta_sub$ID])
  )
  marker_pt <- pc_df |>
    select(ID, pseudotime, Cell_Type) |>
    left_join(df_marker, by = "ID") |>
    pivot_longer(c(Krt5, Arg2), names_to = "protein", values_to = "abundance")
  rsq_df <- marker_pt |>
    group_by(protein) |>
    summarise(rsq = summary(
      lm(abundance ~ pseudotime, na.action = na.omit)
    )$r.squared, .groups = "drop") |>
    mutate(label = sprintf("italic(R)^2 == %.3f", rsq))

  curves   <- slingCurves(sce)
  curve_df <- as.data.frame(curves[[1]]$s)
  colnames(curve_df)[1:2] <- c("PC1", "PC2")
  # slingshot returns curve points already ordered along the trajectory

  cache <- list(pc_df = pc_df, curve_df = curve_df,
                marker_pt = marker_pt, rsq_df = rsq_df)
  saveRDS(cache, cache_path)
}

pc_df     <- cache$pc_df
curve_df  <- cache$curve_df
marker_pt <- cache$marker_pt
rsq_df    <- cache$rsq_df

# ---- Palette & theme ----------------------------------------------------
chan_pal     <- c(heavy = "#C0362E", light = "#3B6FAD")
celltype_pal <- c(Basal = "#B50202", Secratory = "#B606C4")
corr_col     <- "#7F77DD"
heat_low <- "#3B6FAD"; heat_mid <- "white"; heat_high <- "#C0362E"

theme_nat <- theme_classic(base_size = 6.5, base_family = "Helvetica") +
  theme(
    axis.line         = element_line(linewidth = 0.32, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.32, colour = "black"),
    axis.text         = element_text(colour = "black", size = 5.8),
    axis.title        = element_text(colour = "black", size = 6.2),
    legend.title      = element_text(size = 5.8),
    legend.text       = element_text(size = 5.5),
    legend.key.size   = unit(2.5, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-4, 0, -2, 0),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 6.2, face = "bold"),
    plot.title        = element_text(size = 6.8, face = "bold", hjust = 0),
    plot.subtitle     = element_text(size = 6, face = "plain", hjust = 0,
                                     colour = "grey25"),
    plot.tag          = element_text(size = 8.5, face = "bold"),
    panel.grid        = element_blank(),
    plot.margin       = margin(2, 4, 2, 2)
  )
theme_set(theme_nat)

# =========================================================================
# Panel A: PCA + principal curve
# =========================================================================
panel_A <- ggplot() +
  geom_point(data = pc_df,
             aes(PC1, PC2, colour = Cell_Type),
             size = 0.5, alpha = 0.65, stroke = 0) +
  geom_path(data = curve_df, aes(PC1, PC2),
            colour = "black", linewidth = 0.55) +
  scale_colour_manual(values = celltype_pal,
                      labels = c(Basal = "Basal", Secratory = "Secretory")) +
  labs(x = "PC1", y = "PC2", colour = NULL,
       title = "PCA with Basal -> Secretory trajectory") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme(legend.position = "top",
        legend.direction = "horizontal")

# =========================================================================
# Panel B: Krt5 + Arg2 marker fits along principal-curve pseudotime
# =========================================================================
panel_B <- ggplot(marker_pt,
                  aes(pseudotime, abundance, colour = Cell_Type)) +
  geom_point(size = 0.35, alpha = 0.4, stroke = 0) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              colour = "black", linewidth = 0.55) +
  geom_text(data = rsq_df,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, parse = TRUE,
            hjust = 1.1, vjust = 1.7, size = 2.1, colour = "grey25") +
  facet_wrap(~ protein, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = celltype_pal) +
  labs(x = "Pseudotime", y = "Protein abundance", colour = NULL,
       title = "Marker proteins (Krt5, Arg2) along pseudotime") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme(legend.position = "top",
        legend.direction = "horizontal")

# =========================================================================
# Panel C: theoretical (top, half height) / observed (bottom)
# =========================================================================
df_target_markers <- read.csv(file.path(sup_dir, "df_target_markers.csv"))
dfB_obs <- df_target_markers |>
  filter(!is.na(centered_val), protein %in% c("Arg2", "Krt5")) |>
  mutate(channel = factor(channel, levels = c("heavy", "light")))

x_lab_B <- "Principal-curve fitted pseudotime"
y_lab_B <- "Centered log abundance"

set.seed(7)
emp_pt_range <- range(dfB_obs$pseudotime, na.rm = TRUE)
emp_ranges <- dfB_obs |>
  group_by(protein) |>
  summarise(ymin = quantile(centered_val, 0.05, na.rm = TRUE),
            ymax = quantile(centered_val, 0.95, na.rm = TRUE),
            .groups = "drop")

n_sim    <- 200
sim_pt   <- seq(emp_pt_range[1], emp_pt_range[2], length.out = n_sim)
pt_span  <- diff(emp_pt_range)
k_scaled <- 14 / pt_span
sigmoid  <- function(x, x0, k) 1 / (1 + exp(-k * (x - x0)))

wavy <- function(x, amp_frac, span_y, n_waves = 3, seed) {
  set.seed(seed)
  span_x <- diff(range(x))
  freqs  <- runif(n_waves, 1.2, 3.0) / span_x
  phases <- runif(n_waves, 0, 2 * pi)
  amps   <- runif(n_waves, 0.5, 1.0)
  amps   <- amps / sum(amps) * amp_frac * span_y
  rowSums(sapply(seq_len(n_waves), function(i)
    amps[i] * sin(2 * pi * freqs[i] * x + phases[i])
  ))
}

make_sim <- function(name, ymin, ymax, x0_lead_f, x0_lag_f,
                     rising, lead, sh, sl) {
  base_lead <- sigmoid(sim_pt, emp_pt_range[1] + pt_span * x0_lead_f, k_scaled)
  base_lag  <- sigmoid(sim_pt, emp_pt_range[1] + pt_span * x0_lag_f,  k_scaled)
  if (!rising) { base_lead <- 1 - base_lead; base_lag <- 1 - base_lag }
  span <- ymax - ymin
  if (lead == "heavy") { h <- base_lead; l <- base_lag }
  else                 { h <- base_lag;  l <- base_lead }
  tibble(
    pseudotime = sim_pt,
    heavy = ymin + span * h + wavy(sim_pt, 0.05, span, seed = sh),
    light = ymin + span * l + wavy(sim_pt, 0.05, span, seed = sl),
    protein = name
  )
}
arg2_row <- emp_ranges |> filter(protein == "Arg2")
krt5_row <- emp_ranges |> filter(protein == "Krt5")
sim_arg2 <- make_sim("Arg2", arg2_row$ymin, arg2_row$ymax,
                     0.40, 0.55, TRUE,  "heavy", 11, 12)
sim_krt5 <- make_sim("Krt5", krt5_row$ymin, krt5_row$ymax,
                     0.40, 0.55, FALSE, "light", 21, 22)
dfB_theo <- bind_rows(sim_arg2, sim_krt5) |>
  pivot_longer(c(heavy, light), names_to = "channel", values_to = "value") |>
  mutate(channel = factor(channel, levels = c("heavy", "light")))

pC_theo <- ggplot(dfB_theo, aes(pseudotime, value, colour = channel)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_line(linewidth = 0.55) +
  facet_wrap(~ protein, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = chan_pal) +
  labs(x = NULL, y = y_lab_B, colour = NULL,
       title = "Theoretical model",
       subtitle = "heavy leads in synthesis; light leads in degradation") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pC_obs <- ggplot(dfB_obs, aes(pseudotime, centered_val, colour = channel)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_point(size = 0.3, alpha = 0.15, stroke = 0) +
  geom_smooth(method = "loess", span = 0.7, se = FALSE, linewidth = 0.7) +
  facet_wrap(~ protein, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = chan_pal) +
  labs(x = x_lab_B, y = y_lab_B, colour = NULL,
       title = "Observed heavy/light abundance") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.4))) +
  theme(legend.position = "top",
        legend.direction = "horizontal")

panel_C <- pC_theo / pC_obs + plot_layout(heights = c(1, 1))

# =========================================================================
# Panel D: heavy heatmap (top) / light heatmap (bottom, same row order)
# =========================================================================
heavy_z <- as.matrix(read.csv(file.path(sup_dir, "heavy_z.csv"), row.names = 1))
light_z <- as.matrix(read.csv(file.path(sup_dir, "light_z.csv"), row.names = 1))
rownames(heavy_z) <- make.unique(rownames(heavy_z))
rownames(light_z) <- make.unique(rownames(light_z))

# Row order: hclust on heavy (NA -> 0 for distance computation only)
heavy_imp <- heavy_z
heavy_imp[is.na(heavy_imp)] <- 0
hc <- hclust(dist(heavy_imp), method = "complete")
ordered_rows <- rownames(heavy_z)[hc$order]
common_cols <- intersect(colnames(heavy_z), colnames(light_z))

bin_levels <- common_cols  # bin1 ... bin30 in file order
to_long <- function(mat) {
  as.data.frame(mat[ordered_rows, common_cols, drop = FALSE]) |>
    rownames_to_column("protein") |>
    pivot_longer(-protein, names_to = "bin", values_to = "z") |>
    mutate(protein = factor(protein, levels = rev(ordered_rows)),
           bin = factor(bin, levels = bin_levels),
           bin_idx = as.integer(bin))
}
heavy_long <- to_long(heavy_z)
light_long <- to_long(light_z)

z_lim <- c(-2.5, 2.5)
heat_scale <- scale_fill_gradient2(
  low = heat_low, mid = heat_mid, high = heat_high,
  midpoint = 0, limits = z_lim, oob = scales::squish,
  name = "z-score",
  breaks = c(-2, 0, 2),
  guide = guide_colourbar(barheight = unit(10, "mm"),
                          barwidth = unit(1.6, "mm"),
                          ticks.colour = "black",
                          frame.colour = "black",
                          frame.linewidth = 0.3)
)

heat_theme <- theme(
  axis.line   = element_blank(),
  axis.ticks  = element_blank(),
  axis.text.y = element_text(size = 4.6, face = "italic", colour = "black"),
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.32)
)

pD_heavy <- ggplot(heavy_long, aes(bin_idx, protein, fill = z)) +
  geom_raster() +
  heat_scale +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(1, round(length(bin_levels) / 2),
                                length(bin_levels))) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Heavy abundance (z-score)") +
  heat_theme +
  theme(axis.text.x = element_blank())

pD_light <- ggplot(light_long, aes(bin_idx, protein, fill = z)) +
  geom_raster() +
  heat_scale +
  scale_x_continuous(
    expand = c(0, 0),
    name = "Pseudotime bin (Basal -> Secretory)",
    breaks = c(1, round(length(bin_levels) / 2), length(bin_levels))
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = NULL, title = "Light abundance (z-score)") +
  heat_theme

panel_D <- (pD_heavy / pD_light) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

# =========================================================================
# Panel E: static reference window correlations
#   top = center, bottom = mid_sec
# =========================================================================
all_cor_abs <- read.csv(file.path(sup_dir, "all_cor.csv"))
xlim_E <- range(all_cor_abs$window_position)
ylim_E <- range(all_cor_abs$correlation)

make_E <- function(ref_label, title_text, show_x = FALSE) {
  d <- all_cor_abs |> filter(reference == ref_label)
  vx <- unique(d$ref_position)[1]
  p <- ggplot(d, aes(window_position, correlation)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
    geom_vline(xintercept = vx, linetype = "dashed",
               linewidth = 0.32, colour = "grey40") +
    geom_line(colour = corr_col, linewidth = 0.7) +
    coord_cartesian(xlim = xlim_E, ylim = ylim_E) +
    labs(y = "Spearman r", title = title_text) +
    theme(plot.title = element_text(size = 6.4, face = "plain"))
  if (show_x) {
    p + labs(x = "Window position (cell index)")
  } else {
    p + labs(x = NULL) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
}

panel_E <- (make_E("center",  "Reference window: center") /
            make_E("mid_sec", "Reference window: mid_sec", show_x = TRUE)) +
  plot_layout(heights = c(1, 1))

# =========================================================================
# Panel F: sliding-window heavy - light differential (single line)
# =========================================================================
ts_dir            <- "/Users/christina/Desktop/ProteinResearch/Demo/pseudo_time/time_series/5day_old"
ordered_cells     <- read.csv(file.path(sup_dir, "ordered_cells.csv"))
window3_light_cor <- read.csv(file.path(ts_dir,  "window3_light_sliding.csv"))
window3_heavy_cor <- read.csv(file.path(ts_dir,  "window3_heavy_sliding.csv"))
n_cells <- nrow(ordered_cells)
refs <- c(early = 0.10, mid_basal = 0.25, center = 0.50,
          mid_sec = 0.75, late = 0.90) * n_cells

dfF <- window3_heavy_cor |>
  select(ref_position, heavy_diff = differential) |>
  inner_join(window3_light_cor |>
               select(ref_position, light_diff = differential),
             by = "ref_position") |>
  mutate(hml = heavy_diff - light_diff)

ref_df <- tibble(ref_position = unname(refs), label = names(refs))
y_top <- max(dfF$hml, na.rm = TRUE)
y_bot <- min(dfF$hml, na.rm = TRUE)
y_pad <- (y_top - y_bot) * 0.12

panel_F <- ggplot(dfF, aes(ref_position, hml)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_vline(data = ref_df, aes(xintercept = ref_position),
             linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  geom_line(linewidth = 0.7, colour = "#2C3E50") +
  geom_text(data = ref_df,
            aes(x = ref_position, y = y_top + y_pad, label = label),
            inherit.aes = FALSE,
            size = 1.7, colour = "grey25", hjust = 0.5) +
  coord_cartesian(ylim = c(y_bot, y_top + y_pad * 1.4), clip = "off") +
  labs(x = "Reference position (cell index)",
       y = "Heavy - Light differential",
       title = "Sliding window (3 windows) heavy - light differential")

# =========================================================================
# Assemble 2 rows x 3 cols -- wrap_elements freezes each compound cell so
# the outer tag sequence sees 6 cells (a-f), not the inner subplots.
# Row 1: A (PCA)         | B (markers)    | C (theo/obs stacked)
# Row 2: D (heatmap stk) | E (corr stk)   | F (sliding diff)
# =========================================================================
A <- wrap_elements(full = panel_A)
B <- wrap_elements(full = panel_B)
C <- wrap_elements(full = panel_C)
D <- wrap_elements(full = panel_D)
E <- wrap_elements(full = panel_E)
F <- wrap_elements(full = panel_F)

fig <- (A | B | C) / (D | E | F) +
  plot_layout(heights = c(0.65, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.5, face = "bold"))

# ---- Save ---------------------------------------------------------------
w_mm <- 183; h_mm <- 150
w_in <- w_mm / 25.4; h_in <- h_mm / 25.4

base <- file.path(out_dir, "fig_extra")

svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in)
print(fig); invisible(dev.off())

grDevices::pdf(paste0(base, ".pdf"), width = w_in, height = h_in,
               useDingbats = FALSE, family = "Helvetica")
print(fig); invisible(dev.off())

ragg::agg_png(paste0(base, ".png"),
              width = w_in, height = h_in, units = "in", res = 600)
print(fig); invisible(dev.off())

ragg::agg_tiff(paste0(base, ".tiff"),
               width = w_in, height = h_in, units = "in", res = 600,
               compression = "lzw")
print(fig); invisible(dev.off())

cat("Saved to:", out_dir, "\n")
