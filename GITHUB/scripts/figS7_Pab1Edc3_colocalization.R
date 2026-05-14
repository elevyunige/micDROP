## figS7_Pab1Edc3_colocalization.R
## Reproduces Fig S7 panels B and C:
##   B — % of cells with colocalization (Pab1 & Edc3 vs 10 IDRs)
##   C — Estimate of dilute-phase concentration per cell, binned by BFP
##       intensity (Pab1, Edc3, and WT vs 10 IDRs)
##
## Source: Cleaned_Atar_FigS7code.R
##
## NOTE: This figure uses BFP + RFP channels (not GFP + RFP), so we
## bypass load.figure.data() which applies GFP-specific corrections.
## The RDS is loaded directly.

here.repo.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
source(file.path(here.repo.root, "R", "00_setup.R"))

# ---- Load data directly (no GFP corrections) -----------------------

fig.dir <- file.path(data.root, "FigS7")
rds     <- file.path(fig.dir, "results", "raw_list.RDS")
pdef    <- file.path(fig.dir, "platedef_idr.csv")

if (!file.exists(rds))  stop("Missing raw_list.RDS: ", rds,
                             "\nRun convert_csv_to_rds.R first.")
if (!file.exists(pdef)) stop("Missing platedef_idr.csv: ", pdef)

data.list <- readRDS(rds)
# Accept wrapped or flat RDS
if (is.list(data.list) && length(data.list) == 1 &&
    is.list(data.list[[1]]) && is.data.frame(data.list[[1]][[1]])) {
  data.list <- data.list[[1]]
}
# Drop non-data.frame elements
data.list <- data.list[vapply(data.list, is.data.frame, logical(1))]

platedef <- read.csv(pdef)
colnames(platedef)[tolower(colnames(platedef)) == "well"] <- "WELL"

cat("[figS7] Loaded", length(data.list), "wells\n")

# ---- Settings -------------------------------------------------------

well_labels <- c(
  "01" = "DDX4",   "02" = "DYRK3",  "03" = "ER-alpha", "04" = "FUS",
  "05" = "HNRNPA1","06" = "HSPB8",  "07" = "RBM14",    "08" = "TAF15",
  "09" = "TDP43",  "10" = "UBQ2"
)

protein_colors <- c("Pab1" = "purple", "Edc3" = "turquoise", "WT" = "gray")

# ---- Helpers --------------------------------------------------------

process_data <- function(well_prefixes, protein_name) {
  temp_list <- list()
  for (prefix in well_prefixes) {
    for (num in sprintf("%02d", 1:10)) {
      W <- paste0(prefix, num)
      df <- data.list[[W]]
      if (is.null(df) || nrow(df) == 0) next

      df_filtered <- df %>%
        filter(inRFPnfoci > 0, BFP_int_mean > 0) %>%
        mutate(
          Well             = W,
          WellNum          = num,
          Protein          = protein_name,
          WellLabel        = factor(WellNum, levels = names(well_labels),
                                    labels = well_labels),
          log_BFP_int_mean = log10(BFP_int_mean)
        )
      if (nrow(df_filtered) > 0) temp_list[[length(temp_list) + 1]] <- df_filtered
    }
  }
  bind_rows(temp_list)
}

calculate_coloc_fraction <- function(well_prefixes, protein_name) {
  results <- list()
  for (prefix in well_prefixes) {
    for (num in sprintf("%02d", 1:10)) {
      W <- paste0(prefix, num)
      df <- data.list[[W]]
      if (is.null(df) || nrow(df) == 0) next

      df_filt <- df %>%
        filter(inRFPnfoci > 0,
               f1_inRFP_toRFPmean > 150,
               BFP_int_b0 > 5,
               RFP_int_b5 > 30) %>%
        mutate(ratio = f1_inRFP_toBFPmean / BFP_int_b0)

      if (nrow(df_filt) == 0) next

      results[[length(results) + 1]] <- data.frame(
        Well              = W,
        WellNum           = num,
        WellLabel         = well_labels[num],
        Protein           = protein_name,
        Percent_Above_0.6 = mean(df_filt$ratio > 0.6) * 100,
        N                 = nrow(df_filt)
      )
    }
  }
  bind_rows(results)
}

# ---- Panel B: Colocalization boxplot --------------------------------

# Pab1 replicates for colocalization: B, E, F
# Edc3 replicates for colocalization: I, J, K
summary_pab1 <- calculate_coloc_fraction(c("B", "E", "F"), "Pab1")
summary_edc3 <- calculate_coloc_fraction(c("I", "J", "K"), "Edc3")
# WT: no colocalization expected, but shown for reference
summary_wt   <- calculate_coloc_fraction(c("M", "N"),       "WT")

summary_all <- bind_rows(summary_pab1, summary_edc3, summary_wt) %>%
  mutate(
    Protein   = factor(Protein,   levels = c("Pab1", "Edc3", "WT")),
    WellLabel = factor(WellLabel, levels = well_labels)
  )

p_panel_b <- ggplot(summary_all,
                    aes(x = WellLabel, y = Percent_Above_0.6, fill = Protein)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  facet_wrap(~ Protein, scales = "free_x") +
  geom_hline(yintercept = 0, color = "grey", linetype = "solid", linewidth = 0.3) +
  scale_fill_manual(values = protein_colors) +
  labs(y = "% of Cells with Colocalization", x = "IDR-protein") +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold", size = 12),
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none"
  )

# ---- Panel C: Csat boxplots (BFP-binned) ----------------------------

# Pab1 replicates for Csat: A, B, C, E
# Edc3 replicates for Csat: G, H, I, J, K
# WT: M, N
plot_data <- bind_rows(
  process_data(c("A", "B", "C", "E"), "Pab1"),
  process_data(c("G", "H", "I", "J", "K"), "Edc3"),
  process_data(c("M", "N"), "WT")
) %>%
  mutate(
    Protein   = factor(Protein,   levels = c("Pab1", "Edc3", "WT")),
    WellLabel = factor(WellLabel, levels = well_labels)
  )

# Filter to punctate cells with valid measurements
punctate_box_data <- plot_data %>%
  filter(
    inRFPnfoci > 0,
    f1_inRFP_toRFPmean > 150,
    !is.na(RFP_int_b5), !is.na(BFP_int_b5),
    RFP_int_b5 > 0, BFP_int_b5 > 0
  )

# Bin Pab1/Edc3 by BFP intensity thirds; WT is unbinable
binnable_data <- punctate_box_data %>% filter(Protein %in% c("Pab1", "Edc3"))
wt_data       <- punctate_box_data %>% filter(Protein == "WT") %>%
  mutate(BFP_bin = "WT")

# Equalise n across protein x well to the minimum count
min_count <- binnable_data %>%
  group_by(Protein, WellLabel) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pull(Count) %>% min()

set.seed(123)
binnable_data <- binnable_data %>%
  group_by(Protein, WellLabel) %>%
  slice_sample(n = min_count) %>%
  ungroup()

# Tercile bins by log10(BFP_int_b5) within each protein x well
binned_data <- binnable_data %>%
  mutate(log_BFP = log10(BFP_int_b5)) %>%
  group_by(Protein, WellLabel) %>%
  mutate(
    log_min   = min(log_BFP, na.rm = TRUE),
    log_max   = max(log_BFP, na.rm = TRUE),
    thresh_33 = log_min + 0.33 * (log_max - log_min),
    thresh_67 = log_min + 0.67 * (log_max - log_min),
    BFP_bin   = case_when(
      log_BFP <  thresh_33 ~ "Bottom third",
      log_BFP <  thresh_67 ~ "Middle third",
      TRUE                 ~ "Top third"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-log_BFP, -log_min, -log_max, -thresh_33, -thresh_67)

bin_levels     <- c("Bottom third", "Middle third", "Top third", "WT")
protein_levels <- c("Pab1", "Edc3", "WT")

final_data <- bind_rows(binned_data, wt_data) %>%
  mutate(
    BFP_bin   = factor(BFP_bin,   levels = bin_levels),
    Protein   = factor(Protein,   levels = protein_levels),
    WellLabel = factor(WellLabel, levels = well_labels)
  )

# Cap at 500 cells per group
set.seed(123)
final_data <- final_data %>%
  group_by(Protein, WellLabel, BFP_bin) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 500))) %>%
  ungroup()

# Composite x-axis: Well × Protein × BFP_bin
well_levels <- levels(droplevels(final_data$WellLabel))
x_levels <- c()
for (well in well_levels) {
  for (prot in protein_levels) {
    bins     <- if (prot == "WT") "WT" else bin_levels[1:3]
    x_levels <- c(x_levels, paste(well, prot, bins, sep = "_"))
  }
}

final_data <- final_data %>%
  mutate(
    x_group = paste(WellLabel, Protein, BFP_bin, sep = "_"),
    x_group = factor(x_group, levels = x_levels)
  )

# Annotate groups with < 500 cells
label_data <- final_data %>%
  group_by(x_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count < 500) %>%
  mutate(y_pos = max(final_data$RFP_int_b5, na.rm = TRUE) * 1.05)

p_panel_c <- ggplot(final_data,
                    aes(x = x_group, y = RFP_int_b5, fill = Protein)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_text(data = label_data,
            aes(x = x_group, y = y_pos, label = count),
            inherit.aes = FALSE, size = 3, color = "black") +
  scale_y_log10() +
  scale_fill_manual(values = protein_colors) +
  labs(
    x    = "Well / Protein / BFP Bin",
    y    = "Estimate of dilute phase concentration per cell (arb. units)",
    fill = "Protein"
  ) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 90, vjust = 0.5, size = 6),
    legend.position = "top"
  )

# ---- Write PDFs -----------------------------------------------------

pdf(file.path(out.dir, "figS7B_Pab1Edc3_colocalization.pdf"),
    width = 14, height = 6)
print(p_panel_b)
dev.off()

pdf(file.path(out.dir, "figS7C_Pab1Edc3_Csat_boxplots.pdf"),
    width = 20, height = 8)
print(p_panel_c)
dev.off()

cat("Wrote:", file.path(out.dir, "figS7B_Pab1Edc3_colocalization.pdf"), "\n")
cat("Wrote:", file.path(out.dir, "figS7C_Pab1Edc3_Csat_boxplots.pdf"), "\n")
