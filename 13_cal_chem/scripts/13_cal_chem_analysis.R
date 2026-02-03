#!/usr/bin/env Rscript
# =============================================================================
# 13_cal_chem_analysis.R
# Carbonate chemistry and calcification analysis for M. capitata RNA-seq project
# Author: David Armstrong / Claude
# Date: 2025-02-03
# =============================================================================

# Load required libraries
library(tidyverse)
library(seacarb)
library(gt)
library(ggplot2)
library(car)        # for Type II/III ANOVA

# Set working directory paths
input_dir <- "input"
fig_dir <- "figures"
output_dir <- "output"

# Create output directories if they don't exist
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

cat("Loading data...\n")

# Load seawater chemistry data
sw_chem <- read_csv(file.path(input_dir, "sw_chem.csv"), show_col_types = FALSE)

# Load calcification data
calc_data <- read_csv(file.path(input_dir, "calc_data.csv"), show_col_types = FALSE)

# Add treatment and season labels to ALL data
sw_chem_labeled <- sw_chem %>%
  mutate(
    # Convert treatment codes: B = OA, D = Ambient
    treatment_label = case_when(
      treatment == "B" ~ "OA",
      treatment == "D" ~ "Amb.",
      TRUE ~ treatment
    ),
    # Convert season codes
    season_label = case_when(
      season == "S" ~ "Summer",
      season == "W" ~ "Winter",
      TRUE ~ season
    )
  )

# Subset for seacarb calculations (need pH_total, TA, temp, sal)
sw_chem_for_seacarb <- sw_chem_labeled %>%
  filter(!is.na(pH_total) & !is.na(TA_seacarb) & !is.na(temp) & !is.na(sal_cor)) %>%
  mutate(
    # Convert TA from µmol/kg to mol/kg for seacarb
    TA_mol_kg = TA_seacarb / 1e6
  )

cat("Data loaded successfully.\n")
cat(sprintf("  - SW chemistry records: %d total\n", nrow(sw_chem)))
cat(sprintf("  - Records with complete data for seacarb: %d\n", nrow(sw_chem_for_seacarb)))
cat(sprintf("  - Calcification records: %d\n", nrow(calc_data)))

# =============================================================================
# 2. CALCULATE CARBONATE CHEMISTRY USING SEACARB
# =============================================================================

cat("\nCalculating carbonate chemistry with seacarb...\n")

# Apply seacarb carb() function to each row with complete data
carb_results <- sw_chem_for_seacarb %>%
  rowwise() %>%
  mutate(
    carb_calc = list(
      tryCatch(
        carb(
          flag = 8,
          var1 = pH_total,
          var2 = TA_mol_kg,
          S = sal_cor,
          T = temp,
          P = 0,           # Standard pressure (surface)
          Pt = 0,          # Phosphate (assume 0 if not measured)
          Sit = 0,         # Silicate (assume 0 if not measured)
          k1k2 = "l",      # Lueker et al. 2000
          kf = "pf",       # Perez and Fraga 1987
          ks = "d",        # Dickson 1990
          pHscale = "T"    # Total scale
        ),
        error = function(e) NULL
      )
    )
  ) %>%
  ungroup()

# Extract relevant parameters from seacarb output
carb_extracted <- carb_results %>%
  mutate(
    pCO2 = map_dbl(carb_calc, ~ if(!is.null(.x)) .x$pCO2 else NA_real_),
    OmegaAragonite = map_dbl(carb_calc, ~ if(!is.null(.x)) .x$OmegaAragonite else NA_real_)
  ) %>%
  select(-carb_calc)

cat("Seacarb calculations complete.\n")

# =============================================================================
# 3. CREATE SUMMARY TABLE OF CARBONATE CHEMISTRY
# =============================================================================

cat("\nGenerating summary statistics table...\n")

# Helper function for mean ± SD formatting (no n)
format_mean_sd <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(NA_character_)
  m <- mean(x)
  s <- sd(x)
  sprintf("%.*f ± %.*f", digits, m, digits, s)
}

# Function to count n
count_n <- function(x) {
  sum(!is.na(x))
}

# --- Calculate measured parameters from ALL available data ---
# Temperature (use all rows with temp data)
temp_summary <- sw_chem_labeled %>%
  filter(!is.na(temp)) %>%
  group_by(season_label, mesocosm, treatment_label) %>%
  summarise(Temp = format_mean_sd(temp, 1), n_temp = count_n(temp), .groups = "drop")

# Salinity (use all rows with sal_cor data)
sal_summary <- sw_chem_labeled %>%
  filter(!is.na(sal_cor)) %>%
  group_by(season_label, mesocosm, treatment_label) %>%
  summarise(Salinity = format_mean_sd(sal_cor, 2), n_sal = count_n(sal_cor), .groups = "drop")

# pH(T) (use all rows with pH_total data)
ph_summary <- sw_chem_labeled %>%
  filter(!is.na(pH_total)) %>%
  group_by(season_label, mesocosm, treatment_label) %>%
  summarise(`pH(T)` = format_mean_sd(pH_total, 2), n_ph = count_n(pH_total), .groups = "drop")

# --- Calculate AT and derived parameters from seacarb subset ---
at_summary <- carb_extracted %>%
  group_by(season_label, mesocosm, treatment_label) %>%
  summarise(
    AT = format_mean_sd(TA_seacarb, 1),
    pCO2 = format_mean_sd(pCO2, 0),
    OmegaArag = format_mean_sd(OmegaAragonite, 2),
    n_at = count_n(TA_seacarb),
    .groups = "drop"
  )

# Merge all summaries
chem_summary <- temp_summary %>%
  left_join(sal_summary, by = c("season_label", "mesocosm", "treatment_label")) %>%
  left_join(ph_summary, by = c("season_label", "mesocosm", "treatment_label")) %>%
  left_join(at_summary, by = c("season_label", "mesocosm", "treatment_label")) %>%
  rename(
    Season = season_label,
    Sample = mesocosm,
    Treatment = treatment_label
  ) %>%
  # Reorder columns as requested
  select(Season, Treatment, Sample, Temp, Salinity, `pH(T)`, AT, pCO2, OmegaArag, 
         n_temp, n_sal, n_ph, n_at) %>%
  # Sort: Season (Summer first), then Treatment (Amb. first), then Sample
  mutate(
    season_order = ifelse(Season == "Summer", 1, 2),
    treatment_order = ifelse(Treatment == "Amb.", 1, 2)
  ) %>%
  arrange(season_order, treatment_order, Sample) %>%
  select(-season_order, -treatment_order)

# Calculate n per season for footnote
n_per_season <- chem_summary %>%
  group_by(Season) %>%
  summarise(
    n_temp = sum(n_temp),
    n_sal = sum(n_sal),
    n_ph = sum(n_ph),
    n_at = sum(n_at),
    .groups = "drop"
  )

# Create footnote text
footnote_text <- sprintf(
"Carbonate chemistry calculated using seacarb (flag=8: pH and TA). OA = Ocean Acidification treatment. Summer n: Temp=%d, Sal=%d, pH(T)=%d, AT/pCO2/Ωarag=%d. Winter n: Temp=%d, Sal=%d, pH(T)=%d, AT/pCO2/Ωarag=%d.",
  n_per_season$n_temp[n_per_season$Season == "Summer"],
  n_per_season$n_sal[n_per_season$Season == "Summer"],
  n_per_season$n_ph[n_per_season$Season == "Summer"],
  n_per_season$n_at[n_per_season$Season == "Summer"],
  n_per_season$n_temp[n_per_season$Season == "Winter"],
  n_per_season$n_sal[n_per_season$Season == "Winter"],
  n_per_season$n_ph[n_per_season$Season == "Winter"],
  n_per_season$n_at[n_per_season$Season == "Winter"]
)

# Remove n columns before making table
chem_summary_table <- chem_summary %>%
  select(-n_temp, -n_sal, -n_ph, -n_at)

# Create GT table (clean format matching user's style)
chem_gt <- chem_summary_table %>%
  gt() %>%
  cols_label(
    Season = "Season",
    Treatment = "Treatment",
    Sample = "Sample",
    Temp = "Temp (°C)",
    Salinity = "Salinity",
    `pH(T)` = "pH(T)",
    AT = html("AT (µmol kg<sup>−1</sup>)"),
    pCO2 = html("pCO<sub>2</sub> (µatm)"),
    OmegaArag = "Ωarag"
  ) %>%
  tab_footnote(
    footnote = footnote_text,
    locations = cells_column_labels(columns = Season)
  ) %>%
  tab_options(
    table.font.size = px(11),
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.top.width = px(2),
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(2),
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = px(1)
  )

# Save GT table to Word document
gtsave(chem_gt, file.path(output_dir, "carbonate_chemistry_table.docx"))

# Also save as CSV for reference
write_csv(chem_summary_table, file.path(output_dir, "carbonate_chemistry_summary.csv"))

cat("Carbonate chemistry table saved to output/carbonate_chemistry_table.docx\n")

# =============================================================================
# 4. ANOVA ON pH(T)
# =============================================================================

cat("\nRunning two-way ANOVA on pH(T)...\n")

# Prepare data for ANOVA (use ALL pH data)
ph_anova_data <- sw_chem_labeled %>%
  select(season_label, treatment_label, pH_total) %>%
  filter(!is.na(pH_total)) %>%
  mutate(
    season = factor(season_label, levels = c("Summer", "Winter")),
    treatment = factor(treatment_label, levels = c("Amb.", "OA"))
  )

# Two-way ANOVA (additive model: no interaction)
ph_anova <- aov(pH_total ~ treatment + season, data = ph_anova_data)

cat("\n--- ANOVA Results: pH(T) ~ Treatment + Season ---\n")
print(summary(ph_anova))

# Type II ANOVA for unbalanced designs
cat("\n--- Type II ANOVA (car::Anova) ---\n")
ph_anova_typeII <- Anova(ph_anova, type = "II")
print(ph_anova_typeII)

# Save ANOVA results
sink(file.path(output_dir, "pH_anova_results.txt"))
cat("Two-way ANOVA: pH(T) ~ Treatment + Season (Additive Model)\n")
cat("=" , rep("=", 60), "\n\n", sep = "")
cat("Type I ANOVA (Sequential):\n")
print(summary(ph_anova))
cat("\nType II ANOVA (Marginal):\n")
print(ph_anova_typeII)
cat("\n\nGroup Means:\n")
print(ph_anova_data %>%
        group_by(treatment, season) %>%
        summarise(
          mean_pH = mean(pH_total),
          sd_pH = sd(pH_total),
          n = n(),
          .groups = "drop"
        ))
sink()

cat("pH ANOVA results saved to output/pH_anova_results.txt\n")

# =============================================================================
# 5. PROCESS CALCIFICATION DATA
# =============================================================================

cat("\nProcessing calcification data...\n")

# Parse sample names to extract treatment, season, mesocosm
calc_processed <- calc_data %>%
  mutate(
    # Extract components from sample name (e.g., "1BS" = mesocosm 1, treatment B, season S)
    mesocosm = as.numeric(str_extract(sample, "^[0-9]+")),
    treatment_code = str_extract(sample, "[BD]"),
    season_code = str_extract(sample, "[SW]$"),
    # Convert to labels
    treatment_label = case_when(
      treatment_code == "B" ~ "OA",
      treatment_code == "D" ~ "Ambient"
    ),
    season_label = case_when(
      season_code == "S" ~ "Summer",
      season_code == "W" ~ "Winter"
    ),
    # Calculate mass-normalized calcification
    # g_day is in g CaCO3 d-1, tp1 is initial weight in g
    # g_day / tp1 = g CaCO3 g-1 d-1 
    calc_normalized = g_day / tp1
  )

# Verify parsing
cat("Parsed calcification data:\n")
print(calc_processed %>% select(sample, mesocosm, treatment_label, season_label, calc_normalized))

# Calculate means per treatment per season
calc_summary <- calc_processed %>%
  group_by(season_label, treatment_label) %>%
  summarise(
    mean_calc = mean(calc_normalized, na.rm = TRUE),
    se_calc = sd(calc_normalized, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    season = factor(season_label, levels = c("Summer", "Winter")),
    treatment = factor(treatment_label, levels = c("Ambient", "OA"))
  )

cat("\nCalcification summary:\n")
print(calc_summary)

# =============================================================================
# 6. CREATE CALCIFICATION BAR PLOT
# =============================================================================

cat("\nCreating calcification bar plot...\n")

# Define colors (user's preferred colors)
treatment_colors <- c("Ambient" = "#343F3E", "OA" = "#8F91A2")

# Create bar plot
calc_plot <- ggplot(calc_summary, aes(x = season, y = mean_calc, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_calc - se_calc, ymax = mean_calc + se_calc),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = treatment_colors, name = "Treatment") +
  labs(
    x = NULL,
    y = expression("Calcification (g CaCO"[3]*" g"^{-1}*" d"^{-1}*")"),
    title = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Save plot
ggsave(file.path(fig_dir, "calcification_barplot.png"), calc_plot,
       width = 5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "calcification_barplot.pdf"), calc_plot,
       width = 5, height = 4.5)

cat("Calcification bar plot saved to figures/calcification_barplot.png and .pdf\n")

# =============================================================================
# 7. ANOVA ON CALCIFICATION
# =============================================================================

cat("\nRunning two-way ANOVA on calcification...\n")

# Prepare data for ANOVA
calc_anova_data <- calc_processed %>%
  mutate(
    season = factor(season_label, levels = c("Summer", "Winter")),
    treatment = factor(treatment_label, levels = c("Ambient", "OA"))
  )

# Two-way ANOVA (additive model)
calc_anova <- aov(calc_normalized ~ treatment + season, data = calc_anova_data)

cat("\n--- ANOVA Results: Calcification ~ Treatment + Season ---\n")
print(summary(calc_anova))

# Type II ANOVA
cat("\n--- Type II ANOVA (car::Anova) ---\n")
calc_anova_typeII <- Anova(calc_anova, type = "II")
print(calc_anova_typeII)

# Save ANOVA results
sink(file.path(output_dir, "calcification_anova_results.txt"))
cat("Two-way ANOVA: Calcification ~ Treatment + Season (Additive Model)\n")
cat("=", rep("=", 60), "\n\n", sep = "")
cat("Type I ANOVA (Sequential):\n")
print(summary(calc_anova))
cat("\nType II ANOVA (Marginal):\n")
print(calc_anova_typeII)
cat("\n\nGroup Means and Summary:\n")
print(calc_summary)
cat("\n\nIndividual Sample Data:\n")
print(calc_processed %>%
        select(sample, mesocosm, treatment_label, season_label,
               tp1, g_day, calc_normalized) %>%
        arrange(season_label, treatment_label, mesocosm))
sink()

cat("Calcification ANOVA results saved to output/calcification_anova_results.txt\n")

# =============================================================================
# 8. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("\nOutput files generated:\n")
cat("  Tables:\n")
cat("    - output/carbonate_chemistry_table.docx\n")
cat("    - output/carbonate_chemistry_summary.csv\n")
cat("  Figures:\n")
cat("    - figures/calcification_barplot.png\n")
cat("    - figures/calcification_barplot.pdf\n")
cat("  Statistics:\n")
cat("    - output/pH_anova_results.txt\n")
cat("    - output/calcification_anova_results.txt\n")
cat("\n")
