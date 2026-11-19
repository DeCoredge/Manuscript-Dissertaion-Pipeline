rm(list = ls())

setwd("C:/Users/DeCorey Bolton Jr/Documents/GitHub/Manuscript-Dissertaion-Figures-Graphs")

library(ggplot2)
library(readr)
library(cowplot) # A reliable alternative to avoid the layout "+" operator error

# ==============================================================================
# 2. LOAD AND PREPARE DATASETS
# ==============================================================================
ceph_data   <- read_csv("relative_reads_cephalopod_trawl_species_average_pooled.csv")
fish_data   <- read_csv("relative_reads_fish_trawl_species_average_pooled.csv")
invert_data <- read_csv("relative_reads_invertebrate_trawl_species_average_pooled.csv")


# ==============================================================================
# 3. GENERATE THE INDIVIDUAL PLOTS
# ==============================================================================

# Plot 1: Cephalopods
p1 <- ggplot(ceph_data, aes(x = mean_edna_rel_read_count, y = mean_trawl_weight)) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_smooth(aes(color = sample_type, fill = sample_type), 
              method = "lm", linetype = "dashed", alpha = 0.35) +
  scale_y_log10(labels = scales::label_log()) +
  labs(x = NULL,
    y = NULL,
    color = "eDNA Collection Method Across All Trawls (Pooled)", fill = "eDNA Collection Method Across all Trawls (Pooled)"
  ) +
  theme_bw() + theme( axis.line = element_line(color = "black", linewidth = 0.9, linetype = "solid")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold")) + theme(legend.position = "none") # Handled globally by cowplot below

# Plot 2: Invertebrates
p2 <- ggplot(invert_data, aes(x = mean_edna_rel_read_count, y = mean_trawl_weight)) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_smooth(aes(color = sample_type, fill = sample_type), 
              method = "lm", linetype = "dashed", alpha = 0.35) +
  scale_y_log10(labels = scales::label_log()) +
  labs(x = NULL,
    y = NULL
  ) +
  theme_bw() + theme( axis.line = element_line(color = "black", linewidth = 0.9, linetype = "solid")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold")) + theme(legend.position = "none")

# Plot 3: Fish (Keep the legend on this one to extract it clean)
p3 <- ggplot(fish_data, aes(x = mean_edna_rel_read_count, y = mean_trawl_weight)) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_smooth(aes(color = sample_type, fill = sample_type), 
              method = "lm", linetype = "dashed", alpha = 0.35) +
  scale_y_log10(labels = scales::label_log()) +
  scale_color_manual(
    name = "eDNA Collection Method Across All Trawls (Pooled)",
    values = c("Bottom_metaprobe" = "#F8766D", "Slush" = "#00BA38", "Top_metaprobe" = "#619CFF"),
    labels = c("Bottom_metaprobe" = "Bottom Metaprobe", "Slush" = "Slush", "Top_metaprobe" = "Top Metaprobe")) +
  scale_fill_manual(name = "eDNA Collection Method Across All Trawls (Wide)", values = c("Bottom_metaprobe" = "#F8766D", "Slush" = "#00BA38", "Top_metaprobe" = "#619CFF"),
                    labels = c("Bottom_metaprobe" = "Bottom Metaprobe", "Slush" = "Slush", "Top_metaprobe" = "Top Metaprobe")) +
  labs(x = NULL,
       y = NULL,
       color = "eDNA Collection Method Across All Trawls (Pooled)", fill = "eDNA Collection Method Across All Trawls (Pooled)") +
  theme_bw() + theme( axis.line = element_line(color = "black", linewidth = 0.9, linetype = "solid")) +
theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold")) + theme(legend.position = "right")


# ==============================================================================
# 4. REGULAR LINEAR MODELING (lm) WITH TRANSFORMS
# ==============================================================================

# We apply an arcsine square-root transformation to the response variable proportions 
# to better satisfy ordinary least squares assumptions (normality and variance).
# If your proportions contain exact 0s, a tiny offset can be added, e.g., log10(x + 0.0001).

# --- FISH LINEAR MODELS ---
fish_lm_reduced <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight), 
                      data = fish_data)
fish_lm_full    <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight) + sample_type, 
                      data = fish_data)
# Compare nested Gaussian models using an F-test
lrt_fish_lm     <- anova(fish_lm_reduced, fish_lm_full)

# --- INVERTEBRATE LINEAR MODELS ---
invert_lm_reduced <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight), 
                        data = invert_data)
invert_lm_full    <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight) + sample_type, 
                        data = invert_data)
lrt_invert_lm     <- anova(invert_lm_reduced, invert_lm_full)

# --- CEPHALOPOD LINEAR MODELS ---
ceph_lm_reduced <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight), 
                      data = ceph_data)
ceph_lm_full    <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight) + sample_type, 
                      data = ceph_data)
lrt_ceph_lm     <- anova(ceph_lm_reduced, ceph_lm_full)


# ==============================================================================
# 5. PRINT LINEAR ANOVA OUTPUTS
# ==============================================================================
cat("\n==================================================\n")
cat("NESTED F-TEST: FISH (Proportion Scale via lm)\n")
cat("==================================================\n")
print(lrt_fish_lm)

cat("\n==================================================\n")
cat("NESTED F-TEST: INVERTEBRATES (Proportion Scale via lm)\n")
cat("==================================================\n")
print(lrt_invert_lm)

cat("\n==================================================\n")
cat("NESTED F-TEST: CEPHALOPODS (Proportion Scale via lm)\n")
cat("==================================================\n")
print(lrt_ceph_lm)

# Interaction Model: Tests if the predictive power changes by collection type
fish_interaction <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight) * sample_type, data = fish_data)

summary(fish_interaction)


# 1. Fit the reduced linear model (without collection method)
fish_lm_reduced <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight), 
                      data = fish_data)

# 2. Fit the full interaction linear model (with collection method and its interaction)
fish_interaction <- lm(asin(sqrt(mean_edna_rel_read_count)) ~ log10(mean_trawl_weight) * sample_type, 
                       data = fish_data)

# 3. Perform the formal nested ANOVA F-test
nested_f_test <- anova(fish_lm_reduced, fish_interaction)
print(nested_f_test)

# Rebuild Plot 3 to guarantee the legend is present before extraction
p3 <- ggplot(fish_data, aes(x = mean_edna_rel_read_count, y = mean_trawl_weight)) +
  geom_point(aes(color = sample_type), size = 4) +
  geom_smooth(aes(color = sample_type, fill = sample_type), 
              method = "lm", linetype = "dashed", alpha = 0.35) +
  scale_y_log10(labels = scales::label_log()) +
  scale_color_manual(
    name = "eDNA Collection Method Across All Trawls",
    values = c("Bottom_metaprobe" = "#F8766D", "Slush" = "#00BA38", "Top_metaprobe" = "#619CFF"),
    labels = c("Bottom_metaprobe" = "Bottom Metaprobe", "Slush" = "Slush", "Top_metaprobe" = "Top Metaprobe")) +
  scale_fill_manual(
    name = "eDNA Collection Method Across All Trawls", 
    values = c("Bottom_metaprobe" = "#F8766D", "Slush" = "#00BA38", "Top_metaprobe" = "#619CFF"),
    labels = c("Bottom_metaprobe" = "Bottom Metaprobe", "Slush" = "Slush", "Top_metaprobe" = "Top Metaprobe")) +
  labs(x = NULL, 
    y = NULL) +
  theme_bw() + theme( axis.line = element_line(color = "black", linewidth = 0.9, linetype = "solid")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 16, face = "bold")) + theme(legend.position = "right")

# ==============================================================================
# 6. EXTRACT LEGEND & ASSEMBLE SIDE-BY-SIDE WITHOUT THE "+" ERROR
# ==============================================================================

# Extract the shared legend structure from plot 3
shared_legend <- get_legend(p3)

# Extract the shared legend (compatible with ggplot2 >= 3.5.0)
shared_legend <- get_plot_component(p3, "guide-box-right", return_all = TRUE)

# Fallback in case an older ggplot2 version is installed
if (is.null(shared_legend) || length(shared_legend) == 0) {
  shared_legend <- get_legend(p3)}

# Remove the legend from plot 3 so it matches the other two plots visually
p3 <- p3 + theme(legend.position = "none")

# Arrange the three clean charts side-by-side cleanly
plots_row <- plot_grid(p3, p2, p1, ncol = 3, align = 'h', axis = 'b',
                       labels = c("A", "B", "C"),  # Adds panel labels
                       label_size = 20,             # Adjust font size of A, B, C
                       label_fontface = "bold")

# 2. Add the single shared X-axis title directly underneath the 3-plot row
plots_with_x <- add_sub(
  plots_row, 
  "Mean Relative Read Counts", 
  vpadding = grid::unit(1, "lines"),
  y = 0.5, 
  x = 0.5, 
  size = 20, 
  fontface = "bold")

# 3. Create a dedicated vertical label object for the shared Y-axis title
y_axis_label <- ggdraw() + 
  draw_label(
    "Mean Biomass (Log10 Scale)", 
    angle = 90, 
    size = 20, 
    fontface = "bold")

# 4. Combine the Y-axis label, the 3-plot grid (with X-axis), and the shared legend
final_composite_plot <- plot_grid(
  y_axis_label, 
  plots_with_x, 
  shared_legend, 
  ncol = 3, 
  rel_widths = c(0.08, 3, 0.45)) # Allocates dedicated space for Y-axis, plots, and legend

# Output display and save operations
print(final_composite_plot)
ggsave("biomass_read_counts_side_by_side.png", final_composite_plot, width = 22, height = 6, dpi = 300)


rm(list = ls())

#Perform a Pearson & Spearman Statistic test for each sample type
edna_data<- read.csv("relative_reads_fish_trawl_species_average_pooled.csv", header = TRUE)
edna_data$Biomass <- round(edna_data$mean_trawl_weight)
edna_data$MiFish <- edna_data$mean_edna_rel_read_count
edna_data$species <- as.factor(edna_data$species)
edna_data$sample_type <- as.factor(edna_data$sample_type)


# 1. Pearson Correlation for each sample_type

print("--- Pearson Correlation Results ---")
by(edna_data, edna_data$sample_type, function(subsample) {
  cor.test(subsample$Biomass, subsample$MiFish, method = "pearson", na.action = na.omit)
})

# 2. Spearman Correlation for each sample_type
print("--- Spearman Correlation Results ---")
by(edna_data, edna_data$sample_type, function(subsample) {
  cor.test(subsample$Biomass, subsample$MiFish, method = "spearman", na.action = na.omit)
})





rm(list = ls())


#Perform a Poisson Regression GLM using the abundance proportion of Fish species to their read counts by each station csv.
edna_data<- read.csv("SP23096_relative_reads_fish_trawl_species_average_by_station.csv", header = TRUE)
edna_data$Biomass <- round(edna_data$trawl_weight) 
edna_data$MiFish <- edna_data$best_rel_reads
edna_data$species <- as.factor(edna_data$species)

# Correlation Testing (Biomass PER STATION vs MiFish eDNA Read Counts)
# ---------------------------------------------------------

# 1. Pearson Correlation Test (Evaluates linear relationship strength)
pearson_result <- cor.test(edna_data$Biomass, edna_data$MiFish, 
                           method = "pearson")

print("--- PEARSON CORRELATION RESULTS ---")
print(pearson_result)


# 2. Spearman Rank Correlation Test (Evaluates non-linear/monotonic relationship)
# Recommended for skewed biomass counts and proportional eDNA reads
spearman_result <- cor.test(edna_data$Biomass, edna_data$MiFish, 
                            method = "spearman",
                            exact = FALSE) # ADD THIS LINE to silence the ties warning

print("--- SPEARMAN CORRELATION RESULTS ---")
print(spearman_result)


rm(list = ls())


#Perform a Pearson & Spearman Statistic test for each sample type
edna_data<- read.csv("relative_reads_invertebrate_trawl_species_average_pooled.csv", header = TRUE)
edna_data$Biomass <- round(edna_data$mean_trawl_weight)
edna_data$Leray <- edna_data$mean_edna_rel_read_count
edna_data$species <- as.factor(edna_data$species)
edna_data$sample_type <- as.factor(edna_data$sample_type)


# 1. Pearson Correlation for each sample_type

print("--- Pearson Correlation Results ---")
by(edna_data, edna_data$sample_type, function(subsample) {
  cor.test(subsample$Biomass, subsample$Leray, method = "pearson", na.action = na.omit)
})

# 2. Spearman Correlation for each sample_type
print("--- Spearman Correlation Results ---")
by(edna_data, edna_data$sample_type, function(subsample) {
  cor.test(subsample$Biomass, subsample$Leray, method = "spearman", na.action = na.omit)
})


rm(list = ls())

#Perform a Poisson Regression GLM using the abundance proportion of Invertebrate species to their read counts by each station in csv.
edna_data<- read.csv("SP23096_relative_reads_invertebrate_trawl_species_average_by_station.csv", header = TRUE)
edna_data$Biomass <- round(edna_data$trawl_weight) 
edna_data$Leray <- edna_data$best_rel_reads
edna_data$species <- as.factor(edna_data$species)


# Correlation Testing (Biomass across all Trawls vs MiFish eDNA Read Counts)
# ---------------------------------------------------------

# 1. Pearson Correlation Test (Evaluates linear relationship strength)
pearson_result <- cor.test(edna_data$Biomass, edna_data$Leray, 
                           method = "pearson")

print("--- PEARSON CORRELATION RESULTS ---")
print(pearson_result)


# 2. Spearman Rank Correlation Test (Evaluates non-linear/monotonic relationship)
# Recommended for skewed biomass counts and proportional eDNA reads
spearman_result <- cor.test(edna_data$Biomass, edna_data$Leray, 
                            method = "spearman",
                            exact = FALSE) # ADD THIS LINE to silence the ties warning

print("--- SPEARMAN CORRELATION RESULTS ---")
print(spearman_result)


rm(list = ls())


#Perform a Pearson & Spearman Statistic test for each sample type
edna_data<- read.csv("relative_reads_cephalopod_trawl_species_average_pooled.csv", header = TRUE)
edna_data$Biomass <- round(edna_data$mean_trawl_weight)
edna_data$Ceph18s <- edna_data$mean_edna_rel_read_count
edna_data$species <- as.factor(edna_data$species)
edna_data$sample_type <- as.factor(edna_data$sample_type)



rm(list = ls())
