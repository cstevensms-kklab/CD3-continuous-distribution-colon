################################################################################
# CD3+ CELL LOCATION ALONG THE COLONIC CRYPT AXIS
# Boundary-corrected KDE + functional group comparison by age
#
# PURPOSE
#   This script compares the continuous spatial distribution of CD3+ cells
#   along a normalized crypt axis (0 = crypt base, 1 = luminal surface).
#
#   It does NOT test absolute CD3+ abundance. Each animal's KDE is normalized
#   to integrate to 1, so the analysis asks:
#     "Given that a cell is CD3+, does its relative location differ by age?"
#
# PRIMARY STATISTICAL APPROACH
#   1. The animal is the independent unit.
#   2. All animals are evaluated on the same 0-1 grid.
#   3. A single, age-blind bandwidth is used within each tissue analysis.
#   4. Reflection corrects KDE boundary bias at 0 and 1.
#   5. Density curves are centered-log-ratio (CLR) transformed before testing.
#   6. A studentized pointwise statistic (Welch t^2) is calculated.
#   7. Animal age labels are permuted to obtain:
#        - a global max-statistic p-value;
#        - a simultaneous 95% threshold across the whole crypt axis;
#        - an exploratory pointwise 95% threshold.
#
# INTERPRETATION OF THE F-STATISTIC PANEL
#   The black curve is the observed Welch-type F statistic (Welch t^2).
#   The dotted curve is an unadjusted pointwise 95% permutation threshold.
#   The dashed horizontal line is the simultaneous max-statistic threshold.
#   Only regions above the dashed line control family-wise error across the axis.
################################################################################

# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

# Input 2 CSV files: (1) a cell file containing y-coordinates, (2) a covariate file containing animal information
CELL_FILE <- file.path(ANALYSIS_DIR, "colon_cd3_ycoords.csv")
COVARIATE_FILE <- file.path(ANALYSIS_DIR, "demographics_colon.csv")

# The raw y-coordinate must increase from crypt base toward luminal surface.
# Set to TRUE only if your raw coordinate increases in the opposite direction.
REVERSE_AXIS <- FALSE

# Common grid. 301 points includes exactly 0, 1/3, 2/3, and 1.
GRID_N <- 301
GRID <- seq(0, 1, length.out = GRID_N)

# Minimum number of CD3+ cells required for an animal to contribute a KDE.
# There is no universal cutoff. This value should be chosen before examining
# age-group differences and reported in the manuscript.
MIN_CD3_CELLS <- 10

# Final permutation count. Use 999 during code testing; restore 9999 for results.
B_MAIN <- 9999
ALPHA <- 0.05
RANDOM_SEED <- 20260710

# Primary analysis uses a CLR transformation because each KDE is a probability
# density constrained to integrate to 1.
USE_CLR_PRIMARY <- TRUE

# Recommended sensitivity analyses. These take additional computation time.
RUN_SENSITIVITY <- TRUE
B_SENSITIVITY <- 1999
BANDWIDTH_MULTIPLIERS <- c(0.75, 1.25)

# Optional global Fmaxb analysis from fdANOVA. This is a sensitivity analysis,
# not the primary localized test, and it returns a global result only.
RUN_FDANOVA_SENSITIVITY <- FALSE

# Plot colors (colorblind-friendly orange and blue).
GROUP_COLORS <- c(Young = "#D55E00", Old = "#0072B2")


# ==============================================================================
# 2. REQUIRED PACKAGES
# ==============================================================================

required_packages <- c("tidyverse", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(tidyverse)
library(gridExtra)


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# Stop with a clear message when required columns are missing.
check_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(data_name, " is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
}

# Numerical trapezoidal integration.
trapz_numeric <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2) {
    stop("x and y must have the same length and contain at least two values.")
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# Convert age_cat to an explicit Young/Old factor.
# Accepted encodings include 0/1 and Young/Old.
recode_age_group <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- dplyr::case_when(
    z %in% c("0", "young", "y") ~ "Young",
    z %in% c("1", "old", "o")   ~ "Old",
    TRUE                          ~ NA_character_
  )
  factor(out, levels = c("Young", "Old"))
}

# Standardize the compartment labels used by the input file.
recode_compartment <- function(x) {
  z <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    z %in% c("epi", "epithelium", "epithelial") ~ "epi",
    z %in% c("lp", "lamina propria", "lamina_propria") ~ "lp",
    TRUE ~ z
  )
}

# Select a subject-specific bandwidth safely.
# Sheather-Jones is attempted first. The robust nrd0 rule is used as a fallback.
safe_subject_bandwidth <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || length(unique(x)) < 2) return(NA_real_)

  bw <- tryCatch(
    stats::bw.SJ(x, method = "ste"),
    error = function(e) NA_real_,
    warning = function(w) suppressWarnings(
      tryCatch(stats::bw.SJ(x, method = "ste"), error = function(e) NA_real_)
    )
  )

  if (!is.finite(bw) || bw <= 0) bw <- stats::bw.nrd0(x)
  if (!is.finite(bw) || bw <= 0) return(NA_real_)
  bw
}

# Boundary-corrected Gaussian KDE on [0, 1] using reflection.
# Reflected observations are placed at -x and 2-x. The resulting curve is
# renormalized numerically so its area on [0,1] is exactly 1 (up to rounding).
reflection_kde <- function(x, grid, bandwidth, chunk_size = 2000) {
  x <- x[is.finite(x) & x >= 0 & x <= 1]
  if (length(x) < 2) stop("At least two valid positions are required for KDE.")
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("Bandwidth must be positive.")

  kernel_sum <- numeric(length(grid))
  starts <- seq(1, length(x), by = chunk_size)

  for (start in starts) {
    end <- min(start + chunk_size - 1, length(x))
    xx <- x[start:end]

    original <- outer(grid, xx, FUN = "-") / bandwidth
    left_reflection <- outer(grid, -xx, FUN = "-") / bandwidth
    right_reflection <- outer(grid, 2 - xx, FUN = "-") / bandwidth

    kernel_sum <- kernel_sum +
      rowSums(stats::dnorm(original) +
              stats::dnorm(left_reflection) +
              stats::dnorm(right_reflection))
  }

  density_values <- kernel_sum / (length(x) * bandwidth)
  area <- trapz_numeric(grid, density_values)

  if (!is.finite(area) || area <= 0) stop("KDE normalization failed.")
  density_values / area
}

# Transform positive density curves to centered-log-ratio (CLR) curves.
# For each animal, the integral of the transformed curve is approximately zero.
clr_transform <- function(curve_matrix, grid) {
  apply(curve_matrix, 2, function(f) {
    # Gaussian KDEs should be positive, but a very small numerical floor avoids
    # log(0) if values underflow to zero far from all observed cells.
    numerical_floor <- max(f, na.rm = TRUE) * 1e-12
    log_f <- log(pmax(f, numerical_floor))
    mean_log_f <- trapz_numeric(grid, log_f) /
      (max(grid) - min(grid))
    log_f - mean_log_f
  })
}

# Compute the pointwise Welch-type F statistic, which equals Welch t^2 for two
# groups. Studentization makes the statistic less sensitive to unequal group
# variances than the ordinary equal-variance ANOVA F statistic.
welch_t2_curve <- function(curve_matrix, group) {
  group <- factor(group, levels = c("Young", "Old"))
  young <- group == "Young"
  old <- group == "Old"

  n_young <- sum(young)
  n_old <- sum(old)
  if (n_young < 2 || n_old < 2) {
    stop("At least two animals are required in each age group.")
  }

  mean_young <- rowMeans(curve_matrix[, young, drop = FALSE])
  mean_old <- rowMeans(curve_matrix[, old, drop = FALSE])
  var_young <- apply(curve_matrix[, young, drop = FALSE], 1, stats::var)
  var_old <- apply(curve_matrix[, old, drop = FALSE], 1, stats::var)

  denominator <- var_young / n_young + var_old / n_old
  statistic <- (mean_old - mean_young)^2 / denominator
  statistic[!is.finite(statistic) | denominator <= 0] <- 0

  list(
    statistic = statistic,
    mean_young = mean_young,
    mean_old = mean_old,
    difference_old_minus_young = mean_old - mean_young
  )
}

# Animal-level permutation max-statistic test.
# The same permuted age label is used for the entire curve of each animal.
permutation_max_test <- function(curve_matrix, group, grid, B, alpha, seed) {
  observed <- welch_t2_curve(curve_matrix, group)
  observed_stat <- observed$statistic

  permutation_stats <- matrix(NA_real_, nrow = B, ncol = nrow(curve_matrix))
  max_stats <- numeric(B)

  set.seed(seed)
  for (b in seq_len(B)) {
    permuted_group <- sample(group, replace = FALSE)
    stat_b <- welch_t2_curve(curve_matrix, permuted_group)$statistic
    permutation_stats[b, ] <- stat_b
    max_stats[b] <- max(stat_b, na.rm = TRUE)

    if (B >= 5000 && b %% 1000 == 0) {
      message("  Completed ", b, " of ", B, " permutations")
    }
  }

  pointwise_critical <- apply(
    permutation_stats, 2, stats::quantile,
    probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE
  )
  simultaneous_critical <- unname(stats::quantile(
    max_stats, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE
  ))

  # +1 correction prevents a permutation p-value of exactly zero.
  global_p <- (1 + sum(max_stats >= max(observed_stat, na.rm = TRUE))) / (B + 1)

  pointwise_p <- (1 + colSums(
    sweep(permutation_stats, 2, observed_stat, FUN = ">="),
    na.rm = TRUE
  )) / (B + 1)

  simultaneous_p <- vapply(
    observed_stat,
    function(value) (1 + sum(max_stats >= value, na.rm = TRUE)) / (B + 1),
    numeric(1)
  )

  list(
    statistic = observed_stat,
    pointwise_critical = pointwise_critical,
    simultaneous_critical = simultaneous_critical,
    global_p = global_p,
    pointwise_p = pointwise_p,
    simultaneous_p = simultaneous_p,
    significant = observed_stat > simultaneous_critical
  )
}

# Convert adjacent significant grid points into approximate crypt-axis intervals.
extract_significant_intervals <- function(grid, significant, statistic,
                                          original_difference) {
  if (!any(significant)) {
    return(tibble(
      start = numeric(), end = numeric(), peak_position = numeric(),
      peak_statistic = numeric(), density_difference_old_minus_young = numeric(),
      direction_at_peak = character()
    ))
  }

  runs <- rle(significant)
  run_ends <- cumsum(runs$lengths)
  run_starts <- c(1, head(run_ends, -1) + 1)
  keep <- which(runs$values)

  purrr::map_dfr(keep, function(k) {
    idx <- run_starts[k]:run_ends[k]
    peak_idx <- idx[which.max(statistic[idx])]
    difference <- original_difference[peak_idx]

    tibble(
      start = grid[min(idx)],
      end = grid[max(idx)],
      peak_position = grid[peak_idx],
      peak_statistic = statistic[peak_idx],
      density_difference_old_minus_young = difference,
      direction_at_peak = if_else(
        difference > 0,
        "Higher relative density in old",
        "Higher relative density in young"
      )
    )
  })
}

# Prepare data, generate subject-level KDE curves, and perform the functional test.
run_spatial_analysis <- function(
    cd3_data,
    analysis_name,
    compartment = c("total", "epi", "lp"),
    grid = GRID,
    min_cells = MIN_CD3_CELLS,
    bandwidth = NULL,
    use_clr = TRUE,
    B = B_MAIN,
    alpha = ALPHA,
    seed = RANDOM_SEED,
    save_outputs = TRUE,
    make_plot = TRUE) {

  compartment <- match.arg(compartment)
  message("\nRunning analysis: ", analysis_name)

  analysis_data <- cd3_data
  if (compartment != "total") {
    analysis_data <- analysis_data %>% filter(compartment_std == compartment)
  }

  counts <- analysis_data %>%
    count(animal_id, age_group, name = "n_cd3_cells") %>%
    arrange(age_group, animal_id) %>%
    mutate(included = n_cd3_cells >= min_cells)

  excluded <- counts %>% filter(!included)
  if (nrow(excluded) > 0) {
    warning(
      analysis_name, ": excluding ", nrow(excluded),
      " animal(s) with fewer than ", min_cells, " CD3+ cells: ",
      paste(excluded$animal_id, collapse = ", ")
    )
  }

  group_map <- counts %>%
    filter(included) %>%
    transmute(animal_id = as.character(animal_id), age_group) %>%
    arrange(age_group, animal_id)

  analysis_data <- analysis_data %>%
    mutate(animal_id = as.character(animal_id)) %>%
    semi_join(group_map, by = "animal_id")

  group_map <- group_map %>%
    mutate(age_group = factor(age_group, levels = c("Young", "Old")))

  if (any(table(group_map$age_group) < 2)) {
    stop(analysis_name, ": fewer than two included animals in one age group.")
  }

  # Select one common bandwidth without using age labels.
  subject_bandwidths <- analysis_data %>%
    group_by(animal_id) %>%
    summarise(
      n_cd3_cells = n(),
      subject_bandwidth = safe_subject_bandwidth(y_scaled),
      .groups = "drop"
    )

  if (is.null(bandwidth)) {
    bandwidth <- stats::median(
      subject_bandwidths$subject_bandwidth,
      na.rm = TRUE
    )
  }
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop(analysis_name, ": unable to obtain a valid common bandwidth.")
  }

  message("  Included animals: ", nrow(group_map),
          " (Young = ", sum(group_map$age_group == "Young"),
          ", Old = ", sum(group_map$age_group == "Old"), ")")
  message("  Common bandwidth: ", format(round(bandwidth, 5), nsmall = 5))

  # Build a grid-by-animal matrix. Column names and age labels are linked through
  # the explicit group_map, preventing animal_id/group-order mismatches.
  curve_matrix <- vapply(group_map$animal_id, function(id) {
    positions <- analysis_data %>%
      filter(animal_id == id) %>%
      pull(y_scaled)
    reflection_kde(positions, grid, bandwidth)
  }, FUN.VALUE = numeric(length(grid)))

  colnames(curve_matrix) <- group_map$animal_id
  stopifnot(identical(colnames(curve_matrix), group_map$animal_id))

  # Confirm that every density integrates to approximately 1.
  curve_areas <- apply(curve_matrix, 2, function(f) trapz_numeric(grid, f))
  if (any(abs(curve_areas - 1) > 1e-6)) {
    warning(analysis_name, ": one or more KDE curves did not integrate to 1.")
  }

  test_matrix <- if (use_clr) clr_transform(curve_matrix, grid) else curve_matrix

  test_result <- permutation_max_test(
    curve_matrix = test_matrix,
    group = group_map$age_group,
    grid = grid,
    B = B,
    alpha = alpha,
    seed = seed
  )

  # Direction is reported on the original KDE scale, not the CLR scale.
  original_means <- welch_t2_curve(curve_matrix, group_map$age_group)
  intervals <- extract_significant_intervals(
    grid = grid,
    significant = test_result$significant,
    statistic = test_result$statistic,
    original_difference = original_means$difference_old_minus_young
  )

  pointwise_results <- tibble(
    crypt_position = grid,
    mean_density_young = original_means$mean_young,
    mean_density_old = original_means$mean_old,
    density_difference_old_minus_young =
      original_means$difference_old_minus_young,
    welch_t2 = test_result$statistic,
    pointwise_critical_95 = test_result$pointwise_critical,
    simultaneous_critical_95 = test_result$simultaneous_critical,
    pointwise_p_unadjusted = test_result$pointwise_p,
    simultaneous_p_fwer = test_result$simultaneous_p,
    significant_simultaneous = test_result$significant
  )

  curve_long <- as.data.frame(curve_matrix) %>%
    mutate(crypt_position = grid) %>%
    pivot_longer(
      cols = -crypt_position,
      names_to = "animal_id",
      values_to = "relative_cd3_location_density"
    ) %>%
    left_join(group_map, by = "animal_id")

  output_stub <- stringr::str_replace_all(tolower(analysis_name), "[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
  analysis_dir <- file.path(OUTPUT_DIR, output_stub)

  if (save_outputs) {
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(counts, file.path(analysis_dir, "animal_cell_counts_and_inclusion.csv"))
    readr::write_csv(subject_bandwidths,
                     file.path(analysis_dir, "subject_bandwidths.csv"))
    readr::write_csv(curve_long, file.path(analysis_dir, "animal_kde_curves.csv"))
    readr::write_csv(pointwise_results,
                     file.path(analysis_dir, "functional_pointwise_results.csv"))
    readr::write_csv(intervals,
                     file.path(analysis_dir, "simultaneously_significant_intervals.csv"))
  }

  # --------------------------------------------------------------------------
  # Plot: original KDE curves above; test statistic and thresholds below.
  # --------------------------------------------------------------------------
  plot_object <- NULL
  if (make_plot) {
    mean_curves <- curve_long %>%
      group_by(crypt_position, age_group) %>%
      summarise(
        relative_cd3_location_density =
          mean(relative_cd3_location_density),
        .groups = "drop"
      )

    rectangle_data <- intervals %>%
      transmute(xmin = start, xmax = end, ymin = -Inf, ymax = Inf)

    p_curves <- ggplot() +
      geom_rect(
        data = rectangle_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = "grey50", alpha = 0.12
      ) +
      geom_line(
        data = curve_long,
        aes(
          x = crypt_position,
          y = relative_cd3_location_density,
          group = animal_id,
          color = age_group
        ),
        alpha = 0.28, linewidth = 0.45
      ) +
      geom_line(
        data = mean_curves,
        aes(
          x = crypt_position,
          y = relative_cd3_location_density,
          color = age_group
        ),
        linewidth = 1.45
      ) +
      scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        title = analysis_name,
        subtitle = paste0(
          "Boundary-corrected KDE; common bandwidth = ",
          format(round(bandwidth, 4), nsmall = 4)
        ),
        x = NULL,
        y = "Relative CD3+ location density",
        color = "Age group"
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold")
      )

    p_stat <- ggplot(pointwise_results, aes(x = crypt_position)) +
      geom_rect(
        data = rectangle_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = "grey50", alpha = 0.12
      ) +
      geom_line(
        aes(y = welch_t2),
        color = "black", linewidth = 1.0
      ) +
      geom_line(
        aes(y = pointwise_critical_95),
        color = "grey35", linetype = "dotted", linewidth = 1.35
      ) +
      geom_hline(
        yintercept = test_result$simultaneous_critical,
        color = "grey35", linetype = "dashed", linewidth = 1.35
      ) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        subtitle = paste0(
          "Global max-statistic permutation p = ",
          format.pval(test_result$global_p, digits = 3, eps = 1 / (B + 1)),
          "; dotted = pointwise 95%; dashed = simultaneous 95%"
        ),
        x = "Normalized crypt axis (0 = base, 1 = luminal surface)",
        y = if (use_clr) "Welch-type F (t²), CLR scale" else "Welch-type F (t²)"
      ) +
      theme_classic(base_size = 11)

    plot_object <- gridExtra::arrangeGrob(
      p_curves, p_stat, ncol = 1, heights = c(2, 1.15)
    )

    if (save_outputs) {
      ggsave(
        file.path(analysis_dir, "kde_and_functional_test.png"),
        plot = plot_object, width = 8.2, height = 8.5, dpi = 300
      )
      ggsave(
        file.path(analysis_dir, "kde_and_functional_test.pdf"),
        plot = plot_object, width = 8.2, height = 8.5
      )
    }
  }

  list(
    analysis_name = analysis_name,
    compartment = compartment,
    counts = counts,
    group_map = group_map,
    subject_bandwidths = subject_bandwidths,
    bandwidth = bandwidth,
    use_clr = use_clr,
    curve_matrix = curve_matrix,
    test_matrix = test_matrix,
    pointwise_results = pointwise_results,
    intervals = intervals,
    global_p = test_result$global_p,
    max_statistic = max(test_result$statistic, na.rm = TRUE),
    simultaneous_critical = test_result$simultaneous_critical,
    plot = plot_object
  )
}


# ==============================================================================
# 4. READ, CLEAN, MERGE, AND SCALE THE DATA
# ==============================================================================

cell_data <- readr::read_csv("//medctr.ad.wfubmc.edu/dfs/wfusmpath_research$/kklab/Christina Stevens/Project - Colon IHC/Christina IHC analysis/xy coordinate analysis/colon_cd3_ycoords.csv", show_col_types = FALSE)
covariate_data <- readr::read_csv("//medctr.ad.wfubmc.edu/dfs/wfusmpath_research$/kklab/Christina Stevens/Project - Colon IHC/Christina IHC analysis/xy coordinate analysis/demographics_colon.csv", show_col_types = FALSE)

check_columns(
  cell_data,
  c("animal_id", "nuc_label", "lp_epi", "yprime"),
  "Cell-level data"
)
check_columns(
  covariate_data,
  c("animal_id", "age_cat"),
  "Covariate data"
)

# Only age is needed for this analysis. Restricting the merge to needed columns
# avoids accidental duplication from unrelated covariates.
age_map <- covariate_data %>%
  transmute(
    animal_id = as.character(animal_id),
    age_group = recode_age_group(age_cat)
  ) %>%
  distinct()

if (anyNA(age_map$age_group)) {
  bad_ids <- age_map %>% filter(is.na(age_group)) %>% pull(animal_id)
  stop("Unrecognized age_cat values for animal(s): ", paste(bad_ids, collapse = ", "))
}

conflicting_age <- age_map %>%
  count(animal_id, name = "n_age_labels") %>%
  filter(n_age_labels > 1)
if (nrow(conflicting_age) > 0) {
  stop("Conflicting age labels for animal(s): ",
       paste(conflicting_age$animal_id, collapse = ", "))
}

cell_animal_ids <- unique(as.character(cell_data$animal_id))
missing_age_ids <- setdiff(cell_animal_ids, age_map$animal_id)
if (length(missing_age_ids) > 0) {
  stop(
    "The following animals occur in the cell file but have no usable age label: ",
    paste(missing_age_ids, collapse = ", ")
  )
}

combined_data <- cell_data %>%
  transmute(
    animal_id = as.character(animal_id),
    nuc_label = trimws(tolower(as.character(nuc_label))),
    compartment_std = recode_compartment(lp_epi),
    yprime = as.numeric(yprime)
  ) %>%
  filter(!is.na(animal_id), !is.na(yprime)) %>%
  inner_join(age_map, by = "animal_id")

unknown_compartments <- combined_data %>%
  filter(!is.na(compartment_std), !compartment_std %in% c("epi", "lp")) %>%
  distinct(compartment_std) %>%
  pull(compartment_std)
if (length(unknown_compartments) > 0) {
  warning(
    "Unrecognized lp_epi value(s) will contribute to total mucosa but not to ",
    "the epithelium or lamina propria analyses: ",
    paste(unknown_compartments, collapse = ", ")
  )
}

# Direct min-max scaling is equivalent to the earlier midpoint-alignment followed
# by min-max scaling. All nuclei, not only CD3+ nuclei, define each animal's axis.
combined_data <- combined_data %>%
  group_by(animal_id) %>%
  mutate(
    y_min = min(yprime, na.rm = TRUE),
    y_max = max(yprime, na.rm = TRUE),
    y_range = y_max - y_min,
    y_scaled = (yprime - y_min) / y_range
  ) %>%
  ungroup()

zero_range_ids <- combined_data %>%
  distinct(animal_id, y_range) %>%
  filter(!is.finite(y_range) | y_range <= 0) %>%
  pull(animal_id)
if (length(zero_range_ids) > 0) {
  stop("Zero or invalid y-coordinate range for animal(s): ",
       paste(zero_range_ids, collapse = ", "))
}

if (REVERSE_AXIS) {
  combined_data <- combined_data %>% mutate(y_scaled = 1 - y_scaled)
}

if (any(combined_data$y_scaled < -1e-10 | combined_data$y_scaled > 1 + 1e-10)) {
  stop("Scaled y-coordinates fall outside [0,1]. Check the input coordinates.")
}

# Retain only CD3+ nuclei for the conditional spatial-distribution analysis.
cd3_positive <- combined_data %>%
  filter(nuc_label == "positive")

if (nrow(cd3_positive) == 0) {
  stop("No rows with nuc_label == 'positive' were found.")
}

# Basic quality-control output.
qc_by_animal <- combined_data %>%
  group_by(animal_id, age_group) %>%
  summarise(
    n_all_nuclei = n(),
    n_cd3_positive = sum(nuc_label == "positive"),
    n_cd3_epi = sum(nuc_label == "positive" & compartment_std == "epi"),
    n_cd3_lp = sum(nuc_label == "positive" & compartment_std == "lp"),
    raw_y_min = first(y_min),
    raw_y_max = first(y_max),
    raw_y_range = first(y_range),
    .groups = "drop"
  )

readr::write_csv(qc_by_animal, file.path(OUTPUT_DIR, "qc_by_animal.csv"))
readr::write_csv(combined_data,
                 file.path(OUTPUT_DIR, "combined_data_with_y_scaled.csv"))


# ==============================================================================
# 5. PRIMARY FUNCTIONAL ANALYSES
# ==============================================================================

results <- list(
  total = run_spatial_analysis(
    cd3_data = cd3_positive,
    analysis_name = "Total mucosal CD3+ cell distribution",
    compartment = "total",
    use_clr = USE_CLR_PRIMARY,
    B = B_MAIN,
    seed = RANDOM_SEED + 1
  ),
  epithelium = run_spatial_analysis(
    cd3_data = cd3_positive,
    analysis_name = "Epithelial CD3+ cell distribution",
    compartment = "epi",
    use_clr = USE_CLR_PRIMARY,
    B = B_MAIN,
    seed = RANDOM_SEED + 2
  ),
  lamina_propria = run_spatial_analysis(
    cd3_data = cd3_positive,
    analysis_name = "Lamina propria CD3+ cell distribution",
    compartment = "lp",
    use_clr = USE_CLR_PRIMARY,
    B = B_MAIN,
    seed = RANDOM_SEED + 3
  )
)

# Total mucosa can be treated as the prespecified primary analysis. Epithelium
# and lamina propria are secondary compartment analyses. Both all-three and
# secondary-only Holm adjustments are provided for transparency.
global_summary <- purrr::imap_dfr(results, function(res, key) {
  tibble(
    analysis_key = key,
    analysis = res$analysis_name,
    n_young = sum(res$group_map$age_group == "Young"),
    n_old = sum(res$group_map$age_group == "Old"),
    common_bandwidth = res$bandwidth,
    transformation = if_else(res$use_clr, "CLR", "None"),
    max_statistic = res$max_statistic,
    simultaneous_critical_95 = res$simultaneous_critical,
    global_permutation_p = res$global_p,
    n_significant_intervals = nrow(res$intervals)
  )
})

global_summary <- global_summary %>%
  mutate(p_holm_all_three = p.adjust(global_permutation_p, method = "holm"))

global_summary$p_holm_secondary_only <- NA_real_
secondary_rows <- which(global_summary$analysis_key != "total")
global_summary$p_holm_secondary_only[secondary_rows] <- p.adjust(
  global_summary$global_permutation_p[secondary_rows], method = "holm"
)

readr::write_csv(global_summary,
                 file.path(OUTPUT_DIR, "global_functional_test_summary.csv"))


# ==============================================================================
# 6. SENSITIVITY ANALYSES
# ==============================================================================

sensitivity_summary <- tibble()

if (RUN_SENSITIVITY) {
  message("\nRunning sensitivity analyses...")

  analysis_specs <- tribble(
    ~key, ~analysis_name, ~compartment, ~seed_offset,
    "total", "Total mucosal CD3+ cell distribution", "total", 100,
    "epithelium", "Epithelial CD3+ cell distribution", "epi", 200,
    "lamina_propria", "Lamina propria CD3+ cell distribution", "lp", 300
  )

  sensitivity_summary <- purrr::pmap_dfr(
    analysis_specs,
    function(key, analysis_name, compartment, seed_offset) {
      base_bw <- results[[key]]$bandwidth

      bandwidth_results <- purrr::map_dfr(
        BANDWIDTH_MULTIPLIERS,
        function(multiplier) {
          res <- run_spatial_analysis(
            cd3_data = cd3_positive,
            analysis_name = paste0(analysis_name, " [sensitivity]"),
            compartment = compartment,
            bandwidth = base_bw * multiplier,
            use_clr = USE_CLR_PRIMARY,
            B = B_SENSITIVITY,
            seed = RANDOM_SEED + seed_offset + round(multiplier * 10),
            save_outputs = FALSE,
            make_plot = FALSE
          )

          tibble(
            analysis_key = key,
            sensitivity = paste0("Bandwidth x ", multiplier),
            bandwidth = res$bandwidth,
            transformation = if_else(res$use_clr, "CLR", "None"),
            global_permutation_p = res$global_p,
            max_statistic = res$max_statistic,
            simultaneous_critical_95 = res$simultaneous_critical,
            n_significant_intervals = nrow(res$intervals)
          )
        }
      )

      # Also test the untransformed density curves at the primary bandwidth.
      untransformed_result <- run_spatial_analysis(
        cd3_data = cd3_positive,
        analysis_name = paste0(analysis_name, " [untransformed sensitivity]"),
        compartment = compartment,
        bandwidth = base_bw,
        use_clr = FALSE,
        B = B_SENSITIVITY,
        seed = RANDOM_SEED + seed_offset + 50,
        save_outputs = FALSE,
        make_plot = FALSE
      )

      bind_rows(
        tibble(
          analysis_key = key,
          sensitivity = "Primary specification",
          bandwidth = results[[key]]$bandwidth,
          transformation = if_else(results[[key]]$use_clr, "CLR", "None"),
          global_permutation_p = results[[key]]$global_p,
          max_statistic = results[[key]]$max_statistic,
          simultaneous_critical_95 = results[[key]]$simultaneous_critical,
          n_significant_intervals = nrow(results[[key]]$intervals)
        ),
        bandwidth_results,
        tibble(
          analysis_key = key,
          sensitivity = "No CLR transformation",
          bandwidth = untransformed_result$bandwidth,
          transformation = "None",
          global_permutation_p = untransformed_result$global_p,
          max_statistic = untransformed_result$max_statistic,
          simultaneous_critical_95 = untransformed_result$simultaneous_critical,
          n_significant_intervals = nrow(untransformed_result$intervals)
        )
      )
    }
  )

  readr::write_csv(
    sensitivity_summary,
    file.path(OUTPUT_DIR, "sensitivity_analysis_summary.csv")
  )
}


# ==============================================================================
# 7. OPTIONAL fdANOVA Fmaxb GLOBAL SENSITIVITY TEST
# ==============================================================================

if (RUN_FDANOVA_SENSITIVITY) {
  if (!requireNamespace("fdANOVA", quietly = TRUE)) {
    warning("fdANOVA is not installed; skipping optional Fmaxb sensitivity tests.")
  } else {
    for (key in names(results)) {
      res <- results[[key]]
      fmaxb_result <- fdANOVA::fanova.tests(
        x = res$test_matrix,
        group.label = res$group_map$age_group,
        test = "Fmaxb"
      )
      capture.output(
        print(fmaxb_result),
        file = file.path(OUTPUT_DIR, paste0(key, "_fdANOVA_Fmaxb_output.txt"))
      )
    }
  }
}


# ==============================================================================
# 8. SAVE REPRODUCIBILITY INFORMATION
# ==============================================================================

saveRDS(results, file.path(OUTPUT_DIR, "primary_analysis_objects.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUTPUT_DIR, "R_sessionInfo.txt"))

message("\nAnalysis complete. Outputs were saved to:\n", OUTPUT_DIR)
message("\nImportant interpretation reminder:")
message("  This analysis tests relative CD3+ cell location, conditional on CD3 positivity.")
message("  It does not test absolute CD3+ abundance or CD3+ cells per tissue area.")



# ==============================================================================
# CHECK WHETHER THE FALLBACK BANDWIDTH WAS USED
# This does not rerun the functional analyses or permutations.
# ==============================================================================

library(tidyverse)

diagnose_bandwidth <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) < 2 || length(unique(x)) < 2) {
    return(
      tibble(
        sj_bandwidth = NA_real_,
        nrd0_bandwidth = NA_real_,
        bandwidth_used = NA_real_,
        method_used = "Unable to estimate",
        sj_warning = NA_character_,
        sj_error = NA_character_
      )
    )
  }
  
  sj_warning <- NA_character_
  sj_error <- NA_character_
  
  sj_bandwidth <- tryCatch(
    withCallingHandlers(
      stats::bw.SJ(x, method = "ste"),
      warning = function(w) {
        sj_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      sj_error <<- conditionMessage(e)
      NA_real_
    }
  )
  
  nrd0_bandwidth <- stats::bw.nrd0(x)
  
  if (is.finite(sj_bandwidth) && sj_bandwidth > 0) {
    bandwidth_used <- sj_bandwidth
    method_used <- "Sheather-Jones"
  } else if (is.finite(nrd0_bandwidth) && nrd0_bandwidth > 0) {
    bandwidth_used <- nrd0_bandwidth
    method_used <- "nrd0 fallback"
  } else {
    bandwidth_used <- NA_real_
    method_used <- "Unable to estimate"
  }
  
  tibble(
    sj_bandwidth = sj_bandwidth,
    nrd0_bandwidth = nrd0_bandwidth,
    bandwidth_used = bandwidth_used,
    method_used = method_used,
    sj_warning = sj_warning,
    sj_error = sj_error
  )
}


check_bandwidths <- function(data, compartment = "total") {
  
  analysis_data <- data
  
  if (compartment != "total") {
    analysis_data <- analysis_data %>%
      filter(compartment_std == compartment)
  }
  
  analysis_data %>%
    group_by(animal_id) %>%
    group_modify(
      ~ diagnose_bandwidth(.x$y_scaled)
    ) %>%
    ungroup()
}

bandwidth_check_total <- check_bandwidths(
  cd3_positive,
  compartment = "total"
)

bandwidth_check_epi <- check_bandwidths(
  cd3_positive,
  compartment = "epi"
)

bandwidth_check_lp <- check_bandwidths(
  cd3_positive,
  compartment = "lp"
)



bandwidth_check_total
bandwidth_check_epi
Bandwidth_check_lp



table(bandwidth_check_total$method_used)
table(bandwidth_check_epi$method_used)
table(bandwidth_check_lp$method_used)




write.csv(
  bandwidth_check_total,
  file.path(OUTPUT_DIR, "bandwidth_method_check_total.csv"),
  row.names = FALSE
)

write.csv(
  bandwidth_check_epi,
  file.path(OUTPUT_DIR, "bandwidth_method_check_epi.csv"),
  row.names = FALSE
)

write.csv(
  bandwidth_check_lp,
  file.path(OUTPUT_DIR, "bandwidth_method_check_lp.csv"),
  row.names = FALSE
)

sessionInfo()
