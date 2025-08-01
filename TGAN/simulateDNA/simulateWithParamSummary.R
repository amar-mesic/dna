library(simDNAmixtures)
library(dplyr)

# Output directory
out_dir <- "generated_alleles_fixed_params"
dir.create(out_dir, showWarnings = FALSE)

# Load allele frequencies
allele_freqs_file <- system.file("extdata","FBI_extended_Cauc_022024.csv", package = "simDNAmixtures")
allele_freqs <- read_allele_freqs(allele_freqs_file)
gf <- gf_configuration()

# # Filter loci to only those supported by the size regression
# supported_loci <- names(gf$size_regression)
# gf$log_normal_settings$locus_names <- intersect(gf$log_normal_settings$locus_names, supported_loci)

dye_map <- kits$GlobalFiler[, c("Marker", "Color")] %>% distinct(Marker, .keep_all = TRUE)
size_standard_sizes <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                         220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380,
                         400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560,
                         580, 600)
size_standard_df <- data.frame(
  Locus = "LIZ",
  Allele = NA,
  Height = 3000,
  Size = size_standard_sizes,
  Color = "orange"
)

# Template ratio functions
template_ratio_even <- function(n) rep(1/n, n)
template_ratio_increasing <- function(n) { ratios <- 1:n; ratios / sum(ratios) }
template_ratio_last_dominates <- function(n) { ratios <- rep(1/20, n); ratios[n] <- 1 - sum(ratios[-n]); ratios }
template_ratio_first_two_10 <- function(n) {
  ratios <- rep(0, n)
  if (n == 1) {
    ratios[1] <- 1
  } else if (n == 2) {
    ratios[1] <- 0.1
    ratios[2] <- 0.9
  } else {
    ratios[1:2] <- 0.1
    ratios[3:n] <- (1 - 0.2) / (n - 2)
  }
  ratios
}
template_ratio_functions <- list(
  template_ratio_even,
  template_ratio_increasing,
  template_ratio_last_dominates,
  template_ratio_first_two_10
)

# Varying parameters
template_amounts <- c(100, 1000)  # total template amount per sample
degradation_amounts <- c(0.0025, 0.007)  # total degradation amount per sample

# Constants
c2 <- 14.754
k2BackStutter <- 14.48
k2ForwardStutter <- 9.67
k22bpBackStutter <- 3.13
k2DoubleBackStutter <- 6.97

# Randomly sampled parameters will be generated inside the loop



# Number of genotypes per configuration
n_samples <- 1
replicates <- 2
contributors_list <- c(2, 3, 4, 5)


# Must loop over the following elements:
#   - contributors_list
#   - template_amounts
#   - template_ratio_functions
#   - degradation_amounts
#   - n_samples
#   - replicates

# Main simulation loop - Contributors
for (n_contributors in contributors_list) {
    contributors <- paste0("U", seq_len(n_contributors))

  # Template amounts and ratios
  for (template_amount in template_amounts) {
    for (ratio_idx in seq_along(template_ratio_functions)) {
      ratio_func <- template_ratio_functions[[ratio_idx]]
      
      template_ratios <- ratio_func(n_contributors)
      templates <- template_amount * template_ratios

      # Degradation amounts
      for (degradation_idx in seq_along(degradation_amounts)) {
        degradation_value <- degradation_amounts[degradation_idx]
        # Per sample
        for (sample_idx in 1:n_samples) {

          # TODO: What if we generate genotypes in outermost loop? Will the models trained on it be worse?
          genotypes <- sample_contributor_genotypes(
            contributors = contributors,
            loci = gf$log_normal_settings$locus_names,
            freqs = allele_freqs
          )
          
          # Sample LSAE for this sample
          # TODO: Inside or outside replicate loop?
          LSAE <- sample_LSAE(LSAE_variance = gf$log_normal_settings$LSAE_variance_prior,
                             locus_names = gf$log_normal_settings$locus_names)

          # Per replicate
          for (rep in 1:replicates) {
            set.seed(1e7 * n_contributors + 1e5 * template_amount + 1e4 * ratio_idx + 1e3 * sample_idx + 1e2 * degradation_idx + rep)

            # Create parameter summary for this sample (single row, all contributors as columns)
            parameter_summary <- data.frame(
              SampleName = paste0(
                "nC", n_contributors, "_tmp", template_amount, "_ratio", ratio_idx,
                "_deg", degradation_idx, "_s", sample_idx, "_rep", rep
              ),
              # Contributors
              setNames(as.list(contributors), paste0("contributors", seq_len(n_contributors))),
              model = "log_normal_model",
              # Templates
              setNames(as.list(templates), paste0("template", seq_len(n_contributors))),
              # Degradation
              setNames(as.list(rep(degradation_value, n_contributors)), paste0("degradation", seq_len(n_contributors))),
              # c2, k2BackStutter, k2ForwardStutter, etc. (single values)
              c2 = c2,
              k2BackStutter = k2BackStutter,
              k2ForwardStutter = k2ForwardStutter,
              k22bpBackStutter = k22bpBackStutter,
              k2DoubleBackStutter = k2DoubleBackStutter
            )

            # Add LSAE columns
            for (locus in gf$log_normal_settings$locus_names) {
              parameter_summary[[locus]] <- LSAE[locus]
            }

            # Generated mixture
            mixture <- sample_mixtures_fixed_parameters(
              genotypes = genotypes,
              parameter_summary = parameter_summary,
              model_settings = gf$log_normal_settings
            )

            sim_peaks <- mixture$samples[[1]]$annotated_mixture
            sim_peaks_with_dye <- left_join(sim_peaks, dye_map, by = c("Locus" = "Marker"))
            combined_peaks_df <- rbind(sim_peaks_with_dye, size_standard_df)
            combined_peaks_df <- combined_peaks_df %>%
              mutate(Scan = round(Size * 11.2 + 3500))

            out_file <- sprintf(
              "%s/epg_nC%d_tmp%d_ratio%d_deg%d_s%d_rep%d.csv",
              out_dir, n_contributors, template_amount, ratio_idx, degradation_idx, sample_idx, rep
            )
            write.csv(combined_peaks_df, out_file, row.names = FALSE)
          }
        }
      }
    }
  }
}
