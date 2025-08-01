library(simDNAmixtures)
library(dplyr)

# Output directory
dir.create("generated_alleles_mass", showWarnings = FALSE)

# Load allele frequencies
allele_freqs_file <- system.file("extdata","FBI_extended_Cauc_022024.csv", package = "simDNAmixtures")
allele_freqs <- read_allele_freqs(allele_freqs_file)
gf <- gf_configuration()
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

# Parameter grids
template_ranges <- list(
  c(50, 500),
  c(5000, 10000),
  c(50, 10000)
)
degradation_settings <- list(
  list(shape = 2.5, scale = 1e-3),
  list(shape = 3.5, scale = 2e-3)
)
contributors_list <- c(2, 3, 4, 5)
replicates <- 2

# Calculate number of genotypes per configuration
n_per_config <- 2
config_id <- 1

for (n_contributors in contributors_list) {
  for (deg_idx in seq_along(degradation_settings)) {
    for (tmpl_idx in seq_along(template_ranges)) {
      for (genotype in 1:n_per_config) {
        contributors <- paste0("U", seq_len(n_contributors))
        for (rep in 1:replicates) {
          set.seed(100000*config_id + 1000*genotype + rep) # Unique seed per config/genotype/rep

          # Set parameters for this config
          sampling_parameters <- list(
            min_template = template_ranges[[tmpl_idx]][1],
            max_template = template_ranges[[tmpl_idx]][2],
            degradation_shape = degradation_settings[[deg_idx]]$shape,
            degradation_scale = degradation_settings[[deg_idx]]$scale
          )

          mixtures <- sample_mixtures(
            n = n_contributors,
            contributors = contributors,
            freqs = allele_freqs,
            sampling_parameters = sampling_parameters,
            model_settings = gf$log_normal_settings,
            sample_model = sample_log_normal_model,
            number_of_replicates = replicates
          )

          sim_peaks <- mixtures$samples[[1]]$mixture
          sim_peaks_with_dye <- left_join(sim_peaks, dye_map, by = c("Locus" = "Marker"))
          combined_peaks_df <- rbind(sim_peaks_with_dye, size_standard_df)
          combined_peaks_df <- combined_peaks_df %>%
            mutate(Scan = round(Size * 11.2 + 3500))

          out_file <- sprintf(
            "generated_alleles_mass/epg_nC%d_deg%d_tmp%d_geno%04d_rep%d.csv",
            n_contributors, deg_idx, tmpl_idx, genotype, rep
          )
          write.csv(combined_peaks_df, out_file, row.names = FALSE)
        }
      }
      config_id <- config_id + 1
    }
  }
}