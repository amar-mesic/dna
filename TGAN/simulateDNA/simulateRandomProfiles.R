# -----------------------------------------
# Required libraries
# -----------------------------------------
library(simDNAmixtures)
library(dplyr)
# library(readr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
n_contributors <- ifelse(length(args) >= 1, as.integer(args[[1]]), 2)
n_mixtures <- ifelse(length(args) >= 2, as.integer(args[[2]]), 3)


# dir.create("generated_alleles", showWarnings = FALSE)

# Load example allele frequency file that comes with the package
allele_freqs_file <- system.file("extdata","FBI_extended_Cauc_022024.csv",
                                       package = "simDNAmixtures")

# Read allele frequencies to create an empirical distribution
allele_freqs <- read_allele_freqs(allele_freqs_file)

# Not too sure what the point of this is
# seems to be a configuration for a specific kit
gf <- gf_configuration()

# When sampling a genotype, what kind of EPG do we want?
# We can randomzie the template ratio, as well as the degradation
sampling_parameters <- list(min_template = 50., max_template = 10000.,
                            degradation_shape = 2.5, degradation_scale = 1e-3)

# Find the corresponding dyes for each marker in the globalfiler kit
dye_map <- kits$GlobalFiler[, c("Marker", "Color")]
dye_map <- dye_map %>%
  distinct(Marker, .keep_all = TRUE)

# add the internal lane standard
# Define the fragment sizes
size_standard_sizes <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                         220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380,
                         400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560,
                         580, 600)

# Create a DataFrame for the size standard peaks
size_standard_df <- data.frame(
  Locus = "LIZ",
  Allele = NA,
  # Assign a consistent peak height; adjust later if randomness needed
  Height = 3000,
  Size = size_standard_sizes,
  Color = "orange"
)

for (i in 1:n_mixtures) {
  set.seed(i)
  contributors <- paste0("U", seq_len(n_contributors))

  # Generate mixture
  # With this method, we assume that all heights above 75 are true alleles, and those below are false.
  # That means we can expect artefacts to be included in the genotypes.
  # This is not good enough, but serves as a starting point.
  # Bonus with this approach: AMEL (sex) peaks are already included
  mixtures <- sample_mixtures(
    n = n_contributors,
    contributors = contributors,
    freqs = allele_freqs,
    sampling_parameters = sampling_parameters,
    model_settings = gf$log_normal_settings,
    sample_model = sample_log_normal_model,
    results_directory = "."
  )

  # Extract peaks
  # Contains the data of the alleles
  sim_peaks <- mixtures$samples[[1]]$mixture

  # Join dye info
  # Join the table of the peaks with the dye map
  # This will add the color of the dye to the peaks
  sim_peaks_with_dye <- left_join(sim_peaks, dye_map, by = c("Locus" = "Marker"))

  # Add size standard
  # Combine with the existing simulated peaks DataFrame
  combined_peaks_df <- rbind(sim_peaks_with_dye, size_standard_df)

  # Add scan column
  # Add a column representing scan point, which is deduced using formula and current size
  # The values of 11.2 and 3500 are the linear scaling factors that are hard coded from the paper
  combined_peaks_df <- combined_peaks_df %>%
    mutate(Scan = round(Size * 11.2 + 3500))

  # Write to file
  out_file <- sprintf("generated_alleles/simulated_epg_with_dyes_seed_%d.csv", i)
  write.csv(combined_peaks_df, out_file, row.names = FALSE)
}
