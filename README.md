# DNA Electropherogram Simulation and Analysis

This repository provides a comprehensive pipeline for simulating, generating, validating, and visualizing synthetic DNA electropherograms (EPGs) using deep learning and statistical methods. It is designed for forensic and bioinformatics research, enabling the creation of realistic synthetic EPGs and their comparison to real data.

## Repository Structure

### `generated_epgs/`
Contains all generated synthetic EPGs and their ground truth data.
- **epgs/**: Numpy arrays representing the simulated EPG signal values.
- **reference_genotypes/**: Genetic information (genotypes) of contributors to each EPG.
- **alleles_to_genotypes_mapping.csv**: Maps each genotype to the EPG(s) it contributes to.

### `preprocessing/`
Notebooks and scripts for testing and developing preprocessing techniques for EPG data.

### `SynthGenDataValidation/`
Scripts and notebooks for validating synthetic data, including calculation of Frechet Inception Distance (FID) between real and synthetic EPGs.

### `TGAN/`
Codebase for simulating DNA profiles and generating synthetic EPGs.
- **simulateDNA/**: Scripts to generate simulated DNA profiles with both genotypic and artefactual peaks.
- **Generator/**: Converts simulated profiles into synthetic EPGs using a conditional GAN (cGAN). Output EPGs are stored in `generated_epgs/`.

### `TGAN_model/`
Pre-trained model weights for the conditional GAN used to generate synthetic EPGs.

### `visualization/`
Jupyter notebooks for visualizing EPGs and exploring their properties.

## Usage

1. **Simulate DNA Profiles**: Use scripts in `TGAN/simulateDNA` to generate DNA profiles.
2. **Generate Synthetic EPGs**: Use the generator in `TGAN/Generator` to convert profiles into EPGs.
3. **Store and Organize Data**: Generated EPGs and their ground truth are saved in `generated_epgs/`.
4. **Preprocess and Validate**: Use notebooks in `preprocessing` and `SynthGenDataValidation` to preprocess and validate synthetic data.
5. **Visualize Results**: Explore and visualize EPGs using notebooks in `visualization`.

## Requirements

- Python 3.8+
- TensorFlow, NumPy, Pandas, Matplotlib, Plotly, and other scientific libraries

## Citation

If you use this repository for your research, please cite appropriately.

## License

See the `LICENSE` file for details.