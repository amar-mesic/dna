import argparse
from tensorflow.keras.models import load_model
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import sys

# Set working directory to the script's directory
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

def generate_random_latent(batch_size=1, n_lanes=6, length=500):
    # Normal distribution noise per lane, per position
    # return np.random.normal(loc=0.0, scale=0.0, size=(batch_size, n_lanes, length, 1))
    # generate uniformly distributed noise
    return np.random.uniform(low=0, high=1.0, size=(batch_size, n_lanes, length, 1))

# Generate a batch of clean EPGs
def epg_from_df(df, n_lanes=6, epg_length=5000, scan_min=4000, std_dev=4) -> np.ndarray:
    scan_max = scan_min + epg_length
    color_map = {"blue": 0, "green": 1, "yellow": 2, "red": 3, "purple": 4, "orange": 5}
    epg = np.zeros((1, n_lanes, epg_length, 1))
    for _, row in df.iterrows():
        channel = color_map.get(row["Color"])
        center = int(round(row["Scan"]))
        height = row["Height"]
        start = max(center - 4 * std_dev, scan_min)
        end = min(center + 4 * std_dev, scan_max)
        x = np.arange(start, end)
        gaussian = height * np.exp(-0.5 * ((x - center) / std_dev) ** 2)
        epg[0, channel, x - scan_min, 0] += gaussian
    return epg

# Plotting function for EPG
# def plot_epg(epg: np.ndarray, no_dyes=6, title="EPG", save_path='epg', include_lane_standard=True):
#     lane_colors = ["magenta", "blue", "cyan", "green", "yellow", "orange"]
#     lanes = no_dyes if include_lane_standard else no_dyes - 1
#     fig, axes = plt.subplots(lanes, 1, figsize=(14, 2.5 * lanes), sharex=True)
#     if lanes == 1:
#         axes = [axes]
#     for i in range(lanes):
#         axes[i].plot(epg[0, i, :, 0], color=lane_colors[i])
#         axes[i].set_ylabel("Intensity")
#         axes[i].grid(True)
#     axes[-1].set_xlabel("Scan Points")
#     fig.suptitle(title, fontsize=16)
#     plt.tight_layout()
#     plt.savefig(f"{save_path}.pdf")
#     plt.close(fig)


def plot_epg_plotly(epg: np.ndarray, no_dyes=6, title="EPG", save_path="epg_plot",
                    include_lane_standard=True, save_static_pdf=True):
    """
    Plot electropherogram using Plotly with optional saving as HTML and static image.

    Parameters:
        epg (np.ndarray): EPG array of shape (1, no_dyes, scan_points, 1)
        no_dyes (int): Number of dye channels (default: 6)
        title (str): Title of the plot
        save_path (str): Path (without extension) for saving output
        include_lane_standard (bool): Whether to include internal lane standard channel
        save_static_pdf (bool): Whether to save a static PDF version (requires kaleido)
    """
    lane_colors = ["magenta", "blue", "cyan", "green", "yellow", "orange"]
    lanes = no_dyes if include_lane_standard else no_dyes - 1

    fig = make_subplots(rows=lanes, cols=1, shared_xaxes=True,
                        vertical_spacing=0.02, subplot_titles=[f"Dye Channel {i+1}" for i in range(lanes)])
    
    for i in range(lanes):
        signal = epg[0, i, :, 0]
        fig.add_trace(
            go.Scatter(
                x=np.arange(len(signal)),
                y=signal,
                mode='lines',
                line=dict(color=lane_colors[i], width=1),
                name=f"Dye {i+1}"
            ),
            row=i + 1,
            col=1
        )
        fig.update_yaxes(title_text="Intensity (RFU)", row=i + 1, col=1)

    fig.update_xaxes(title_text="Scan Points", row=lanes, col=1)
    fig.update_layout(height=250 * lanes, width=1000,
                      title_text=title,
                      hovermode="x unified",
                      showlegend=False)
    
    # Show the interactive figure
    # fig.show()

    # Save as interactive HTML
    fig.write_html(f"{save_path}.html")

# Main function to generate EPGs from CSV files and save them as PDFs
def main(csv_dir, output_dir, batch_size=1, epg_length=5000, scan_min=4000):
    generator = load_model('../../TGAN_model/Profile_generator_GAN_trained_model_random_24-09-2023.h5')
    csv_files = sorted(glob.glob(os.path.join(csv_dir, "*.csv")))
    n_files = len(csv_files)
    if n_files == 0:
        print(f"No CSV files found in {csv_dir}")
        return
    dfs = [pd.read_csv(f) for f in csv_files]
    clean_epgs = np.concatenate([
        epg_from_df(df, epg_length=epg_length, scan_min=scan_min) for df in dfs
    ], axis=0)
    latent = generate_random_latent(batch_size=n_files, n_lanes=6, length=500)
    for start in range(0, n_files, batch_size):
        end = min(start + batch_size, n_files)
        batch_clean_epgs = clean_epgs[start:end]
        batch_latent = latent[start:end]
        batch_csv_files = csv_files[start:end]
        batch_generated_epgs = generator.predict([batch_clean_epgs, batch_latent], batch_size=(end-start))
        for _, (csv_file, gen_epg) in enumerate(zip(batch_csv_files, batch_generated_epgs)):
            base = os.path.splitext(os.path.basename(csv_file))[0]
            plot_epg_plotly(
                gen_epg[np.newaxis, ...],
                no_dyes=6,
                title=base,
                save_path=os.path.join(output_dir, base),
                include_lane_standard=True
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate EPG PDFs from CSV files.")
    parser.add_argument("--csv_dir", type=str, help="Directory containing input CSV files.")
    parser.add_argument("--output_dir", type=str, help="Directory to save output EPG PDFs.")
    parser.add_argument("--batch_size", type=int, default=8, help="Batch size for generator model.")
    parser.add_argument("--epg_shape", type=int, nargs=2, metavar=('SCAN_MIN', 'EPG_LENGTH'),
                        default=[0, 9641],
                        # default=[4000, 5000],
                        help="EPG length and scan_min as two integers, e.g. --epg_shape 5000 4000")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    scan_min, epg_length = args.epg_shape
    main(csv_dir=args.csv_dir, output_dir=args.output_dir, batch_size=args.batch_size, epg_length=epg_length, scan_min=scan_min)

    # example command:
    # python Generator/generatePDFs.py --csv_dir='../generated/generated_alleles_fixed_ratios_base_template_500_5000_improved/epgs' --output_dir='../../generated_epgs/fixed_ratios_base_template_500_5000_improved_pdfs' --batch_size=64 --epg_shape 4000 5000