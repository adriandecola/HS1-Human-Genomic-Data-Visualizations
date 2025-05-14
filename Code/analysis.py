"""
Author: Adrian deCola
Date Created: May 12, 2025
Title: analysis.py 
Purpose: To visualise 6‑mer counts for HS1 vs HG38
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


# ────────────────────────── Constants and Configs ────────────────────────── #
### Constants 
MITO_WEIGHT = 100  # each mt-DNA hit counts as 100 nuclear hits
RELATIVE_CHANGE_TO_BE_OUTLIAR = 1.7 # Only considers increases
# number of entries for most and least frequent sixmers tables/txt file
X_MOST_AND_LEAST_COMMON = 15

### File Paths
# Path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up one level (from Code/ to Sixmers Code/)
PARENT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))

# Now point into the Sixmer Data folder
DATA_DIR = os.path.join(PARENT_DIR, "Sixmer Data")

# Directory to save the ouputed figured to
FIG_DIR = os.path.join(PARENT_DIR, "Figures")
os.makedirs(FIG_DIR, exist_ok=True) # makes figures directory if required

FILES = {
    "HG38_raw":          os.path.join(DATA_DIR, "HG38.txt"),
    "HS1_raw":           os.path.join(DATA_DIR, "HS1.txt"),
    "HG38_mito_raw":     os.path.join(DATA_DIR, "HG38_Mitochondria.txt"),
    "HS1_mito_raw":      os.path.join(DATA_DIR, "HS1_Mitochondria.txt"),
}



# ──────────────────────────── Helper Functions ────────────────────────────── #
def read_kmer_file(path: str) -> pd.DataFrame:
    """Load AAAAAA,27 format to DataFrame[kmer,count]."""
    data = []
    with open(path) as fileHandle:
        for line in fileHandle:
            kmer, count = line.rstrip().split(',') # remove whitespace and split into list
            data.append((kmer, int(count))) 
    return pd.DataFrame(data, columns=["kmer", "count"])

def revcomp(kmer: str) -> str:
    """Return the reverse complement of a DNA 6-mer."""
    return kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ──────────────────────────── main ───────────────────────────────── #
def main() -> None:
    #####################################
    #      Loading the Data             #
    #####################################
    
    ### Raw Genome with no correction for mitochondrial DNA frequencies ###
    hs1_raw_df  = read_kmer_file(FILES["HS1_raw"])
    hg38_raw_df = read_kmer_file(FILES["HG38_raw"])
    # Merging the dataframes
    genome_df = pd.merge(
        hs1_raw_df,
        hg38_raw_df,
        on="kmer",
        suffixes=("_HS1_raw", "_HG38_raw"),
    )
    
    ### Mitochondria Genome ###
    hs1_mito_raw_df  = read_kmer_file(FILES["HS1_mito_raw"])
    hg38_mito_raw_df = read_kmer_file(FILES["HG38_mito_raw"])
    # Merging the dataframes
    mito_df = pd.merge(
        hs1_mito_raw_df,
        hg38_mito_raw_df,
        on="kmer",
        suffixes=("_HS1_raw", "_HG38_raw"),
    )

    ####################################################################
    #     Creating apply mito-DNA correction for a “blood-like” mix    #
    ####################################################################
    corrected_df = pd.DataFrame({
    "kmer": genome_df["kmer"],
    "count_HS1_corr": genome_df["count_HS1_raw"]
                       + (MITO_WEIGHT-1) * mito_df["count_HS1_raw"],
    "count_HG38_corr": genome_df["count_HG38_raw"]
                        + (MITO_WEIGHT-1) * mito_df["count_HG38_raw"],
    })

    #######################################################################
    #   Getting the 10 most and least popular sixmers (HS1 corrected)     #
    #######################################################################
    # Prepare outfile
    out_path = os.path.join(FIG_DIR, "Most and Least Frequent Sixmers.txt")
    with open(out_path, "w") as f:
        # sort by corrected HS1 count
        sorted_hs1 = corrected_df.sort_values("count_HS1_corr")
    
        least_common = sorted_hs1.head(X_MOST_AND_LEAST_COMMON)
        most_common  = sorted_hs1.tail(X_MOST_AND_LEAST_COMMON)

        ### print bottom {X_MOST_AND_LEAST_COMMON} ###
        # print title
        f.write(f"\n{X_MOST_AND_LEAST_COMMON} least abundant 6-mers (HS1 corrected):\n")
        # bold header
        f.write("  " + "K-mer         Frequency\n")
        # dashed line
        f.write("  " + "-" * 30 + "\n")
        for kmer, cnt in zip(least_common["kmer"], least_common["count_HS1_corr"]):
            f.write(f"  {kmer:<6s}  {cnt:>12,d}\n")

        ### print top 10 (reverse so biggest first) ###
        # print title
        f.write(f"\n{X_MOST_AND_LEAST_COMMON} most abundant 6-mers (HS1 corrected):\n")
        # bold header
        f.write("  " + "K-mer         Frequency\n")
        # dashed line
        f.write("  " + "-" * 30 + "\n")
        for kmer, cnt in zip(most_common["kmer"].iloc[::-1], most_common["count_HS1_corr"].iloc[::-1]):
            f.write(f"  {kmer:<6s}  {cnt:>12,d}\n")

    print(f"✓   Top/bottom sixmers written to: {out_path}")


    ################################################
    #   Plotting the raw frequency differences     #
    ################################################
    fig, ax = plt.subplots()

    # Axis Labels
    ax.set_xlabel("HS1_raw 6‑mer count")
    ax.set_ylabel("HG38_raw 6‑mer count")

    # Title
    ax.set_title("Raw 6‑mer frequencies: HS1 vs HG38 (nuclear genome)")

    # Compute the data range
    min_val = genome_df[["count_HS1_raw","count_HG38_raw"]].values.min()  # smallest 6-mer count
    max_val = genome_df[["count_HS1_raw","count_HG38_raw"]].values.max()  # largest 6-mer count

    # Scatter raw HS1 vs HG38 points
    ax.scatter(genome_df["count_HS1_raw"], genome_df["count_HG38_raw"],
               s=6, alpha=0.5)  # small, semi-transparent dots

    # Log–log scaling
    ax.set_xscale("log")  # x axis from min_val to max_val on log scale
    ax.set_yscale("log")  # y axis from min_val to max_val on log scale

    # Force both axes to start at the same minimum
    ax.set_xlim(min_val*.8, max_val*1.2)  # x from min to max with some buffer
    ax.set_ylim(min_val*.8, max_val*1.2)  # y from min to max with some buffer

    # Draw y = x reference line over the full plot
    ax.plot([min_val*.8, max_val*1.2], [min_val*.8, max_val*1.2],
            ls="--", lw=1)  # diagonal

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "Raw Genome Comparison Plot.png"), dpi=1000)
    plt.close(fig)

    print("✓   Raw Frequency Plot Saved")

    ################################################
    #   Plotting the corrected frequencies         #
    ################################################
    fig, ax = plt.subplots()
    
    # Axis Labels
    ax.set_xlabel("HS1 corrected count")
    ax.set_ylabel("HG38 corrected count")

    # Title
    ax.set_title("6-mer frequencies with 100× mito correction")

    # scatter corrected HS1 vs HG38
    ax.scatter(corrected_df["count_HS1_corr"],
           corrected_df["count_HG38_corr"],
           s=6, alpha=0.5)  

    # Compute the data range
    min_val = corrected_df[["count_HS1_corr","count_HG38_corr"]].values.min()  # smallest 6-mer count
    max_val = corrected_df[["count_HS1_corr","count_HG38_corr"]].values.max()  # largest 6-mer count

    # Log–log scaling
    ax.set_xscale("log")  # x axis from min_val to max_val on log scale
    ax.set_yscale("log")  # y axis from min_val to max_val on log scale

    # Force both axes to start at the same minimum
    ax.set_xlim(min_val*.8, max_val*1.2)  # x from min to max with some buffer
    ax.set_ylim(min_val*.8, max_val*1.2)  # y from min to max with some buffer

    # Draw y = x reference line over the full plot
    ax.plot([min_val*.8, max_val*1.2], [min_val*.8, max_val*1.2],
            ls="--", lw=1)  # diagonal

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "Corrected Genome Comparison Plot.png"), dpi=1000)
    plt.close(fig)

    print("✓   Mitochondria Corrected Frequency Plot Saved")

    ################################################
    #   Plotting GC-content colored scatter        #
    ################################################
    # Compute GC count (0–6) for each 6-mer
    corrected_df["gc_count"] = corrected_df["kmer"].apply(
        lambda k: k.count("C") + k.count("G")
    )

    fig, ax = plt.subplots()

    # Plotting points with color by GC count
    sc = ax.scatter(
        corrected_df["count_HS1_corr"],      # x = HS1 corrected
        corrected_df["count_HG38_corr"],     # y = HG38 corrected
        c=corrected_df["gc_count"],          # color by GC count
        cmap="viridis",                      # perceptually uniform colormap
        s=7,                                 # still small points
        alpha=0.7,                            # good visibility
        edgecolors='none',                   # edges same as face
        linewidths=0.2                       # thin outline
    )

    # log–log axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # labels & title
    ax.set_xlabel("HS1 corrected count")
    ax.set_ylabel("HG38 corrected count")
    ax.set_title("Corrected counts coloured by GC content (0–6)")

    # diag reference from min→max
    min_val = corrected_df[["count_HS1_corr","count_HG38_corr"]].values.min()
    max_val = corrected_df[["count_HS1_corr","count_HG38_corr"]].values.max()

    # Force both axes to start at the same minimum
    ax.set_xlim(min_val*.8, max_val*1.2)  # x from min to max with some buffer
    ax.set_ylim(min_val*.8, max_val*1.2)  # y from min to max with some buffer


    # Draw y = x reference line over the full plot
    ax.plot(
        [min_val*0.8, max_val*1.2],
        [min_val*0.8, max_val*1.2],
        color="C1",
        linestyle="--",
        linewidth=1,
        label="x = y"
    )
    
    ### Draw y = {RELATIVE_CHANGE_TO_BE_OUTLIAR}* x trend line in a second color
    # Build a smooth line on log axes for x = 1.2 * y
    y_vals = np.logspace(
        np.log10(min_val * 0.8),
        np.log10(max_val * 1.2),
        200
    )
    x_vals = RELATIVE_CHANGE_TO_BE_OUTLIAR * y_vals

    ax.plot(
        x_vals,
        y_vals,
        color="C3",                # contrasting default color
        linestyle="-.",            # dotted-dash for distinction
        linewidth=1,
        label=f"x = {RELATIVE_CHANGE_TO_BE_OUTLIAR}y"
    )

    # Add a legend to distinguish the two lines
    ax.legend(loc="upper left", frameon=False, fontsize="small")
    
    # add colourbar legend
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Number of G or C in 6-mer", rotation=270, labelpad=15)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "GC-Color Plot Comparison.png"), dpi=300)
    plt.close(fig)

    print("✓   GC-color plot saved")


    ##############################################
    #   Write sixmers with large relative change #
    ##############################################
    # Add column for ratio of frequency in HS1 to frequency in HG38
    corrected_df["ratio"]   = corrected_df["count_HS1_corr"] / corrected_df["count_HG38_corr"]
    # Add column for percentage increase from HG38 to HS1
    corrected_df["pct_inc"] = (corrected_df["ratio"] - 1.0) * 100

    # select outliers where ratio ≥ RELATIVE_CHANGE_TO_BE_OUTLIAR
    outliers = corrected_df[corrected_df["ratio"] >= RELATIVE_CHANGE_TO_BE_OUTLIAR] \
                   .sort_values("ratio", ascending=False)

    # prepare output file name
    fname = f"Sixmers that change by {round((RELATIVE_CHANGE_TO_BE_OUTLIAR-1)*100)} Percent.txt"
    out_path = os.path.join(FIG_DIR, fname)

    # write table
    with open(out_path, "w") as f:
        # title
        f.write(f"Sixmers with HS1/HG38 ≥ {RELATIVE_CHANGE_TO_BE_OUTLIAR}"
                f"  (≥ {(RELATIVE_CHANGE_TO_BE_OUTLIAR-1)*100:.1f}% increase)\n")
        f.write("-" * 60 + "\n")
        # header row
        f.write(f"{'6-mer':<10s}{'HG38 Frequency':>18s}{'HS1 Frequency':>18s}{'% increase':>14s}\n")
        f.write("-" * 60 + "\n")
        # each outlier
        for _, row in outliers.iterrows():
            k = row["kmer"]
            h = int(row["count_HG38_corr"])
            s = int(row["count_HS1_corr"])
            p = row["pct_inc"]
            f.write(f"{k:<10s}{h:>18,d}{s:>18,d}{p:>14.1f}%\n")

    print(f"✓   Outliers table written to: {out_path}")

    #########################################################
    #   Plotting sixmer and its reverse complement together #
    #########################################################
    # add revcomplement column
    corrected_df["revcomp"]    = corrected_df["kmer"].apply(revcomp)

    # Prepare reverse counts
    rev_df = corrected_df[["kmer", "count_HS1_corr", "count_HG38_corr"]].rename(
        columns={
            "kmer": "revcomp",
            "count_HS1_corr": "revcomp_HS1_corr",
            "count_HG38_corr": "revcomp_HG38_corr"
        }
    )


    # Merge to bring in revcomp counts
    corrected_df = corrected_df.merge(rev_df, on="revcomp", how="left")

    corrected_df["sum_HS1_corr"] = corrected_df["count_HS1_corr"] + corrected_df["revcomp_HS1_corr"]
    corrected_df["sum_HG38_corr"] = corrected_df["count_HG38_corr"] + corrected_df["revcomp_HG38_corr"]

    # get cananical identifier
    corrected_df["canon"] = corrected_df.apply(
        lambda r: min(r["kmer"], r["revcomp"]), axis=1
    )

    # drop duplicates, keeping just one of each canon pair
    uniq_df = corrected_df.drop_duplicates("canon").copy()

    # compute GC count for each canonical 6-mer
    uniq_df["gc_count"] = uniq_df["kmer"].apply(
        lambda k: k.count("C") + k.count("G")
    )


    # plotting
    fig, ax = plt.subplots()

    sc = ax.scatter(
        uniq_df["sum_HS1_corr"],  
        uniq_df["sum_HG38_corr"],
        c=uniq_df["gc_count"],
        cmap="viridis",
        s=7,
        alpha=0.7,
        edgecolors='none',
        linewidths=0.2
    )

    # log-log axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # labels & title
    ax.set_xlabel("HS1 corrected count")
    ax.set_ylabel("HG38 corrected count")
    ax.set_title("Reduced (2048) 6-mer counts colored by GC content")

    # limits and trend lines
    min_val = uniq_df[["sum_HS1_corr", "sum_HG38_corr"]].values.min()
    max_val = uniq_df[["sum_HS1_corr", "sum_HG38_corr"]].values.max()

    ax.set_xlim(min_val * 0.8, max_val * 1.2)
    ax.set_ylim(min_val * 0.8, max_val * 1.2)

    ax.plot(
        [min_val * 0.8, max_val * 1.2],
        [min_val * 0.8, max_val * 1.2],
        color="C1", linestyle="--", linewidth=1, label="x = y"
    )

    y_vals = np.logspace(np.log10(min_val * 0.8), np.log10(max_val * 1.2), 200)
    x_vals = RELATIVE_CHANGE_TO_BE_OUTLIAR * y_vals
    ax.plot(
        x_vals,
        y_vals,
        color="C3", linestyle="-.", linewidth=1,
        label=f"x = {RELATIVE_CHANGE_TO_BE_OUTLIAR}y"
    )

    ax.legend(loc="upper left", frameon=False, fontsize="small")

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Number of G or C in 6-mer", rotation=270, labelpad=15)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "Canonical GC-Color Plot Comparison.png"), dpi=300)
    plt.close(fig)

    print("✓   Conanical GC-color plot saved")

    ########################################################
    # Interactive version of canonical GC-color plot       #
    ########################################################
    print("✓   Creating interactive Plotly GC scatter")

    df = uniq_df.copy()
    df["log_sum_HS1"] = np.log10(df["sum_HS1_corr"])
    df["log_sum_HG38"] = np.log10(df["sum_HG38_corr"])

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df["sum_HS1_corr"],
        y=df["sum_HG38_corr"],
        mode='markers',
        marker=dict(
            color=df["gc_count"],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="GC Count"),
            size=5,
            opacity=0.7,
        ),
        text=df["canon"],
        hovertemplate="K-mer: %{text}<br>HS1: %{x}<br>HG38: %{y}<extra></extra>"
    ))

    fig.update_layout(
        title="Interactive Canonical 6-mer Plot (log scale)",
        xaxis_title="HS1 corrected count",
        yaxis_title="HG38 corrected count",
        xaxis_type="log",
        yaxis_type="log",
        height=700,
    )

    # Save as interactive HTML
    html_path = os.path.join(FIG_DIR, "Interactive Canonical Scatter.html")
    fig.write_html(html_path)
    print(f"✓   Interactive plot saved: {html_path}")

    ##############################################
    #   Plotly Interactive Scatter for Canonical  #
    ##############################################
    print("✓   Creating interactive draggable matplotlib plot")

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)

    # Scatter plot
    sc = ax.scatter(
        uniq_df["sum_HS1_corr"], uniq_df["sum_HG38_corr"],
        c=uniq_df["gc_count"], cmap="viridis", s=7, alpha=0.7,
        edgecolors='none', linewidths=0.2
    )

    # Log–log axes
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("HS1 corrected count")
    ax.set_ylabel("HG38 corrected count")
    ax.set_title("Canonical 6-mer GC plot (interactive)")

    # Reference lines
    min_val = uniq_df[["sum_HS1_corr", "sum_HG38_corr"]].values.min() * 0.8
    max_val = uniq_df[["sum_HS1_corr", "sum_HG38_corr"]].values.max() * 1.2
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)

    # x = y line
    ax.plot([min_val, max_val], [min_val, max_val], ls="--", lw=1, color="gray", label="x = y")
    ax.legend(loc="lower right", fontsize="small")

    # Draggable threshold line
    threshold = uniq_df["sum_HS1_corr"].median()
    vline = ax.axvline(threshold, color='red', linestyle='--')

    # Percent annotation
    total = len(uniq_df)
    left = (uniq_df["sum_HS1_corr"] < threshold).sum()
    right = total - left
    text = ax.text(
        0.05, 0.95,
        f"< {threshold:.0f}: {left/total:.1%}\n≥ {threshold:.0f}: {right/total:.1%}",
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray")
    )

    # Dragging logic
    # Initial threshold tracker
    dragging = {"active": False}
    current_thresh = {"val": threshold}

    def on_press(event):
        if event.inaxes != ax or event.xdata is None:
            return
        if abs(np.log10(event.xdata) - np.log10(current_thresh["val"])) < 0.05:
            dragging["active"] = True

    def on_release(event):
        dragging["active"] = False

    def on_motion(event):
        if dragging["active"] and event.inaxes == ax and event.xdata:
            thresh = event.xdata
            current_thresh["val"] = thresh
            vline.set_xdata([thresh, thresh])
            left = (uniq_df["sum_HS1_corr"] < thresh).sum()
            right = total - left
            text.set_text(f"< {thresh:.0f}: {left/total:.1%}\n≥ {thresh:.0f}: {right/total:.1%}")
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("button_release_event", on_release)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)

    plt.show() 


    ##################################################
    #   Write table of all 2048 canonical 6-mers     #
    ##################################################
    out_path = os.path.join(
        FIG_DIR, f"Canonical Most and Least Frequent Sixmers.txt"
    )

    with open(out_path, "w") as f:
        f.write(f"{X_MOST_AND_LEAST_COMMON} least and most abundant canonical 6-mers (HS1 mitochondria weighted)\n")
        f.write("-" * 60 + "\n")

        # Sort by HS1 corrected count
        sorted_uniq = uniq_df.sort_values("sum_HS1_corr")

        # Get bottom N and top N (reversed)
        least = sorted_uniq.head(X_MOST_AND_LEAST_COMMON)
        most  = sorted_uniq.tail(X_MOST_AND_LEAST_COMMON).iloc[::-1]

        # Bottom section
        f.write(f"\nBottom {X_MOST_AND_LEAST_COMMON} canonical 6-mers:\n")
        f.write("\n")
        f.write(f"{'CANONICAL / REVERSE_COMPLEMENT':<36s}{'Frequency (HS1)':>16s}\n")
        f.write("-" * 52 + "\n")
        for _, row in least.iterrows():
            pair = f"{row['kmer']}/{row['revcomp']}"
            f.write(f"{pair:<36s}{row['sum_HS1_corr']:>16,d}\n")

        # Top section
        f.write(f"\nTop {X_MOST_AND_LEAST_COMMON} canonical 6-mers:\n")
        f.write("\n")
        f.write(f"{'CANONICAL / REVERSE_COMPLEMENT':<36s}{'Frequency (HS1)':>16s}\n")
        f.write("-" * 52 + "\n")
        for _, row in most.iterrows():
            pair = f"{row['kmer']}/{row['revcomp']}"
            f.write(f"{pair:<36s}{row['sum_HS1_corr']:>16,d}\n")

    print(f"✓   Canonical top/bottom sixmers table saved: {out_path}")

    ##############################################
    # Write canonical sixmers with large increase #
    ##############################################
    uniq_df["ratio"] = uniq_df["sum_HS1_corr"] / uniq_df["sum_HG38_corr"]
    uniq_df["pct_inc"] = (uniq_df["ratio"] - 1.0) * 100

    # select canonical outliers
    canon_outliers = uniq_df[uniq_df["ratio"] >= RELATIVE_CHANGE_TO_BE_OUTLIAR] \
                         .sort_values("ratio", ascending=False)

    # write canonical outlier table
    fname = f"Canonical Sixmers that change by {round((RELATIVE_CHANGE_TO_BE_OUTLIAR-1)*100)} Percent.txt"
    out_path = os.path.join(FIG_DIR, fname)

    with open(out_path, "w") as f:
        f.write(f"Canonical sixmers with HS1/HG38 ≥ {RELATIVE_CHANGE_TO_BE_OUTLIAR}"
                f"  (≥ {(RELATIVE_CHANGE_TO_BE_OUTLIAR-1)*100:.1f}% increase)\n")
        f.write("-" * 75 + "\n")
        f.write(f"{'CANONICAL / REVERSE_COMPLEMENT':<36s}"
                f"{'HG38 Freq':>12s}{'HS1 Freq':>12s}{'% Increase':>14s}\n")
        f.write("-" * 75 + "\n")
        for _, row in canon_outliers.iterrows():
            pair = f"{row['kmer']}/{row['revcomp']}"
            h = int(row["sum_HG38_corr"])
            s = int(row["sum_HS1_corr"])
            p = row["pct_inc"]
            f.write(f"{pair:<36s}{h:>12,d}{s:>12,d}{p:>14.1f}%\n")

    print(f"✓   Canonical outliers table written to: {out_path}")

    ##############################################
    #   Write canonical sixmers txt file         #
    ##############################################

    out_path = os.path.join(FIG_DIR, "Canonical Sixmer Frequencies (HS1 Corrected).txt")
    with open(out_path, "w") as f:
        for _, row in uniq_df.iterrows():
            pair = f"{row['kmer']}/{row['revcomp']}"
            count = row["sum_HS1_corr"]
            f.write(f"{pair},{count}\n")

    print(f"✓   Canonical frequency list written to: {out_path}")




# only run if executed directly
if __name__ == "__main__":
    main()