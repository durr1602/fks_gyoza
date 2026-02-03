"""Plot correlation between replicates.

This script is adapted from gyoza v1.2.0 scripts.plot_scores.py
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations

plt.rcParams["svg.fonttype"] = "none"


def concatenate_df(df_files):
    r"""Opens and concatenates multiple dataframes.

    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.

    Returns
    -------
    pandas.DataFrame
        Single concatenated dataframe.
    """
    list_df = []
    for f in df_files:
        list_df.append(pd.read_csv(f))
    df = pd.concat(list_df, ignore_index=True)
    return df


def plot_replicate_scatter(df, spearman_df, outpath):
    r"""Plot correlation between first two replicates, for each sample group.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores. Should contain columns:

        * ``Sample attributes`` (**str**, sample group identifier)
        * ``Compared timepoints`` (**str**)
        * ``r1`` (**str**, first replicate)
        * ``r2`` (**str**, second replicate)

    spearman_df : pandas.DataFrame
        Dataframe of Spearman correlation coeffs and p-values, with:

        * ``Sample attributes`` (**str**, sample group identifier)
        * ``Compared timepoints`` (**str**)
        * ``rho`` (**float**)
        * ``pval`` (**float**)

    outpath : str
        Path to save scatter plot as SVG (should end with ``.svg``).
    """
    sns.set_theme(
        rc={
            "font.family": "Arial",
            "font.size": 8,
            "legend.title_fontsize": 8,
            "legend.fontsize": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.major.pad": 2,
            "ytick.major.pad": 2,
            "xtick.bottom": True,
            "ytick.left": True,
            "xtick.major.size": 2,
            "ytick.major.size": 2,
        },
        style="ticks",
    )

    # Identify unique pairs to set up consistent color palette
    unique_pairs = sorted(df["pair"].unique())
    palette_colors = sns.color_palette("mako", len(unique_pairs))
    color_map = dict(zip(unique_pairs, palette_colors))

    # Identify unique row/col values for consistent order
    l_order = sorted(df["locus"].unique())
    c_order = sorted(df["compound"].unique())

    g = sns.FacetGrid(
        df,
        row="locus",
        row_order=l_order,
        col="compound",
        col_order=c_order,
        hue="pair",
        palette=color_map,
        margin_titles=True,
        height=1.8,
    )

    # Map scatter
    g.map_dataframe(sns.scatterplot, x="x", y="y", s=8, alpha=0.2, edgecolor="none")

    # Annotation loop
    for (row_val, col_val), ax in g.axes_dict.items():
        sub_stats = spearman_df[
            (spearman_df["locus"] == row_val) & (spearman_df["compound"] == col_val)
        ]

        for i, (_, row) in enumerate(sub_stats.iterrows()):
            pair_name = row["pair"]
            ax.text(
                0.05,
                0.92 - (i * 0.08),
                f"{pair_name}: ρ={row['rho']:.2f}",
                transform=ax.transAxes,
                fontsize=8,
                color=color_map[pair_name],  # Matches the scatter color
                fontweight="bold",
            )

    g.set_axis_labels("s (Rep. X)", "s (Rep. Y)")
    g.set_titles(row_template="{row_name}", col_template="{col_name}")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.savefig(outpath.replace(".svg", ".png"), dpi=300)

    return


def main(df_files, outpath):
    r"""Aggregate data and plot correlation between replicates.

    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.
    outpath : str
        Path to save scatter plots as SVG (should end with ``.svg``).
    """
    df = concatenate_df(df_files)
    df["Replicate"] = df["Replicate"].astype(str)

    # Reshape dataframe
    repwide = df.pivot(
        index=["Sample attributes", "Nham_aa", "aa_seq"],
        columns="Replicate",
        values="s",
    ).reset_index()

    # Extract metadata from sample attributes
    split_data = repwide["Sample attributes"].str.split("__")
    repwide["locus"] = split_data.str[-2]
    repwide["compound"] = split_data.str[-1].str.title()

    # Restrict to non-control conditions
    repwide = repwide[
        repwide["compound"].isin(["Anidulafungin", "Caspofungin", "Micafungin"])
    ]

    # Get combinations of replicates
    replicates = ["1", "2", "3"]
    pairs = list(combinations(replicates, 2))

    plot_data_list = []
    stats_list = []

    # Process each pair
    for r_x, r_y in pairs:
        pair_label = f"{r_x} vs {r_y}"
        temp = repwide[["locus", "compound", r_x, r_y]].dropna(subset=[r_x, r_y])

        for (locus, compound), group in temp.groupby(["locus", "compound"]):
            if len(group) > 1:
                rho, pval = stats.spearmanr(group[r_x], group[r_y])
                stats_list.append(
                    {
                        "locus": locus,
                        "compound": compound,
                        "pair": pair_label,
                        "rho": rho,
                        "pval": pval,
                    }
                )
        temp = temp.rename(columns={r_x: "x", r_y: "y"})
        temp["pair"] = pair_label
        plot_data_list.append(temp)

    full_plot_df = pd.concat(plot_data_list)
    spearman_df = pd.DataFrame(stats_list)

    plot_replicate_scatter(full_plot_df, spearman_df, outpath)

    return


if __name__ == "__main__":
    import os

    filedir = "../../results/df/agg_aa/"
    main(
        [os.path.join(filedir, f) for f in os.listdir(filedir) if "R1158__FKS1__" not in f],
        "../graphs/replicates.svg",
    )
