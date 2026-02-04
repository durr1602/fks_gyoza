"""Plot correlation between different backgrounds."""

import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"


# Specify paths and constants
s_path = "../../results/df/avg_scores.csv"
outpath = "../graphs/"
loci = ["FKS1-HS1", "FKS1-HS2"]
compounds = ["None", "Anidulafungin", "Caspofungin", "Micafungin"]

# Import data
df = pd.read_csv(s_path)
df["compound"] = df["compound"].astype(str).str.title()
wide = (
    df[(df.Mutated_seq.isin(loci)) & (df.compound.isin(compounds)) & (df.Nham_aa == 1)]
    .pivot(
        index=["Mutated_seq", "compound", "aa_pos", "alt_aa"],
        columns="strain",
        values="fitness_T2",
    )
    .reset_index()
)

# Plot
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

g = sns.FacetGrid(
    wide,
    row="Mutated_seq",
    row_order=loci,
    col="compound",
    col_order=compounds,
    margin_titles=True,
    height=2,
)

# Map scatter
g.map_dataframe(
    sns.scatterplot,
    x="BY4741",
    y="R1158",
    s=8,
    alpha=0.2,
    color=sns.color_palette("mako", 1),
    edgecolor="none",
)

# Spearman correlation
for i, l in enumerate(loci):
    for j, c in enumerate(compounds):
        subset = wide[(wide.Mutated_seq == l) & (wide.compound == c)].dropna(
            subset=["BY4741", "R1158"]
        )
        sr, sp = stats.spearmanr(subset.BY4741, subset.R1158)
        g.axes[i][j].text(
            0.5,
            0.05,
            rf"$\rho$ = {sr:.2f}" + "\n$\it{p}$-val = " + f"{sp:.1e}",
            transform=g.axes[i][j].transAxes,
            fontsize=8,
        )

g.set_axis_labels("s ($\it{FKS2}$)", "s ($\it{fks2}$)")
g.set_titles(row_template="{row_name}", col_template="{col_name}")

plt.tight_layout()
plt.savefig(f"{outpath}/BYvR11.svg", format="svg", dpi=300)
plt.savefig(f"{outpath}/BYvR11.png", format="png", dpi=300)
