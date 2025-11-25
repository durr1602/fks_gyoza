import pandas as pd

df = pd.read_csv("../../results/df/avg_scores_HS1_ortho.csv").rename(
    columns={
        "Mutated_seq": "locus",
        "fitness_T2": "s",
        "lower_err_T2": "min_s",
        "upper_err_T2": "max_s",
    }
)
GMM = (
    df[(df.strain == "BY4741") & (df.locus == "FKS1-HS1")]
    .groupby(["compound", "aa_seq"])[["s", "min_s", "max_s"]]
    .agg(s=("s", "first"), min_s=("min_s", "first"), max_s=("max_s", "first"))
    .reset_index()
)
GMM.to_csv("../classified/BY4741_FKS1-HS1/refined_classification_ortho.csv")
