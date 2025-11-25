import pandas as pd
df = pd.read_csv("../FKS1-HS1_ortholog_nt.csv")
df["nt_seq"] = df["nt_seq"].astype(str).str.upper()
df["Mutated_seq"] = "FKS1-HS1"
df["WT_seq"] = "TTTTTAGTTTTATCTTTGAGAGATCCA"
df.to_csv("../FKS1-HS1.csv.gz", index=False)