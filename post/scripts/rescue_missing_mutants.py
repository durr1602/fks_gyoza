"""
Add missing mutants and correct biased estimates.
"""

import pandas as pd

rescued = pd.read_csv("../growth_data/validation_DMS_missing_estimates.csv")
fitness = pd.read_csv("../classified/BY4741_FKS1-HS1/refined_classification.csv")

# Replacing the two DMS datapoints which were overestimated to replace by validation data
# Note: even though values obtained with caspofungin and anidulafungin were very well correlated, we still replace with the inferred score + that way we don't have duplicate values in the df
fitness.drop(fitness[fitness.aa_seq.isin(rescued.aa_seq.unique())].index, inplace=True)
DMS_with_missing = pd.concat([fitness, rescued], ignore_index=True)
DMS_with_missing.to_csv("../classified/BY4741_FKS1-HS1/refined_classification_with_missing.csv", index=False)