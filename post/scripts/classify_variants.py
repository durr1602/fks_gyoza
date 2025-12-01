"""
Classify variants based on score meeting thresholds.
"""

import pathlib
import argparse
import sys
import pandas as pd
from custom_functions import refine_class


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", type=str, required=True, help="Path to input data basename"
    )
    parser.add_argument("-o", type=str, required=True, help="Outpath basename")
    parser.add_argument("--strain", type=str, required=True, help="Strain identifier")
    parser.add_argument("--locus", type=str, required=True, help="Locus identifier")
    args = parser.parse_args()

    # Input
    input_path = pathlib.Path(args.i)
    if not input_path.is_absolute():
        data_path = f'../classified/{"_".join([args.strain,args.locus])}/{args.i}'
    else:
        data_path = args.i

    thresholds = f'../classified/{"_".join([args.strain,args.locus])}/thresholds.csv'

    # Output
    classified_outpath = f'../classified/{"_".join([args.strain,args.locus])}/{args.o}'

    # Import data
    df = pd.read_csv(data_path).rename(columns={"fitness_T2": "s"})
    stddf = pd.read_csv(thresholds)

    # Apply thresholds
    df["refined_class"] = df.apply(
        lambda row: refine_class(row.s, row.compound, stddf), axis=1
    )

    # Collapse classification
    df["sensres"] = df.refined_class.replace(
        {
            "intermediary": "resistant",
            "WT-like": "sensitive",
            "slightly deleterious": "sensitive",
            "deleterious": "deleterious",
        }
    )

    # Export classified data
    df.to_csv(f"{classified_outpath}", index=False)

    return


if __name__ == "__main__":
    # When running a script with %run, it might sometimes try to parse Jupyter's
    # internal arguments. The next two lines handle this more gracefully.
    if "__file__" not in globals() and "ipykernel" in sys.modules:
        # We are likely running in Jupyter via %run
        # The sys.argv will contain the arguments after the script name
        # We just rely on argparse to handle the rest.
        pass

    main()
