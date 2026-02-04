#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 13:35:03 2024
Major update 13-01-2026
@author: aliciapageau
"""

# %% Load packages
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# %% Paths and parameters
wkdir = "/home/rodur28/dev/durr1602/fks_gyoza/pre/fungamr/"
os.chdir(wkdir)
date = "160126"  # DDMMYY


# %% Functions
def mutate_sequence(sequence, position, alt_aa, mut_type):
    # Convert position to integer
    position = int(position)

    # Create a Biopython Seq object
    seq = Seq(sequence)

    # Mutate the sequence at the specified position
    if mut_type == "SNPs":
        mutated_seq = seq[: position - 1] + alt_aa + seq[position:]
    elif mut_type == "Indels":
        if alt_aa == "del":
            mutated_seq = seq[: position - 1] + seq[position:]
        elif alt_aa.startswith("ins"):
            insert = alt_aa.replace("ins", "")
            mutated_seq = seq[:position] + insert + seq[position:]
    # elif mut_type == 'stop':
    #     mutated_seq = seq[:position-1] + alt_aa

    return str(mutated_seq)


def dataframe_to_multifasta(df, output_file):
    counter = 1
    fasta_records = []
    for idx, row in df.iterrows():
        numerical_id = counter
        species = row["species"]
        gene_name = row["ortho_homolog"]
        uniprot_id = row["uniprot"]
        mutation = row["mutation"]
        sequence = str(row["mutated_seq"])

        # Creating a Bio.SeqRecord object
        seq = Seq(sequence)
        record = SeqRecord(
            seq,
            id=f"{numerical_id}",
            description=f"| gene_name={gene_name} | species={species} | uniprot={uniprot_id} | mutation={mutation}",
        )

        fasta_records.append(record)
        counter += 1

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")


def get_positions(alignment_file):
    """
    Reads an alignment file, retrieves positions and corresponding indices
    for each sequence, and returns the result.
    My positions are NOT zero-based (first position = 1)

    Args:
        alignment_file (str): The path to the input alignment file.

    Returns:
        dict: A dictionary containing positions and indices for each sequence in the alignment.
    """
    alignment = AlignIO.read(alignment_file, "fasta")

    positions = {}
    for record in alignment:
        name = record.description
        sequence = str(record.seq)
        positions[name] = {}

        ungapped_index = 1
        for alignment_index, char in enumerate(sequence):
            positions[name][alignment_index + 1] = (char, "NA")
            if char != "-":
                positions[name][alignment_index + 1] = (char, ungapped_index)
                ungapped_index += 1

    return positions


def position_to_df(alignment_file):
    """
    Converts the position dictionary obtained from an alignment file into a long format dataframe.

    Args:
        alignment_file (str): The path to the input alignment file.

    Returns:
        pandas.DataFrame: A dataframe representing the positions and corresponding information
        for each sequence in the alignment.
    """
    positions = get_positions(alignment_file)

    df_positions = pd.DataFrame.from_dict(
        {
            (name, alignment_index): positions[name][alignment_index]
            for name in positions.keys()
            for alignment_index in positions[name].keys()
        },
        orient="index",
        columns=["aa", "Ungapped position"],
    )

    # Resetting index and splitting 'info' column into separate columns
    df_positions = df_positions.reset_index(names=["info"])
    df_positions[["info", "Alignment position"]] = df_positions["info"].apply(
        lambda x: pd.Series([x[0], x[1]])
    )
    df_positions[["ID", "ortho_homolog", "species", "uniprot", "mutation"]] = (
        df_positions["info"].str.split("|", expand=True)
    )

    # Extracting values after '=' in each column and filling missing values
    df_positions["ortho_homolog"] = (
        df_positions["ortho_homolog"].str.split("=").str[1].str.strip().fillna("NA")
    )
    df_positions["species"] = (
        df_positions["species"].str.split("=").str[1].str.strip().fillna("NA")
    )
    df_positions["uniprot"] = (
        df_positions["uniprot"].str.split("=").str[1].str.strip().fillna("NA")
    )
    df_positions["mutation"] = (
        df_positions["mutation"].str.split("=").str[1].str.strip().fillna("NA")
    )

    # Dropping the original 'info' column
    df_positions.drop(columns=["info", "ID"], inplace=True)

    return df_positions


def get_wild_type_hotspots(row, df):
    """
    Function to retrieve the wild-type hotspot sequences based on species
    """
    wild_type_row = df[
        (df["species"] == row["species"])
        & (df["mutation"] == "Wild-type")
        & (df["ortho_homolog"] == row["ortho_homolog"])
        & (df["Hotspot"] == row["Hotspot"])
    ]
    if not wild_type_row.empty:
        wild_type_hotspot = wild_type_row["aa"].iloc[0]
        return wild_type_hotspot
    else:
        return None


def update_position(row, alignment_positions):
    alignment_pos = row["alignment_pos"]
    gene = row["gene or protein name"]
    specie = row["conversion_species"]  # Retrieve species from added column

    if pd.notna(alignment_pos):
        convert_pos = alignment_positions.loc[
            (alignment_positions["gene or protein name"] == gene)
            & (alignment_positions["alignment_pos"] == alignment_pos)
            & (alignment_positions["species"] == specie),
            "position",
        ]

        if not convert_pos.empty and not pd.isna(convert_pos.iloc[0]):
            return convert_pos.iloc[0]
        else:
            return None

    else:
        return None


def update_wt(row, alignment_positions):
    alignment_pos = row["alignment_pos"]
    gene = row["gene or protein name"]
    specie = row["conversion_species"]  # Retrieve species from added column

    if pd.notna(alignment_pos):
        convert_wt = alignment_positions.loc[
            (alignment_positions["gene or protein name"] == gene)
            & (alignment_positions["alignment_pos"] == alignment_pos)
            & (alignment_positions["species"] == specie),
            "wt_AA",
        ]

        if not convert_wt.empty and not pd.isna(convert_wt.iloc[0]):
            return convert_wt.iloc[0]
        else:
            return None
    else:
        return None


# %% Import files
sequences = pd.read_csv(f"{wkdir}/ref_seq_030425.csv", index_col=0)
database = pd.read_csv(f"{wkdir}FungAMR_070425.csv", index_col=0)

# Extract only Fks1 and Fks2 single mutation entries in fungAMR database
fks = database[database["gene or protein name"].str.contains("Fks1|Fks2")]
fks = fks[fks["mutation_composition"] == "single"]
fks = fks[fks["first author name"] != "Durand"]  # Remove mut from Romain's DMS

# Remove frameshift mutations and wild-type
fks = fks[~fks["mutation"].str.contains("fs|Wild-type")]
fks_mut = fks[
    ["species", "gene or protein name", "ortho_homolog", "mutation", "mutation_type"]
].drop_duplicates()

# Split mutation in wt AA, position and alt AA
fks_mut[["wt_AA", "position", "alt_AA"]] = fks_mut["mutation"].str.extract(
    r"(\D+)(\d+)(\D+)"
)

# Fix insertion cases: extract first number as position, AA after "ins" as alt_AA
ins_mask = fks_mut["mutation"].str.contains("ins")
fks_mut.loc[ins_mask, "position"] = fks_mut.loc[ins_mask, "mutation"].str.extract(
    r"(\d+)"
)[0]
fks_mut.loc[ins_mask, "alt_AA"] = fks_mut.loc[ins_mask, "mutation"].str.extract(
    r"(ins[A-Za-z]+)"
)[0]


# %% Mutate FKS sequences with mutation from fungAMR database
fks_mut = pd.merge(
    fks_mut,
    sequences,
    on=["species", "ortho_homolog", "gene or protein name"],
    how="left",
)
fks_mut = fks_mut.dropna(subset=["sequence"])
fks_mut["mutated_seq"] = fks_mut.apply(
    lambda row: mutate_sequence(
        row["sequence"], row["position"], row["alt_AA"], row["mutation_type"]
    ),
    axis=1,
)

# Add wild-type orthologues sequences
fks_wt = sequences[
    (sequences["gene or protein name"] == "Fks1")
    | (sequences["gene or protein name"] == "Fks2")
].drop_duplicates()
fks_wt = fks_wt.dropna(subset=["sequence"])  # No ref seq for S.cer x S.par hybrid
fks_wt["mutation"] = "Wild-type"

fks_mut = pd.concat([fks_mut, fks_wt]).reset_index(drop=True)
fks_mut.loc[fks_mut["mutation"] == "Wild-type", "mutated_seq"] = fks_mut["sequence"]

# Split Fks1 and Fks2 for alignment
fks1 = fks_mut[fks_mut["gene or protein name"] == "Fks1"]
fks2 = fks_mut[fks_mut["gene or protein name"] == "Fks2"]

# %% Align mutated sequence
"""
I turn my mutations df into a fasta file and align sequences with muscle5.
$ muscle5 -align Fks1_mutated_{date}.fasta -output Fks1_mutated_aligned_{date}.fasta
I take resulting MSA files and convert it back to a df while extracting positions.
"""

fks1_mutated_msa = f"{wkdir}/Fks1_mutated_{date}.fasta"
fks2_mutated_msa = f"{wkdir}/Fks2_mutated_{date}.fasta"
dataframe_to_multifasta(df=fks1, output_file=fks1_mutated_msa)
dataframe_to_multifasta(df=fks2, output_file=fks2_mutated_msa)

fks1_mutated_align = f"{wkdir}/Fks1_mutated_aligned_{date}.fasta"
fks2_mutated_align = f"{wkdir}/Fks2_mutated_aligned_{date}.fasta"
fks1_mutated_positions = position_to_df(fks1_mutated_align)
fks2_mutated_positions = position_to_df(fks2_mutated_align)

# %% Extract FKS hotspots
Sc_fks1_hs1_start = 639
Sc_fks1_hs1_end = 647
Sc_fks1_hs2_start = 1353
Sc_fks1_hs2_end = 1360
Sc_fks1_hs3_start = 690
Sc_fks1_hs3_end = 700

Sc_fks2_hs1_start = 658
Sc_fks2_hs1_end = 666
Sc_fks2_hs2_start = 1372
Sc_fks2_hs2_end = 1379
Sc_fks2_hs3_start = 709
Sc_fks2_hs3_end = 719


def extract_hotspot_sequences(
    df, species, ortho, start_pos, end_pos, hotspot_label, ortho_label
):
    # Get alignment positions for start and end
    start = df[
        (df["species"] == species)
        & (df["mutation"] == "Wild-type")
        & (df["ortho_homolog"] == "Fks")
        & (df["Ungapped position"] == start_pos)
    ]["Alignment position"].values[0]

    end = df[
        (df["species"] == species)
        & (df["mutation"] == "Wild-type")
        & (df["ortho_homolog"] == "Fks")
        & (df["Ungapped position"] == end_pos)
    ]["Alignment position"].values[0]

    # Filter and aggregate
    hotspot_df = df[df["Alignment position"].between(start, end)].drop_duplicates()
    hotspot_df = (
        hotspot_df.groupby(["species", "ortho_homolog", "uniprot", "mutation"])
        .agg({"aa": "".join, "Ungapped position": lambda x: x.tolist()})
        .reset_index()
    )

    # Rename and add metadata
    hotspot_df = hotspot_df.rename(columns={"Ungapped position": "positions"})
    hotspot_df["Hotspot"] = hotspot_label
    hotspot_df["ortho_homolog"] = ortho_label

    return hotspot_df


# Fks1
Fks1_HS1_mutated_seq = extract_hotspot_sequences(
    fks1_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks1_hs1_start,
    Sc_fks1_hs1_end,
    "HS1",
    "Fks1",
)
Fks1_HS2_mutated_seq = extract_hotspot_sequences(
    fks1_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks1_hs2_start,
    Sc_fks1_hs2_end,
    "HS2",
    "Fks1",
)
Fks1_HS3_mutated_seq = extract_hotspot_sequences(
    fks1_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks1_hs3_start,
    Sc_fks1_hs3_end,
    "HS3",
    "Fks1",
)

# Fks2
Fks2_HS1_mutated_seq = extract_hotspot_sequences(
    fks2_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks2_hs1_start,
    Sc_fks2_hs1_end,
    "HS1",
    "Fks2",
)
Fks2_HS2_mutated_seq = extract_hotspot_sequences(
    fks2_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks2_hs2_start,
    Sc_fks2_hs2_end,
    "HS2",
    "Fks2",
)
Fks2_HS3_mutated_seq = extract_hotspot_sequences(
    fks2_mutated_positions,
    "Saccharomyces cerevisiae",
    "Fks",
    Sc_fks2_hs3_start,
    Sc_fks2_hs3_end,
    "HS3",
    "Fks2",
)

# Combine all mutated seq
Fks1_mutated_seq = pd.concat(
    [Fks1_HS1_mutated_seq, Fks1_HS2_mutated_seq, Fks1_HS3_mutated_seq]
)
Fks2_mutated_seq = pd.concat(
    [Fks2_HS1_mutated_seq, Fks2_HS2_mutated_seq, Fks2_HS3_mutated_seq]
)

Fks_mutated_seq = pd.concat([Fks1_mutated_seq, Fks2_mutated_seq])

# %% Filter out mutation not in hotspot sequences (same as WT)
# Get the wild-type hotspot sequences for each species
Fks_mutated_seq["wild_type_hotspot"] = Fks_mutated_seq.apply(
    lambda row: get_wild_type_hotspots(row, Fks_mutated_seq), axis=1
)

# Filter out rows where both hotspot sequences are the same as the wild-type hotspot sequences
Fks_mutated_seq_filtered = Fks_mutated_seq[
    Fks_mutated_seq["aa"] != Fks_mutated_seq["wild_type_hotspot"]
]
Fks_mutated_seq_filtered = Fks_mutated_seq_filtered.drop(columns="wild_type_hotspot")

# %% Found mutation resitance profil to echinochandin
fks_metadata = fks[
    [
        "species",
        "gene or protein name",
        "drug",
        "mutation",
        "pubmedid",
        "confidence score",
        "ortho_mut",
        "MIC",
        "Strongest resistance evidence reported",
        "Strongest sensitivity evidence reported",
    ]
]

fks_mutated_final = pd.merge(
    Fks_mutated_seq_filtered,
    fks_metadata,
    left_on=["species", "ortho_homolog", "mutation"],
    right_on=["species", "gene or protein name", "mutation"],
    how="left",
)

# %% Convert column "mutation" to Scer numbering
alignment_positions = pd.read_csv(
    f"{wkdir}/master_positions_df_030425.csv", index_col=0
).drop(columns={"uniprot"})

fks_mutated_final["ortho_homolog"] = "Fks"
fks_mutated_final[["wt_AA", "position", "alt_AA"]] = fks_mutated_final[
    "mutation"
].str.extract(r"(\D+)(\d+)(\D+)")

fks_mutated_final["position"] = pd.to_numeric(fks_mutated_final["position"])
fks_mutated_final_aligned = pd.merge(
    fks_mutated_final,
    alignment_positions,
    on=["species", "ortho_homolog", "gene or protein name", "position", "wt_AA"],
    how="left",
)

fks_mutated_final_aligned["conversion_species"] = "Saccharomyces cerevisiae"

fks_mutated_final_aligned["convert_pos"] = fks_mutated_final_aligned.apply(
    update_position, axis=1, alignment_positions=alignment_positions
)
fks_mutated_final_aligned["convert_wt"] = fks_mutated_final_aligned.apply(
    update_wt, axis=1, alignment_positions=alignment_positions
)

fks_mutated_final_aligned["Scer_mutation"] = (
    fks_mutated_final_aligned["convert_wt"]
    + fks_mutated_final_aligned["convert_pos"].astype("Int64").astype(str)
    + fks_mutated_final_aligned["alt_AA"]
)
fks_mutated_final_aligned.drop(
    columns=[
        "wt_AA",
        "position",
        "alt_AA",
        "conversion_species",
        "convert_pos",
        "convert_wt",
        "alignment_pos",
        "accession",
    ],
    inplace=True,
)
fks_mutated_final_aligned.to_csv(
    f"{wkdir}Fks_mutated_hotspot_seq_{date}.csv"
)  # S1 Data after a few adjustments

# %% Aggregate
fks_mutated_final_aligned = pd.read_csv(
    "/home/rodur28/dev/durr1602/fks_gyoza/pre/fungamr/Fks_mutated_hotspot_seq_160126.csv",
    index_col=0,
)
gby = (
    fks_mutated_final_aligned.groupby(["Hotspot", "Scer_mutation", "drug"])[
        [
            "pubmedid",
            "species",
            "aa",
            "mutation",
            "Strongest resistance evidence reported",
            "Strongest sensitivity evidence reported",
        ]
    ]
    .agg(
        pubmed_ids=("pubmedid", "unique"),
        species=("species", "unique"),
        mutations=("mutation", "unique"),
        aa=("aa", "unique"),
        best_res=("Strongest resistance evidence reported", "min"),
        best_sens=("Strongest sensitivity evidence reported", "max"),
    )
    .reset_index()
)


# %% Label resistant/sensitive
def annotate_phenotype(res_score, sens_score):
    if np.isnan(res_score):
        return "sensitive"
    elif np.isnan(sens_score):
        return "resistant"
    else:
        if res_score < -1 * sens_score:
            return "resistant"
        else:
            return "sensitive"


gby["phenotype"] = gby.apply(
    lambda row: annotate_phenotype(row.best_res, row.best_sens), axis=1
)

gby[
    gby.drug.isin(["Caspofungin", "Micafungin", "Anidulafungin", "Echinocandin"])
].reset_index(drop=True).to_csv(f"{wkdir}fungamrmut_df.csv", index=None)
