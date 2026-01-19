"""Generate gyoza input for FKS1-HS1 orthologs.

This script is composed of a quick import of orthologous sequences
followed by part of the gyoza annotation module to quickly retrieve
all silent mutants.

"""

import pandas as pd


def get_mutations(seq, wt, codon_dic):
    r"""Collect differences between a mutated DNA sequence and the wild-type.

    Parameters
    ----------
    seq : str
        DNA sequence of the mutant (with bases either ``A``, ``C``, ``G`` or ``T``).
        Length should be the same as `wt` and be a multiple of 3.
    wt : str
        Wild-type DNA sequence (with bases either ``A``, ``C``, ``G`` or ``T``).
        Length should be the same as `seq` and be a multiple of 3.
    codon_dic : dict
        Codon table associating codons to amino acid residues.

    Returns
    -------
    is_wt : bool
        ``True`` if `seq` is `wt`
    aa_seq : str
        Protein sequence translated from `seq`
    Nham_codons : int
        Number of codon changes
    Nham_nt : int
        Number of nucleotide changes
    Nham_aa : int
        Number of amino acid changes
    nonsilent_pos : list
        Positions in the protein sequence of each non-silent mutation
    nonsilent_alt_aa : list
        Alternative residue in the protein sequence for each non-silent mutation
    nonsilent_wt_aa : list
        Wild-type residue in the protein sequence at the position of
        each non-silent mutation
    mutated_codon : list
        1-based indexes for each mutated codon (nth mutated codon)
    mutation_pos : list
        Positions in `seq` of each mutation
    mutation_alt_codon : list
        Alternative codon in `seq` for each mutation
    mutation_alt_aa : list
        Alternative residue translated from `mutation_alt_codon`
    mutation_type : list
        Either ``silent``, ``missense`` or ``nonsense`` based on ``mutation_alt_aa``

    Raises
    ------
    ValueError
        If the length of `seq` is not a multiple of 3.
    ValueError
        If the lengths of `seq` and `wt` are not equal.
    ValueError
        If `seq` contains unrecognized characters.

    Notes
    -----
    Mutations are formatted as # mutated codon / position / alternative codon /
    alternative amino acid, in lists with matching indexes to be able
    to quickly convert to 1 row per mutation per mutated codon.

    The alternative and corresponding wild-type codons are translated into
    their corresponding amino acid using the `codon_dic`.

    Sequence-level attributes include the Hamming distances (``Nham``),
    i.e. the number of codon, nucleotide and amino acid changes.
    """
    if len(seq) % 3 != 0:
        raise ValueError(
            f"Error.. the length of the DNA sequence is not a multiple of 3."
        )

    # Note: the two following checks are validated early on
    # but let's keep them just in case

    if len(seq) != len(wt):
        raise ValueError(
            f"Error.. Cannot annotate expected mutants because at least one sequence is of different length than wild-type."
        )

    if not set(seq).issubset({"A", "C", "G", "T"}):
        raise ValueError(
            f"Error.. one of the provided nucleotide sequences contains an unrecognized character."
        )

    is_wt = seq == wt

    mutation_pos, mutation_alt_codon, mutation_alt_aa, mutation_type = [], [], [], []
    full_aa_seq, full_wt_aa = [], []
    Nham_nt = 0
    Nham_aa = 0

    wt_codons = [
        wt[i : i + 3]
        for i in range(
            0, len(wt), 3
        )  # Converting WT nucleotide sequence to list of codons
    ]
    seq_codons = [
        seq[i : i + 3]
        for i in range(
            0, len(seq), 3
        )  # Converting nucleotide sequence of variant to list of codons
    ]

    for i, (wtc, c) in enumerate(zip(wt_codons, seq_codons)):  # Loop through codons
        wt_aa = codon_dic.get(wtc)
        alt_aa = codon_dic.get(c)
        full_aa_seq.append(alt_aa)
        full_wt_aa.append(wt_aa)

        if c != wtc:
            mutation_pos.append(i)
            mutation_alt_codon.append(c)
            mutation_alt_aa.append(alt_aa)
            Nham_nt += sum(x != y for x, y in zip(wtc, c))
            if alt_aa != wt_aa:
                Nham_aa += 1
                mutation_type.append("nonsense" if alt_aa == "*" else "missense")
            else:
                mutation_type.append("silent")

    # Protein-level comparison
    aa_seq = "".join(full_aa_seq)
    wt_aa_seq = "".join(full_wt_aa)

    nonsilent_pos, nonsilent_alt_aa, nonsilent_wt_aa = [], [], []

    for i, (wt_aa, alt_aa) in enumerate(zip(wt_aa_seq, aa_seq)):
        if wt_aa != alt_aa:
            nonsilent_pos.append(int(i)),
            nonsilent_alt_aa.append(alt_aa)
            nonsilent_wt_aa.append(wt_aa)

    Nham_codons = len(mutation_pos)

    if Nham_codons > 0:
        mutated_codon = list(
            range(1, Nham_codons + 1)
        )  # 1-based index of mutated codons
        if Nham_aa == 0:
            nonsilent_alt_aa = ["not-applicable"]
            nonsilent_wt_aa = ["not-applicable"]
    else:  # Handle WT (no mutations)
        mutated_codon = [0]
        mutation_pos = ["not-applicable"]
        mutation_alt_codon = ["not-applicable"]
        mutation_alt_aa = ["not-applicable"]
        mutation_type = ["wt"]
        nonsilent_alt_aa = ["not-applicable"]
        nonsilent_wt_aa = ["not-applicable"]

    return (
        is_wt,
        aa_seq,
        Nham_codons,
        Nham_nt,
        Nham_aa,
        nonsilent_pos,
        nonsilent_alt_aa,
        nonsilent_wt_aa,
        mutated_codon,
        mutation_pos,
        mutation_alt_codon,
        mutation_alt_aa,
        mutation_type,
    )


def annotate_mutants(df, codon_dic):
    r"""Annotate a dataframe of mutated DNA sequences with mutations.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns: ``nt_seq`` and ``WT_seq``
    codon_dic : dict
        Codon table associating codons to amino acid residues.

    Returns
    -------
    pandas.DataFrame
        `df` with additional columns to describe mutations.

    Notes
    -----
    Uses custom function to collect mutations by comparing each sequence
    to the corresponding wild-type on the matching row.
    """
    per_seq_cols = [
        "WT",
        "aa_seq",
        "Nham_codons",
        "Nham_nt",
        "Nham_aa",
        "pos",
        "alt_aa",
        "wt_aa",
    ]
    per_mut_cols = [
        "mutated_codon",
        "mutation_pos",
        "mutation_alt_codons",
        "mutation_alt_aa",
        "mutation_type",
    ]
    new_cols = per_seq_cols + per_mut_cols

    # Making sure sequences are capitalized
    df["nt_seq"] = df["nt_seq"].str.upper()
    df["WT_seq"] = df["WT_seq"].str.upper()

    collected_mutations = [
        get_mutations(seq, wt, codon_dic) for seq, wt in zip(df["nt_seq"], df["WT_seq"])
    ]

    mutations_dict = dict(zip(new_cols, zip(*collected_mutations)))
    df = df.assign(**mutations_dict)
    df = df.explode(per_mut_cols).reset_index(drop=True)

    return df


def generate_expected_mutants(locus):
    """Generate list of expected mutants to be used as gyoza input.

    Parameters
    ----------
    locus : str
        Unique identifier for the mutated locus.

    Returns
    -------
    pandas.DataFrame
    """
    # Import codon dictionary
    codon_table = pd.read_csv("../../config/project_files/codon_table.csv", header=0)
    codon_table["aminoacid"] = codon_table["aminoacid"].str.upper()
    codon_table["codon"] = codon_table["codon"].str.upper()
    codon_dic = dict(zip(codon_table["codon"], codon_table["aminoacid"]))

    # Import WT sequence
    # Make sure HS3 WT features in the file (quick hack)
    wt_df = pd.read_csv("../../config/project_files/wt_seq.csv")
    wt = wt_df.loc[wt_df.Mutated_seq == locus, "WT_seq"].values[0].upper()

    # Import orthologous sequences
    df = pd.read_csv(f"../ortholog_data/{locus}_ortholog_nt.csv")
    df["nt_seq"] = df["nt_seq"].astype(str).str.upper()
    df["Mutated_seq"] = locus
    df["WT_seq"] = wt

    # Import and annotate silent mutants
    # Note: the file is generated by running gyoza on the dataset of single mutants
    # with --no-temp --until generate_mutants
    # also needs --forcerun generate_mutants if gyoza has already run
    # (quick hack for HS3 - copy file in correct folder)
    all_mut = pd.read_csv(f"../../config/project_files/expected_mut/{locus}.csv.gz")[
        ["Mutated_seq", "WT_seq", "nt_seq"]
    ]
    annotated = annotate_mutants(all_mut, codon_dic)
    silent = annotated[annotated.Nham_aa == 0][["Mutated_seq", "WT_seq", "nt_seq"]]

    return pd.concat([df, silent])


# Export in proper format for gyoza
for l in ["FKS1-HS1", "FKS1-HS2", "FKS1-HS3", "FKS2-HS1", "FKS2-HS2"]:
    df = generate_expected_mutants(l)
    df.to_csv(f"../ortholog_data/{l}.csv.gz", index=False)

# The files are then moved in the appropriate folder to start gyoza with the 'provided' mode
