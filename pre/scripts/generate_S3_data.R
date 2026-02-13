#### Title: generate_S3_data
####
#### Author: Mathieu Giguere
####
#### Last updated: 2026-02-13
####
#### Pseudocode:
#### 1. Extract Hotspots.
#### 2. Add taxonomy.
#### 3. Combine dataframes.
#### 4. Check if Hotspots are in the oPools.
####
#### Requirements: R4.4.3, dplyr package
#### We recommend using RStudio to reproduce the analysis.
#### Make sure to open a new R project and select the upstream directory
#### as working directory (e.g. {repo-name}/pre/),
#### then open scripts/generate_S3_data.R
#### Input files should be located in {repo-name}/pre/ortholog_data

################################################################################
# Step 1. Extract Hotspots from alignments.
################################################################################

###
### FUNCTION: ExtractionHotspot(aln_path)
###
##
## Find and extract hotspots 1, HS2 and HS3 for each sequences of the FKS genes.
## Returns a dataframe with the hotspot sequences, the species and
## their identifier and other informations.
##
##
## Arguments:
##
## aln_path: path to a .aln file of sequences of 1 Fks gene.
##
##
## Value:
##
## A dataframe of 8 columns describing the species, their identifier,
## the hotspots sequences of these species, the presence of gaps in the hotspots,
## and the presence of the hotspot in the oPools.
##
##
## Example
##
## ExtractionHotspot("FKS1_filtered.aln")
##

ExtractionHotspot <- function(aln_path)
{
  # Validation of the input file.
  stopifnot("Invalid file for analysis. Require a *.aln file" = grepl("\\.aln$", aln_path))

  # Read the file of filtered sequences.
  sequences <- readLines(aln_path)
  
  # Extract and Format identifiers using des REGEX.
  id <- gsub("\\|.*", "", grep(">", sequences, value = TRUE))
  id <- gsub(" .*", "", id)
  id <- gsub(">", "", id)
  
  # Extract and Format the species using REGEX.
  raw_species <- gsub("^.*\\[", "", grep("\\[*\\]", sequences, value = TRUE))
  formated_species <- gsub("\\]", "", raw_species)
  
  # Create data frame.
  df <- data.frame(Identifier = id, Species = formated_species)
  
  # Create hotspot lists.
  hs1_list <- list()
  hs2_list <- list()
  hs3_list <- list()
  
  # Find the sequences corresponding to the Hotspots using their fixed positions.
  # For each header line (so for each sequence)
  for (i in grep(">", sequences))
  {
    # if FKS1
    if (grepl("FKS1", aln_path))
    {
      # Add to the Hotspot list the 8 or 9 amino acids at these fixed positions.
      hs1_list <- append(hs1_list, substr(sequences[i + 17], start = 52, stop = 60))
      hs2_list <- append(hs2_list, substr(sequences[i + 38], start = 1, stop = 8))
      hs3_list <- append(hs3_list, paste(
        substr(sequences[i + 18], start = 49, stop = 60),
        substr(sequences[i + 19], start = 0, stop = 12),
        sep = ""
      ))
    }
    # if FKS2
    else
    {
      # Add to the Hotspot list the 8 or 9 amino acids at these fixed positions.
      # FKS2 Hotspot1 is 1 row farther and 1 aa smaller than the FKS1 Hotspot1.
      # FKS2 Hotspot2 is farther on its line than the FKS1 Hotspot2.
      hs1_list <- append(hs1_list, substr(sequences[i + 19], start = 12, stop = 19))
      hs2_list <- append(hs2_list, paste(
        substr(sequences[i + 38], start = 56, stop = 60),
        substr(sequences[i + 39], start = 0, stop = 3),
        sep = ""
      ))
    }
  }
  
  # if FKS1
  if (grepl("FKS1", aln_path))
  {
    # Add Hotspots columns to the data frame.
    df <- cbind(df, Hotspot1 = unlist(hs1_list))
    df <- cbind(df, Hotspot2 = unlist(hs2_list))
    df <- cbind(df, Hotspot3 = unlist(hs3_list))
    
    ## Adding to the data frame columns indicating gaps in the hotspots.
    GapsHotspot1 <- grepl("-", df$Hotspot1)
    GapsHotspot2 <- grepl("-", df$Hotspot2)
    GapsHotspot3 <- grepl("-", df$Hotspot3)
    df <- cbind(df, GapsHotspot1)
    df <- cbind(df, GapsHotspot2)
    df <- cbind(df, GapsHotspot3)
  }
  # if FKS2
  else
  {
    # Add Hotspots columns to the data frame.
    df <- cbind(df, Hotspot1 = unlist(hs1_list))
    df <- cbind(df, Hotspot2 = unlist(hs2_list))
    
    ## Adding to the data frame columns indicating gaps in the hotspots.
    GapsHotspot1 <- grepl("-", df$Hotspot1)
    GapsHotspot2 <- grepl("-", df$Hotspot2)
    df <- cbind(df, GapsHotspot1)
    df <- cbind(df, GapsHotspot2)
  }
  
  df
}

### Make dataframes with hotspot data.
df_FKS1 <- ExtractionHotspot("ortholog_data/FKS1_filtered.aln")
df_FKS2 <- ExtractionHotspot("ortholog_data/FKS2_filtered.aln")

################################################################################
# Step 2. Add taxonomy to the dataframes.
################################################################################

# Create file for NCBI Common Taxonomy Tree tool. One time only.
#write(c(df_FKS1$Species, df_FKS2$Species), file = "allspecies.txt")

# Read NCBI Common Taxonomy Tree output.
MyTree <- readLines("ortholog_data/commontree-2026.txt")

# Format the file so it can be properly used.
MyTreeFormate <- gsub("[|+-]", "", MyTree)
MyTreeFormate <- gsub("\\", "", MyTreeFormate, fixed = TRUE)
MyTreeFormate <- gsub("^ +", "", MyTreeFormate)
MyTreeFormate <- gsub("\\[", "", MyTreeFormate)
MyTreeFormate <- gsub("\\]", "", MyTreeFormate)

# Create list and add the species to it.
TaxoTree <- list(
  Fungi = list(
    Blastocladiomycota = MyTreeFormate[3],
    Basidiomycota = MyTreeFormate[4:33],
    Ascomycota = MyTreeFormate[34:471],
    Mucoromycota = MyTreeFormate[472:477]
  ),
  Viridiplantae = list(Streptophyta = MyTreeFormate[481], Chlorophyta = MyTreeFormate[483])
)

# Add Specific cases.
TaxoTree$Fungi$Mucoromycota <- append(TaxoTree$Fungi$Mucoromycota, "Rhizopus oryzae")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Phanerochaete chrysosporium")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Melampsora larici-populina")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Melampsora laricipopulina")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Chaetomium thermophilum")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida glabrata")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida castellii")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida bracarensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida nivariensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Arthrobotrys oligospora")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Rutstroemia sp. NJR-2017a BBW")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Rutstroemia sp. NJR-2017a WRK4")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Marssonina brunnea")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Byssochlamys spectabilis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Penicillium sp. 'occitanis'")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Talaromyces cellulolyticus")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Leptosphaeria maculans")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Pyrenophora tritici-repentis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Acidomyces richmondensis")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Torrubiella hemipterigena")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota,
                                    "Saccharomycetaceae sp. 'Ashbya aceri'")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Saccharomyces sp. 'boulardii'")
TaxoTree$Viridiplantae$Chlorophyta <- append(TaxoTree$Viridiplantae$Chlorophyta,
                                             "Chlamydomonas reinhardtii")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Ustilago maydis")
TaxoTree$Fungi$Basidiomycota <- append(TaxoTree$Fungi$Basidiomycota, "Cryptococcus gattii VGI")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Rhynchosporium commune")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Fonsecaea erecta")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Grosmannia clavigera")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Acremonium chrysogenum")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Podospora anserina")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida auris")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida intermedia")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Candida inconspicua")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Tetrapisispora blattae")
TaxoTree$Fungi$Ascomycota <- append(TaxoTree$Fungi$Ascomycota, "Kazachstania naganishii")

AddTaxonomy <- function(Vector_of_species,
                        dataframe,
                        TaxonomyTree = TaxoTree)
{
  # Create vectors
  Kingdom <- vector()
  Phylum <- vector()
  
  # Find the kingdom and phylum of every specie in the data.
  for (i in Vector_of_species)
  {
    if (i %in% TaxonomyTree$Fungi$Blastocladiomycota)
    {
      Kingdom <- append(Kingdom, "Fungi")
      Phylum <- append(Phylum, "Blastocladiomycota")
      next
    }
    if (i %in% TaxonomyTree$Fungi$Basidiomycota)
    {
      Kingdom <- append(Kingdom, "Fungi")
      Phylum <- append(Phylum, "Basidiomycota")
      next
    }
    if (i %in% TaxonomyTree$Fungi$Ascomycota)
    {
      Kingdom <- append(Kingdom, "Fungi")
      Phylum <- append(Phylum, "Ascomycota")
      next
    }
    if (i %in% TaxonomyTree$Fungi$Mucoromycota)
    {
      Kingdom <- append(Kingdom, "Fungi")
      Phylum <- append(Phylum, "Mucoromycota")
      next
    }
    if (i %in% TaxonomyTree$Viridiplantae$Streptophyta)
    {
      Kingdom <- append(Kingdom, "Viridiplantae")
      Phylum <- append(Phylum, "Streptophyta")
      next
    }
    if (i %in% TaxonomyTree$Viridiplantae$Chlorophyta)
    {
      Kingdom <- append(Kingdom, "Viridiplantae")
      Phylum <- append(Phylum, "Chlorophyta")
      next
    }
    else
    {
      Kingdom <- append(Kingdom, "-")
      Phylum <- append(Phylum, "-")
      next
    }
  }
  
  # Add to Kingdom and Phylum to FKS data frame.
  dataframe <- cbind(dataframe, Kingdom)
  dataframe <- cbind(dataframe, Phylum)
  
  dataframe
}


# Is human pathogen ? Based on WHO fungal priority pathogens list to guide
# research, development and public health action.
IsPathogen <- function(Species)
{
  return (
    Species %in% c(
      'Cryptococcus neoformans',
      'Candida auris',
      'Aspergillus fumigatus',
      'Candida albicans',
      'Nakaseomyces glabrata',
      'Candida glabrata',
      'Histoplasma spp.',
      'Mucorales',
      'Fusarium spp.',
      'Candida tropicalis',
      'Candida parapsilosis',
      'Scedosporium spp.',
      'Cryptococcus gattii',
      'Lomentospora prolificans',
      'Talaromyces marneffei',
      'Coccidiodes spp.',
      'Pneumocystis jirovecii',
      'Pichia kudriavzeveii',
      'Candida krusei',
      'Paracoccidioides spp.'
    )
  )
}


# Execute functions for taxonomy.
df_FKS1 <- AddTaxonomy(df_FKS1$Species, df_FKS1)
df_FKS2 <- AddTaxonomy(df_FKS2$Species, df_FKS2)

CriticalHumanPathogen1 <- unlist(lapply(df_FKS1$Species, IsPathogen))
CriticalHumanPathogen2 <- unlist(lapply(df_FKS2$Species, IsPathogen))

df_FKS1 <- cbind(df_FKS1, Is_Human_Pathogen = CriticalHumanPathogen1)
df_FKS2 <- cbind(df_FKS2, Is_Human_Pathogen = CriticalHumanPathogen2)

################################################################################
# Step 3. Combine dataframes, reorder and rename columns.
################################################################################

library(dplyr)

# make new columns.
df_FKS1$Hotspot3_gaps <- df_FKS1$Hotspot3  # hotspot3 with gaps
df_FKS1$Hotspot3 <- gsub("-", "", df_FKS1$Hotspot3)  # hotspot3 without gaps
df_FKS1$Paralog <- "Fks1"
df_FKS2$Paralog <- "Fks2"

# combine rows.
s3_data <- bind_rows(df_FKS1, df_FKS2)

# Reorder columns.
s3_data_reordered <- s3_data[, c(
  "Identifier",
  "Species",
  "Kingdom",
  "Phylum",
  "Paralog",
  "Hotspot1",
  "Hotspot2",
  "Hotspot3_gaps",
  "Hotspot3",
  "GapsHotspot1",
  "GapsHotspot2",
  "GapsHotspot3",
  "Is_Human_Pathogen"
)]

# Rename columns.
colnames(s3_data_reordered)[1] <- "Metaphors_ID"
colnames(s3_data_reordered)[5] <- "Hit_when_using_query"

################################################################################
# Step 4. Check if hotspot sequences are in the oPools.
################################################################################

# Read the files carrying the sequences representing the oligos.
FKS1_HS1_uniques <- (readLines("ortholog_data/FKS1_HS1_orthologs_unique.fa"))[-1]
FKS1_HS2_uniques <- (readLines("ortholog_data/FKS1_HS2_orthologs_unique.fa"))[-1]
FKS1_HS3_uniques <- (readLines("ortholog_data/FKS1_HS3_orthologs_unique.fa"))[-1]
FKS2_HS1_uniques <- (readLines("ortholog_data/FKS2_HS1_orthologs_unique.fa"))[-1]
FKS2_HS2_uniques <- (readLines("ortholog_data/FKS2_HS2_orthologs_unique.fa"))[-1]

# Add columns checking if hotspot sequences are in a file.
s3_data_reordered$HS1_in_oPool_FKS1_HS1 <- with(s3_data_reordered,
                                                ifelse(Hotspot1 %in% FKS1_HS1_uniques, TRUE, FALSE))
s3_data_reordered$HS1_in_oPool_FKS2_HS1 <- with(s3_data_reordered,
                                                ifelse(Hotspot1 %in% FKS2_HS1_uniques, TRUE, FALSE))

s3_data_reordered$HS2_in_oPool_FKS1_HS2 <- with(s3_data_reordered,
                                                ifelse(Hotspot2 %in% FKS1_HS2_uniques, TRUE, FALSE))
s3_data_reordered$HS2_in_oPool_FKS2_HS2 <- with(s3_data_reordered,
                                                ifelse(Hotspot2 %in% FKS2_HS2_uniques, TRUE, FALSE))

s3_data_reordered$HS3_in_oPool_FKS1_HS3 <- with(s3_data_reordered,
                                                ifelse(Hotspot3 %in% FKS1_HS3_uniques, TRUE, FALSE))

# change names of Candida glagrata and Candida auris.
s3_data_reordered$Species <- gsub('Candida glabrata', 'Nakaseomyces glabratus', s3_data_reordered$Species)
s3_data_reordered$Species <- gsub('Candida auris', 'Candidozyma auris', s3_data_reordered$Species)



write.csv(file = "ortholog_data/S3_data.csv",
          x = s3_data_reordered,
          row.names = FALSE,
          quote = FALSE)