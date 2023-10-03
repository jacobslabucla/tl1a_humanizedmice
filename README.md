# tl1a_humanizedmice

PROJECT TITLE = INVESTIGATING THE ROLE OF TL1A-RECEPTOR IN MEDIATING IBD SYMPTOMS

PROJECT GOALS = TO IDENTIFY CANDIDATE METABOLITE/METABOLIC PATHWAYS THAT MEDIATES IBD

ANALYSIS PIPELINE =
  A. METABOLOMICS
    Input = Peak height tables from UCSF group = 1) primary metabolites - GCMS, 2) lipidomics, 3) biogenic amines; Data are normalized with internal standard TIC; Some metabolites have already been identified based on spectra matching
    1. Data preparation:
      1.1. Annotation
      1.2. Imputation
      1.3. Log transformation
        QC =  Check normal distribution per feature
      1.4. Further normalization
        QC = Check normal distribution per feature
      1.5. Summary table/pivot table: N of samples for each donor group; n of features per sample, split known and unknowns
        
    2. Data analysis:
      2.1. Global/High-level analysis: beta diversity (euclidean) -> PCoA -> PERMANOVA
      2.2. Identification of differential metabolies based on donor groups (OTU-type and IBD status) -> generalized linear model
      2.3. Identification of metabolites associated with IBD severity (correlate metabolite amount with inflammation scoring)
      2.4. Identification of metabolites associated with immune cell differentiation (correlate metabolite amount with cell sorting data)

  B. SHOTGUN METAGENOMICS
