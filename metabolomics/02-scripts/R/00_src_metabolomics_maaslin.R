
library(Hmisc)
library(tidyverse)
library(ggplot2)
library(vegan)
library(htmlwidgets)
library(readxl)

setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metabolomics")

##reading the files

lipid_pp <- read.csv("03-data/intermediate_data/00_lipidall_preprocessed.csv")
amines_pp <- read.csv("03-data/intermediate_data/00_aminesall_preprocessed.csv")
gcms_pp <- read.csv("03-data/intermediate_data/00_gcmsall_preprocessed.csv")

lipid_pp_nan <- read.csv("03-data/intermediate_data/00_lipidall_preprocessed_nan.csv")
amines_pp_nan <- read.csv("03-data/intermediate_data/00_aminesall_preprocessed_nan.csv")
gcms_pp_nan <- read.csv("03-data/intermediate_data/00_gcmsall_preprocessed_nan.csv")

fixDuplicatedNames <- function(df){
  d <- df
  
  duplicated_names <- d %>% select(one_of(c("name", "feat_id"))) %>% distinct() %>%
    group_by(name) %>% filter(n()>1)
  null_names <- d %>% select(one_of(c("name", "feat_id"))) %>% distinct() %>%
    filter(is.na(name))
  
  d.fixed <- d %>% 
    mutate(nameFix = ifelse(feat_id %in% null_names$feat_id, feat_id, paste(name, strsplit(feat_id, "_") %>% sapply(`[[`, 3), sep = "_"))) %>%
    mutate(nameFix = ifelse(feat_id %in% duplicated_names$feat_id, nameFix, name))
  
  return(d.fixed)
}

lipid_pp <- fixDuplicatedNames(lipid_pp)
amines_pp <- fixDuplicatedNames(amines_pp)
gcms_pp <- fixDuplicatedNames(gcms_pp)

# lipid_pp_nan <- fixDuplicatedNames(lipid_pp_nan)
# amines_pp_nan <- fixDuplicatedNames(amines_pp_nan)
# gcms_pp_nan <- fixDuplicatedNames(gcms_pp_nan)

lipid_pp_tl1a <- lipid_pp %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
lipid_pp_il10 <- lipid_pp %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)

amines_pp_tl1a <- amines_pp %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
amines_pp_il10 <- amines_pp %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)

gcms_pp_tl1a <- gcms_pp %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
gcms_pp_il10 <- gcms_pp %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)

# lipid_pp_tl1a_nan <- lipid_pp_nan %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# lipid_pp_il10_nan <- lipid_pp_nan %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# 
# amines_pp_tl1a_nan <- amines_pp_nan %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# amines_pp_il10_nan <- amines_pp_nan %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# 
# gcms_pp_tl1a_nan <- gcms_pp_nan %>% filter(gene_background == "TL1A_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# gcms_pp_il10_nan <- gcms_pp_nan %>% filter(gene_background == "IL10_mice") %>% mutate(interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)
# 

metadata <- read_excel("03-data/raw_data/metadata.xlsx", sheet = 'clean', na = "NA") %>%
  mutate(d_colonize = factor(d_colonize), d_collect = factor(d_collect), 
         donor_type = paste0(OTU_type, "_", IBD_stat),
         interactIBD_otu2 = as.data.frame(model.matrix(~IBD_stat*OTU_type, data = .))$`IBD_statIBD:OTU_typeotu_2`)


#===========001 UNSUPERVISED CLUSTERING===========#

pca_result <- function(imputedf, metadata) {
  #Create PCA plot of raw data, colored according to potential confounding factors
  #PCA should be with or without blanks. Samples that are similar to blanks are removed and deemed defective.
  #Aim is to provide initial diagnostic plot
  
  d <- imputedf
  
  d.compact <- d %>% 
    select(one_of(c("sample_id", "feat_id", "log_pkh_imputed"))) %>%
    pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>% 
    column_to_rownames(var = 'sample_id')
  
  factors <- c("mouse_sex", "IBD_stat", "OTU_type")
  
  m <- metadata %>% select(one_of(append(factors, "sample_id")))
  
  gene_bg <<- unique(d$gene_background)
  assay <<- unique(d$assay)
  
  d.dist <- d.compact %>% vegdist(method = "euclidean", na.rm = TRUE)
  
  pca <- wcmdscale(d.dist, k = 2) %>% as.data.frame() %>%
    rownames_to_column(var = "sample_id") %>% 
    left_join(m, by = "sample_id") %>%
    left_join(distinct(select(d, one_of(c("sample_id", "gene_background")))), by = "sample_id")
  
  for (x in factors){
    p <- ggplot(pca, aes(x = V1, y = V2, text = sample_id, color = eval(parse(text = x)))) +
      geom_point() +
      theme_bw() +
      labs(x = "PC1", y = "PC2") +
      guides(color = guide_legend(title = x))+
      ggtitle(paste0("PCA \n", assay, "_", gene_bg, "_", x))
    
    gp <- ggplotly(p)
    
    saveWidgetFix(gp, file = paste0("04-analysis/01_pcaresult/", assay, "_", gene_bg, "_", x, "_pca.html"))
  }
}

pca_result(lipid_pp_tl1a, metadata)
pca_result(lipid_pp_il10, metadata)
pca_result(amines_pp_tl1a, metadata)
pca_result(amines_pp_il10, metadata)
pca_result(gcms_pp_tl1a, metadata)
pca_result(gcms_pp_il10, metadata)

#looking at correlation between different features
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

corr_features <- function(imputedf){
  d <- imputedf
  
  d.compact <- d %>% 
    select(one_of(c("sample_id", "feat_id", "log_pkh_imputed"))) %>%
    pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>% 
    column_to_rownames(var = 'sample_id')
  
  d.corr <- rcorr(as.matrix(d.compact))
  
  d.flat <- flattenCorrMatrix(d.corr$r, d.corr$P)
  
  ggplot(d.flat, aes(x =  d$column, y, fill= Z)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum()
}


#===============002 LINEAR MODELLING=================#

##==================MAASLIN2======================##
library(Maaslin2)

maaslin.dplyr <- function(df, factors){
  d <- df
  
  assay <- unique(d$assay)
  gene_bg <- unique(d$gene_background)
  names <- d %>% select(one_of(c("nameFix", "feat_id", "assay", 'gene_background'))) %>% distinct()
  
  input <- d %>% select(one_of(c("sample_id", "log_pkh_imputed", "feat_id"))) %>%
    pivot_wider(names_from = feat_id, values_from = log_pkh_imputed) %>%
    column_to_rownames("sample_id")
  
  m <- d %>% select(one_of(append(factors, "sample_id"))) %>%
    distinct() %>%
    column_to_rownames("sample_id")
  
  if(gene_bg == "TL1A_mice"){
    factors = factors
  }
  else{
    factors = c("OTU_type", "IBD_stat", "mouse_sex", "interactIBD_otu2")
  }
  
  fit_data <- Maaslin2(input_data=input, 
                       input_metadata=m, 
                       output = paste0("04-analysis/01_maaslin2/",assay, "_", gene_bg), 
                       fixed_effects = factors,
                       reference = c("mouse_sex,F", "OTU_type,otu_1"),
                       normalization="NONE", 
                       transform ="NONE")
  
  fit_data_df <- fit_data[["results"]] %>% as.data.frame() %>%
    left_join(names, by = c("feature" = "feat_id"))
  
  return(fit_data_df)
}

factors <- c("OTU_type", "IBD_stat", "mouse_sex", "infl_score", "interactIBD_otu2")

factors_woibd <- c("OTU_type", "mouse_sex", "infl_score")

factors_wointeract <- c("OTU_type", "IBD_stat", "mouse_sex", "infl_score")

factors_wointeractinfl <- c("OTU_type", "IBD_stat", "mouse_sex")

factors_onlyinfl <- c("mouse_sex", "infl_score")

lipid_tl1a_maaslin_wointeractinfl <- maaslin.dplyr(lipid_pp_tl1a, factors_wointeractinfl)
amines_tl1a_maaslin_wointeractinfl <- maaslin.dplyr(amines_pp_tl1a, factors_wointeractinfl)
gcms_tl1a_maaslin_wointeractinfl <- maaslin.dplyr(gcms_pp_tl1a, factors_wointeractinfl)

lipid_il10_maaslin_wointeractinfl <- maaslin.dplyr(lipid_pp_il10, factors_wointeractinfl)
amines_il10_maaslin_wointeractinfl <- maaslin.dplyr(amines_pp_il10, factors_wointeractinfl)
gcms_il10_maaslin_wointeractinfl <- maaslin.dplyr(gcms_pp_il10, factors_wointeractinfl)

lipid_tl1a_maaslin_onlyinfl <- maaslin.dplyr(lipid_pp_tl1a, factors_onlyinfl)
amines_tl1a_maaslin_onlyinfl <- maaslin.dplyr(amines_pp_tl1a, factors_onlyinfl)
gcms_tl1a_maaslin_onlyinfl <- maaslin.dplyr(gcms_pp_tl1a, factors_onlyinfl)

lipid_il10_maaslin_onlyinfl <- maaslin.dplyr(lipid_pp_il10, factors_onlyinfl)
amines_il10_maaslin_onlyinfl<- maaslin.dplyr(amines_pp_il10, factors_onlyinfl)
gcms_il10_maaslin_onlyinfl <- maaslin.dplyr(gcms_pp_il10, factors_onlyinfl)

# write.csv(lipid_tl1a_maaslin_wointeractinfl, "03-data/intermediate_data/01_lipid_tl1a_maaslin_wointeractinfl.csv")
# write.csv(amines_tl1a_maaslin_wointeractinfl, "03-data/intermediate_data/01_amines_tl1a_maaslin_wointeractinfl.csv")
# write.csv(gcms_tl1a_maaslin_wointeractinfl, "03-data/intermediate_data/01_gcms_tl1a_maaslin_wointeractinfl.csv")
# 
# write.csv(lipid_il10_maaslin_wointeractinfl, "03-data/intermediate_data/01_lipid_il10_maaslin_wointeractinfl.csv")
# write.csv(amines_il10_maaslin_wointeractinfl, "03-data/intermediate_data/01_amines_il10_maaslin_wointeractinfl.csv")
# write.csv(gcms_il10_maaslin_wointeractinfl, "03-data/intermediate_data/01_gcms_il10_maaslin_wointeractinfl.csv")

lipid_tl1a_maaslin_noibd <- maaslin.dplyr(lipid_pp_tl1a, factors_woibd)
amines_tl1a_maaslin_noibd <- maaslin.dplyr(amines_pp_tl1a, factors_woibd)
gcms_tl1a_maaslin_noibd <- maaslin.dplyr(gcms_pp_tl1a, factors_woibd)

lipid_il10_maaslin_noibd <- maaslin.dplyr(lipid_pp_il10, factors_woibd)
amines_il10_maaslin_noibd <- maaslin.dplyr(amines_pp_il10, factors_woibd)
gcms_il10_maaslin_noibd <- maaslin.dplyr(gcms_pp_il10, factors_woibd)

# write.csv(lipid_tl1a_maaslin_noibd, "03-data/intermediate_data/01_lipid_tl1a_maaslin_noibd.csv")
# write.csv(amines_tl1a_maaslin_noibd, "03-data/intermediate_data/01_amines_tl1a_maaslin_noibd.csv")
# write.csv(gcms_tl1a_maaslin_noibd, "03-data/intermediate_data/01_gcms_tl1a_maaslin_noibd.csv")
# 
# write.csv(lipid_il10_maaslin_noibd, "03-data/intermediate_data/01_lipid_il10_maaslin_noibd.csv")
# write.csv(amines_il10_maaslin_noibd, "03-data/intermediate_data/01_amines_il10_maaslin_noibd.csv")
# write.csv(gcms_il10_maaslin_noibd, "03-data/intermediate_data/01_gcms_il10_maaslin_noibd.csv")

lipid_tl1a_maaslin_woint <- maaslin.dplyr(lipid_pp_tl1a, factors_wointeract)
amines_tl1a_maaslin_woint <- maaslin.dplyr(amines_pp_tl1a, factors_wointeract)
gcms_tl1a_maaslin_woint <- maaslin.dplyr(gcms_pp_tl1a, factors_wointeract)

lipid_il10_maaslin_woint <- maaslin.dplyr(lipid_pp_il10, factors_wointeract)
amines_il10_maaslin_woint <- maaslin.dplyr(amines_pp_il10, factors_wointeract)
gcms_il10_maaslin_woint <- maaslin.dplyr(gcms_pp_il10, factors_wointeract)

# write.csv(lipid_tl1a_maaslin_woint, "03-data/intermediate_data/01_lipid_tl1a_maaslin_woint.csv")
# write.csv(amines_tl1a_maaslin_woint, "03-data/intermediate_data/01_amines_tl1a_maaslin_woint.csv")
# write.csv(gcms_tl1a_maaslin_woint, "03-data/intermediate_data/01_gcms_tl1a_maaslin_woint.csv")
# 
# write.csv(lipid_il10_maaslin_woint, "03-data/intermediate_data/01_lipid_il10_maaslin_woint.csv")
# write.csv(amines_il10_maaslin_woint, "03-data/intermediate_data/01_amines_il10_maaslin_woint.csv")
# write.csv(gcms_il10_maaslin_woint, "03-data/intermediate_data/01_gcms_il10_maaslin_woint.csv")

lipid_tl1a_maaslin_nan <- maaslin.dplyr(lipid_pp_tl1a_nan, factors)
amines_tl1a_maaslin_nan <- maaslin.dplyr(amines_pp_tl1a_nan, factors)
gcms_tl1a_maaslin_nan <- maaslin.dplyr(gcms_pp_tl1a_nan, factors)

lipid_il10_maaslin_nan <- maaslin.dplyr(lipid_pp_il10_nan, factors)
amines_il10_maaslin_nan <- maaslin.dplyr(amines_pp_il10_nan, factors)
gcms_il10_maaslin_nan <- maaslin.dplyr(gcms_pp_il10_nan, factors)

write.csv(lipid_tl1a_maaslin_nan, "03-data/intermediate_data/01_lipid_tl1a_maaslin_nan.csv")
write.csv(amines_tl1a_maaslin_nan, "03-data/intermediate_data/01_amines_tl1a_maaslin_nan.csv")
write.csv(gcms_tl1a_maaslin_nan, "03-data/intermediate_data/01_gcms_tl1a_maaslin_nan.csv")

write.csv(lipid_il10_maaslin_nan, "03-data/intermediate_data/01_lipid_il10_maaslin_nan.csv")
write.csv(amines_il10_maaslin_nan, "03-data/intermediate_data/01_amines_il10_maaslin_nan.csv")
write.csv(gcms_il10_maaslin_nan, "03-data/intermediate_data/01_gcms_il10_maaslin_nan.csv")