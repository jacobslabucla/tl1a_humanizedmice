# #===========README===========#
# input:
#   1. csv containing all datapoints in long format
# 
# processes:
#   1. standardize data
#   2. create distance matrix for each sample, for each method, lipidomics + prim metabolite; amine + prim metabolite; amine + lipidomics; all 3 combined
#   3. create pcoa plot for each distance matrix
# 
# output:
#   1. figures:
#     1.1. from processes 2.

library(tidyverse)
library(ggplot2)
library(vegan)

setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metabolomics")

#===========005 UNSUPERVISED CLUSTERING===========#

#all_metab <- read_csv("04-analysis/compiled_metabolites_longformat.csv")
lipid <- all_metab %>% filter(assay == "Lipidomics") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean", 'gene_background')))
amines <- all_metab %>% filter(assay == "BiogenicAmines") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean", 'gene_background')))
gcms <- all_metab %>% filter(assay == "PrimaryMetabolite") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean", 'gene_background')))

lipid.compact <- all_metab %>% filter(assay == "Lipidomics") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
                  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
amines.compact <- all_metab %>% filter(assay == "BiogenicAmines") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
                  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
gcms.compact <- all_metab %>% filter(assay == "PrimaryMetabolite") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
                  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')

lipid.dist <- lipid.compact %>% decostand("standardize") %>% vegdist(method = "euclidean")
lipid_pcoa <- wcmdscale(d = lipid.dist, eig = TRUE)
lipid_pcoa_df <- data.frame(lipid_pcoa$points) %>% rownames_to_column('sample_name') %>% left_join(lipid, by = 'sample_name')

amines.dist <- amines.compact %>% decostand("standardize") %>% vegdist(method = "euclidean")
amines_pcoa <- wcmdscale(d = amines.dist, eig = TRUE)
amines_pcoa_df <- data.frame(amines_pcoa$points) %>% rownames_to_column('sample_name') %>% left_join(lipid, by = 'sample_name')

gcms.dist <- gcms.compact %>% decostand("standardize") %>% vegdist(method = "euclidean")
gcms_pcoa <- wcmdscale(d = gcms.dist, eig = TRUE)
gcms_pcoa_df <- data.frame(gcms_pcoa$points) %>% rownames_to_column('sample_name') %>% left_join(lipid, by = 'sample_name')

ggplot(data = lipid_pcoa_df) +
  geom_point(aes(x = Dim1, y = Dim2, color = gene_background)) + theme_bw()

ggsave("04-analysis/01_summarystatistics/pcoa__lipidomics.png", width = 6, height = 6, limitsize = FALSE)

ggplot(data = amines_pcoa_df) +
  geom_point(aes(x = Dim1, y = Dim2, color = gene_background)) + theme_bw()

ggsave("04-analysis/01_summarystatistics/pcoa__biogenicamines.png", width = 6, height = 6, limitsize = FALSE)

ggplot(data = gcms_pcoa_df) +
  geom_point(aes(x = Dim1, y = Dim2, color = gene_background)) + theme_bw()

ggsave("04-analysis/01_summarystatistics/pcoa__primarymetabolite.png", width = 6, height = 6, limitsize = FALSE)

#test for normality of each feature
feature.normtest <- all_metab %>% group_by(assay, UID) %>% summarise(feat_norm_shapiro_pval = shapiro.test(peak_height_log)$p.value)

feature.test <- all_metab %>% filter(assay == "BiogenicAmines" & UID == "0.20_832.22_ESI pos") %>% select(one_of(c('assay', 'UID', 'sample_name', 'peak_height_impmean', 'peak_height_log')))

ggplot(feature.test, aes(x = peak_height_log)) + geom_histogram()
ggplot(feature.test, aes(x = peak_height_impmean)) + geom_histogram()
