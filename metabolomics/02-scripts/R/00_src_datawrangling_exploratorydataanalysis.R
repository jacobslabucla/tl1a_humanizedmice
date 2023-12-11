#===========README===========#
# input:
#   1. raw excel sheets from service providers
#   2. metadata sheet linking sample name with experimental treatments
# 
# processes:
#   1. data wrangling: 
#     1.1. read the data and metadata sheets + check if all contains same amount of samples
#     1.2. append columns to determine the methodology (i.e. lipidomics, prim_metab, or biogenicamines)
#     1.3. convert tables to 'long' format
#     1.4. convert all values 0 or below to NA
#   2. first data filtering:
#     2.1. for each gene background, stool donor, sample loc, and sex (i.e. biological replicates), filter features according to presence of feature in % of samples. do not cut out more than 10% of features found.
#       2.1.F. histogram of presence of feature in % of samples (00_31F_(lipid,gcms,amines)_hist__featurepercentpresence)
#     2.2. for each gene background, stool donor, sample loc, and sex (i.e. biological replicates), filter features according to coef of variance (CV) of feature peak height (that is not NA) that is more than 50%.
#       2.2.F. boxplot, x = stool donor, y = feature variance (00_32F_(lipid,gcms,amines)_boxp__featurevariance_stooldonor)
#   3. imputation:
#     3.1. determine the missing value distribution after first data filtering
#     3.2. for each gene background, stool donor, and sex: impute NA values with kNN method
#   4. second data filtering:
#     4.1. for each gene background, stool donor, sample loc, and sex (i.e. biological replicates), filter features according to percentile peak height value (try 5th and 10th percentile). do not cut out more than 5% of features found.
#       4.1.F. boxplot, x = percentile peak height value cutoff, y = alpha_diversity (00_33F_(lipid,gcms,amines)_boxp__peakheightcutoff_adiversity). to determine the effects of absolute value cutoff towards a-diversity
#   5. exploratory data analysis (EDA) for each method and sample:
#     2.1. determine the mean and median; SD and IQR; min and max; skewness and kurtosis of peak height
#     2.2. determine the counts of metabolite found and percent identified/assigned
#     2.3. determine their alpha diversity
#     2.4. compile data into a table and export csv
#     2.5. visualize 2.1 with boxplot (with 3 rows pertaining to the 3 methods, x = sample_name, y = peak height)
#     2.6. visualize 2.2 with stacked histogram (with x = sample name, y = counts of metabolite (identified and not))
#     2.7. visualize 2.3 with barplot
# 
# output:
#   1. tables:
#     1.1 compiled long table of all metabolites
#     1.2 summary table of the various EDA variables
#   2. figures:
#     2.1 from 1.6
#     2.2 from 2.5-2.7

library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(plotly)
library(forcats)
library(stats)
library(vegan)
library(lme4)
library(nlme)

setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metabolomics")

#===========001 READING FILE AND INITIAL DATA WRANGLING===========#

#read the raw excel files from service provider
lipid <- read_excel("03-data/raw_data/mx 738840_Jacob_Lipidomics_humanstool_mousefeces_Submit_07-2023.xlsx", sheet = 'clean_1', na = c("na", ""))
amines <- read_excel("03-data/raw_data/mx 739116_Jacob_HILIC_biogenicamines_humanstool_09-2023 submit.xlsx", sheet = 'clean_1', na = c("na", ""))
gcms <- read_excel("03-data/raw_data/mx 738540_Jacob_GCTOF_humanstool_mousestool_08-2023 submit.xlsx", sheet = 'clean_1', na = c("na", ""))

#read the metadata file
metadata <- read_excel("03-data/raw_data/metadata.xlsx", sheet = 'clean', na = c("NA", "#N/A")) %>%
            mutate(d_colonize = factor(d_colonize), d_collect = factor(d_collect))

#create unique feat_id column, for lipid and amines as both can have same mz, ret time but different ESI mode
lipid <- lipid %>% mutate(feat_id = paste("lip", esi_mode, mz, r.time, sep = "_"))
amines <- amines %>% mutate(feat_id = paste("amin", esi_mode, mz, r.time, sep = "_"))
gcms <- gcms %>% mutate(feat_id = paste("gcms", mz, r.time, sep = "_"))

#check if all mass are unique
which(duplicated(lipid$feat_id) == TRUE)
which(duplicated(amines$feat_id) == TRUE)
which(duplicated(gcms$feat_id) == TRUE)

#generate 'long' dataframe for each
name_start <- c("BSD", "TL1A", "IL", "blank")
samp_names <- colnames(lipid[,grepl(paste(name_start, collapse="|"), names(lipid))])

name_start_gcms <- c("BSD", "TL1A", "IL")
samp_names_gcms <- colnames(gcms[,grepl(paste(name_start_gcms, collapse="|"), names(gcms))])

lipid_td <- lipid %>% 
  select(one_of(c('name', 'feat_id', samp_names))) %>%
  pivot_longer(cols = all_of(samp_names), 
               names_to='sample_id', 
               values_to='pkh')

amines_td <- amines %>% 
  select(one_of(c('name', 'feat_id', samp_names))) %>%
  pivot_longer(cols = all_of(samp_names), 
               names_to='sample_id', 
               values_to='pkh')

gcms_td <- gcms %>% 
  select(one_of(c('name', 'feat_id', samp_names_gcms))) %>%
  pivot_longer(cols = all_of(samp_names_gcms), 
               names_to='sample_id', 
               values_to='pkh')

#check if each row corresponds to each observation in data
nrow(lipid_td) == length(unique(lipid$feat_id)) * length(samp_names)
nrow(amines_td) == length(unique(amines$feat_id)) * length(samp_names)
nrow(gcms_td) == length(unique(gcms$feat_id)) * length(samp_names_gcms)

#check if all samples are present in all methods
setdiff(unique(lipid_td$sample_id), unique(amines_td$sample_id))
setdiff(unique(lipid_td$sample_id), unique(gcms_td$sample_id))

#add assay and gene_background column
lipid_td$assay <- "Lipidomics"
amines_td$assay <- "BiogenicAmines"
gcms_td$assay <- "PrimaryMetabolite"

lipid_td$gene_background <- with(lipid_td, 
                                 ifelse(grepl('BSD', sample_id),'Human',
                                        ifelse(grepl('TL1A', sample_id),'TL1A_mice', 
                                               ifelse(grepl('IL', sample_id),'IL10_mice', 
                                                      'Blank'))))
amines_td$gene_background <- with(amines_td,  
                                  ifelse(grepl('BSD', sample_id),'Human',
                                         ifelse(grepl('TL1A', sample_id),'TL1A_mice', 
                                                ifelse(grepl('IL', sample_id),'IL10_mice', 
                                                       'Blank'))))
gcms_td$gene_background <- with(gcms_td,  
                                ifelse(grepl('BSD', sample_id),'Human',
                                       ifelse(grepl('TL1A', sample_id),'TL1A_mice', 
                                              ifelse(grepl('IL', sample_id),'IL10_mice', 
                                                     'Blank'))))

#create a dataframe containing the amount of unique features for each step
lipid_unique_feat <- lipid_td %>% group_by(sample_id) %>% summarize(raw_unique_features = sum(!is.na(pkh)))
amines_unique_feat <- amines_td %>% group_by(sample_id) %>% summarize(raw_unique_features = sum(!is.na(pkh)))
gcms_unique_feat <- gcms_td %>% group_by(sample_id) %>% summarize(raw_unique_features = sum(!is.na(pkh)))

#filter out human samples, split data according to TL1A and IL10 mice
lipid_tl1a_td <- lipid_td %>% filter(gene_background == "TL1A_mice" | gene_background == "Blank")
amines_tl1a_td <- amines_td %>% filter(gene_background == "TL1A_mice" | gene_background == "Blank")
gcms_tl1a_td <- gcms_td %>% filter(gene_background == "TL1A_mice" | gene_background == "Blank")

lipid_il10_td <- lipid_td %>% filter(gene_background == "IL10_mice" | gene_background == "Blank")
amines_il10_td <- amines_td %>% filter(gene_background == "IL10_mice" | gene_background == "Blank")
gcms_il10_td <- gcms_td %>% filter(gene_background == "IL10_mice" | gene_background == "Blank")

#====== DIAGNOSTIC PLOTS FOR RAW DATA ======#
raw_summary <- function(rawtidydf, pltitle, gplot = TRUE){
  d <- rawtidydf
  
  d.sum <- d %>% group_by(gene_background, sample_id) %>%
            summarise(total_pkh = sum(pkh, na.rm = TRUE), found_features = sum(!is.na(pkh)))
  
  if(gplot == TRUE){
    
    pltitle = pltitle
    
    p <- ggplot(data = d.sum, aes(x = reorder(sample_id, total_pkh, max, na.rm = T, decreasing = TRUE), y = total_pkh)) +
          geom_segment(aes(xend = sample_id, yend = 0)) +
          geom_point(size = 2, color = "red") +
          theme_bw()  +
          facet_grid(rows = "gene_background") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4)) +
          labs(x = "Sample ID", y = "Sum of feature peak heights") +
          ggtitle(pltitle)
          
  }
  
  return(p)
}

p.sumpkh.lipid <- raw_summary(lipid_td, "Lipidomics Peak Sum per Sample (Raw Data)")
p.sumpkh.amines <- raw_summary(amines_td, "Amines Peak Sum per Sample (Raw Data)")
p.sumpkh.gcms <- raw_summary(gcms_td, "GCMS Peak Sum per Sample (Raw Data)")

ggsave("04-analysis/00_summarystatistics/raw_sumpeakheight_lipid.png", plot = p.sumpkh.lipid, width = 30, height = 10, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/raw_sumpeakheight_amines.png", plot = p.sumpkh.amines, width = 30, height = 10, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/raw_sumpeakheight_gcms.png", plot = p.sumpkh.gcms, width = 30, height = 10, limitsize = FALSE)

##======VIOLIN PLOT ========##
library(ggpubr)

fibrosis <- metadata %>% 
  mutate(donor_id = paste0(IBD_stat, "_", OTU_type)) %>%
  mutate(donor_id = factor(donor_id, levels=c("Healthy_otu_1", "Healthy_otu_2", "IBD_otu_1", "IBD_otu_2"))) %>%
  select(one_of(c("donor_id", "fibrosis_duo", "fibrosis_jej", "fibrosis_ile", "fibrosis_colon"))) %>%
  pivot_longer(-donor_id, names_to = "location", values_to = "fibrosis_score")

loc_name <- c('fibrosis_colon' = "Colon",
              'fibrosis_duo' = "Duodenum",
              'fibrosis_jej' = "Jejunum",
              'fibrosis_ile' = "Ileum")

compar <- list(c("Healthy_otu_1", "Healthy_otu_2"),
               c("Healthy_otu_1", "IBD_otu_1"),
               c("Healthy_otu_2", "IBD_otu_1"),
               c("IBD_otu_1", "IBD_otu_2"))

ggplot(fibrosis, aes(x = donor_id, y = fibrosis_score)) +
  geom_violin(aes(color = donor_id)) + 
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               stackdir = "center",
               dotsize = 2,
               alpha = .7) +
  scale_x_discrete(breaks = unique(fibrosis$donor_id), labels = as.character(c("Healthy - Non Dysbiosis",
                                                                               "IBD - Non Dysbiosis",
                                                                               "IBD - Dysbiosis",
                                                                               "Healthy - Dysbiosis"
                                                                               ))) +
  facet_grid(~location, labeller = as_labeller(loc_name)) +
  theme_bw() +
  scale_color_manual(values = c("#B0E0E6","#0000FF","#F0CCB0","#D2691E")) +
  labs(x = "", y = "Histology Score") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  stat_compare_means(method = "anova", label.y = 53) +
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = compar, label.y = c(42, 44, 48, 51)) +
  

#=====PCA PLOT RAW DATA=====#  
pca_diag <- function(rawdf, metadata, confound_fac, blank_incl = TRUE) {
  #Create PCA plot of raw data, colored according to potential confounding factors
  #PCA should be with or without blanks. Samples that are similar to blanks are removed and deemed defective.
  #Aim is to provide initial diagnostic plot
  
  d <- rawdf
  
  if(blank_incl){
      d.compact <- d %>% 
                mutate(log_pkh = log10(pkh+1)) %>%
                mutate(log_pkh = ifelse(is.na(log_pkh),0, log_pkh)) %>%
                select(one_of(c("sample_id", "feat_id", "log_pkh"))) %>%
                pivot_wider(names_from = "feat_id", values_from = "log_pkh") %>% 
                column_to_rownames(var = 'sample_id')
  }
  else
  {
    d.compact <- d %>% 
      filter(gene_background != "Blank") %>%
      mutate(log_pkh = log10(pkh+1)) %>%
      mutate(log_pkh = ifelse(is.na(log_pkh),0, log_pkh)) %>%
      select(one_of(c("sample_id", "feat_id", "log_pkh"))) %>%
      pivot_wider(names_from = "feat_id", values_from = "log_pkh") %>% 
      column_to_rownames(var = 'sample_id')
  }
  
  m <- metadata %>% select(one_of(append(confound_fac, "sample_id")))
  
  gene_bg <<- unique(d$gene_background)[1]
  assay <<- unique(d$assay)
  
  d.dist <- d.compact %>% vegdist(method = "euclidean", na.rm = TRUE)
  
  pca <- wcmdscale(d.dist, k = 2) %>% as.data.frame() %>%
         rownames_to_column(var = "sample_id") %>% 
         left_join(m, by = "sample_id") %>%
         left_join(distinct(select(d, one_of(c("sample_id", "gene_background")))), by = "sample_id")
  
  confound_fac <- append(confound_fac, "gene_background")
  
  for (x in confound_fac){
    p <- ggplot(pca, aes(x = V1, y = V2, text = sample_id, color = eval(parse(text = x)))) +
        geom_point() +
        theme_bw() +
        labs(x = "PC1", y = "PC2") +
        guides(color = guide_legend(title = x))
    
    gp <- ggplotly(p)
    
    htmlwidgets::saveWidget(gp, file = paste0("04-analysis/00_pcadiag/", assay, "_", gene_bg, "_", x, "_raw_data_pca.html"))
  }
}

confound_fac = c("d_colonize", "d_collect", "d_interval", "cage_id", "isolator_id", "batch")

pca_diag(lipid_tl1a_td, metadata, confound_fac)
pca_diag(amines_tl1a_td, metadata, confound_fac)
pca_diag(gcms_tl1a_td, metadata, confound_fac)
pca_diag(lipid_il10_td, metadata, confound_fac)
pca_diag(amines_il10_td, metadata, confound_fac)
pca_diag(gcms_il10_td, metadata, confound_fac)

#===========002 BLANK CORRECTION===========#

blank_corr <- function(df_tidy, blank_missing_cutoff, blank_log_diff, gplot = FALSE){
  
  #Function: to do blank_correction given a tidy-formatted dataframe.
  #Workflow: take dataframe > split into blank and non-blank sample > log10 transform peak height
  # > calculate blank average per feature > filter out (turn to 0) features that are not present in 75% of samples
  # > calculate log difference between blank and sample log peak height > filter out (turn to NA) features that are
  # below difference of 1
  
  d <- df_tidy
  d.names <- names(d)
  
  #separate out the samples and the blank
  d.sample <- d %>% filter(gene_background != "Blank") %>% mutate(log_pkh = log10(pkh))
  d.blank <- d %>% filter(gene_background == "Blank") %>% mutate(log_pkh = log10(pkh))
  
  d.blank.avg <- d.blank %>% group_by(feat_id) %>% 
                  summarize(log_blank_avg = mean(log_pkh, na.rm = T), percent_presence = 1 - sum(is.na(log_pkh))/n()) %>%
                  mutate(filt_blank_avg = ifelse(percent_presence >= blank_missing_cutoff, log_blank_avg, 0))  #reason for 75% is because 21 out of 28 measurements. We want our blanks to be more certain than 50%
  
  d.corrected <- d.sample %>% left_join(d.blank.avg, by = "feat_id") %>%
                  mutate(blank_diff = log_pkh - filt_blank_avg) %>%
                  mutate(log_pkh_corrected = ifelse(blank_diff >= blank_log_diff, log_pkh, NA)) %>%
                  select(one_of(c(d.names, "log_pkh_corrected", "filt_blank_avg")))
  
  return(d.corrected)
  
  if (gplot == TRUE){
    
    
  }
}

lipid_tl1a_bc <- blank_corr(lipid_tl1a_td, 0.8, 1) %>% 
                left_join(metadata, by = "sample_id") %>% 
                mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_"))
amines_tl1a_bc <- blank_corr(amines_tl1a_td, 0.8, 1) %>% 
                  left_join(metadata, by = "sample_id") %>% 
                  mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_"))
gcms_tl1a_bc <- gcms_tl1a_td %>% 
                mutate(log_pkh_corrected = log10(pkh + 1)) %>% 
                left_join(metadata, by = "sample_id") %>% 
                mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_")) #gcms no blanks were run or it was already subtracted by vendor

lipid_il10_bc <- blank_corr(lipid_il10_td, 0.8, 1) %>% 
                left_join(metadata, by = "sample_id") %>% 
                mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_"))
amines_il10_bc <- blank_corr(amines_il10_td, 0.8, 1) %>% 
                  left_join(metadata, by = "sample_id") %>% 
                  mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_"))
gcms_il10_bc <- gcms_il10_td %>% 
                mutate(log_pkh_corrected = log10(pkh + 1)) %>% 
                left_join(metadata, by = "sample_id") %>% 
                mutate(donor_type = paste(OTU_type, IBD_stat, sep = "_")) #gcms no blanks were run or it was already subtracted by vendor

lipid_unique_feat_bl <- bind_rows(lipid_tl1a_bc, lipid_il10_bc) %>% 
                      group_by(sample_id)%>%
                      summarize(ufeat_post_blank = sum(!is.na(log_pkh_corrected))) %>%
                      right_join(lipid_unique_feat, by="sample_id")

amines_unique_feat_bl <- bind_rows(amines_tl1a_bc, amines_il10_bc) %>% 
                          group_by(sample_id)%>%
                          summarize(ufeat_post_blank = sum(!is.na(log_pkh_corrected))) %>%
                          right_join(amines_unique_feat, by="sample_id")

gcms_unique_feat_bl <- bind_rows(gcms_tl1a_bc, gcms_il10_bc) %>% 
                        group_by(sample_id)%>%
                        summarize(ufeat_post_blank = sum(!is.na(log_pkh_corrected))) %>%
                        right_join(gcms_unique_feat, by="sample_id")

rm(lipid, lipid_td, lipid_il10_td, lipid_tl1a_td, lipid_unique_feat)
rm(amines, amines_td, amines_il10_td, amines_tl1a_td, amines_unique_feat)
rm(gcms, gcms_td, gcms_il10_td, gcms_tl1a_td, gcms_unique_feat)

#===========002 DATA FILTERING AND QUALITY CONTROL===========#
## filtering based on percent 'missingness' of the data, by their ICC, and whether or not the feature is annotated
fisher.test.dplyr <- function(df){
  
  cont_table <- df %>% 
                  select(one_of(c('donor_type', 'n_missing', 'n_non_missing'))) %>%
                  column_to_rownames("donor_type")
  
  fisher_pval <- fisher.test(cont_table)$p.value
  
  return(fisher_pval)
}

icc.dpylr <- function(df){
  
  data <- df %>%
          select(one_of("log_pkh_corrected", "donor_type"))
  
  mm <- lme(log_pkh_corrected~ 1, random = ~1|donor_type, data=data, na.action = na.omit)
  
  varcorr <- as.numeric(VarCorr(mm)[1:2])
  
  
  return(varcorr[1]/sum(varcorr))
}

feature_filtering<- function(df_blanked, max_missingness, missing_cutoff, icc_cutoff, gplot = FALSE, filtname = TRUE){
  
  d <- df_blanked
  
  d.fisher <- d %>%
              select(one_of(c("sample_id","feat_id", "donor_type", "log_pkh_corrected"))) %>%
              group_by(donor_type, feat_id) %>%
              summarize(n_missing = sum(is.na(log_pkh_corrected)), n_non_missing = sum(!is.na(log_pkh_corrected))) %>%
              ungroup()%>%
              group_by(feat_id) %>%
              do(fisher_pval = fisher.test.dplyr(.)) %>% 
              unnest(fisher_pval) %>% ungroup() 
  
  #the percent missingness below is across all samples, and our filtering will remove the feature across all of our samples
  
  d.sum.miss <- d %>%
            group_by(feat_id) %>%
            summarize(percent_missing = sum(is.na(log_pkh_corrected))/ (sum(is.na(log_pkh_corrected)) + sum(!is.na(log_pkh_corrected)))) %>%
            left_join(d.fisher, by = "feat_id")
  
  d.filt.miss.clean <-  d %>%
                    left_join(d.sum.miss, by = "feat_id") %>%
                    filter(percent_missing <= missing_cutoff)
  
  d.filt.miss.unclean <- d %>%
                   left_join(d.sum.miss, by = "feat_id") %>% 
                   filter(percent_missing < max_missingness & percent_missing > missing_cutoff & fisher_pval <= 0.05) %>% #0.9 is to just remove features that are heavily missing
                   mutate(log_pkh_corrected = ifelse(is.na(log_pkh_corrected), 0, log_pkh_corrected)) #this is 0 imputation for values missing above the cutoff (as KNN imputation will be spurious for these features)
  
  d_zero_imputed <- d.filt.miss.unclean %>% 
                    group_by(sample_id) %>% 
                    summarize(ufeat_zero_imputed = sum(log_pkh_corrected == 0))
  
  d.filt.miss <- bind_rows(d.filt.miss.clean, d.filt.miss.unclean)
                  
  # d.filt.miss <- d %>%
  #           left_join(d.sum.miss, by = "feat_id") %>%
  #           filter(percent_missing <= 0.75) %>% #at least 1 donor_type out of 4 has values
  #           filter(percent_missing <= missing_cutoff | fisher_pval <= 0.05) #recommended less than 20% missing as cutoff
  
  d_unique_feat_miss <- d.filt.miss %>%
                   group_by(sample_id)%>%
                   summarize(ufeat_post_missing_filt = sum(!is.na(log_pkh_corrected))) %>%
                   right_join(d_zero_imputed, by="sample_id")
  
  #ICC is used to determine the variability between different treatment (or traditionally test giver in clinical setting) 
  # is larger or smaller than the variability within each treatment. In our case of metabolomics is for each feature,
  # we ask the question whether the variability across 4 donor types is larger or smaller than variability within donor type
  # if variability across 4 donor types is larger than within donor type, then we will keep the feature
  
  d.icc <- d.filt.miss %>%
            group_by(feat_id) %>%
            do(icc = icc.dpylr(.)) %>%
            unnest(icc) %>% ungroup()
  
  d.filt.icc <- d.filt.miss %>%
                  left_join(d.icc, by = "feat_id") %>%
                  filter(icc >= icc_cutoff)
  
  d_unique_feat_icc <- d.filt.icc %>%
                        group_by(sample_id)%>%
                        summarize(ufeat_post_icc_filt = sum(!is.na(log_pkh_corrected)))%>%
                        right_join(d_unique_feat_miss, by="sample_id")
  
  if (filtname == TRUE){
    d.filt.name <- d.filt.icc %>%
      filter(!is.na(name))
  }
  
  else{
    d.filt.name <- d.filt.icc 
  }
  
  d_unique_feat_name <- d.filt.name %>%
                        group_by(sample_id)%>%
                        summarize(ufeat_post_annotated_filt = sum(!is.na(log_pkh_corrected)))%>%
                        right_join(d_unique_feat_icc, by="sample_id")
  
  if (gplot == TRUE){
    
    # p <- ggplot(data = filter(d.filt.miss, feat_id == "lip_esi_neg_392.24619999999999_0.864"), aes(x = donor_type, y = log_pkh_corrected)) +
    #   geom_violin() +
    #   geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5)
    # p
    
  }
  
  returnlist <- list("filt_df" = d.filt.name, "filt_summary" = d_unique_feat_name)
  
  return(returnlist)
}

filt_lipid_tl1a_stringent <- feature_filtering(lipid_tl1a_bc, 0.9, 0.2, 0.1)
filt_amines_tl1a_stringent <- feature_filtering(amines_tl1a_bc, 0.9, 0.2, 0.1)
filt_gcms_tl1a_stringent <- feature_filtering(gcms_tl1a_bc, 0.9, 0.2, 0.1)

filt_lipid_il10_stringent <- feature_filtering(lipid_il10_bc, 0.9, 0.2, 0.1)
filt_amines_il10_stringent <- feature_filtering(amines_il10_bc, 0.9, 0.2, 0.1)
filt_gcms_il10_stringent <- feature_filtering(gcms_il10_bc, 0.9, 0.2, 0.1)

filt_lipid_tl1a_nonstringent <- feature_filtering(lipid_tl1a_bc, 0.5, 0.5, 0)
filt_amines_tl1a_nonstringent <- feature_filtering(amines_tl1a_bc, 0.5, 0.5, 0)
filt_gcms_tl1a_nonstringent <- feature_filtering(gcms_tl1a_bc, 0.5, 0.5, 0)

filt_lipid_il10_nonstringent <- feature_filtering(lipid_il10_bc, 0.5, 0.5, 0)
filt_amines_il10_nonstringent <- feature_filtering(amines_il10_bc, 0.5, 0.5, 0)
filt_gcms_il10_nonstringent <- feature_filtering(gcms_il10_bc, 0.5, 0.5, 0)

# filt_lipid_tl1a <- feature_filtering(lipid_tl1a_bc, 0.2, 0.1, filtname = FALSE)
# filt_amines_tl1a <- feature_filtering(amines_tl1a_bc, 0.2, 0.1, filtname = FALSE)
# filt_gcms_tl1a <- feature_filtering(gcms_tl1a_bc, 0.2, 0.1, filtname = FALSE)
# 
# filt_lipid_il10 <- feature_filtering(lipid_il10_bc, 0.2, 0.1, filtname = FALSE)
# filt_amines_il10 <- feature_filtering(amines_il10_bc, 0.2, 0.1, filtname = FALSE)
# filt_gcms_il10 <- feature_filtering(gcms_il10_bc, 0.2, 0.1, filtname = FALSE)

lipid_tl1a_fl_s <- filt_lipid_tl1a_stringent$filt_df
amines_tl1a_fl_s <- filt_amines_tl1a_stringent$filt_df
gcms_tl1a_fl_s <- filt_gcms_tl1a_stringent$filt_df

lipid_il10_fl_s <- filt_lipid_il10_stringent$filt_df
amines_il10_fl_s <- filt_amines_il10_stringent$filt_df
gcms_il10_fl_s <- filt_gcms_il10_stringent$filt_df

lipid_tl1a_fl_ns <- filt_lipid_tl1a_nonstringent$filt_df
amines_tl1a_fl_ns <- filt_amines_tl1a_nonstringent$filt_df
gcms_tl1a_fl_ns <- filt_gcms_tl1a_nonstringent$filt_df

lipid_il10_fl_ns <- filt_lipid_il10_nonstringent$filt_df
amines_il10_fl_ns <- filt_amines_il10_nonstringent$filt_df
gcms_il10_fl_ns <- filt_gcms_il10_nonstringent$filt_df

lipid_n_of_feat <- bind_rows(filt_lipid_tl1a$filt_summary, filt_lipid_il10$filt_summary) %>% 
                    right_join(lipid_unique_feat_bl, by="sample_id") %>%
                    mutate(ufeat_post_blank_0imput = ufeat_zero_imputed + ufeat_post_blank) %>%
                    relocate(ufeat_post_blank_0imput, .after = ufeat_post_missing_filt) %>%
                    rev() %>%
                    relocate(sample_id)
amines_n_of_feat <- bind_rows(filt_amines_tl1a$filt_summary, filt_amines_il10$filt_summary) %>% 
                    right_join(amines_unique_feat_bl, by="sample_id") %>%
                    mutate(ufeat_post_blank_0imput = ufeat_zero_imputed + ufeat_post_blank) %>%
                    relocate(ufeat_post_blank_0imput, .after = ufeat_post_missing_filt) %>%
                    rev() %>%
                    relocate(sample_id)
gcms_n_of_feat <- bind_rows(filt_gcms_tl1a$filt_summary, filt_gcms_il10$filt_summary) %>% 
                    right_join(gcms_unique_feat_bl, by="sample_id") %>%
                    mutate(ufeat_post_blank_0imput = ufeat_zero_imputed + ufeat_post_blank) %>%
                    relocate(ufeat_post_blank_0imput, .after = ufeat_post_missing_filt) %>%
                    rev() %>%
                    relocate(sample_id)

rm(lipid_il10_bc, lipid_tl1a_bc, filt_lipid_il10, filt_lipid_tl1a, p.sumpkh.lipid, lipid_unique_feat_bl)
rm(amines_il10_bc, amines_tl1a_bc, filt_amines_il10, filt_amines_tl1a, p.sumpkh.amines, amines_unique_feat_bl)
rm(gcms_il10_bc, gcms_tl1a_bc, filt_gcms_il10, filt_gcms_tl1a, p.sumpkh.gcms, gcms_unique_feat_bl)

#===========003 IMPUTATION#===========#
##KNN Imputation for features present in >75% of samples in a given donor_type

impute.knn.obs.sel <- function(dat, K=10) { #rownames are samples, columns are features
  
  results <- list()
  cor.cutoff <- 0.2     # use only variables with cor>0.2 for distance computation
  
  da1 <- dat 
  da1list <- da2list <- rep(list(dat),length(K)) 
  
  incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
  incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
  
  Cor <- cor(da1,use="p")
  
  D2list <- lapply(incom.vars, function(j) {
    varsel <- which(abs(Cor[j,])>cor.cutoff)  
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
    D2 <- as.matrix(dist(scale(da1[,varsel])),upper=T,diag=T) 
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist(scale(da1)),upper=T,diag=T)*sqrt(length(varsel)/ncol(da1)) 
      D2[is.na(D2)] <- D2a[is.na(D2)] 
    }
    diag(D2) <- NA
    D2})
  names(D2list) <- incom.vars
  
  for (i in incom.obs){
    comvars <-  complete.cases(as.numeric(da1[i,]))
    for (j in which(!comvars)) {
      D2 <- D2list[[as.character(j)]]                                 
      if(any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last=NA)
        KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da1[ii,j])))] 
      } else {
        KNNids  <- NULL
      }
      da1list <- lapply(1:length(da1list),function(ii) {
        k <- K[ii] 
        da <-  da1list[[ii]]
        if(!is.null(KNNids)) {
          KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
        }
        if(length(KNNids_sel)<1) {
          KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]
        } else if (length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] 
        if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
          da_sel <- da[KNNids_sel,j]
          da[i,j] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
        da}) 
    }
  }
  da1list <- lapply(da1list, function(da) {
    da <- apply(da,2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
      x}) 
    da})
  
  results <- c(results,list(da1list))
  names(results)[length(results)] <- "knn.sample.euc.sel" 
  rm(da1list,da2list)  
  return(results$knn.sample.euc.sel[[1]])
}

imput.knn.dplyr <- function(inputdf, col_names, feat_name_col, feat_values_col, K= 10){
  datr <- inputdf %>% 
    select(one_of(col_names)) %>%
    pivot_wider(names_from = feat_name_col, values_from = feat_values_col) %>% ## the pivot wider needs to be inside the function; if not, the pivot will force filtered out features as NA into the dataframe
    column_to_rownames("sample_id") %>% 
    select(matches(paste(c("lip", 'amin', 'gcms'), collapse = "|")))
  
  knn.res <- impute.knn.obs.sel(datr, K = 10) %>% as.data.frame()
  
  output <- knn.res %>% rownames_to_column(var="sample_id") %>%
            pivot_longer(-sample_id, names_to = 'feat_id', values_to = "log_pkh_imputed")
            
  
  return(output)
}

lipid_tl1a_knn_s <- lipid_tl1a_fl_s %>% group_by(donor_type) %>% 
                  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                  unnest(imput) %>% ungroup() %>%
                  right_join(lipid_tl1a_fl_s, by = c("feat_id", "sample_id", "donor_type"))

amines_tl1a_knn_s <- amines_tl1a_fl_s %>% group_by(donor_type) %>% 
                  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                  unnest(imput) %>% ungroup() %>%
                  right_join(amines_tl1a_fl_s, by = c("feat_id", "sample_id", "donor_type"))

gcms_tl1a_knn_s <- gcms_tl1a_fl_s %>% group_by(donor_type) %>% 
                  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                  unnest(imput) %>% ungroup() %>%
                  right_join(gcms_tl1a_fl_s, by = c("feat_id", "sample_id", "donor_type"))

lipid_il10_knn_s <- lipid_il10_fl_s %>% group_by(donor_type) %>% 
                  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                  unnest(imput) %>% ungroup() %>%
                  right_join(lipid_il10_fl_s, by = c("feat_id", "sample_id", "donor_type"))

amines_il10_knn_s <- amines_il10_fl_s %>% group_by(donor_type) %>% 
                    do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                    unnest(imput) %>% ungroup() %>%
                    right_join(amines_il10_fl_s, by = c("feat_id", "sample_id", "donor_type"))

gcms_il10_knn_s <- gcms_il10_fl_s %>% group_by(donor_type) %>% 
                do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
                unnest(imput) %>% ungroup() %>%
                right_join(gcms_il10_fl_s, by = c("feat_id", "sample_id", "donor_type"))

lipid_tl1a_knn_ns <- lipid_tl1a_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(lipid_tl1a_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

amines_tl1a_knn_ns <- amines_tl1a_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(amines_tl1a_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

gcms_tl1a_knn_ns <- gcms_tl1a_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(gcms_tl1a_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

lipid_il10_knn_ns <- lipid_il10_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(lipid_il10_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

amines_il10_knn_ns <- amines_il10_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(amines_il10_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

gcms_il10_knn_ns <- gcms_il10_fl_ns %>% group_by(donor_type) %>% 
  do(imput = imput.knn.dplyr(., c("feat_id", "log_pkh_corrected", "sample_id"), "feat_id", "log_pkh_corrected")) %>% 
  unnest(imput) %>% ungroup() %>%
  right_join(gcms_il10_fl_ns, by = c("feat_id", "sample_id", "donor_type"))

##=======Combining the two filtering results for visualization=======##
# kruskal_filt <- function(dffil1, dffil2){
#   d1 <- dffil1
#   d2 <- dffil2
#   
#   d1ktest <- d1 %>%
#             group_by(feat_id) %>%
#             do(kruskald1 = kruskal.test(log_pkh_imputed ~ donor_type, data = .)$p.value)
#   
#   d2ktest <- d2 %>%
#             group_by(feat_id) %>%
#             do(kruskald2 = kruskal.test(log_pkh_imputed ~ donor_type, data = .)$p.value)
#   
#   dktest <- d1ktest %>% left_join(d2ktest, by = "feat_id") %>% 
#     filter(!is.null(kruskald1) & !is.null(kruskald2)) %>%
#     mutate(kruskald1 = as.numeric(kruskald1),
#            kruskald2 = as.numeric(kruskald2)) %>%
#     left_join(distinct(select(d1, one_of("feat_id", "percent_missing"))), by = "feat_id")
#   
#   g <- ggplot(dktest, aes(x = kruskald1, y = kruskald2, text = feat_id, color = percent_missing)) +
#         geom_point() +
#         labs(x = "Kruskal-Wallis pval method 1",
#              y = "Kruskal-Wallis pval method 2") +
#         scale_x_continuous(trans = "log10") +
#         scale_y_continuous(trans = "log10") +
#         theme_bw()
#   
#   return(ggplotly(g))
# }
# 
# lipid_tl1a_filtest <- kruskal_filt(lipid_tl1a_knn_s, lipid_tl1a_knn_ns)
# amines_tl1a_filtest <- kruskal_filt(amines_tl1a_knn_s, amines_tl1a_knn_ns)
# 
# pre1feature_test <- amines_tl1a_knn_s %>% filter(feat_id == "amin_esi_pos_113.0694_1.251") %>% 
#                 mutate(filmethod = "method1") %>% 
#                 select(one_of(c("log_pkh_imputed", "filmethod", "donor_type", "feat_id")))
# 
# pre2feature_test <- amines_tl1a_knn_ns %>% filter(feat_id == "amin_esi_pos_113.0694_1.251") %>%
#                 mutate(filmethod = "method2") %>% 
#                 select(one_of(c("log_pkh_imputed", "filmethod", "donor_type", "feat_id")))
# 
# feature_test <- bind_rows(pre1feature_test, pre2feature_test)
# 
# ggplot(feature_test, aes(x = donor_type, y = log_pkh_imputed)) +
#   geom_violin() +
#   facet_grid(~filmethod)

##=======Preparing Files for Saving========##

lipid_knn <- bind_rows(lipid_tl1a_knn, lipid_il10_knn)
amines_knn <- bind_rows(amines_tl1a_knn, amines_il10_knn)
gcms_knn <- bind_rows(gcms_tl1a_knn, gcms_il10_knn)

csvsave <- function(combinedf){
  d <- combinedf
  
  assay <- unique(d$assay)
  
  for (x in unique(d$gene_background)){
    d.subset <- d %>%
                filter(gene_background == x) %>%
                select(one_of(c("log_pkh_imputed", "sample_id", "feat_id"))) %>%
                pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed")
    
    d.metadata <- d %>%
                  filter(gene_background == x) %>%
                  select(one_of(c("sample_id", "OTU_type","IBD_stat", "mouse_sex"))) %>%
                  distinct()
    
    write.csv(d.subset, 
              paste0("03-data/intermediate_data/", "00_",assay,"_",x,"_preprocesseddata.csv"),
              row.names = FALSE)
    write.csv(d.metadata,
              paste0("03-data/intermediate_data/", "00_",assay,"_",x,"_metadata.csv"),
              row.names = FALSE)
    
  }
  
}

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

lipid_pp <- fixDuplicatedNames(lipid_knn)
amines_pp <- fixDuplicatedNames(amines_knn)
gcms_pp <- fixDuplicatedNames(gcms_knn)

csvsave(lipid_pp)
csvsave(amines_pp)
csvsave(gcms_pp)

# write.csv(lipid_knn, "03-data/intermediate_data/00_lipidall_preprocessed.csv")
# write.csv(amines_knn, "03-data/intermediate_data/00_aminesall_preprocessed.csv")
# write.csv(gcms_knn, "03-data/intermediate_data/00_gcmsall_preprocessed.csv")

# write.csv(lipid_knn, "03-data/intermediate_data/00_lipidall_preprocessed_nan.csv")
# write.csv(amines_knn, "03-data/intermediate_data/00_aminesall_preprocessed_nan.csv")
# write.csv(gcms_knn, "03-data/intermediate_data/00_gcmsall_preprocessed_nan.csv")

rm(lipid_il10_fl, lipid_tl1a_fl, lipid_tl1a_knn, lipid_il10_knn)
rm(amines_il10_fl, amines_tl1a_fl, amines_tl1a_knn, amines_il10_knn)
rm(gcms_il10_fl, gcms_tl1a_fl, gcms_tl1a_knn, gcms_il10_knn)

#===========004 SUMMARY STATISTICS#===========#

write.csv(lipid_n_of_feat, "04-analysis/00_summarystatistics/00_lipid_summary_uniquefeat.csv")
write.csv(amines_n_of_feat, "04-analysis/00_summarystatistics/00_amines_summary_uniquefeat.csv")
write.csv(gcms_n_of_feat, "04-analysis/00_summarystatistics/00_gcms_summary_uniquefeat.csv")

p.feathist.lipid.tl1a <- lipid_knn %>%
                    filter(gene_background == "TL1A_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)
p.feathist.lipid.il10 <- lipid_knn %>%
                    filter(gene_background == "IL10_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)

p.feathist.amines.tl1a <- amines_knn %>%
                    filter(gene_background == "TL1A_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)
p.feathist.amines.il10 <- amines_knn %>%
                    filter(gene_background == "IL10_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)

p.feathist.gcms.tl1a <- gcms_knn %>%
                    filter(gene_background == "TL1A_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)
p.feathist.gcms.il10 <- gcms_knn %>%
                    filter(gene_background == "IL10_mice") %>%
                    ggplot(aes(x = log_pkh_imputed)) + geom_histogram() + facet_wrap(~feat_id)


ggsave("04-analysis/00_summarystatistics/imputed_feathist_lipid_tl1a.png", plot = p.feathist.lipid.tl1a, width = 30, height = 30, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/imputed_feathist_amines_tl1a.png", plot = p.feathist.amines.tl1a, width = 30, height = 30, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/imputed_feathist_gcms_tl1a.png", plot = p.feathist.gcms.tl1a, width = 30, height = 30, limitsize = FALSE)

ggsave("04-analysis/00_summarystatistics/imputed_feathist_lipid_il10.png", plot = p.feathist.lipid.il10, width = 30, height = 30, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/imputed_feathist_amines_il10.png", plot = p.feathist.amines.il10, width = 30, height = 30, limitsize = FALSE)
ggsave("04-analysis/00_summarystatistics/imputed_feathist_gcms_il10.png", plot = p.feathist.gcms.il10, width = 30, height = 30, limitsize = FALSE)





# #===========004 NORMALISATION===========# NOT IMPLEMENTED AS OF 11/01/2023
# vsn.normalization <- function(imputed_df){
#   d <- imputed_df
#   d.vsn <- d %>%
#     select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
#     pivot_wider(names_from = "sample_id", values_from = "log_pkh_imputed") %>%
#     column_to_rownames("feat_id") %>%
#     as.matrix() %>%
#     justvsn()%>%
#     t() %>% as.data.frame() %>%
#     rownames_to_column(var="sample_id")  %>%
#     pivot_longer(cols = -sample_id, names_to = 'feat_id', values_to = "log_pkh_norm")
# }
# 
# library(vsn)
# #using vsn normalization technique
# lipid.vsn.norm <- lipid_filtered_imputed %>% 
#   select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#   pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#   column_to_rownames("UID") %>%
#   as.matrix() %>%
#   justvsn() %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column(var="sample_id") %>%
#   pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#   filter(!is.na(peak_height_imputed_normalized))
# lipid_fin <- lipid_filtered_imputed %>% left_join(lipid.vsn.norm, by = c("UID", "sample_id"))
# 
# amines.vsn.norm <- amines_filtered_imputed %>% 
#   select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#   pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#   column_to_rownames("UID") %>%
#   as.matrix() %>% replace(is.na(.), 0) %>%
#   justvsn() %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column(var="sample_id") %>%
#   pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#   filter(!is.na(peak_height_imputed_normalized))
# amines_fin <- amines_filtered_imputed %>% left_join(amines.vsn.norm, by = c("UID", "sample_id"))
# 
# gcms.vsn.norm <- gcms_filtered_imputed %>% 
#   select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#   pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#   column_to_rownames("UID") %>%
#   as.matrix() %>% replace(is.na(.), 0) %>%
#   justvsn() %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column(var="sample_id") %>%
#   pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#   filter(!is.na(peak_height_imputed_normalized))
# gcms_fin <- gcms_filtered_imputed %>% left_join(gcms.vsn.norm, by = c("UID", "sample_id"))
# 
# 
# 
# #check if features need normalization
# shapiro_norm_check <- function(imputedf) {
#   d <- imputedf
#   
#   shapiro_all <- d %>%
#     group_by(feat_id) %>%
#     summarise(
#       shapiro_p_value_all = shapiro.test(log_pkh_imputed)$p.value,
#       total_samples = n(),
#       zero_samples = sum(log_pkh_imputed == 0)
#     )
#   
#   shapiro_nonzero <- d %>%
#     filter(log_pkh_imputed != 0) %>%
#     group_by(feat_id) %>%
#     summarise(
#       shapiro_p_value_nonzero = shapiro.test(log_pkh_imputed)$p.value
#     )
#   
#   output <- left_join(shapiro_all, shapiro_nonzero, by = "feat_id")
#   
#   return(output)
# }
# 
# lipid_tl1a_norm.check <- shapiro_norm_check(lipid_tl1a_knn)
# amines_tl1a_norm.check <- shapiro_norm_check(amines_tl1a_knn)
# gcms_tl1a_norm.check <- shapiro_norm_check(gcms_tl1a_knn)
# 
# lipid_il10_norm.check <- shapiro_norm_check(lipid_il10_knn)
# amines_il10_norm.check <- shapiro_norm_check(amines_il10_knn)
# gcms_il10_norm.check <- shapiro_norm_check(gcms_il10_knn)

#lookup the donor stool value, incl ibd status and otu type
#check first if all sample_id in data matches sample_id in metadata
# setdiff(unique(metadata$sample_id), unique(lipid_td$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted
# setdiff(unique(metadata$sample_id), unique(amines_td$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted
# setdiff(unique(metadata$sample_id), unique(gcms_td$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted
# 
# lipid_td <- lipid_td %>% left_join(metadata, by = "sample_id") %>% mutate(donor_type = paste(IBD_stat, OTU_type, sep = "_"), days_to_collection = as.numeric(d_collect) - as.numeric(d_colonize))
# amines_td <- amines_td %>% left_join(metadata, by = "sample_id") %>% mutate(donor_type = paste(IBD_stat, OTU_type, sep = "_"), days_to_collection = as.numeric(d_collect) - as.numeric(d_colonize))
# gcms_td <- gcms_td %>% left_join(metadata, by = "sample_id") %>% mutate(donor_type = paste(IBD_stat, OTU_type, sep = "_"), days_to_collection = as.numeric(d_collect) - as.numeric(d_colonize))
# 
# #change all 0s and negative values to NA; negative value is from subtracting with blank ctrl signal;
# #useful to do before doing log transform to prevent error
# lipid_td <- lipid_td %>% 
#             mutate(pk_h = replace(pk_h, which(pk_h <= 0), NA)) %>% 
#             mutate(pk_h_log10 = log10(pk_h))
# amines_td <- amines_td %>% 
#             mutate(pk_h = replace(pk_h, which(pk_h <= 0), NA)) %>% 
#             mutate(pk_h_log10 = log10(pk_h))
# gcms_td <- gcms_td %>% 
#             mutate(pk_h = replace(pk_h, which(pk_h <= 0), NA)) %>% 
#             mutate(pk_h_log10 = log10(pk_h))
# 
# #check the data distribution
# p <- lipid_td %>%
#       ggplot(aes(x = reorder(sample_id, pk_h_log10, mean, na.rm = T, decreasing = TRUE), y = pk_h_log10)) + 
#       geom_boxplot() + facet_wrap(~gene_background, dir = "v", scales = "free_x") +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       labs(x = "Sample ID", y = "Log10 Peak Height")
# 
# p <- amines_td %>%
#       ggplot(aes(x = reorder(sample_id, pk_h_log10, mean, na.rm = T, decreasing = TRUE), y = pk_h_log10)) + 
#       geom_boxplot() + facet_wrap(~gene_background, dir = "v", scales = "free_x") +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       labs(x = "Sample ID", y = "Log10 Peak Height")
# 
# p <- gcms_td %>%
#       ggplot(aes(x = reorder(sample_id, pk_h_log10, mean, na.rm = T, decreasing = TRUE), y = pk_h_log10)) + 
#       geom_boxplot() + facet_wrap(~gene_background, dir = "v", scales = "free_x") +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       labs(x = "Sample ID", y = "Log10 Peak Height")
# 
# #Flagging samples that might be an outlier
# lipid_zscores <- lipid_td %>% 
#                   group_by(sample_id) %>% 
#                   summarize(mean_pk_h_log10 = mean(pk_h_log10, na.rm = T)) %>% 
#                   mutate(z_score = scale(mean_pk_h_log10))
# 
# amines_zscores <- amines_td %>% 
#                   group_by(sample_id) %>% 
#                   summarize(mean_pk_h_log10 = mean(pk_h_log10, na.rm = T)) %>% 
#                   mutate(z_score = scale(mean_pk_h_log10))
# 
# gcms_zscores <- gcms_td %>% 
#                   group_by(sample_id) %>% 
#                   summarize(mean_pk_h_log10 = mean(pk_h_log10, na.rm = T),
#                             shapiro_test = shapiro.test(pk_h_log10)) %>% 
#                   mutate(z_score = scale(mean_pk_h_log10))
# 
# library(ggpubr)
# ggqqplot(lipid_zscores$mean_pk_h_log10)
# 
# shapiro.test(lipid_zscores$mean_pk_h_log10)
# shapiro.test(amines_zscores$mean_pk_h_log10)
# shapiro.test(gcms_zscores$mean_pk_h_log10)

#===========002 DATA FILTERING AND QUALITY CONTROL===========#

#Filtering of samples that has means that are lower than other samples
#look at the mean of each sample and calculate Z-score

# 
# 
# test <- lipid_td %>% filter(sample_id == "TL1A_1229") %>% summarize(median =  median(pk_h, na.rm =T))
# 
# library(moments)
# 
# check_filtered_out <- function(prefiltereddf, postfiltereddf) {
#   feat_prefil = nrow(prefiltereddf) / length(unique(prefiltereddf$sample_id))
#   feat_postfil = nrow(postfiltereddf) / length(unique(prefiltereddf$sample_id))
#   percent_filt = (1 - feat_postfil/feat_prefil)*100
#   
#   cat("All features from data: ", feat_prefil, "\n")
#   cat("Features after removing unidentified ones: ", feat_postfil, "\n")
#   cat("Percent features removed: ", percent_filt, "% \n")
#   
#   return(feat_prefil)
# }
# #turn feature pk_h values to NA if below the x percentile of the entire dataset
# 
# 
# #remove features that are only 
# 
# #filter out features that are only present below 75% of the biological replicate (same fecal donor)
# 
# #filter out features whose value is below the x percentile of a given sample
# lipid_filt <- lipid_td
# 
# #remove features that are not annotated within the database
# lipid_filt <- lipid_td %>% filter(!is.na(name))
# amines_filt <- amines_td %>% filter(!is.na(name))
# gcms_filt <- gcms_td %>% filter(!is.na(name))
# 
# lipid_1st_filter <- check_filtered_out(lipid_td, lipid_filt)
# amines_1st_filter <- check_filtered_out(amines_td, amines_filt)
# gcms_1st_filter <- check_filtered_out(gcms_td, gcms_filt)
# 
# 
# 
# #checking the amount of biological replicates for each gene background, sampling location, donor stool, and sex
# lipid_nreplicate <- lipid_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates
# amines_nreplicate <- amines_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates
# gcms_nreplicate <- gcms_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates
# 
# #filter features according to: 1) presence <50% of the biological replicate samples, 2) CV less than 150, 3) kurtosis less than 4, and 4)peak height mean is less than 10th percentile of 
# lipid_bioreplicate <- lipid_lon %>% group_by(gene_background, location, donor, UID) %>% 
#                           summarize(tot_na = sum(is.na(peak_height_na)), 
#                                     identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
#                                     n_repl = n_distinct(sample_id),
#                                     sdev_bioreplicate = sd(peak_height_na, na.rm=T),
#                                     avg_bioreplicate = mean(peak_height_na, na.rm=T),
#                                     kurt_bioreplicate = kurtosis(peak_height_na, na.rm=T),
#                                     skew_bioreplicate =  skewness(peak_height_na, na.rm=T)) %>%
#                           mutate(percent_presence_across_bioreplicates = 1 - (tot_na / n_repl), CV_bioreplicate = sdev_bioreplicate/avg_bioreplicate*100)
# 
# lipid_percentile10 <- quantile(lipid_bioreplicate$avg_bioreplicate, .10, na.rm = T)
# 
# lipid_filt_stats <- lipid_bioreplicate %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_bioreplicates > 0),
#                                                                              feat_less50 = sum(percent_presence_across_bioreplicates < 0.5)) %>%
#                                                                     mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 50%
# 
# # lipid_lon_filtered <- lipid_lon %>% left_join(lipid_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
# #                         filter(percent_presence_across_bioreplicates>=0.5) %>%
# #                         filter(CV_bioreplicate <= 150) %>%
# #                         filter(kurt_bioreplicate <= 4) %>%
# #                         filter(avg_bioreplicate >= lipid_percentile10 )
# 
# lipid_lon_filtered <- lipid_lon %>% left_join(lipid_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
#   mutate(peak_height_na = ifelse(percent_presence_across_bioreplicates>=0.5&CV_bioreplicate <= 150&kurt_bioreplicate <= 4&avg_bioreplicate >= lipid_percentile10, peak_height_na, 0))
# 
# amines_bioreplicate <- amines_lon %>% group_by(gene_background, location, donor, UID) %>% 
#   summarize(tot_na = sum(is.na(peak_height_na)), 
#             identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
#             n_repl = n_distinct(sample_id),
#             sdev_bioreplicate = sd(peak_height_na, na.rm=T),
#             avg_bioreplicate = mean(peak_height_na, na.rm=T),
#             kurt_bioreplicate = kurtosis(peak_height_na, na.rm=T),
#             skew_bioreplicate =  skewness(peak_height_na, na.rm=T)) %>%
#   mutate(percent_presence_across_bioreplicates = 1 - (tot_na / n_repl), CV_bioreplicate = sdev_bioreplicate/avg_bioreplicate*100)
# 
# amines_percentile10 <- quantile(amines_bioreplicate$avg_bioreplicate, .10, na.rm = T)
# 
# amines_filt_stats <- amines_bioreplicate %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_bioreplicates > 0),
#                      feat_less50 = sum(percent_presence_across_bioreplicates < 0.5)) %>%
#                       mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 50%
# 
# # amines_lon_filtered <- amines_lon %>% left_join(amines_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
# #   filter(percent_presence_across_bioreplicates>=0.5) %>%
# #   filter(CV_bioreplicate <= 150) %>%
# #   filter(kurt_bioreplicate <= 4) %>%
# #   filter(avg_bioreplicate >= lipid_percentile10 )
# 
# amines_lon_filtered <- amines_lon %>% left_join(amines_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
#   mutate(peak_height_na = ifelse(percent_presence_across_bioreplicates>=0.5&CV_bioreplicate <= 150&kurt_bioreplicate <= 4&avg_bioreplicate >= amines_percentile10, peak_height_na, 0))
# 
# 
# gcms_bioreplicate <- gcms_lon %>% group_by(gene_background, location, donor, UID) %>% 
#   summarize(tot_na = sum(is.na(peak_height_na)), 
#             identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
#             n_repl = n_distinct(sample_id),
#             sdev_bioreplicate = sd(peak_height_na, na.rm=T),
#             avg_bioreplicate = mean(peak_height_na, na.rm=T),
#             kurt_bioreplicate = kurtosis(peak_height_na, na.rm=T),
#             skew_bioreplicate =  skewness(peak_height_na, na.rm=T)) %>%
#   mutate(percent_presence_across_bioreplicates = 1 - (tot_na / n_repl), CV_bioreplicate = sdev_bioreplicate/avg_bioreplicate*100)
# 
# gcms_percentile10 <- quantile(gcms_bioreplicate$avg_bioreplicate, .10, na.rm = T)
# 
# gcms_filt_stats <- gcms_bioreplicate %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_bioreplicates > 0),
#                    feat_less50 = sum(percent_presence_across_bioreplicates < 0.5)) %>%
#                    mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 50%
# 
# # gcms_lon_filtered <- gcms_lon %>% left_join(gcms_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
# #   filter(percent_presence_across_bioreplicates>=0.5) %>%
# #   filter(CV_bioreplicate <= 150) %>%
# #   filter(kurt_bioreplicate <= 4) %>%
# #   filter(avg_bioreplicate >= lipid_percentile10 )
# 
# gcms_lon_filtered <- gcms_lon %>% left_join(gcms_bioreplicate, by = c('UID', 'gene_background', 'location', 'donor')) %>% 
#   mutate(peak_height_na = ifelse(percent_presence_across_bioreplicates>=0.5&CV_bioreplicate <= 150&kurt_bioreplicate <= 4&avg_bioreplicate >= gcms_percentile10, peak_height_na, 0))
# 
# #===========003 IMPUTATION#===========#
# ##MEAN Imputation (seems to perform better than zero or half lowest and easy to implement according to Kokla, Virtanen, et al., 2019) (NOT USED)
# # mean_peak_per_sample <- lipid_lon_filtered %>% group_by(gene_background, location, donor, UID) %>% filter(!is.na(peak_height_na)) %>%
# #   summarise(mean_peak = mean(peak_height_na))
# # lipid_lon_filtered <- lipid_lon_filtered %>% left_join(mean_peak_per_sample, by = c("gene_background", "location", "donor", "UID"))  %>%
# #   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))
# # 
# # mean_peak_per_sample <- amines_lon_filtered %>% group_by(sample_id) %>% filter(!is.na(peak_height_na)) %>%
# #   summarise(mean_peak = mean(peak_height_na))
# # amines_lon_filtered <- amines_lon_filtered %>% left_join(mean_peak_per_sample, by = "sample_id")  %>%
# #   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))
# # 
# # mean_peak_per_sample <- gcms_lon_filtered %>% group_by(sample_id) %>% filter(!is.na(peak_height_na)) %>%
# #   summarise(mean_peak = mean(peak_height_na))
# # gcms_lon_filtered <- gcms_lon_filtered %>% left_join(mean_peak_per_sample, by = "sample_id")  %>%
# #   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))
# 
# ##KNN Imputation
# #impute with these columns = gene background, otu type, ibd status, peak height na, mz, r.time, identifier as row name
# library(VIM)
# 
# lipid_knn <- lipid_lon_filtered %>% group_by(gene_background, donor) %>% do(imput = imput.knn.dplyr(.)) %>% unnest(imput) %>% ungroup()
# lipid_filtered_imputed <- lipid_lon_filtered %>% left_join(lipid_knn, by = c("UID", "donor", "sample_id", "gene_background")) %>% mutate(peak_height_na_imputed = (replace_na(peak_height_na_imputed, 0)))
# 
# gcms_knn <- gcms_lon_filtered %>% group_by(gene_background, donor) %>% do(imput = imput.knn.dplyr(.)) %>% unnest(imput) %>% ungroup()
# gcms_filtered_imputed <- gcms_lon_filtered %>% left_join(gcms_knn, by = c("UID", "donor", "sample_id", "gene_background"))%>% mutate(peak_height_na_imputed = (replace_na(peak_height_na_imputed, 0)))
# 
# amines_knn <- amines_lon_filtered %>% group_by(gene_background, donor) %>% do(imput = imput.knn.dplyr(.)) %>% unnest(imput) %>% ungroup()
# amines_filtered_imputed <- amines_lon_filtered %>% left_join(amines_knn, by = c("UID", "donor", "sample_id", "gene_background"))%>% mutate(peak_height_na_imputed = (replace_na(peak_height_na_imputed, 0)))
# 
# # test <- setdiff(test_knn_2[!duplicated(test_knn_2[c("UID", "gene_background", "donor", "sample_id")]),]%>% select(one_of(c("gene_background","donor", "UID", "sample_id"))), 
# #                             lipid_lon_filtered[!duplicated(lipid_lon_filtered[c("UID", "gene_background", "donor", "sample_id")]),] %>% select(one_of(c("gene_background","donor", "UID", "sample_id")))) ##trying to figure out why is my pivot wider causing more observations
# 
# #===========004 NORMALISATION#===========#
# library(vsn)
# #using vsn normalization technique
# lipid.vsn.norm <- lipid_filtered_imputed %>% 
#                     select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#                     pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#                     column_to_rownames("UID") %>%
#                     as.matrix() %>%
#                     justvsn() %>%
#                     t() %>% as.data.frame() %>%
#                     rownames_to_column(var="sample_id") %>%
#                     pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#                     filter(!is.na(peak_height_imputed_normalized))
# lipid_fin <- lipid_filtered_imputed %>% left_join(lipid.vsn.norm, by = c("UID", "sample_id"))
# 
# amines.vsn.norm <- amines_filtered_imputed %>% 
#   select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#   pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#   column_to_rownames("UID") %>%
#   as.matrix() %>% replace(is.na(.), 0) %>%
#   justvsn() %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column(var="sample_id") %>%
#   pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#   filter(!is.na(peak_height_imputed_normalized))
# amines_fin <- amines_filtered_imputed %>% left_join(amines.vsn.norm, by = c("UID", "sample_id"))
# 
# gcms.vsn.norm <- gcms_filtered_imputed %>% 
#   select(one_of(c("UID", "peak_height_na_imputed", "sample_id")))%>%
#   pivot_wider(names_from = "sample_id", values_from = "peak_height_na_imputed") %>% 
#   column_to_rownames("UID") %>%
#   as.matrix() %>% replace(is.na(.), 0) %>%
#   justvsn() %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column(var="sample_id") %>%
#   pivot_longer(cols = starts_with("esi"), names_to = 'UID', values_to = "peak_height_imputed_normalized") %>%
#   filter(!is.na(peak_height_imputed_normalized))
# gcms_fin <- gcms_filtered_imputed %>% left_join(gcms.vsn.norm, by = c("UID", "sample_id"))
#                     
# #===========005 COMBINING PREPROCESSED VALUES#===========#
# 
# # #create log transformed data
# # lipid_filtered_imputed <- lipid_filtered_imputed %>% mutate(peak_height_log = log(peak_height_na_imputed, 10))
# # amines_filtered_imputed <- amines_filtered_imputed %>% mutate(peak_height_log = log(peak_height_na_imputed, 10))
# # gcms_filtered_imputed <- gcms_filtered_imputed %>% mutate(peak_height_log = log(peak_height_na_imputed, 10))
# 
# #combine the different assay output into one big table
# base_cols <- intersect(names(lipid_fin), names(amines_fin))
# 
# gcms_fin$mz <- as.character(gcms_fin$mz)
# 
# all_metab <- lipid_fin %>% select(one_of(base_cols)) %>% 
#               bind_rows(select(amines_fin, one_of(base_cols))) %>% 
#               bind_rows(select(gcms_fin, one_of(base_cols)))
# 
# #write.csv(all_metab, "04-analysis/compiled_metabolites_imputed_normalized_longformat.csv") ##commented so don't accidentally save the file whenever script is run
# 
# #===========003 EXPLORATORY DATA ANALYSIS===========#
# #all_metab <- read_csv("04-analysis/compiled_metabolites_knnimputed_longformat.csv")
# library(vegan)
# summary_stats <- all_metab %>% group_by(assay, sample_name)%>%
#                    summarise(avg = mean(peak_height_impmean), 
#                                         n_metab = sum(!is.na(peak_height_na)), 
#                                         n_metab_identified = sum(!is.na(peak_height_na) & !is.na(name)),
#                                         median = median(peak_height_impmean),
#                                         sd = sd(peak_height_impmean),
#                                         iqr = IQR(peak_height_impmean),
#                                         min_peak_height = min(peak_height_impmean),
#                                         max_peak_height = max(peak_height_impmean),
#                                         alpha_diver = diversity(peak_height_impmean),
#                                         norm_test_shapiro_pval = shapiro.test(peak_height_impmean)$p.value,
#                                         log_norm_test_shapiro_pval = shapiro.test(peak_height_log)$p.value
#                             ) %>%
#                    mutate(percent_identified = n_metab_identified/n_metab *100, avg_log2 = log2(avg))
# 
# #write.csv(summary_stats, "04-analysis/01_summarystatistics/summary_stats.csv")
# 
# #===========004 INITIAL FIGURES TO SEE DATA STRUCTURE===========#
# library(hrbrthemes)
# library(viridis)
# library(forcats)
# library(pheatmap)
# 
# #facetted histogram of the avg log2 intensity of the sample reads. This is to determine if there is a machine error for any of the samples.
# all_metab %>% ggplot(aes(x = sample_name, y = peak_height_na, fill = assay)) + 
#                 geom_boxplot() + 
#                 facet_grid(assay~.) +
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave("04-analysis/01_summarystatistics/box__avg_peakheight.png", width = 30, height = 10, limitsize = FALSE)
# 
# #clustered pearson correlation heatmap. To check if all biological replicates are similar to each other
# lipid.compact <- all_metab %>% filter(assay == "Lipidomics") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
#   pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
# amines.compact <- all_metab %>% filter(assay == "BiogenicAmines") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
#   pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
# gcms.compact <- all_metab %>% filter(assay == "PrimaryMetabolite") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
#   pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
# 
# lipid.corr <- lipid.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")
# amines.corr <- amines.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")
# gcms.corr <- gcms.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")
# 
# pheatmap(lipid.compact, scale = "row")
# ggsave("04-analysis/01_summarystatistics/heatmp__lipidomics.png", width = 20, height = 20, limitsize = FALSE)
#   
# #check if all features are normally distributed 
# 
# #grouped histogram of features that are missing from n% samples; to determine imputation strategy
# ##define that a biological replicate are mice with the same donor group
# ##x would be the features that are missing from n% samples, with the 3 assay grouped, y is the count
# missing_counts <- all_metab %>% group_by(UID, assay) %>% summarise(missing_count = sum(is.na(peak_height_na)))
# total_samples <- all_metab %>% group_by(UID) %>% summarise(total_samples = n_distinct(sample_code))
# missing_percentage <- missing_counts %>%
#   left_join(total_samples, by = "UID") %>%
#   mutate(percentage_missing = (missing_count / total_samples) * 100) %>%
#   select(UID, assay, percentage_missing) #table with columns=identifiers, percent missing, and assay out
# 
# ggplot(missing_percentage, aes(x=percentage_missing, color = assay, fill = assay)) +
#   geom_histogram(alpha = 0.6, binwidth = 5)+
#   scale_fill_viridis(discrete=TRUE) +
#   scale_color_viridis(discrete=TRUE) +
#   theme(
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     strip.text.x = element_text(size = 8)
#   ) + facet_wrap(~assay) + scale_x_continuous(trans = 'log10')
# 
# ggsave("04-analysis/01_summarystatistics/hist__percentagemissing.png", width = 7, height = 5)
# 
# #stripplot of assay and peak_heights_log; to see the data distribution for each
# ggplot(all_metab, aes(x=assay, y=peak_height_log)) + geom_jitter()
# 
# ggsave("04-analysis/01_summarystatistics/strip__assay_peak_heights_log.png", width = 7, height = 5)
# 
# #stripplot of each sample, with height is peak height, and color is the assay
# ggplot(all_metab, aes(x=sample_code, y=peak_height_log, color=assay)) + geom_jitter(position=position_dodge(0.5), alpha = 1/10)
# 
# ggsave("04-analysis/01_summarystatistics/strip__sample_peak_heights_log.png", width = 50, height = 8, limitsize = FALSE)
# 
# 
# #===========003 IMPUTATION#===========#
# ##ZERO
# all_metab$peak_height_impzero <- replace(all_metab$peak_height_na, which(is.na(all_metab$peak_height_na)), 0)
# 
# ##HALF LOWEST
# min_peak_per_sample <- all_metab %>% group_by(sample_code) %>% filter(!is.na(peak_height_na)) %>% 
#                         summarise(min_peak = min(peak_height_na)) %>% mutate(half_min_peak = min_peak/2)
# 
# all_metab <- all_metab %>% left_join(min_peak_per_sample, by = "sample_code") %>%
#               mutate(peak_height_imphalf = ifelse(is.na(peak_height_na), half_min_peak, peak_height_na))
# 
# ###Important columns for classifier-based imputation: m/z, retention time, peak_height, method, donor_type, gene_background
# all_metab_imput <- all_metab %>% select(one_of(c('m/z', 'ret. time', 'peak_height_na', 'assay', 'method', 'gene_background')))
# 
# ##RF (data too big)
# library(missForest)
# all_metab_imput <- all_metab_imput %>% mutate(`m/z` = str_extract(all_metab_imput$`m/z`, "[^_]+$")) 
# all_metab_imput <- all_metab_imput %>% mutate(`m/z` = as.numeric(all_metab_imput$`m/z`))#this is a temporary idea until can understand how to get neutral mass based on adduct
# all_metab_imput$assay <- as.factor(all_metab_imput$assay)
# all_metab_imput$method <- as.factor(all_metab_imput$method)
# all_metab_imput$gene_background <- as.factor(all_metab_imput$gene_background)
# 
# RFimput <- missForest::missForest(xmis = all_metab_imput)
# 
# ##KNN
# library(VIM)
# KNNimput <- kNN(all_metab_imput)
# 
# # Set the seed for reproducibility
# set.seed(42)
# 
# # Create a matrix of random values between 1 and 100
# data_matrix <- matrix(runif(6, min = 1, max = 100), nrow = 1)
# 
# # Create the dataframe with the desired row and column names
# df <- data.frame(data_matrix)
# rownames(df) <- c("feature1")
# colnames(df) <- c("S1", "S3", "S14", "S15", "QC1", "QC3")
# 
# 
# biosamples <- c("S1", "S3", "S14", "S15")
# qcsamples <- c("QC1", "QC3") 
# 
# reps <- factor(c(1:length(biosamples),rep(length(biosamples)+1,length(qcsamples))))
# 
# data <- data.frame(y=c(as.numeric(df[1,biosamples]),as.numeric(df[1,qcsamples])),reps=reps)
# 
# mm <- lme(y~ 1, random = ~1|reps, data=data,na.action = na.omit)
# 
# 
# 
# test <- d.filt.miss %>% filter(feat_id == "lip_esi_neg_392.24619999999999_0.864") %>% select(one_of("log_pkh_corrected", "donor_type")) %>% 
#   mutate(log_pkh_corrected = ifelse(is.na(log_pkh_corrected), 0, log_pkh_corrected))
# 
# mm_test <- lme(log_pkh_corrected~ 1, random = ~1|donor_type, data=test,na.action = na.omit)
# 
# varcorr <- as.numeric(VarCorr(mm_test)[1:2])
# 
# varcorr[1]/sum(varcorr)
# 
# ######
# 
# test_ori <- d.filt.miss %>% filter(feat_id == "lip_esi_neg_392.24619999999999_0.864") %>% select(one_of("log_pkh_corrected", "donor_type"))
# 
# mm_test_ori <- lme(log_pkh_corrected~ 1, random = ~1|donor_type, data=test_ori,na.action = na.omit)
# 
# varcorr_ori <- as.numeric(VarCorr(mm_test_ori)[1:2])
# 
# varcorr_ori[1]/sum(varcorr_ori)
