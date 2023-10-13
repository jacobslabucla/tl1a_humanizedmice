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
#     2.1. for each gene background, stool donor, sample loc, and sex (i.e. biological replicates), filter according to presence of feature in % of samples. do not cut out more than 10% of features found.
#       2.1.F. histogram of presence of feature in % of samples (00_31F_(lipid,gcms,amines)_hist__featurepercentpresence)
#     3.2. for each gene background, stool donor, sample loc, and sex (i.e. biological replicates), filter according to coef of variance (CV) of feature peak height (that is not NA) that is more than 50%.
#       3.2.F. boxplot, x = stool donor, y = feature variance (00_32F_(lipid,gcms,amines)_boxp__featurevariance_stooldonor)
#   3. imputation:
#     2.1. determine the missing value distribution after first data filtering
#     2.2. for each gene background, stool donor, and sex: impute NA values with kNN method
#   2. second data filtering (NOT IMPLEMENTED YET):
#     3.3. for each stool donor (i.e. biological replicates), filter according to percentile peak height value (try 2.5, 5, and 10th percentile). do not cut out more than 5% of features found.
#       3.3.F. boxplot, x = percentile peak height value cutoff, y = alpha_diversity (00_33F_(lipid,gcms,amines)_boxp__peakheightcutoff_adiversity). to determine the effects of absolute value cutoff towards adiversity
#   4. exploratory data analysis (EDA) for each method and sample:
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
library(ggplot2)

setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metabolomics")

#===========001 READING FILE AND DATA WRANGLING===========#

#read the raw excel files from service provider
lipid <- read_excel("03-data/raw_data/mx 738840_Jacob_Lipidomics_humanstool_mousefeces_Submit_07-2023.xlsx", sheet = 'clean')
amines <- read_excel("03-data/raw_data/mx 739116_Jacob_HILIC_biogenicamines_humanstool_09-2023 submit.xlsx", sheet = 'clean')
gcms <- read_excel("03-data/raw_data/mx 738540_Jacob_GCTOF_humanstool_mousestool_08-2023 submit.xlsx", sheet = 'clean')

#read the metadata files
metadata <- read_excel("03-data/raw_data/metadata.xlsx", sheet = 'clean')

#check if all mass are unique
which(duplicated(lipid$identifier) == TRUE)
which(duplicated(amines$identifier) == TRUE)
which(duplicated(gcms$identifier) == TRUE)

#make a new id incorporating ESI mode
lipid <- lipid %>% unite(identifier, esi_mode, sep = "_", col = "UID")
amines <- amines %>% unite(identifier, esi_mode, sep = "_", col = "UID")
gcms <- gcms %>% unite(identifier, esi_mode, sep = "_", col = "UID")

##IDEA: FILTER OUT THE UNIDENTIFIED ONES

#generate 'long' dataframe for each
lipid_lon <- lipid %>% pivot_longer(cols = names(lipid[7:274]), 
                                    names_to='sample_id',
                                    values_to='peak_height')
amines_lon <- amines %>% pivot_longer(cols = names(amines[8:275]), 
                                     names_to='sample_id',
                                     values_to='peak_height')
gcms_lon <- gcms %>% pivot_longer(cols = names(gcms[11:278]), 
                                   names_to='sample_id',
                                   values_to='peak_height')

#check if all samples are present in all methods
setdiff(unique(lipid_lon$sample_id), unique(amines_lon$sample_id)) #0 diff, meaning all same sample set
setdiff(unique(lipid_lon$sample_id), unique(gcms_lon$sample_id)) #0 diff, meaning all same samples set

#add assay and method column
lipid_lon$assay <- "Lipidomics"
amines_lon$assay <- "BiogenicAmines"
gcms_lon$assay <- "PrimaryMetabolite"

lipid_lon$method <- "LCMS"
amines_lon$method <- "LCMS"
gcms_lon$method <- "GCMS"

lipid_lon$gene_background <- with(lipid_lon, ifelse(grepl('BSD', sample_id),'Human',ifelse(grepl('TL1A', sample_id),'TL1A_mice', 'IL10_mice')))
amines_lon$gene_background <- with(amines_lon, ifelse(grepl('BSD', sample_id),'Human',ifelse(grepl('TL1A', sample_id),'TL1A_mice', 'IL10_mice')))
gcms_lon$gene_background <- with(gcms_lon, ifelse(grepl('BSD', sample_id),'Human',ifelse(grepl('TL1A', sample_id),'TL1A_mice', 'IL10_mice')))

#lookup the donor stool value, incl ibd status and otu type
#check first if all sample_id in data matches sample_id in metadata
setdiff(unique(metadata$sample_id), unique(lipid_lon$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted
setdiff(unique(metadata$sample_id), unique(amines_lon$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted
setdiff(unique(metadata$sample_id), unique(gcms_lon$sample_id)) #proxcolon, duod, jej, and ile didn't get submitted

lipid_lon <- lipid_lon %>% left_join(metadata, by = "sample_id")
amines_lon <- amines_lon %>% left_join(metadata, by = "sample_id")
gcms_lon <- gcms_lon %>% left_join(metadata, by = "sample_id")

#change all 0s and negative values to NA; negative value is from subtracting with blank ctrl signal;
#useful to do before doing log transform to prevent error
lipid_lon$peak_height_na <- replace(lipid_lon$peak_height, which(lipid_lon$peak_height <= 0), NA)
amines_lon$peak_height_na <- replace(amines_lon$peak_height, which(amines_lon$peak_height <= 0), NA)
gcms_lon$peak_height_na <- replace(gcms_lon$peak_height, which(gcms_lon$peak_height <= 0), NA)

#===========002 FIRST FILTERING===========#
#checking the amount of biological replicates for each gene background, sampling location, donor stool, and sex
lipid_nreplicate <- lipid_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates
amines_nreplicate <- amines_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates
gcms_nreplicate <- gcms_lon %>% group_by(gene_background, location, donor) %>% summarize(n_repl = n_distinct(sample_id)) #most are 4-5 replicates

#filter features that are only present in <50% of the samples (pretty conservative filter here)
lipid_filtering <- lipid_lon %>% group_by(gene_background, location, donor, UID) %>% 
                          summarize(tot_na = sum(is.na(peak_height_na)), 
                                    identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
                                    n_repl = n_distinct(sample_id)) %>%
                          mutate(percent_presence_across_replicates = 1 - (tot_na / n_repl))
lipid_filt_stats <- lipid_filtering %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_replicates > 0),
                                                                             feat_less50 = sum(percent_presence_across_replicates < 0.5)) %>%
                                                                    mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 10%
lipid_lon_filtered <- lipid_lon %>% left_join(lipid_filtering, by = c('UID', 'gene_background', 'location', 'donor')) %>% filter(percent_presence_across_replicates>=0.5)

amines_filtering <- amines_lon %>% group_by(gene_background, location, donor, UID) %>% 
                                    summarize(tot_na = sum(is.na(peak_height_na)), 
                                              identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
                                              n_repl = n_distinct(sample_id)) %>%
                                    mutate(percent_presence_across_replicates = 1 - (tot_na / n_repl))
amines_filt_stats <- amines_filtering %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_replicates > 0),
                                                                             feat_less50 = sum(percent_presence_across_replicates < 0.5)) %>%
                      mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 10%
amines_lon_filtered <- amines_lon %>% left_join(amines_filtering, by = c('UID', 'gene_background', 'location', 'donor')) %>% filter(percent_presence_across_replicates>=0.5)

gcms_filtering <- gcms_lon %>% group_by(gene_background, location, donor, UID) %>% 
                                summarize(tot_na = sum(is.na(peak_height_na)), 
                                          identified = ifelse(sum(is.na(name)) > 0, FALSE, TRUE), 
                                          n_repl = n_distinct(sample_id)) %>%
                                mutate(percent_presence_across_replicates = 1 - (tot_na / n_repl))
gcms_filt_stats <- gcms_filtering %>% group_by(gene_background, location, donor) %>% summarize(n_feat = sum(percent_presence_across_replicates > 0),
                                                                             feat_less50 = sum(percent_presence_across_replicates < 0.5)) %>%
                                        mutate(percent_less_than_50 = feat_less50 / n_feat)#check the number of features that are below 10%
gcms_lon_filtered <- gcms_lon %>% left_join(gcms_filtering, by = c('UID', 'gene_background', 'location', 'donor')) %>% filter(percent_presence_across_replicates>=0.5)

#===========002 IMPUTATION#===========#
##MEAN Imputation (seems to perform better than zero or half lowest and easy to implement according to Kokla, Virtanen, et al., 2019)
# mean_peak_per_sample <- lipid_lon_filtered %>% group_by(sample_id) %>% filter(!is.na(peak_height_na)) %>%
#   summarise(mean_peak = mean(peak_height_na))
# lipid_lon_filtered <- lipid_lon_filtered %>% left_join(mean_peak_per_sample, by = "sample_id")  %>%
#   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))
# 
# mean_peak_per_sample <- amines_lon_filtered %>% group_by(sample_id) %>% filter(!is.na(peak_height_na)) %>%
#   summarise(mean_peak = mean(peak_height_na))
# amines_lon_filtered <- amines_lon_filtered %>% left_join(mean_peak_per_sample, by = "sample_id")  %>%
#   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))
# 
# mean_peak_per_sample <- gcms_lon_filtered %>% group_by(sample_id) %>% filter(!is.na(peak_height_na)) %>%
#   summarise(mean_peak = mean(peak_height_na))
# gcms_lon_filtered <- gcms_lon_filtered %>% left_join(mean_peak_per_sample, by = "sample_id")  %>%
#   mutate(peak_height_impmean = ifelse(is.na(peak_height_na), mean_peak, peak_height_na))

##KNN Imputation
#impute with these columns = gene background, otu type, ibd status, peak height na, mz, r.time, identifier as row name
library(VIM)

test_knn <- lipid_lon_filtered %>% filter(gene_background == "TL1A_mice")%>% 
              select(one_of(c("UID", "peak_height_na", "donor", "sample_id"))) %>%
              mutate(peak_height_na_log = log(peak_height_na)) %>% 
              left_join(kNN(.,variable = "peak_height_na_log"), by = c("UID","sample_id"))

test_knn2 <- kNN(test_knn,variable = "peak_height_na")


lipid_knn <- lipid_lon_filtered %>% select(one_of(c("UID","gene_background", "donor", "peak_height_na"))) %>% 
              filter(gene_background != "Human") %>% group_by(gene_background, donor) %>% 
              mutate(peak_height_knn = kNN(., variable = "peak_height_na"))


#create log transformed data
lipid_lon <- lipid_lon %>% mutate(peak_height_log = log(peak_height_impmean, 10))
amines_lon <- amines_lon %>% mutate(peak_height_log = log(peak_height_impmean, 10))
gcms_lon <- gcms_lon %>% mutate(peak_height_log = log(peak_height_impmean, 10))

#combine the different assay output into one big table
base_cols <- intersect(names(lipid_lon), names(amines_lon))

all_metab <- merge(gcms_lon, merge(lipid_lon, amines_lon, by = base_cols, all.x = TRUE, all.y = TRUE), 
                   by = base_cols, all.x = TRUE, all.y = TRUE)

#write.csv(all_metab, "04-analysis/compiled_metabolites_longformat.csv")

#===========003 EXPLORATORY DATA ANALYSIS===========#
#all_metab <- read_csv("04-analysis/compiled_metabolites_longformat.csv")
library(vegan)
summary_stats <- all_metab %>% group_by(assay, sample_name)%>%
                   summarise(avg = mean(peak_height_impmean), 
                                        n_metab = sum(!is.na(peak_height_na)), 
                                        n_metab_identified = sum(!is.na(peak_height_na) & !is.na(name)),
                                        median = median(peak_height_impmean),
                                        sd = sd(peak_height_impmean),
                                        iqr = IQR(peak_height_impmean),
                                        min_peak_height = min(peak_height_impmean),
                                        max_peak_height = max(peak_height_impmean),
                                        alpha_diver = diversity(peak_height_impmean),
                                        norm_test_shapiro_pval = shapiro.test(peak_height_impmean)$p.value,
                                        log_norm_test_shapiro_pval = shapiro.test(peak_height_log)$p.value
                            ) %>%+
                   mutate(percent_identified = n_metab_identified/n_metab *100, avg_log2 = log2(avg))

#write.csv(summary_stats, "04-analysis/01_summarystatistics/summary_stats.csv")

#===========004 INITIAL FIGURES TO SEE DATA STRUCTURE===========#
library(hrbrthemes)
library(viridis)
library(forcats)
library(pheatmap)

#facetted histogram of the avg log2 intensity of the sample reads. This is to determine if there is a machine error for any of the samples.
all_metab %>% ggplot(aes(x = sample_name, y = peak_height_na, fill = assay)) + 
                geom_boxplot() + 
                facet_grid(assay~.) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("04-analysis/01_summarystatistics/box__avg_peakheight.png", width = 30, height = 10, limitsize = FALSE)

#clustered pearson correlation heatmap. To check if all biological replicates are similar to each other
lipid.compact <- all_metab %>% filter(assay == "Lipidomics") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
amines.compact <- all_metab %>% filter(assay == "BiogenicAmines") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')
gcms.compact <- all_metab %>% filter(assay == "PrimaryMetabolite") %>% select(one_of(c("sample_name", "UID", "peak_height_impmean"))) %>%
  pivot_wider(names_from = "UID", values_from = "peak_height_impmean") %>% column_to_rownames(var = 'sample_name')

lipid.corr <- lipid.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")
amines.corr <- amines.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")
gcms.corr <- gcms.compact %>% t() %>% cor(use = "pairwise.complete.obs", method = "pearson")

pheatmap(lipid.compact, scale = "row")
ggsave("04-analysis/01_summarystatistics/heatmp__lipidomics.png", width = 20, height = 20, limitsize = FALSE)
  
#check if all features are normally distributed 

#grouped histogram of features that are missing from n% samples; to determine imputation strategy
##define that a biological replicate are mice with the same donor group
##x would be the features that are missing from n% samples, with the 3 assay grouped, y is the count
missing_counts <- all_metab %>% group_by(UID, assay) %>% summarise(missing_count = sum(is.na(peak_height_na)))
total_samples <- all_metab %>% group_by(UID) %>% summarise(total_samples = n_distinct(sample_code))
missing_percentage <- missing_counts %>%
  left_join(total_samples, by = "UID") %>%
  mutate(percentage_missing = (missing_count / total_samples) * 100) %>%
  select(UID, assay, percentage_missing) #table with columns=identifiers, percent missing, and assay out

ggplot(missing_percentage, aes(x=percentage_missing, color = assay, fill = assay)) +
  geom_histogram(alpha = 0.6, binwidth = 5)+
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) + facet_wrap(~assay) + scale_x_continuous(trans = 'log10')

ggsave("04-analysis/01_summarystatistics/hist__percentagemissing.png", width = 7, height = 5)

#stripplot of assay and peak_heights_log; to see the data distribution for each
ggplot(all_metab, aes(x=assay, y=peak_height_log)) + geom_jitter()

ggsave("04-analysis/01_summarystatistics/strip__assay_peak_heights_log.png", width = 7, height = 5)

#stripplot of each sample, with height is peak height, and color is the assay
ggplot(all_metab, aes(x=sample_code, y=peak_height_log, color=assay)) + geom_jitter(position=position_dodge(0.5), alpha = 1/10)

ggsave("04-analysis/01_summarystatistics/strip__sample_peak_heights_log.png", width = 50, height = 8, limitsize = FALSE)


#===========003 IMPUTATION#===========#
##ZERO
all_metab$peak_height_impzero <- replace(all_metab$peak_height_na, which(is.na(all_metab$peak_height_na)), 0)

##HALF LOWEST
min_peak_per_sample <- all_metab %>% group_by(sample_code) %>% filter(!is.na(peak_height_na)) %>% 
                        summarise(min_peak = min(peak_height_na)) %>% mutate(half_min_peak = min_peak/2)

all_metab <- all_metab %>% left_join(min_peak_per_sample, by = "sample_code") %>%
              mutate(peak_height_imphalf = ifelse(is.na(peak_height_na), half_min_peak, peak_height_na))

###Important columns for classifier-based imputation: m/z, retention time, peak_height, method, donor_type, gene_background
all_metab_imput <- all_metab %>% select(one_of(c('m/z', 'ret. time', 'peak_height_na', 'assay', 'method', 'gene_background')))

##RF (data too big)
library(missForest)
all_metab_imput <- all_metab_imput %>% mutate(`m/z` = str_extract(all_metab_imput$`m/z`, "[^_]+$")) 
all_metab_imput <- all_metab_imput %>% mutate(`m/z` = as.numeric(all_metab_imput$`m/z`))#this is a temporary idea until can understand how to get neutral mass based on adduct
all_metab_imput$assay <- as.factor(all_metab_imput$assay)
all_metab_imput$method <- as.factor(all_metab_imput$method)
all_metab_imput$gene_background <- as.factor(all_metab_imput$gene_background)

RFimput <- missForest::missForest(xmis = all_metab_imput)

##KNN
library(VIM)
KNNimput <- kNN(all_metab_imput)




