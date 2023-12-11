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

saveWidgetFix <- function (widget,file,...) {
  ## A wrapper to saveWidget which compensates for arguable BUG in
  ## saveWidget which requires `file` to be in current working
  ## directory.
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  saveWidget(widget,file=file,...)
}

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

##=====DATA VIS=====##

##volcanoplots to visualize the results of differential abundance analysis
# 
# volplot_linmodel <- function(linmodel_output){
#   d <- linmodel_output
#   
#   assay <- unique(d$assay)
#   gene_bg <- unique(d$gene_background)
#   
#   d.filtered <- d %>%
#                 select(one_of(c("coef", "qval", "feature", "nameFix", "metadata", "value"))) %>%
#                 mutate(log10_pval_adj = -log10(qval),
#                        upordown = ifelse(coef > 0, "Up", "Down"),
#                        log2_FC = coef / log10(2)) %>%
#                 mutate(upordown = ifelse(log10_pval_adj <1.3, "Not Significant", upordown))
#   
#   write.csv(d.filtered, 
#             paste0("04-analysis/01_mseainput/", assay, "_", gene_bg, "_sighit_maaslin2.csv"),
#             row.names = FALSE)
#   
#   # top5up <- deg.filtered %>% 
#   #   filter(log10_pval_adj > 1.3 & upordown == "Up") %>%
#   #   top_n(5, log10_pval_adj)
#   # 
#   # top5down <- deg.filtered %>% 
#   #   filter(log10_pval_adj > 1.3 & upordown == "Down") %>%
#   #   top_n(5, log10_pval_adj)
#   # 
#   # topupdown <- bind_rows(top5down, top5up)
#   
#   for (x in unique(d$metadata)){
#     
#     if (x == "OTU_type"){
#       title = "OTU type 1 vs OTU type 2"
#     }
#     else if (x == "IBD_stat") {
#       title = "Healthy vs IBD"
#     }
#     else if (x == "mouse_sex") {
#       title = "Female vs Male"
#     }
#     else {
#       title = "Interaction Term OTU type 2 - IBD"
#     }
#     
#     d.filtered1 <- d.filtered %>% filter(metadata == x)
#     
#     p <- ggplot(d.filtered1, aes(x = log2_FC, y = log10_pval_adj, color = upordown, text = nameFix)) +
#         geom_point() +
#         theme_bw() +
#         labs(x = "log2FC", y = "-log10(p-value)", color = NULL) + 
#         scale_color_manual(values = c("red", "blue", "grey"), 
#                            limits = c("Up", "Down", "Not Significant"))+
#         ggtitle(paste0("Volcano plot \n", assay, "_", gene_bg, "\n", title))
#     
#     gp <- ggplotly(p)
#     
#     saveWidgetFix(gp, file = paste0("04-analysis/01_volplot_maaslin2/", assay, "_", gene_bg, "_", x, "_volcanoplot_maaslin2.html"))
#   }
# }
# 
# volplot_linmodel(lipid_tl1a_maaslin)
# volplot_linmodel(amines_tl1a_maaslin)
# volplot_linmodel(gcms_tl1a_maaslin)
# volplot_linmodel(lipid_il10_maaslin)
# volplot_linmodel(amines_il10_maaslin)
# volplot_linmodel(gcms_il10_maaslin)

metab_tl1a_maaslin_wointeractinfl <- amines_tl1a_maaslin_wointeractinfl %>%
  bind_rows(bind_rows(lipid_tl1a_maaslin_wointeractinfl, gcms_tl1a_maaslin_wointeractinfl)) %>%
  filter(qval <= 0.05) %>%
  pivot_wider(names_from = value, values_from = qval) %>%
  select(one_of(c("feature","nameFix", "coef", "otu_2", "IBD")))

metab_tl1a_maaslin_onlyinfl <- amines_tl1a_maaslin_onlyinfl %>% 
  bind_rows(bind_rows(lipid_tl1a_maaslin_onlyinfl, gcms_tl1a_maaslin_onlyinfl)) %>%
  filter(value == "infl_score", qval <= 0.1) %>%
  mutate(nameFix = reorder(factor(nameFix), coef)) %>%
  select(one_of(c("feature","nameFix", "coef", "qval")))

combined_tl1a <- metab_tl1a_maaslin_wointeractinfl %>% left_join(metab_tl1a_maaslin_onlyinfl, by = c("feature", "nameFix"))

# write.csv(metab_tl1a_maaslin_onlyinfl, "inflamationrelatedmetab2.csv")
metab_tl1a_annot <- read.csv("inflamationrelatedmetab.csv") %>% 
  mutate(nameFix = reorder(factor(nameFix), coef), sig = -log10(qval))

p <- ggplot(data = metab_tl1a_annot, aes(x = coef, y = nameFix, size = -qval, color = taxonomy2)) +
        geom_point(alpha = .65) +
        theme_classic() +
        labs(x = "Correlation with inflammation",
             y = "") +
        scale_size_continuous(name = "Adj. p-value", breaks = c(-0.005, -0.05, -0.1), labels = c(0.005, 0.05, 0.1)) +
        scale_color_discrete(name = "Chemical Taxonomy") +
        geom_vline(xintercept = 0, linetype="dotted", linewidth=.5)

ggsave("inflamationrelatedmetab.pdf", width =10, height = 6)

##==================zeroinflated gaussian======================##
library()


##==================DASEV======================##
#Get prior distribution for sigma using MLE with group means seperated
Getprior <- function(params.prior, lod, covars, yvec) {
  #print(params.prior)
  n.Params <- ncol(covars)
  mu.all <- covars%*%params.prior[1:n.Params]
  sd0 <- exp(params.prior[n.Params+1])
  mu.nonzero <- mu.all[yvec>-Inf]
  nonzero <- yvec[yvec>-Inf]
  A <- pnorm(lod, mu.nonzero, sd0) # control the CDF at detection limit is less than 1
  if (sum(A>0.9999999999)>0){
    mle= -100000000000000
  }else{
    mle <- -sum(log((1 - A)*sqrt(2*pi)*sd0)+(nonzero-mu.nonzero)^2/(2*sd0^2))
  }
  return(-mle)
}



loglikestep1 <- function(Params, lsd, selectmodel, lod, s0, d0, nbeta, cov.matrix, test_cov, yvec){
  n.betas <- 1:nbeta #mean paras
  n.gammas <- (nbeta+1):length(Params)#zero paras
  betas<-Params[n.betas]
  gammas <-Params[n.gammas]
  
  if(selectmodel=="F"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="N"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas,length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  if(selectmodel=="M"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="P"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  mu.nonzero <- mu.all[yvec>-Inf]
  p.nonzero <- p.all[yvec>-Inf]
  mu.zero <- mu.all[yvec==-Inf]
  p.zero <- p.all[yvec==-Inf]
  nonzero <- yvec[yvec>-Inf] # non-zero obs
  zero <- yvec[yvec==-Inf] #zero obs
  sd <- exp(lsd)
  if (length(mu.zero)>0){
    A <- pnorm(lod, mu.zero, sd) # control the CDF at detection limit is less than 1 #A[A=="NaN"]<-1
    #  print(sum(A>0.5)>0)
    if (sum(A>0.99999999999999)>0){
      logl= -100000000000000
    }else{
      loglzero <- sum(log(p.zero + (1 - p.zero)*A))
      loglnonzero <- sum(log((1 - p.nonzero))) - sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
      logl <- loglzero + loglnonzero #+ loglprior
    }
  }else{
    loglnonzero <- sum(log((1 - p.nonzero))) - sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
    logl <- loglnonzero #+ loglprior
  }
  return <- (-logl)
}


loglikestep2 <- function(Params, betagamma, selectmodel, lod, s0, d0,
                         nbeta, cov.matrix, test_cov, yvec){
  n.betas <- 1:nbeta #mean paras
  n.gammas <- (nbeta+1):length(betagamma)#zero paras
  betas<-betagamma[n.betas]
  gammas <-betagamma[n.gammas]
  
  if(selectmodel=="F"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="N"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  if(selectmodel=="M"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="P"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  
  mu.nonzero <- mu.all[yvec>-Inf]
  p.nonzero <- p.all[yvec>-Inf]
  mu.zero <- mu.all[yvec==-Inf]
  p.zero <- p.all[yvec==-Inf]
  nonzero <- yvec[yvec>-Inf] # non-zero obs
  zero <- yvec[yvec==-Inf] #zero obs
  sd <- exp(Params)
  
  if (length(mu.zero)>0){
    A <- pnorm(lod, mu.zero, sd) # control the CDF at detection limit is less than 1
    if (sum(A>0.999999999999)>0){
      logl= -100000000000000
    }else{
      loglzero <- sum(log(p.zero + (1 - p.zero)*A))
      loglnonzero <- sum(log((1 - p.nonzero))) -sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
      loglprior <- (d0/2)*log(d0*s0^2/2) - (d0/2 +1)*log(sd^2) - log(gamma(d0/2)) - d0*s0^2/(2*sd^2)
      logl <- loglzero + loglnonzero + loglprior
    }
  }else{
    loglnonzero <- sum(log((1 - p.nonzero))) -sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
    loglprior <- (d0/2)*log(d0*s0^2/2) - (d0/2 +1)*log(sd^2) - log(gamma(d0/2)) - d0*s0^2/(2*sd^2)
    logl <- loglnonzero + loglprior
  }
  
  return <- (-logl)
}


#' Differential abundance analysis for zero-inflated data.
#'
#' This function processes the two group differential abundance (DA) analysis for zero-inflated data.
#' It has the flexibility to let users to choose various detection limits (DL). The function returns
#' estimates on group means, group biological point mass value (BPMV) proportions, standard deviation,
#' and p-values.
#'
#' @param indata Specifies the input data matrix. Rows are features and columns are subjects.
#' @param cov.matrix A matrix for variables. Categorical data should be represented by dummy variables
#' contain only 0 and 1.
#' @param test_cov A vector corresponding to the coloums in \code{cov.matrix} which specify the
#' variables to be tested.
#' @param test_cov_conti A vector corresponding to the coloums in \code{cov.matrix} which specify the
#' continuous variables.
#' @param min.non0n The minimum number of nonzero observations for features to be included in the analysis.
#' Default and the minimum value is 3. Please input a number equal or larger than 3.
#' @param requiredn The minimum number of nonzero obs while getting prior distribution of variance.
#' Default value is 10.
#' @param requiredn2 The minimum number of features while getting prior distribution of variance.
#' Default value is 30.\cr
#' If requiredn is specified, it will retrun features used in getting prior distribution of variance.
#' If the number of features used is less than requiredn2, the function returns top requiredn2
#' features with the smallest all point mass value (PMV) proportions.
#' @param DL_method Specifies the detection limit method. The options are:\cr
#' Fixed Difference (default): for each feature, DL is the minimun value for all nonzero observation
#' minus a number sepcified by DL_value\cr
#' Fixed Rate: for each feature, DL is the minimun value for all nonzero observation devided by
#' a number sepcified by DL_value\cr
#' Fixed Value: DL is the same value (a number sepcified by DL_value) for all features.
#' @param DL_value Custom specified number. Default is 0.1 for "Fixed Difference".
#' @param maxit_MLE The maximum number of iterations while re-estimating the model parameters.
#' The default is 100.
#' @param maxit The maximum number of iterations applied to the optimization function.
#' @param test_model A vector of models to be tested. The default value is c("Both","Mean","Pzero").\cr
#' The  default is to achieve all three comparsions listed as follow:\cr
#' If test_model contains "Both", the function compares both difference in group means and BPMV proportions.\cr
#' If test_model contains "Mean", the function compares difference in group means.\cr
#' If test_model contains "Pzero", the function compares difference in BPMV proportions.
#'
#' @details
#' The optimization function used here is \code{optim} with the method option as \code{"BFGS"}.
#'
#' If using "Fixed Difference" as the method to get detection limit, a small positive value is
#' recommended for \code{DL_value}.\cr
#' If using "Fixed Rate" as the method to get detection limit, a positive value larger than 1
#' is recommended for \code{DL_value}.\cr
#' If using "Fixed Value" as the method to get detection limit, a value smaller than the log value
#' of the minimum observation in the dataset is recommended for \code{DL_value}.
#'
#' This function requires at least one zero and one non-zero observations in each group.
#'
#' @return \item{feature_names}{A vector of feature names substracted from the input data.}
#' @return \item{pvalue_both}{A vector of estimated p-values for comparing both difference in
#' group means and BPMV proportions. pvalue_both will be NULL if test_model doesn't contain "Both".}
#' @return \item{pvalue_mean}{A vector of estimated p-values for comparing only difference
#' in group means. pvalue_mean will be NULL if test_model doesn't contain "Mean".}
#' @return \item{pvalue_zero}{A vector of estimated p-values for comparing only difference
#' in group BPMV proportions.pvalue_zero will be NULL if test_model doesn't contain "Pzero".}
#' @return \item{DL}{A vector of detection limits used.}
#' @return \item{estimates}{A matrix of estimates on optimization parameters, which can be
#' used to calculate group means, group BPMV proportions, and standard deviation.}
#'
#' @export
#' @seealso \code{\link{Getsample}}
#' @examples
#' #Get simulation samples#
#' data(simpool)
#' sim <- Getsample(numc=1000,numobs=100,pdiff=0.2,lfc=c(log(2),-log(2)),pzerodiff=NULL, simpool)
#'
#' data<- sim$simdata
#' paradata <- sim$Parameters
#' indall <- c()
#' for(i in 1:nrow(data)){
#'   ind <- (!all(data[i,][1:100]!=0)
#'           &!all(data[i,][1:100]==0)
#'           &!all(data[i,][101:200]!=0)
#'           &!all(data[i,][101:200]==0))
#'   indall<- c(indall, ind)
#' }
#' data_used <- data[indall,]
#'
#' example_result <- DASEV(indata=data_used,cov.matrix=cbind(rep(1,200),c(rep(0,100),rep(1,100))),
#' test_cov= 2,test_cov_conti=NULL, min.non0n=3, requiredn=10, requiredn2=30,
#' DL_method= "Fixed Difference",DL_value=0.1,maxit_MLE=100,maxit=1000,
#' test_model=c("Both","Mean","Pzero"))

DASEV <- function(indata, cov.matrix, test_cov, test_cov_conti=NULL, min.non0n=3, requiredn=10,
                  requiredn2=30, DL_method= "Fixed Difference", DL_value=0.1,
                  maxit_MLE=100, maxit=10000, test_model=c("Both", "Mean", "Pzero")){
  if (!(DL_method == "Fixed Difference"|DL_method == "Fixed Rate"|DL_method ==
        "Fixed Value")){
    stop("Please enter the correct method to calculate the detect limit.")
  }
  
  if (is.na(DL_value)){
    stop("Please enter the correct value to calculate the detect limit.")
  }
  
  if (DL_method == "Fixed Rate" & DL_value < 1) {
    message("You are using the fixed rate method to calculate the detect limit,
            your input DL_value is less than 1, if this is not correct, please
            stop the program and check your input.")
  }
  
  if(min.non0n<3){
    stop("The minimum number of nonzero observation should be 3.")
  }
  
  #get data ready based on feature inclusion criteria
  n.ParamsF <- ncol(as.matrix(cov.matrix))
  n.ParamsN <- ncol(as.matrix(cov.matrix[, -test_cov]))
  ldata <- as.matrix(indata)
  idint<- rowSums(ldata> 0) >= min.non0n 
  ldata <- ldata[idint,] #to filter features that have excessive 0s
  
  check_cov <- length(test_cov) - length(test_cov_conti)
  if (length(test_cov_conti)>0){
    test_cov_dummy <- test_cov[-test_cov_conti]
  }else{
    test_cov_dummy <- test_cov
  }
  if (check_cov > 0){
    if(length(test_cov_dummy) > 1){
      test.cov <- cov.matrix[, test_cov_dummy]
      fslc <- function(test.cov){
        f.slc <- function(ldata){
          ind<- (!all(ldata[test.cov==0]!=0)
                 &!all(ldata[test.cov==0]==0)
                 &!all(ldata[test.cov==1]!=0)
                 &!all(ldata[test.cov==1]==0))
        }
        indall <-apply(ldata, 1, f.slc)
      }
      allind <- apply(test.cov, 2, fslc)
      allindsub <- rowSums(allind==TRUE)==ncol(allind)
    }else {
      f.slc <- function(ldata_row){
              ind <- (!all(ldata_row[test.cov==0]!=0) #for feature A in healthy patients, are all not 0?
                     &!all(ldata_row[test.cov==0]==0) #for feature A in healthy patients, are all 0?
                     &!all(ldata_row[test.cov==1]!=0) #for feature A in IBD patients, are all not 0?
                     &!all(ldata_row[test.cov==1]==0)) #for feature A in IBD patients, are all 0?
      }
      test.cov <- cov.matrix[, test_cov]
      allindsub <-apply(ldata, 1, f.slc)
    }
    ldata2 <- ldata[allindsub,]
  }
  
  #ldata <- log(ldata) #######TURN THIS OFF COS WE ALREADY HAVE LOG TRANSFORMED DATA
  feature.names <- rownames(ldata)
  
  varprior <- c()
  for (i in 1:nrow(ldata)) {
    yvec <- ldata[i,]
    n <- (yvec > -Inf)
    covar.nonzero <- cov.matrix[n,]
    estbeta <- lm(yvec[n]~0+covar.nonzero)
    betasF <- coefficients(estbeta)
    betasF[is.na(betasF)] <- 0
    lsd <- log(sd(yvec[n]))
    if (DL_method == "Fixed Difference") {
      lod<-min(yvec[n])-DL_value
    } else if (DL_method == "Fixed Rate") {
      lod<-min(yvec[n])/DL_value
    } else if (DL_method == "Fixed Value") {
      lod <- DL_value
    } else {
      stop("Please enter the correct method to calculate the detect limit.")
    }
    
    mleest <- optim(par= c(betasF,lsd), fn=Getprior, yvec=yvec, lod=lod, covars=cov.matrix,
                    method="BFGS", control = list(maxit = maxit))
    if(mleest$convergence!=0)
      cat("\nGet prior distribution Convergence problems for feature ", i)
    varprior <- c(varprior, (exp(mleest$par[length(mleest$par)]))^2)
  }
  
  idint1<- rowSums(ldata> 0) >= requiredn
  varprior.sub <- varprior[idint1]
  count <- rowSums(ldata> -Inf)
  if(sum(idint1) < requiredn2){
    sort.varprior <- varprior[order(-count)]
    varprior.sub <- sort.varprior[1:requiredn2]
  }
  p.mean1 <- mean(varprior.sub)
  p.var1 <- var(varprior.sub)
  d0 <- 2*p.mean1^2/p.var1 + 4
  s0 <- sqrt(p.mean1*(d0-2)/d0)
  
  
  
  
  ################################################################################
  # Getting post statistics
  Para.Full <- c()
  MLE.Full <- c()
  MLE.Null <- c()
  MLE.Mean <- c()
  MLE.Pzero <- c()
  
  Pvalue <- c() #pvalue for Full model vs null model
  Pmean <-c() #pvalue for Full model vs same mean
  Pzero <-c() #pvalue for Full model vs same zero proportion
  
  DL <- c()
  for (i in 1:nrow(ldata)) {
    #print(i)
    yvec <- ldata[i,]
    n <- (yvec>-Inf)
    covar.nonzero <- cov.matrix[n,]
    if (sum(covar.nonzero[, 1] == covar.nonzero[, 2]) == nrow(covar.nonzero)){
      covar.nonzero[, 1]=1-covar.nonzero[, 1]
    }
    estbeta <- lm(yvec[n]~0+covar.nonzero)
    betasF <- coefficients(estbeta)
    betasF[is.na(betasF)] <- 0
    
    if (DL_method == "Fixed Difference") {
      lod<-min(yvec[n])-DL_value
    } else if (DL_method == "Fixed Rate") {
      lod<-min(yvec[n])/DL_value
    } else if (DL_method == "Fixed Value") {
      lod <- DL_value
    } else {
      stop("Please enter the correct method to calculate the detect limit.")
    }
    DL <- c(DL, lod)
    lsd0 <- log(sqrt(varprior[i]))# initial value for variance parameters
    mu <- mean(yvec[n])
    p <- 1 - length(yvec[n])/length(yvec)/(1-pnorm(lod, mu, exp(lsd0)))#true zero proportion starting value
    if (p>0){
      gamma0 <- log((1-p)/p)
    }else{
      gamma0 <- 1e+200
    }
    gammasF <-  c(gamma0, rep(0, ncol(cov.matrix)-1))
    
    covar.nonzero <- cov.matrix[n,]
    estbetaN <- lm(yvec[n]~0+covar.nonzero[, -test_cov])
    betasN <- coefficients(estbetaN)
    betasN[is.na(betasN)] <- 0
    gammasN <- gammasF[-test_cov]
    
    xF <- c(betasF, gammasF)#Full model params
    xN <- c(betasN, gammasN)#Null model params
    xM <- c(betasN, gammasF)#Params for model with different zero percentages but same means
    xP <- c(betasF, gammasN)#Params for model with different means but same zero percentages
    
    
    #Full Model
    nbeta <- length(betasF)
    Diff <- 1
    time <- 0
    lsd <- lsd0
    xFP <- xF
    while (Diff > 0.01 & time < maxit_MLE){
      time <- time +1
      outF <- optim(xFP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="F",
                    lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                    cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
      xFP <- outF$par
      if(outF$convergence!=0)
        cat("\nStep 1 Get post statistics full model Convergence problems for feature ",i)
      out <- optim(lsd, fn=loglikestep2, selectmodel="F", method = "BFGS", betagamma=xFP,
                   lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                   cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
      Diff <- abs(out$par-lsd)
      lsd <- out$par
      if(out$convergence!=0)
        cat("\nStep 2 Get post statistics full model Convergence problems for feature ",i)
    }
    
    Para.Full <- rbind(Para.Full, c(outF$par, exp(out$par)))
    MLE.Full <- c(MLE.Full, outF$value)
    
    #Null model
    if ("Both" %in% test_model){
      nbeta <- length(betasN)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xNP <- xN
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outN <- optim(xNP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="N",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xNP <- outN$par
        if(outN$convergence!=0)
          cat("\nStep 1 Get post statistics null model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="N", method = "BFGS", betagamma=xNP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics null model Convergence problems for feature ",i)
        
      }
      MLE.Null<-c(MLE.Null, outN$value)
    }else{
      MLE.Null <- c()
    }
    
    
    #Test mean
    if ("Mean" %in% test_model){
      nbeta <- length(betasN)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xMP <- xM
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outM <- optim(xMP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="M",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xMP <- outM$par
        if(outM$convergence!=0)
          cat("\nStep 1 Get post statistics mean model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="M", method = "BFGS", betagamma=xMP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics mean model Convergence problems for feature ",i)
        
      }
      MLE.Mean <-c(MLE.Mean, outM$value)
    }else{
      MLE.Mean <- c()
    }
    
    
    #Test zero proportion
    if ("Pzero" %in% test_model){
      nbeta <- length(betasF)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xPP <- xP
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outP <- optim(xPP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="P",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xPP <- outP$par
        if(outP$convergence!=0)
          cat("\nStep 1 Get post statistics mean model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="P", method = "BFGS", betagamma=xPP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics mean model Convergence problems for feature ",i)
        
      }
      MLE.Pzero <-c(MLE.Pzero, outP$value)
    }else{
      MLE.Pzero <- c()
    }
  }
  n.Params = ncol(cov.matrix)
  for (k in 1:n.Params){
    colnames(Para.Full)[k] <- paste("beta.", k, sep="")
  }
  for (k in (n.Params+1):(2*n.Params)){
    colnames(Para.Full)[k] <- paste("gamma.", k-n.Params, sep="")
  }
  colnames(Para.Full)[2*n.Params+1] <- "sd"
  
  if ("Both" %in% test_model){
    Pvalue <- c(Pvalue, 1-pchisq(2*(MLE.Null-MLE.Full), 2))
  }else{
    Pvalue <- c()
  }#Null vs Full
  if ("Mean" %in% test_model){
    Pmean <- c(Pmean, 1-pchisq(2*(MLE.Mean-MLE.Full), 1))
  }else{
    Pmean <- c()
  }#Mean vs Full
  if ("Pzero" %in% test_model){
    Pzero <- c(Pzero, 1-pchisq(2*(MLE.Pzero-MLE.Full), 1))
  }else{
    Pzero <- c()
  }#Pzero vs Full
  
  result <- list(feature_names =feature.names,
                 pvalue_both=Pvalue,
                 pvalue_mean=Pmean,
                 pvalue_zero=Pzero,
                 DL=DL,
                 estimates=Para.Full)
  
  return(result)
  
}

test.data <- amines_pp_tl1a %>% select(one_of(c("sample_id", "log_pkh_imputed", "nameFix"))) %>%
            pivot_wider(names_from = nameFix, values_from = log_pkh_imputed) %>%
            column_to_rownames("sample_id") %>% t()

indata <- test.data

test.covmatrix.df <- amines_pp_tl1a %>% 
                    select(one_of(c("sample_id", "log_pkh_imputed", "nameFix", "OTU_type", "IBD_stat", "mouse_sex"))) %>% 
                    pivot_wider(names_from = nameFix, values_from = log_pkh_imputed) 
                  
test.covmatrix <- model.matrix( ~ OTU_type + IBD_stat + OTU_type*IBD_stat + mouse_sex, data = test.covmatrix.df)

cov.matrix <-test.covmatrix

result <- as.data.frame(result)

dasev.res <- DASEV(indata, cov.matrix, test_cov = 2) %>% as.data.frame()

##==================SDAMS======================##
library(SDAMS)

factors <- c("OTU_type", "IBD_stat", "mouse_sex", "interactIBD_otu2")

sdams.dplyr <- function(df, factors) {
  d <- df
  
  input <- d %>% select(one_of(c("sample_id", "log_pkh_imputed", "feat_id"))) %>%
    pivot_wider(names_from = feat_id, values_from = log_pkh_imputed) %>%
    column_to_rownames("sample_id") %>% t()
  
  metadata <- d %>% select(one_of(append(factors, "sample_id"))) %>%
    distinct() %>%
    column_to_rownames("sample_id")
  
  mm <- model.matrix(as.formula(paste0("~", paste(factors, collapse = "+"))), data=metadata) %>% 
        .[,2:length(colnames(.))]
  
  sumexp.input <- createSEFromMatrix(feature = input, colData = mm)
  
  fit_data <- SDA(sumexp.input)
  
}



###===================MAPPING TO HMDB====================###
library(S4Vectors)
library(xml2)
getdis = function(file="03-data/published_data/feces_metabolites_hmdb.xml") {
  rr = read_xml(file)
  rrl = as_list(rr)
  acc = unlist(sapply(rrl, "[[", "accession"))
  nm = unlist(sapply(rrl, "[[", "name"))
  dl = sapply(rrl, "[[", "diseases")
  diss = sapply(dl, function(x) try(unlist(sapply(x, "[[", "name"))))
  omim = sapply(dl, function(x) try(unlist(sapply(x, "[[", "omim_id"))))
  diss = sapply(diss, function(x) {
    if (inherits(x, "try-error")) x = NA
    x 
  })
  omim = sapply(omim, function(x) {
    if (inherits(x, "try-error")) x = NA
    x 
  })
  ns = unlist(sapply(diss, length))
  ons = unlist(sapply(omim, length))
  accs = rep(acc, ns)
  nms = rep(nm, ns)
  oaccs = rep(acc, ons)
  disframe = DataFrame(accession=accs, name=nms, disease=unlist(diss))
  omimframe = DataFrame(accession=oaccs, omim=unlist(omim))
  plist = function (pp, prop="uniprot_id") {
    #g = sapply(pp, "[[", "gene_name")
    ps = try(unlist(sapply(pp[[2]], "[[", prop)))
    if (inherits(ps, "try-error")) {
      df = DataFrame(accession=pp[[1]], val=NA)
      names(df)[2] = prop
      return(df)
    }
    n = length(ps)
    df = DataFrame(accession = rep(pp[[1]], n), val=ps)
    names(df)[2] = prop
    df
  }
  ppl = lapply(rrl, function(x) list(x$accession[[1]], x$protein))
  aa = lapply(ppl, plist)
  gg = lapply(ppl, plist, prop="gene_name")
  list(disframe=disframe, omimframe=omimframe,
       prots=do.call(rbind, aa), genes=do.call(rbind, gg))
  #
}
mm = getdis()

library(XML)

xml_process <- function(file) {      
  tryCatch({
    # PARSE XML TO DATA FRAME
    doc <- xmlParse(file)        
    df <- xmlToDataFrame(doc)
    
    unlink(file)# DESTROY TEMP XML
    
    # RETURN XML DF
    return(df)
  }, error = function(e) NA)      
}

df <- xml_process("03-data/published_data/hmdb_metabolites.xml") %>% 
      select(one_of(c("name", "synonyms", "inchikey", "accession", "kegg_id", "biocyc_id")))

write.csv(df, file = "03-data/published-data/hmdb_metabolites_id.csv")
