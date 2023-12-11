###===========INITIALIZATION===========###

##Loading required libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(plotly)
library(forcats)
library(stats)
library(vegan)
library(lme4)
library(nlme)

##Setting working directory
setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metabolomics")

##Setting the sample metadata

data.inflscore <- read_excel("03-data/raw_data/metadata_metab.xlsx", sheet = 'inflamation', na = c("NA", "#N/A"))

data.inflscore.imputed <-data.inflscore %>% group_by(donor_status) %>% summarize(infl_colon_med = median(Colon, na.rm =TRUE),
                                                                                 infl_duo_med = median(Duodenum, na.rm =TRUE),
                                                                                 infl_jej_med = median(Jejunum, na.rm =TRUE),
                                                                                 infl_ile_med = median(Ileum, na.rm =TRUE)) %>%
  right_join(data.inflscore, by = "donor_status") %>% mutate(Colon = ifelse(is.na(Colon), infl_colon_med, Colon),
                                                             Duodenum = ifelse(is.na(Duodenum), infl_duo_med, Duodenum),
                                                             Jejunum = ifelse(is.na(Jejunum), infl_jej_med, Jejunum),
                                                             Ileum = ifelse(is.na(Ileum), infl_ile_med, Ileum)) %>%
  rename_at(vars(c("Colon", "Duodenum", "Jejunum", "Ileum")), ~ c("infl_col", "infl_duo", "infl_jej", "infl_ile")) %>% 
  select(-one_of(c("Cecum", "infl_colon_med", "infl_duo_med", "infl_jej_med", "infl_ile_med")))

data.sampmeta <- read_excel("03-data/raw_data/metadata_metab.xlsx", sheet = 'clean', na = c("NA", "#N/A")) %>%
  mutate(d_colonize = factor(d_colonize), d_collect = factor(d_collect), donor_type = paste0(IBD_stat, "_", OTU_type)) %>%
  left_join(data.inflscore.imputed, by = c("mouse_ID" = "mouse_id"))


###===========DATA UPLOAD and WRANGLING===========###

data_upload_preprocessing <- function(
    folder_dir, sheet,
    assay_type,
    blanking_presence_cutoff = 0.8,
    blanking_log_difference = 1,
    missing_cutoff = 0.5){
  
  data <- read_excel(folder_dir, sheet = sheet, na = c("na", "", "NA", "#N/A")) %>% 
    mutate(feat_id = paste(assay_type, esi_mode, mz, r.time, sep = "_")) %>%
  filter(!is.na(name)) #filter the non identified features #create unique id for each feature; required because same feature might come up in both esi mode
  
  #separate data and feature metadata
  data.raw <- data %>% select(matches(c("feat_id", "BSD", "TL1A", "IL", "blank")))
  data.featmeta <- data %>% select(-matches(c("BSD", "TL1A", "IL", "blank", "pool")))
  
  if (assay_type == "PriMeta"){
    name_start <- c("BSD", "TL1A", "IL") #gcms data don't have any blanks, will give error if we force to select blank
  }
  
  else{
    name_start <- c("BSD", "TL1A", "IL", "blank")
  }
  
  samp_names <- colnames(data.raw[,grepl(paste(name_start, collapse="|"), names(data.raw))]) #get a list of the sample names
  
  data.long <- data.raw %>% 
    select(one_of(c('feat_id', samp_names))) %>%
    pivot_longer(cols = all_of(samp_names), 
                 names_to='sample_id', 
                 values_to='pkh') %>%
    mutate(assay = assay_type,
           gene_background = ifelse(grepl('BSD', sample_id),'Human',
                                    ifelse(grepl('TL1A', sample_id),'TL1A_mice', 
                                           ifelse(grepl('IL', sample_id),'IL10_mice', 
                                                  'Blank'))))
  
  
  ##BLANK CORRECTION
  
  #Separate blank and non-blank
  data.sample <- data.long %>% filter(gene_background != "Blank") %>% mutate(log_pkh = log10(pkh+1))
  data.blank <- data.long %>% filter(gene_background == "Blank") %>% mutate(log_pkh = log10(pkh+1))
  
  d.blank.avg <- data.blank %>% group_by(feat_id) %>% 
    summarize(log_blank_avg = mean(log_pkh, na.rm = T), percent_presence = 1 - sum(is.na(log_pkh))/n()) %>%
    mutate(filt_blank_avg = ifelse(percent_presence >= blanking_presence_cutoff, log_blank_avg, 0))  #summarize the average log of blank and the percent presence for each feature, then assume the features that are missing in >80% of the blanks as 0 
  
  d.corrected <- data.sample %>% left_join(d.blank.avg, by = "feat_id") %>%
    mutate(blank_diff = log_pkh - filt_blank_avg) %>%
    mutate(log_pkh_corrected = ifelse(blank_diff >= blanking_log_difference, log_pkh, NA)) %>%
    select(one_of(c(names(data.sample), "log_pkh_corrected", "filt_blank_avg"))) #if difference between data and blank < 1 log, then we assume the data as NA
  
  if (assay_type == "PriMeta"){
    d.corrected <- data.sample %>%
      mutate(log_pkh_corrected = log_pkh)
  } #because gcms don't have blanks
  
  ##DATA FILTERING
  
  d.filtered <<- d.corrected %>%
    group_by(feat_id) %>%
    summarize(percent_missing = sum(is.na(log_pkh_corrected))/ (sum(is.na(log_pkh_corrected)) + sum(!is.na(log_pkh_corrected)))) %>%
    right_join(d.corrected, by = "feat_id") %>%
    filter(percent_missing <= missing_cutoff) #filter by missingness
  
  ##IMPUTATION
  
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
          } else if (length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) ){
            KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]} 
          
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
  
  d.knninput <- d.filtered %>%
    select(matches(c("feat_id", "log_pkh_corrected", "sample_id"))) %>%
    pivot_wider(names_from = "feat_id", values_from = "log_pkh_corrected") %>%
    column_to_rownames("sample_id") %>%
    select(matches(assay_type))

  d.knn.res <- impute.knn.obs.sel(d.knninput, K = 10) %>% as.data.frame()%>% rownames_to_column(var="sample_id") %>%
    pivot_longer(-sample_id, names_to = 'feat_id', values_to = "log_pkh_imputed") %>%
    right_join(d.filtered, by = c("feat_id","sample_id"))
  
  d.out <- d.knn.res 
  
  return(list(data.long = d.out, data.featmeta = data.featmeta))
}

lipid <- data_upload_preprocessing (
  "03-data/raw_data/mx 738840_Jacob_Lipidomics_humanstool_mousefeces_Submit_07-2023.xlsx", "clean_1",
  "Lipidomics")

lipid.df.tl1a <- lipid$data.long %>% filter(gene_background == "TL1A_mice")
write.csv(lipid.df.tl1a, "lipid.df.tl1a.csv")

amines <- data_upload_preprocessing (
  "03-data/raw_data/mx 739116_Jacob_HILIC_biogenicamines_humanstool_09-2023 submit.xlsx", "clean_1",
  "Amines")

amines.df.tl1a <- amines$data.long %>% filter(gene_background == "TL1A_mice")
write.csv(amines.df.tl1a, "amines.df.tl1a.csv")

primeta <- data_upload_preprocessing (
  "03-data/raw_data/mx 738540_Jacob_GCTOF_humanstool_mousestool_08-2023 submit.xlsx", "clean_1",
  "PriMeta")

primeta.df.tl1a <- primeta$data.long %>% filter(gene_background == "TL1A_mice")
write.csv(primeta.df.tl1a, "primeta.df.tl1a.csv")

##======VSN Normalization======## NOT IMPLEMENTED YET
# library(vsn)
# 
# vsnnorm <- function(df, assay_type){
#   d <- df %>%
#     select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
#     pivot_wider(names_from = "sample_id", values_from = "log_pkh_imputed") %>% ## the pivot wider needs to be inside the function; if not, the pivot will force filtered out features as NA into the dataframe
#     column_to_rownames("feat_id") %>%
#     as.matrix() %>%
#     justvsn() %>%
#     t() %>% as.data.frame() %>%
#     rownames_to_column(var="sample_id") %>%
#     pivot_longer(cols = starts_with(assay_type), names_to = 'UID', values_to = "log_pkh_norm")
#     
#   return(d)
# }
# 
# lipid.norm <- vsnnorm(lipid.df, "Lipidomics")
# 
# d <- lipid.df%>%
#   select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
#   pivot_wider(names_from = "sample_id", values_from = "log_pkh_imputed") %>% ## the pivot wider needs to be inside the function; if not, the pivot will force filtered out features as NA into the dataframe
#   column_to_rownames("feat_id") %>%
#   as.matrix() %>%
#   justvsn() %>%
#   t() %>% as.data.frame()

##=======PCOA PLOT======##
pcoa_result <- function(df, meta, plot_title) {
  
  d <- df
  meta1 <- meta %>% right_join(distinct(select(d, one_of("sample_id"))), by = 'sample_id')
  
  plot_title <<- plot_title
  
  d.sh <- d %>%
        select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
        pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>%
        column_to_rownames("sample_id")
  
  d.dist <- d.sh %>% decostand(method = "normalize") %>% vegdist(method = "euclidean", na.rm = TRUE)
  
  d.adonis.IBD <- adonis2(d.dist ~ IBD_stat + mouse_sex,
                          method = "euclidean",
                          na.rm = TRUE,
                          permutations = 9999,
                          data = meta1)
  d.adonis.OTU <- adonis2(d.dist ~ OTU_type + mouse_sex,
                          method = "euclidean",
                          na.rm = TRUE,
                          permutations = 9999,
                          data = meta1)
  d.adonis.infl <- adonis2(d.dist ~ infl_score + mouse_sex,
                           method = "euclidean",
                           na.rm = TRUE,
                           permutations = 9999,
                           data = meta1)
  
  d.pcoa <- wcmdscale(d.dist, eig = TRUE)
  
  d.pcoa.points <- d.pcoa$points %>% as.data.frame() %>% select(1,2) %>%
    rownames_to_column(var = "sample_id") %>% 
    left_join(meta1, by = "sample_id")
  
  d.pcoa.eig <- d.pcoa$eig %>% as.data.frame %>% mutate(var = ./sum(.)*100)
  
  p1 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = OTU_type)) +
    geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
              y = .95 * max(d.pcoa.points$Dim2),
              hjust = 0,
              label = paste0("P.val = ", signif(d.adonis.OTU$"Pr(>F)"[1], digits = 2),
                             "\n","R² = ", signif(d.adonis.OTU$R2[1], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean() +
    labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
         y = paste0("PC2"," (",signif(d.pcoa.eig$var[2], digits = 2),"%)")) +
    scale_color_discrete(name = "",
                         labels = c("OTU type 1",
                                    "OTU type 2"))
  
  p2 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = IBD_stat)) +
    geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
              y = .95 * max(d.pcoa.points$Dim2),
              hjust = 0,
              label = paste0("P.val = ", signif(d.adonis.IBD$"Pr(>F)"[1], digits = 2),
                             "\n","R² = ", signif(d.adonis.IBD$R2[1], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean() +
    labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
         y = element_blank()) +
    scale_color_manual(values = c("Healthy" = "#b37d8b", "IBD" = "#b3d4ff"),
                       name = "",
                       labels = c("Healthy",
                                  "IBD"))
  
  p3 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = infl_col)) +
    geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
              y = .95 * max(d.pcoa.points$Dim2),
              hjust = 0,
              label = paste0("P.val = ", signif(d.adonis.infl$"Pr(>F)"[1], digits = 2),
                             "\n","R² = ", signif(d.adonis.infl$R2[1], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean() +
    scale_color_continuous(name = "Inflammation Score", type = "viridis") +
    labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
         y = element_blank())
  
  pall <- p1 + p2 + plot_annotation(plot_title, theme=theme(plot.title=element_text(hjust=0.5)))
  
  return(list(plot = pall, df = d.sh, dfpcoa = d.pcoa.points, dfeig = d.pcoa.eig))
}

lipid.df.tl1a <- read.csv("withunknown/lipid.df.tl1a.csv")
amines.df.tl1a <- read.csv("withunknown/amines.df.tl1a.csv")
primeta.df.tl1a <- read.csv("withunknown/primeta.df.tl1a.csv")

# lipid.df.tl1a <- read.csv("withoutunknown/lipid.df.tl1a.csv")
# amines.df.tl1a <- read.csv("withoutunknown/amines.df.tl1a.csv")
# primeta.df.tl1a <- read.csv("withoutunknown/primeta.df.tl1a.csv")

lipid.pcoa <- pcoa_result(lipid.df.tl1a, data.sampmeta, "Lipidomics")
lipid.pcoa$plot
ggsave("metabolomics_pcoa_tl1a_lipid_norm_infl.pdf", height = 5, width = 10)

amines.pcoa <- pcoa_result(amines.df.tl1a, data.sampmeta, "Biogenic Amines")
amines.pcoa$plot
ggsave("metabolomics_pcoa_tl1a_amines_norm_infl.pdf", height = 5, width = 10)

primeta.pcoa <- pcoa_result(primeta.df.tl1a, data.sampmeta, "Primary Metabolite")
primeta.pcoa$plot
ggsave("metabolomics_pcoa_tl1a_primeta_norm_infl.pdf", height = 5, width = 10)

#COMBINING ALL VALUES
d.amines.sh <- amines.df.tl1a %>%
  select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
  pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>%
  column_to_rownames("sample_id") %>% decostand(method = "normalize")%>%
  rownames_to_column(var = "sample_id")

d.lipid.sh <- lipid.df.tl1a %>%
  select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
  pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>%
  column_to_rownames("sample_id") %>% decostand(method = "normalize")%>%
  rownames_to_column(var = "sample_id")

d.primeta.sh <- primeta.df.tl1a %>%
  select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
  pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>%
  column_to_rownames("sample_id") %>% decostand(method = "normalize") %>%
  rownames_to_column(var = "sample_id")

d.sh <- d.amines.sh %>% left_join(d.lipid.sh, by = "sample_id") %>% left_join(d.primeta.sh, by = "sample_id") %>%
  column_to_rownames("sample_id") %>% as.matrix()
  

combined.df.tl1a <- lipid.df.tl1a %>% bind_rows(amines.df.tl1a) %>% bind_rows(primeta.df.tl1a)
combined.pcoa <- pcoa_result(combined.df.tl1a, data.sampmeta, "Metabolomics")
combined.pcoa$plot
ggsave("metabolomics_pcoa_tl1a_combined_norm_infl.pdf", height = 5, width = 10)

#=====MAASLIN2======##
library(Maaslin2)

metab_maaslin <- function(factors, df, metadata, feat_meta, sample_subset_string_name){
  
  d <- df
  meta <- metadata
  f.meta <- feat_meta
  
  d.sh <- d %>%
    select(one_of(c("feat_id", "log_pkh_imputed", "sample_id"))) %>%
    pivot_wider(names_from = "feat_id", values_from = "log_pkh_imputed") %>% 
    column_to_rownames("sample_id")
  
  m <- meta %>% select(one_of(c(factors, "sample_id"))) %>% column_to_rownames("sample_id")
  
  fit_data <- Maaslin2(input_data=d.sh, 
                       input_metadata=m, 
                       output = paste0("04-analysis/01_maaslin2/metabolomics_", sample_subset_string_name), 
                       fixed_effects = factors,
                       reference = c("mouse_sex,F", "OTU_type,otu_1", "IBD_stat,Healthy"),
                       normalization="NONE", 
                       transform ="NONE",
                       )
  
  fit_data.df <- fit_data[["results"]] %>% as.data.frame() %>%
    left_join(f.meta, by = c("feature" = "feat_id"))
  
  fit_data.df.colormap <- fit_data.df %>% filter(metadata != "mouse_sex" & qval <= .05) %>% 
    group_by(nameFix, metadata) %>% summarize(meancoef = mean(coef)) %>%
    pivot_wider(names_from = metadata, values_from = meancoef) %>%
    mutate(same_dir_signif = ifelse(is.na(IBD_stat)|is.na(OTU_type), "ivory4", 
                                    ifelse(IBD_stat * OTU_type > 0, "cornflowerblue", "coral1"))
           
    )
  
  otu.colmap <- fit_data.df.colormap %>% filter(!is.na(OTU_type)) %>% 
    arrange(OTU_type)
  IBD.colmap <- fit_data.df.colormap %>% filter(!is.na(IBD_stat)) %>% 
    arrange(IBD_stat)
  
  p1 <- fit_data.df %>% filter(metadata == "OTU_type" & qval <= .05) %>%
    ggplot(aes(x = coef, y = reorder(nameFix, coef))) +
    geom_point(alpha = .65) +
    theme_classic() +
    labs(x = "Log2FC OTU type 1 vs OTU type 2",
         y = "") +
    geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
    theme(
      axis.text.y = element_text(size = 8, color = otu.colmap$same_dir_signif)
    )
  
  p2 <- fit_data.df %>% filter(metadata == "IBD_stat" & qval <= .05) %>%
    ggplot(aes(x = coef, y = reorder(nameFix, coef))) +
    geom_point(alpha = .65) +
    theme_classic() +
    labs(x = "Log2FC Healthy vs IBD",
         y = "") +
    geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
    theme(
      axis.text.y = element_text(size = 8, color = IBD.colmap$same_dir_signif)
    )
  
  pall <- p1 + p2 + 
    plot_layout(guides = "collect") + 
    plot_annotation(sample_subset_string_name, theme=theme(plot.title=element_text(hjust=0.5)))
  
  return(list(fit_df = fit_data.df, plot = pall))
  
}

filterNonNames <- function(df, rawfeatmeta){
  d <- df
  fm <- rawfeatmeta
  
  d.filt <- d %>% left_join(fm, by = "feat_id") %>% filter(!is.na(name))
  
  return(d.filt)
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

factors <- c("OTU_type", "IBD_stat", "mouse_sex")

lipid.df.metafeat <- fixDuplicatedNames(lipid$data.featmeta)
lipid.df.tl1a.Named <- filterNonNames(lipid.df.tl1a, lipid$data.featmeta)
maaslin2_lipid <- metab_maaslin(factors, lipid.df.tl1a.Named, data.sampmeta, lipid.df.metafeat, "Lipidomics")
maaslin2_lipid$plot
ggsave("metabolomics_maaslin_tl1a_lipidomics.pdf", height = 11, width = 14)

amines.df.metafeat <- fixDuplicatedNames(amines$data.featmeta)
amines.df.tl1a.Named <- filterNonNames(amines.df.tl1a, amines$data.featmeta)
maaslin2_amines <- metab_maaslin(factors, amines.df.tl1a.Named, data.sampmeta, amines.df.metafeat, "Biogenic Amines")
maaslin2_amines$plot
ggsave("metabolomics_maaslin_tl1a_biogenicamines.pdf", height = 11, width = 14, limitsize = FALSE)

primeta.df.metafeat <- fixDuplicatedNames(primeta$data.featmeta)
primeta.df.tl1a.Named <- filterNonNames(primeta.df.tl1a, primeta$data.featmeta)
maaslin2_primeta <- metab_maaslin(factors, primeta.df.tl1a.Named, data.sampmeta, primeta.df.metafeat, "Primary Metabolite")
maaslin2_primeta$plot
ggsave("metabolomics_maaslin_tl1a_primeta.pdf", height = 11, width = 14)


#========= reading the json table of HMDB ===========#
library(fuzzyjoin)
hmdb <-  rjson::fromJSON(file = "03-data/databases/hmdb_full.json")

lookup_hmdb <- function(compound) {
  # Check if compound matches any 'name' or 'synonyms' in HMDB entries
  
  compound <<- compound
  
  if (nchar(compound) > 55){
    compound <- substr(compound, nchar(compound)-55, nchar(compound))
  }
  
  matching_entries <- lapply(hmdb, function(entry) {
    name_match <- grepl(compound, entry$name, ignore.case = TRUE)
    synonyms_match <- any(grepl(compound, entry$synonyms, ignore.case = TRUE))
    name_match | synonyms_match
  })
  
  # Extract the names of the matching entries
  matching_entry_names <- names(hmdb)[unlist(matching_entries)]
  
  # Return a list containing entry names and taxonomy information
  
  taxonomy<- lapply(matching_entry_names, function(entry_name) {
    list(input_name = compound,
         hmdb_id = entry_name, 
         name = c(hmdb[[entry_name]]$name, hmdb[[entry_name]]$synonyms),
         tax_kingdom = hmdb[[entry_name]]$tax_kingdom,
         tax_superclass = hmdb[[entry_name]]$tax_superclass,
         tax_class = hmdb[[entry_name]]$tax_class,
         tax_subclass = hmdb[[entry_name]]$tax_subclass,
         tax_direct_parent = hmdb[[entry_name]]$tax_direct_parent)
  })
  
  taxonomy <<- taxonomy
  
  if (length(taxonomy) == 0){
    taxonomy[[1]] <- list(input_name = compound,
                          hmdb_id = NA, 
                          name = NA,
                          tax_kingdom = NA,
                          tax_superclass = NA,
                          tax_class = NA,
                          tax_subclass = NA,
                          tax_direct_parent = NA)
    
    output <- data.table::rbindlist(taxonomy, fill = TRUE) %>% 
      mutate(string_dist = adist(name, compound))
  }
  
  else{
    output <- data.table::rbindlist(taxonomy, fill = TRUE) %>% 
      mutate(string_dist = adist(name, compound)) %>%
      filter(string_dist == min(string_dist)) %>%
      distinct()
  }
  
  return(output)
}

# test <- lookup_hmdb("(2-Methylpropyl)-3,6-dioxopiperazin-2-yl]propanoic acid")

# Apply the lookup function to each chemical compound
chem_compound.lipid <- maaslin2_lipid$fit_df %>% filter(qval < 0.05) %>% select(name.y) %>% distinct() %>% pull(name.y)
chem_compound.amines <- maaslin2_amines$fit_df %>% filter(qval < 0.05) %>% select(name.y) %>% distinct()%>% pull(name.y)
chem_compound.primeta <- maaslin2_primeta$fit_df %>% filter(qval < 0.05) %>% select(name.y) %>% distinct()%>% pull(name.y)

result.lipid <- lapply(chem_compound.lipid, lookup_hmdb)
result.amines <- lapply(chem_compound.amines, lookup_hmdb)
result.primeta <- lapply(chem_compound.primeta, lookup_hmdb)

##creating figure for inflammation scoring
library(rstatix)
colmap1=c("Healthy OTU 1"="#B0E0E6","Healthy OTU 2"="#0000FF","IBD OTU 1"="#F0CCB0","IBD OTU 2"="#D2691E")

inflamation <- read_excel("03-data/raw_data/metadata_metab.xlsx", sheet = 'inflamation', na = c("NA", "#N/A")) %>%
  pivot_longer(cols = c("Colon", "Duodenum", "Jejunum", "Ileum", "Cecum"), names_to = "location", values_to = "infl_score") %>%
  mutate(location = factor(location, levels = c("Duodenum", "Jejunum", "Ileum","Cecum", "Colon"))) %>%
  filter(location != "Cecum")

inflamation.stat <- inflamation %>% group_by(location) %>% dunn_test(infl_score ~ donor_status, p.adjust.method = "BH") %>% 
  filter(p.adj <= 0.05)

ggplot(inflamation, aes(x = donor_status, y = infl_score)) +
  geom_violin(aes(color = donor_status), position = position_dodge(0.9)) + 
  geom_boxplot(position = position_dodge(0.9), width = .1, aes(color = donor_status)) +
  scale_color_manual(values=colmap1, name = "Donor Status\n", 
                     labels = c("Healthy - OTU type 1",
                                "Healthy - OTU type 2",
                                "IBD - OTU type 1",
                                "IBD - OTU type 2")) +
  ylim(c(0, 1.2*max(inflamation$infl_score))) +
  facet_grid(.~location , labeller = labeller(location = c(Cecum = "Cecum", Duodenum = "Duodenum", Jejunum = "Jejunum", Ileum = "Ileum", Colon = "Colon"))) +
  stat_compare_means(method = "kruskal.test", label.y = 10, size = 3, aes(group = donor_status)) +
  labs(y = "Histologic Inflammation Score") +
  theme_pubclean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(colour = guide_legend(nrow = 2)) +
  stat_pvalue_manual(inflamation.stat, label = "p.adj.signif", y.position = c(8, 8.5, 9, 9.5), size = 2.5,
                     tip.length = 0.01)
  
  # stat_compare_means(comparisons = list(c("Healthy OTU 1", "Healthy OTU 2"), c("Healthy OTU 1", "IBD OTU 1"),c("Healthy OTU 1", "IBD OTU 2")),
  #                    method = "wilcox.test", label = "p.signif", size = 2.5, tip.length = 0)

# ggsave("inflammation_score_violinplot_nolines.pdf", width = 4.5, height = 9.5)

fibrosis <- data.sampmeta %>% 
  mutate(donor_id = paste0(IBD_stat, "_", OTU_type)) %>%
  mutate(donor_id = factor(donor_id, levels=c("Healthy_otu_1", "Healthy_otu_2", "IBD_otu_1", "IBD_otu_2"))) %>%
  select(one_of(c("donor_id", "fibrosis_duo", "fibrosis_jej", "fibrosis_ile", "fibrosis_colon"))) %>%
  pivot_longer(-donor_id, names_to = "location", values_to = "fibrosis_score") %>%
  mutate(location = factor(location, levels = c("fibrosis_duo", "fibrosis_jej", "fibrosis_ile","fibrosis_colon")))

fibrosis.stat <- fibrosis %>% group_by(location) %>% dunn_test(fibrosis_score ~ donor_id, p.adjust.method = "BH") %>% 
  filter(p.adj <= 0.05)

loc_name <- c('fibrosis_colon' = "Colon",
              'fibrosis_duo' = "Duodenum",
              'fibrosis_jej' = "Jejunum",
              'fibrosis_ile' = "Ileum")

colmap2=c("Healthy_otu_1"="#B0E0E6","Healthy_otu_2"="#0000FF","IBD_otu_1"="#F0CCB0","IBD_otu_2"="#D2691E")

ggplot(fibrosis, aes(x = donor_id, y = fibrosis_score)) +
  geom_violin(aes(color = donor_id), position = position_dodge(0.9)) + 
  geom_boxplot(position = position_dodge(0.9), width = .1, aes(color = donor_id)) +
  scale_color_manual(values=colmap2, name = "Donor Status\n", 
                     labels = c("Healthy - OTU type 1",
                                "Healthy - OTU type 2",
                                "IBD - OTU type 1",
                                "IBD - OTU type 2")) +
  ylim(c(0, 1.25*max(fibrosis$fibrosis_score, na.rm = TRUE))) +
  facet_grid(.~location , labeller = labeller(location = c(Cecum = "Cecum", fibrosis_duo = "Duodenum", 
                                                           fibrosis_jej = "Jejunum", fibrosis_ile = "Ileum", fibrosis_colon = "Colon"))) +
  stat_compare_means(method = "kruskal.test", label.y = 1.19*max(fibrosis$fibrosis_score, na.rm = TRUE), size = 3, aes(group = donor_id)) +
  labs(y = "Histologic Fibrosis Score") +
  theme_pubclean()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(colour = guide_legend(nrow = 2)) +
  stat_pvalue_manual(fibrosis.stat, label = "p.adj.signif", 
                     y.position = c(1.02*max(fibrosis$fibrosis_score, na.rm = TRUE),
                                    1.08*max(fibrosis$fibrosis_score, na.rm = TRUE), 
                                    1.14*max(fibrosis$fibrosis_score, na.rm = TRUE)), 
                     size = 2.5,
                     tip.length = 0.01)
  