library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(vegan)
library(plotly)
setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/metagenomics")

data.inflscore <- read_excel("03-data/raw_data/metadata_metag.xlsx", sheet = 'inflamation', na = c("NA", "#N/A"))

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


dat <- read_tsv("03-data/raw_data/ASV_count_table_TL1A_gnotobiotic.txt")
metaraw <- read_xlsx("03-data/raw_data/metadata_metag.xlsx", sheet = "clean") %>%
  mutate(samp_subset = ifelse(grepl("-DTT", sample_id, fixed = TRUE), "mucosal", "lumen")) %>%
  left_join(data.inflscore.imputed, by = "mouse_id")

##FIXING WEIRD MAPPINGS
samp_names_mucosal <- colnames(dat[,grepl("-DTT", names(dat))])
samp_names_lumen <- colnames(dat[,grepl("-FP", names(dat))])
samp_names_others <- colnames(dat[,!grepl("-FP", names(dat)) & !grepl("-DTT", names(dat))]) %>% str_sort() %>% as.data.frame() %>%
  rename(., unmapped_names = .) %>%
  mutate(unmapped_names2 = unmapped_names)

mucosal_prefix <- samp_names_mucosal %>%
  str_sort() %>%
  as.data.frame() %>%
  separate_wider_delim( data = ., cols = ., delim = "-", names = c("prefix", 'muc'), 
                        too_many = "debug", cols_remove = FALSE)

lumen_prefix <- samp_names_lumen %>%
  str_sort() %>%
  as.data.frame() %>%
  separate_wider_delim( data = ., cols = ., delim = "-", names = c("prefix", 'lum'), 
                        too_many = "debug", cols_remove = FALSE)

compare_lumen_mucosal_mapping <- mucosal_prefix %>% full_join(lumen_prefix, by = "prefix") %>% select(one_of(c("prefix", "..x", "..y"))) %>%
  separate_wider_delim(data = ., cols = prefix, cols_remove = FALSE, delim = "TL1A", names = c("number", "num_prefix")) %>%
  select(-number) %>% full_join(samp_names_others, by = c("num_prefix" = "unmapped_names"))

# compare_lumen_mucosal_mapping_knownmap <- compare_lumen_mucosal_mapping %>% filter(!is.na(..x) & !is.na(..y))

compare_lumen_mucosal_mapping_unknownnumber <- compare_lumen_mucosal_mapping %>% 
  filter(is.na(..x) | is.na(..y)) %>% filter((!is.na(..x) | !is.na(..y)) & !is.na(unmapped_names2)) %>%
  mutate(probable_mapping = ifelse(is.na(..y), paste0("TL1A",unmapped_names2,"-FP"), paste0("TL1A",unmapped_names2,"-DTT"))) %>% 
  select(one_of("unmapped_names2", "probable_mapping"))

compare_lumen_mucosal_mapping_weird <- compare_lumen_mucosal_mapping %>% filter(is.na(..x) | is.na(..y)) %>% filter(is.na(..x) & is.na(..y)) %>%
  filter(grepl("TL1a", unmapped_names2)) %>%
  separate_wider_delim(data = ., cols = unmapped_names2, cols_remove = FALSE, delim = "_", names = c("numberee", "watt", "hello")) %>%
  mutate(probable_mapping = str_replace(toupper(numberee), "(TL1A)([A-Z]+)([0-9]+)", "\\1\\3-\\2"),
         platform = "novaseq") %>% 
  select(one_of("unmapped_names2", "probable_mapping", "platform"))
           
corrected_mapping <- bind_rows(compare_lumen_mucosal_mapping_unknownnumber, compare_lumen_mucosal_mapping_weird)

###

meta <- metaraw %>% left_join(corrected_mapping, by = c("sample_id" = "probable_mapping")) %>% 
  mutate(seq_id = ifelse(is.na(unmapped_names2), sample_id, unmapped_names2),
         platform = ifelse(is.na(platform), "miseq", platform)) %>%
  select(one_of("seq_id", 'sample_id', "platform", "OTU_type", "IBD_stat", "location", "mouse_sex", "infl_col"))

# samp_names_mucosalfix <- colnames(dat[,grepl("-DTT", names(dat))]) %>% c(., pull(filter(testcombi_final, grepl("-DTT", ggfix)), wat))
# samp_names_lumenfix <- colnames(dat[,grepl("-FP", names(dat))]) %>% c(., pull(filter(testcombi_final, grepl("-FP", ggfix)), wat))
# 
# df.lumen <- dat %>% 
#   select(c(samp_names_lumenfix, "taxonomy")) %>% 
#   pivot_longer(cols = samp_names_lumenfix, names_to = "sample_id", values_to = "count") %>%
#   left_join(testcombi_final, by = c("sample_id" = "wat")) %>%
#   mutate(ggfix = ifelse(is.na(ggfix), sample_id, ggfix))
#   
# df.mucosal <- dat %>% 
#   select(one_of(c(samp_names_mucosalfix, "taxonomy"))) %>%  
#   pivot_longer(cols = samp_names_mucosalfix, names_to = "sample_id", values_to = "count")%>%
#   left_join(testcombi_final, by = c("sample_id" = "wat")) %>%
#   mutate(ggfix = ifelse(is.na(ggfix), sample_id, ggfix))


pca_result <- function(rawdf, metadata, sample_names_suffix, sample_subset_string_name, sample_loc = ".all") {
  #Create PCA plot of raw data, colored according to potential confounding factors
  #PCA should be with or without blanks. Samples that are similar to blanks are removed and deemed defective.
  #Aim is to provide initial diagnostic plot
  
  d <- rawdf
  
  if (sample_loc != ".all"){
    seq_names <- meta %>% filter(grepl(sample_names_suffix, sample_id) & grepl(sample_loc, location)) %>% pull(seq_id)
  }
  else{
    seq_names <- meta %>% filter(grepl(sample_names_suffix, sample_id)) %>% pull(seq_id)
  }
  
  d.filt <- d %>% select(-1) %>% select(one_of("taxonomy", seq_names))

  d.sh <- d.filt %>% group_by(taxonomy) %>%
    summarize(noftax = n()) %>%
    right_join(d.filt, by = 'taxonomy') %>%
    mutate(row_num = row_number()) %>%
    mutate(tax_fixed = ifelse(noftax>1, paste0(taxonomy, "_otu", row_num), taxonomy)) %>%
    column_to_rownames("tax_fixed") %>% select(one_of(seq_names)) %>% t()
  
  totcount <- d.sh %>% rowSums() %>% as.data.frame() %>% rownames_to_column(var = "seq_id")
  colnames(totcount)[colnames(totcount) == "."] <- "countsum" 
  specnum <- specnumber(d.sh) %>% as.data.frame() %>% rownames_to_column(var = "seq_id")
  colnames(specnum)[colnames(specnum) == "."] <- "totspec" 
  qctable <- left_join(totcount, specnum, by = "seq_id") %>% mutate(logcountsum = log2(countsum+1))
  
  d.dist <<- d.sh %>% vegdist(method = "bray", na.rm = TRUE)
  
  metaadonis <<- meta %>% 
    right_join(as.data.frame(row.names(d.sh)), by = c("seq_id" = "row.names(d.sh)")) %>%
    left_join(qctable, by = "seq_id")
  
  d.adonis.IBD <- adonis2(d.dist ~ platform + mouse_sex + IBD_stat,
                          method = "bray",
                          na.rm = TRUE,
                          permutations = 9999,
                          data = metaadonis)
  d.adonis.OTU <- adonis2(d.dist ~ platform + mouse_sex + OTU_type,
                          method = "bray",
                          na.rm = TRUE,
                          permutations = 9999,
                          data = metaadonis)
  d.adonis.infl <- adonis2(d.dist ~ platform + mouse_sex + OTU_type + IBD_stat + infl_col,
                           method = "bray",
                           na.rm = TRUE,
                           permutations = 9999,
                           data = metaadonis)
  # d.adonis.loc <- adonis2(d.dist ~ location + mouse_sex,
  #                          method = "bray",
  #                          na.rm = TRUE,
  #                          permutations = 9999,
  #                          data = metaadonis)
  
  d.pcoa <- wcmdscale(d.dist, eig = TRUE)
  
  d.pcoa.points <- d.pcoa$points %>% as.data.frame() %>% select(1,2) %>%
    rownames_to_column(var = "seq_id") %>% 
    left_join(metaadonis, by = "seq_id")
  
  d.pcoa.eig <- d.pcoa$eig %>% as.data.frame %>% mutate(var = ./sum(.)*100)
  
  p1 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2, text = sample_id)) +
    geom_point(aes(color = OTU_type)) + 
    geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
              y = 1.1 * max(d.pcoa.points$Dim2),
              hjust = 0,
              label = paste0("P.val = ", signif(d.adonis.OTU["OTU_type", "Pr(>F)"], digits = 2),
                             "\n","R² = ", signif(d.adonis.OTU["OTU_type", "R2"], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean()+
    ylim(1.05 * min(d.pcoa.points$Dim2), 1.15 * max(d.pcoa.points$Dim2)) +
    labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
         y = paste0("PC2"," (",signif(d.pcoa.eig$var[2], digits = 2),"%)")) +
    scale_color_discrete(name = "",
                         labels = c("OTU type 1",
                                    "OTU type 2"))
  
  p2 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2, text = sample_id)) +
    geom_point(aes(color = IBD_stat)) +
    geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
              y = 1.1 * max(d.pcoa.points$Dim2),
              hjust = 0,
              label = paste0("P.val = ", signif(d.adonis.IBD["IBD_stat", "Pr(>F)"], digits = 2),
                             "\n","R² = ", signif(d.adonis.IBD["IBD_stat", "R2"], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean() +
    ylim(1.05 * min(d.pcoa.points$Dim2), 1.15 * max(d.pcoa.points$Dim2)) +
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
              label = paste0("P.val = ", signif(d.adonis.infl["infl_col", "Pr(>F)"], digits = 2),
                             "\n","R² = ", signif(d.adonis.infl["infl_col", "R2"], digits = 2)),
              show.legend = FALSE,
              size = 2.5) +
    theme_pubclean() +
    scale_color_distiller(name = "Inflammation Score \n", palette = "YlOrRd", direction = 1) +
    labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"),
         y = element_blank())
  # 
  # p4 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2)) +
  #   geom_point(aes(color = location)) +
  #   geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
  #             y = .95 * max(d.pcoa.points$Dim2),
  #             hjust = 0,
  #             label = paste0("P.val = ", signif(d.adonis.loc$"Pr(>F)"[1], digits = 2),
  #                            "\n","R² = ", signif(d.adonis.loc$R2[1], digits = 2)),
  #             show.legend = FALSE,
  #             size = 2.5) +
  #   theme_pubclean() +
  #   scale_color_brewer(name = "", palette = "Set2") +
  #   labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
  #        y = element_blank())
  
  p5 <- ggplot(d.pcoa.points, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = logcountsum)) +
    scale_color_distiller(name = "Log2 Sum of all OTU counts \n", palette = "Spectral", direction = 1) +
      labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"), 
           y = paste0("PC2"," (",signif(d.pcoa.eig$var[2], digits = 2),"%)")) + coord_fixed()
  
  # +
  #   geom_text(x = min(d.pcoa.points$Dim1) + 0.01*(max(d.pcoa.points$Dim1) - min(d.pcoa.points$Dim1)),
  #             y = .95 * max(d.pcoa.points$Dim2),
  #             hjust = 0,
  #             label = paste0("P.val = ", signif(d.adonis.infl["infl_col", "Pr(>F)"], digits = 2),
  #                            "\n","R² = ", signif(d.adonis.infl["infl_col", "R2"], digits = 2)),
  #             show.legend = FALSE,
  #             size = 2.5) +
  #   theme_pubclean() +
  #   scale_color_distiller(name = "Inflammation Score \n", palette = "YlOrRd", direction = 1) +
  #   labs(x = paste0("PC1"," (",signif(d.pcoa.eig$var[1], digits = 2),"%)"),
  #        y = element_blank())
  
  pall <- p1 + p2 + p3 + plot_annotation(paste0(sample_subset_string_name), theme=theme(plot.title=element_text(hjust=0.5)))
  
  return(list(plot = pall, df = d.sh, dfpcoa = d.pcoa.points, dfeig = d.pcoa.eig))
}

for (x in unique(meta$location)){
  pcoa_dtt <- pca_result(dat, meta, sample_names_suffix = "-DTT", "Mucosal", sample_loc = x)
  pcoa_dtt$plot
  ggsave(paste0("metagenomics_pcoa_mucosal_",x,".pdf"), height = 4, width = 6)
}

pcoa_dtt <- pca_result(dat, meta, sample_names_suffix = "-DTT", "Mucosal", sample_loc = "colon")
pcoa_dtt$plot
ggsave(paste0("metagenomics_pcoa_mucosal_colon_infl.pdf"), height = 5, width = 11)

pcoa_dtt <- pca_result(dat, meta, sample_names_suffix = "-FP", "Lumen", sample_loc = "colon")
pcoa_dtt$plot
ggsave(paste0("metagenomics_pcoa_lumen_colon_infl.pdf"), height = 5, width = 11)

# 
# pcoa_dtt_proxcol <- pca_result(dat, meta, "-DTT", "Mucosal", sample_loc = "proximal_colon", "Proximal Colon")
# pcoa_dtt_proxcol$plot
# ggsave("metagenomics_pcoa_mucosal_new_proxcol.pdf", height = 6, width = 12)
# 
# pcoa_dtt_discol <- pca_result(dat, meta, "-DTT", "Mucosal", sample_loc = "distal_colon", "Distal Colon")
# pcoa_dtt_discol$plot
# ggsave("metagenomics_pcoa_mucosal_new_discol.pdf", height = 6, width = 12)
# 
# pcoa_dtt_jej <- pca_result(dat, meta, "-DTT", "Mucosal", sample_loc = "Jejunum", "Jejunum")
# pcoa_dtt_jej$plot
# ggsave("metagenomics_pcoa_mucosal_new_jej.pdf", height = 6, width = 12)
# 
# pcoa_dtt_ile <- pca_result(dat, meta, "-DTT", "Mucosal", sample_loc = "Ileum", "Ileum")
# pcoa_dtt_ile$plot
# ggsave("metagenomics_pcoa_mucosal_new_ile.pdf", height = 6, width = 12)
# 
# pcoa_dtt_duo <- pca_result(dat, meta, "-DTT", "Mucosal", sample_loc = "Duodenum", "Duodenum")
# pcoa_dtt_duo$plot
# ggsave("metagenomics_pcoa_mucosal_new_duo.pdf", height = 6, width = 12)
# 
# pcoa_dtt <- pca_result(dat, meta, "-DTT", "Mucosal", location_title = "Metagenomics")
# pcoa_dtt$plot
# ggsave("metagenomics_pcoa_mucosal_new_loc.pdf", height = 6, width = 12)

pcoa_fp <- pca_result(dat, meta, "-FP", "Lumen")
pcoa_fp$plot
ggsave("metagenomics_pcoa_lumen_new_colon.pdf", height = 6, width = 12)

checkdtt <-pcoa_dtt$dfpcoa %>% group_by(OTU_type, IBD_stat) %>% summarize (n = n())
checkfp <- pcoa_fp$dfpcoa %>% group_by(OTU_type, IBD_stat) %>% summarize (n = n())

savemeta <- meta %>% select(one_of(c("wat", "mouse_id", "sample_id", "IBD_stat", "OTU_type"))) %>% distinct()
write.csv(savemeta, "new_mapping_16s_tl1a.csv")


####MAASLIN######
library(Maaslin2)

# factors <- c("OTU_type", "IBD_stat", "mouse_sex")
# 
# d <- pcoa_dtt$df
# m <- pcoa_dtt$dfpcoa %>% select(one_of(c(factors, "sample_id"))) %>% column_to_rownames("sample_id")
# m2 <- meta %>% select(one_of(c(factors, "sample_id"))) %>% column_to_rownames("sample_id")
# 
# fit_data <- Maaslin2(input_data=d, 
#                      input_metadata=m2, 
#                      output = "04-analysis/01_maaslin2/metagenomics_dtt", 
#                      fixed_effects = factors,
#                      reference = c("mouse_sex,F", "OTU_type,otu_1", "IBD_stat,Healthy"),
#                      normalization="TSS", 
#                      transform ="LOG",
#                      plot_heatmap = FALSE,
#                      plot_scatter = FALSE)
# 
# fit_data.df <- fit_data[["results"]] %>% as.data.frame() %>%
#   separate_wider_delim(feature, delim = regex("(k__)|(..p__)|(..c__)|(..o__)|(..g__)|(..f__)|(..s__)"), 
#                        names = c("rm1","tax_king", "tax_phy", "tax_class", 
#                                   "tax_ord", "tax_fam", "tax_gen", "tax_sp")) %>%
#   select(-rm1) %>%
#   mutate(otu_name = ifelse(tax_gen != "", paste0(tax_gen, "_", tax_sp), ifelse(
#                                 tax_fam != "", paste0(tax_fam, "_", tax_sp), paste0(tax_ord, "_", tax_sp))),
#          tax_genfixed = ifelse(tax_gen == "", ifelse(tax_fam == "", ifelse(tax_ord == "", tax_class, tax_ord), tax_fam), tax_gen)
#          )
# 
# fit_data.df.colormap <- fit_data.df %>% filter(metadata != "mouse_sex" & qval <= .05) %>% 
#   group_by(tax_genfixed, metadata) %>% summarize(meancoef = mean(coef)) %>%
#   pivot_wider(names_from = metadata, values_from = meancoef) %>%
#   mutate(same_dir_signif = ifelse(is.na(IBD_stat)|is.na(OTU_type), "ivory4", 
#                                   ifelse(IBD_stat * OTU_type > 0, "cornflowerblue", "coral1"))
#     
#     )
# 
# otu.colmap <- fit_data.df.colormap %>% filter(!is.na(OTU_type)) %>% 
#   arrange(OTU_type)
# IBD.colmap <- fit_data.df.colormap %>% filter(!is.na(IBD_stat)) %>% 
#   arrange(IBD_stat)
# 
# 
# p1 <- fit_data.df %>% filter(metadata == "OTU_type" & qval <= .05) %>%
#   ggplot(aes(x = coef, y = reorder(tax_genfixed, coef), color = tax_phy)) +
#   geom_point(alpha = .65) +
#   theme_classic() +
#   labs(x = "Log2FC OTU type 1 vs OTU type 2",
#        y = "") +
#   scale_color_discrete(name = "Phylum") +
#   geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
#   theme(
#     axis.text.y = element_text(size = 8, color = otu.colmap$same_dir_signif)
#   )
# 
# p2 <- fit_data.df %>% filter(metadata == "IBD_stat" & qval <= .05) %>%
#   ggplot(aes(x = coef, y = reorder(tax_genfixed, coef), color = tax_phy)) +
#   geom_point(alpha = .65) +
#   theme_classic() +
#   labs(x = "Log2FC Healthy vs IBD",
#        y = "") +
#   scale_color_discrete(name = "Phylum") +
#   geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
#   theme(
#     axis.text.y = element_text(size = 8, color = IBD.colmap$same_dir_signif)
#   )
# 
# p1 + p2 + plot_layout(guides = "collect")

metag_maaslin <- function(factors, df, metadata, sample_names_suffix, sample_subset_string_name, sample_loc){
  
  d <- df
  meta <- metadata
  
  if (sample_loc != ".all"){
    samp_names <- meta %>% filter(grepl(sample_names_suffix, sample_id) & grepl(sample_loc, location)) %>% pull(wat)
  }
  else{
    samp_names <- meta %>% filter(grepl(sample_names_suffix, sample_id)) %>% pull(wat)
  }
  
  d.sh <- d %>% select(-1) %>%
    group_by(taxonomy) %>%
    summarize(noftax = n()) %>%
    right_join(d, by = "taxonomy") %>%
    select(one_of(c("taxonomy", "noftax", samp_names))) %>%
    mutate(row_num = row_number()) %>%
    mutate(tax_fixed = ifelse(noftax>1, paste0(taxonomy, "_otu", row_num), taxonomy)) %>%
    column_to_rownames("tax_fixed") %>%
    select(one_of(samp_names)) %>%
    t()
  
  m <- meta %>% select(one_of(c(factors, "sample_id"))) %>% column_to_rownames("sample_id")
  
  fit_data <- Maaslin2(input_data=d.sh, 
                       input_metadata=m, 
                       output = paste0("04-analysis/01_maaslin2/metagenomics_", sample_subset_string_name), 
                       fixed_effects = factors,
                       reference = c("mouse_sex,F", "OTU_type,otu_1", "IBD_stat,Healthy", "infl_score,0"),
                       normalization="TSS", 
                       transform ="LOG",
                       plot_scatter = FALSE,
                       cores = 1)
  
  fit_data.df <- fit_data[["results"]] %>% as.data.frame() %>%
    separate_wider_delim(feature, delim = regex("(k__)|(..p__)|(..c__)|(..o__)|(..g__)|(..f__)|(..s__)"), 
                         names = c("rm1","tax_king", "tax_phy", "tax_class", 
                                   "tax_ord", "tax_fam", "tax_gen", "tax_sp")) %>%
    select(-rm1) %>%
    mutate(otu_name = ifelse(tax_gen != "", paste0(tax_gen, "_", tax_sp), ifelse(
      tax_fam != "", paste0(tax_fam, "_", tax_sp), paste0(tax_ord, "_", tax_sp))),
      tax_genfixed = ifelse(tax_gen == "", ifelse(tax_fam == "", ifelse(tax_ord == "", tax_class, tax_ord), tax_fam), tax_gen)
    )
  
  fit_data.df.colormap <- fit_data.df %>% filter(metadata != "mouse_sex" & qval <= .05) %>% 
    group_by(tax_genfixed, metadata) %>% summarize(meancoef = mean(coef)) %>%
    pivot_wider(names_from = metadata, values_from = meancoef) %>%
    mutate(same_dir_signif = ifelse(is.na(IBD_stat)|is.na(OTU_type), "ivory4", 
                                    ifelse(IBD_stat * OTU_type > 0, "cornflowerblue", "coral1"))
           
    )
  
  otu.colmap <- fit_data.df.colormap %>% filter(!is.na(OTU_type)) %>% 
    arrange(OTU_type)
  IBD.colmap <- fit_data.df.colormap %>% filter(!is.na(IBD_stat)) %>% 
    arrange(IBD_stat)
  
  p1 <- fit_data.df %>% filter(metadata == factors[1] & qval <= .05) %>%
    ggplot(aes(x = coef, y = reorder(tax_genfixed, coef), color = tax_phy)) +
    geom_point(alpha = .65) +
    theme_classic() +
    labs(x = "Log2FC OTU type 1 vs OTU type 2",
         y = "") +
    scale_color_discrete(name = "Phylum") +
    geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
    theme(
      axis.text.y = element_text(size = 8, color = otu.colmap$same_dir_signif)
    )
  
  p2 <- fit_data.df %>% filter(metadata == factors[2] & qval <= .05) %>%
    ggplot(aes(x = coef, y = reorder(tax_genfixed, coef), color = tax_phy)) +
    geom_point(alpha = .65) +
    theme_classic() +
    labs(x = "Log2FC Healthy vs IBD",
         y = "") +
    scale_color_discrete(name = "Phylum") +
    geom_vline(xintercept = 0, linetype="dotted", linewidth=.5) +
    theme(
      axis.text.y = element_text(size = 8, color = IBD.colmap$same_dir_signif)
    )
  
  pall <- p1 + p2  + 
    plot_annotation(sample_subset_string_name, theme=theme(plot.title=element_text(hjust=0.5)))
  
  return(list(fit_df = fit_data.df, plot = pall))
  
}

factors <- c("OTU_type", "IBD_stat", "mouse_sex")

for (x in unique(meta$location)){
  maaslin2_dtt <- metag_maaslin(factors, dat, meta, sample_names_suffix = "-DTT", "Mucosal", sample_loc = x)
  maaslin2_dtt$plot
  ggsave(paste0("metagenomics_maaslin_mucosal_",x,".pdf"), height = 11, width = 14)
}

maaslin2_dtt <- metag_maaslin(factors, dat, meta, "-DTT", "Mucosal", sample_loc = "colon")
maaslin2_dtt$plot
ggsave("metagenomics_maaslin_mucosal_colon.pdf", height = 9, width = 11)

maaslin2_fp <- metag_maaslin(factors, dat, meta, "-FP", "Lumen", sample_loc = "colon")
maaslin2_fp$plot
ggsave("metagenomics_maaslin_lumen_colon.pdf", height = 9, width = 11) 

factors1 <- c("infl_score", "mouse_sex")

maaslin2_dtt_infl <- metag_maaslin(factors1, dat, meta, "-DTT", "Mucosal", sample_loc = "colon")
maaslin2_dtt_infl$plot
ggsave("metagenomics_maaslin_mucosal_colon_infl.pdf", height = 9, width = 11)

maaslin2_fp_infl <- metag_maaslin(factors1, dat, meta, "-FP", "Lumen", sample_loc = "colon")
maaslin2_fp_infl$plot
ggsave("metagenomics_maaslin_lumen_colon_infl.pdf", height = 9, width = 11) 
  

###ALPHA DIVERSITY###  

seq_names <- meta %>% filter(grepl("-DTT", sample_id)) %>% pull(seq_id)

dat.filt <- dat %>% select(-1) %>% select(one_of("taxonomy", seq_names))

dat.clean <- dat.filt %>% group_by(taxonomy) %>%
  summarize(noftax = n()) %>%
  right_join(dat.filt, by = 'taxonomy') %>%
  mutate(row_num = row_number()) %>%
  mutate(tax_fixed = ifelse(noftax>1, paste0(taxonomy, "_otu", row_num), taxonomy)) %>%
  column_to_rownames("tax_fixed") %>% select(-one_of(c("taxonomy", "noftax", "row_num"))) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>%
  filter(!grepl("Blank", rowname)) %>%
  column_to_rownames()


totcount <- rowSums(dat.clean) %>% as.data.frame() %>% rownames_to_column(var = "sample_id")
colnames(totcount)[colnames(totcount) == "."] <- "countsum" 
specnum <- specnumber(dat.clean) %>% as.data.frame() %>% rownames_to_column(var = "sample_id")
colnames(specnum)[colnames(specnum) == "."] <- "totspec" 
qctable <- left_join(totcount, specnum, by = "sample_id")

pqc1 <- qctable %>% arrange(countsum)%>%
  group_by((row_number()-1) %/% (n()/2)) %>%
  mutate(group = sample(1:1000, 1),
         logrowsum = log10(countsum)) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, countsum), y = countsum)) + geom_bar(stat = "identity") +
  #facet_wrap(~group, scales = "free_x", nrow = 2) +
  ylab("Sum of all OTU counts") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2),
        axis.title.x = element_blank(),
        strip.text = element_blank())

#ggsave("totalcountmetagenomics_lumen.pdf", width = 25, height = 10)

pqc2 <- qctable %>% arrange(countsum)%>%
  group_by((row_number()-1) %/% (n()/2)) %>%
  mutate(group = sample(1:1000, 1),
         logrowsum = log10(totspec)) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, countsum), y = totspec)) + geom_bar(stat = "identity") +
  #facet_wrap(~group, scales = "free_x", nrow = 2) +
  ylab("Total unique OTU count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2),
        axis.title.x = element_blank(),
        strip.text = element_blank())

#ggsave("totalspeciesmetagenomics.pdf", width = 25, height = 10)

pqc3 <- qctable %>%
  ggplot(aes(x = countsum, y = totspec)) + geom_point() +
  #facet_wrap(~group, scales = "free_x", nrow = 2) +
  ylab("Total unique OTU") +
  xlab("Sum of all OTU counts")

design <- "AAAC
           BBBD"

pqc1  + pqc3 + pqc2 + p5
  plot_layout(design = design)

ggsave("metagenomicstotalcount_mucosal.pdf", width = 25, height = 10)
  
##rarefy species count

  
  














####=============== CODE GRAVEYARD ============###


# dat_ln_DTT <- dat %>% select(-1) %>%
#               select(one_of(c("taxonomy", samp_names_mucosal))) %>%
#               pivot_longer(cols = all_of(samp_names),
#                                names_to = "sample_id",
#                                values_to = "count")
# dat_sh_DTT <- dat %>% select(-1) %>%
#   group_by(taxonomy) %>%
#   summarize(noftax = n()) %>%
#   right_join(dat, by = "taxonomy") %>%
#   select(one_of(c("taxonomy", "noftax", samp_names_mucosal))) %>%
#   mutate(row_num = row_number()) %>%
#   mutate(tax_fixed = ifelse(noftax>1, paste0(taxonomy, row_num), taxonomy)) %>%
#   column_to_rownames("tax_fixed") %>%
#   select(one_of(samp_names_mucosal)) %>%
#   t()
# 
# dat_sh_DTT.dist <- dat_sh_DTT %>% vegdist(method = "bray", na.rm = TRUE)
# 
# dat_sh_DTT.pca <- wcmdscale(dat_sh_DTT.dist, eig = TRUE)$points %>% as.data.frame() %>% select(1,2) %>%
#   rownames_to_column(var = "sample_id") %>% 
#   left_join(meta, by = "sample_id")
# 
# dat_sh_DTT.eig <- wcmdscale(dat_sh_DTT.dist, eig = TRUE)$eig
# 
# p <- ggplot(dat_sh_DTT.pca, aes(x = Dim1, y = Dim2, text = sample_id, color = OTU_type)) +
#   geom_point() +
#   theme_pubclean() +
#   labs(x = paste0("PC1"," (",signif(dat_sh_DTT.eig[1], digits = 2),"%)"), 
#        y = paste0("PC2"," (",signif(dat_sh_DTT.eig[2], digits = 2),"%)"),
#        title = "Mucosal") +
#   scale_color_discrete(name = "",
#                        labels = c("OTU type 1",
#                                 "OTU type 2"))+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p2 <- ggplot(dat_sh_DTT.pca, aes(x = Dim1, y = Dim2, text = sample_id, color = IBD_stat)) +
#   geom_point() +
#   theme_pubclean() +
#   labs(x = paste0("PC1"," (",signif(dat_sh_DTT.eig[1], digits = 2),"%)"), 
#        y = element_blank()) +
#   scale_color_manual(values = c("Healthy" = "#b37d8b", "IBD" = "#b3d4ff"),
#                      name = "",
#                      labels = c("Healthy",
#                                 "IBD"))

