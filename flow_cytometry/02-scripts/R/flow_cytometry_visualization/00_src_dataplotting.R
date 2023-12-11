# library(vioplot)
# 
# library(FSA)
# 
# library(ggplot2)
# 
# library(ggpubr)
# 
# library(lmerTest)
# 
# library(tidyr)
# install.packages(c("vioplot","FSA","ggplot2","ggpubr","lmerTest"))
# 
# getwd()
# 
# setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Jacob's Lab/TL1A/TL1A")
# 
# 
# 
# data<-read.csv("Fibrosis DUO.csv",header=T)
# 
# names(data)
# 
# data$OTU.type<-factor(data$OTU.type,levels=c("Healthy OTU 1","Healthy OTU 2","IBD OTU 1","IBD OTU 2"))
# 
# data$Area.Fibrosis.DUO<-as.numeric(data$Area.Fibrosis.DUO)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###Stats
# 
# kruskal.test(Area.Fibrosis.DUO~OTU.type,data=data)
# 
# dunnTest(Area.Fibrosis.DUO~OTU.type, data=data,method="bh")
# colors=c("Healthy OTU 1"="#B0E0E6","Healthy OTU 2"="#0000FF","IBD OTU 1"="#F0CCB0","IBD OTU 2"="#D2691E")
# ggplot(data=data,aes(x=OTU.type,y=Area.Fibrosis.DUO,fill=OTU.type))+
#   scale_fill_manual(values=colors)+
#   geom_violin(alpha=0.25,size=1,color="black",draw_quantiles=c(0.5))+
#   geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.6,aes(fill=OTU.type))+
#   theme_pubr()
# 
# 
# 
# ###Plots
# 
# Compare_vector <- list(c("Healthy OTU 1", "Healthy OTU 2"),
#                        
#                        c("Healthy OTU 1", "IBD OTU 1"),
#                        
#                        c("Healthy OTU 1", "IBD OTU 2"),
#                        
#                        c("Healthy OTU 2", "IBD OTU 1"),
#                        
#                        c("Healthy OTU 2", "IBD OTU 2"),
#                        
#                        c("IBD OTU 1", "IBD OTU 2"))
# 
# 
# 
# ggplot(data=data,aes(x=Group,y=Score_1_inflammatory_infiltration,fill=Group))+scale_fill_viridis_d(option="D")+
#   geom_violin(alpha=0.25,size=1,color="black",draw_quantiles=c(0.5))+geom_dotplot(binaxis = "y", stackdir = "center",
#                                                                                   dotsize=0.6,aes(fill=Group))+
#   theme_pubr()+
#   
#   theme(plot.title = element_text(hjust = 0.5)) +
#   
#   ggtitle("Subscore inflammation")+
#   
#   labs(y="Score_1_inflammatory_infiltration", x="")+
#   
#   stat_compare_means(comparisons = Compare_vector,
#                      
#                      method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(rstatix)
setwd("~/UCLA Documents/Rotations/jjacobs_fall2023/tl1a_humanizedmice/flow_cytometry")

dc <- read_excel("03-data/raw_data/dc_percenttotal.xlsx", sheet = "clean", na = c("na", "", "NA", "#N/A")) %>% 
      pivot_longer(colnames(select(., matches("CD"))), names_to = "marker", values_to = "percent_total") %>%
      group_by(mouse_id, location) %>%
      mutate(countperE6 = percent_total * 1e4, logcountperE6 = log10(countperE6 + 1),
             location = factor(location, levels = c("duodenum", "jejunum", "ileum", "colon")),
             percentcd11c = if_else(marker != "CD11c", countperE6 / countperE6[marker == "CD11c"], 1) * 100,
             percentcd11c_corr = if_else(percentcd11c > 100, NA, percentcd11c))

tcel <- read_excel("03-data/raw_data/cd3_percenttotal.xlsx", sheet = "clean", na = c("na", "", "NA", "#N/A"))%>% 
      pivot_longer(colnames(select(., matches(c("Foxp", "ROR", "CD")))), names_to = "marker", values_to = "percent_total") %>%
      group_by(mouse_id, location) %>%    
      mutate(countperE6 = percent_total * 1e4, logcountperE6 = log10(countperE6 + 1),
             location = factor(location, levels = c("duodenum", "jejunum", "ileum", "colon")),
             percentcd3 = if_else(marker != "CD3+", percent_total / percent_total[marker == "CD3+"], 1) * 100)

colmap=c("Healthy_otu_1"="#B0E0E6","Healthy_otu_2"="#0000FF","IBD_otu_1"="#F0CCB0","IBD_otu_2"="#D2691E")

dc.stat <- dc %>% filter(marker != "CD11c" & marker != "CD11c+|CD11b+|CD103-") %>% 
  group_by(location, marker) %>% dunn_test(percentcd11c ~ donor_type, p.adjust.method = "BH") %>% filter(p.adj <= 0.04)

p1 <- dc %>% filter(marker != "CD11c" & marker != "CD11c+|CD11b+|CD103-") %>%
  ggplot(aes(x = donor_type, y = percentcd11c)) +
  geom_violin(aes(color = donor_type), position = position_dodge(0.9), scale = "width") + 
  geom_boxplot(position = position_dodge(0.9), width = .1, aes(color = donor_type)) +
  scale_color_manual(values=colmap, name = "Donor Status", 
                     labels = c("Healthy - OTU type 1",
                                "Healthy - OTU type 2",
                                "IBD - OTU type 1",
                                "IBD - OTU type 2")) + 
  ylim(c(-1, 85)) +
  facet_grid(location~marker, labeller = labeller(location = c(duodenum = "Duodenum", jejunum = "Jejunum", ileum = "Ileum", colon = "Colon")))  +
  labs(y = "% of CD11c cells") +
  theme_pubclean() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  stat_compare_means(method = "kruskal.test", label.y = 79, size = 3, aes(group = donor_type)) +
  stat_pvalue_manual(dc.stat, label = "p.adj.signif", 
                     y.position = c(62, 66, 58, 70),
                     tip.length = 0,
                     size = 2.5)

labeling <- dc %>% filter(marker != "CD11c" & marker != "CD11c+|CD11b+|CD103-") %>% 
            separate_wider_delim(marker, "|", names = c("CD11c","CD11b", "CD103"), cols_remove = FALSE) %>%
            mutate_all(~gsub("CD11c|CD11b|CD103", "", .))

p2 <- labeling %>%
      pivot_longer(c("CD11c","CD11b", "CD103")) %>%
      ggplot(aes(x = marker, y = name, label = value)) +
      geom_text(size = 6) +
      labs(x = NULL, y = NULL) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_text(size = 12)) +
      coord_fixed(ratio = .2)

p1 / p2
# ggsave("DC_percentofCD11c_v2.pdf", width = 7, height = 12)

dc_a <- dc %>% filter(marker == "CD11c") %>% distinct() %>%
        mutate(cd11cpercenttotal = percent_total) %>%
        select(one_of("mouse_id", "location", "cd11cpercenttotal")) %>%
        right_join(dc, by = c("mouse_id", "location")) %>% filter(marker != "CD11c") %>%
        group_by(location, donor_type, marker) %>%
        summarize(meanpercentcd11c = mean(percentcd11c, na.rm = TRUE))

dcx <- dc_a %>% group_by(location, donor_type) %>%
       summarize(sumpercentcd11c = sum(meanpercentcd11c)) %>%
       mutate("CD11c+|CD11b-|CD103-" = 100 - sumpercentcd11c) %>%
       pivot_longer("CD11c+|CD11b-|CD103-", names_to = "marker", values_to = "meanpercentcd11c") %>%
       select(-sumpercentcd11c) %>% bind_rows(dc_a)

colmap2=c("CD11c+|CD11b-|CD103-"="lightgrey","CD11c+|CD11b-|CD103+"="#e07a5f","CD11c+|CD11b+|CD103-"="#fde4cf","CD11c+|CD11b+|CD103+"="#bee9e8")

p1a <- dcx %>%
  ggplot(aes(fill = marker, y = meanpercentcd11c, x = donor_type)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(location~., labeller = labeller(location = c(duodenum = "Duodenum", jejunum = "Jejunum", ileum = "Ileum", colon = "Colon"))) +
  scale_x_discrete(name = "", 
                   labels = c("Healthy - OTU type 1",
                              "Healthy - OTU type 2",
                              "IBD - OTU type 1",
                              "IBD - OTU type 2"),
                   guide = guide_axis(angle = 45)) +
  scale_fill_manual(name = "Cell Marker", values = colmap2) +
  labs(y = "Mean % of cell type w.r.t CD11C+ cells") +
  theme_pubclean()+
  guides(fill = guide_legend(nrow = 4))

plotdesign <- "
12
32

"

p1 + p1a + p2 + plot_layout(design = plotdesign, widths = c(2,1))

# ggsave("DC_meanpercentcelltype.pdf", width = 3.5, height = 8)

##TCELL

tcell.stat <- tcel %>% filter(marker != "RORγt+/Foxp3+" & marker != "CD3+") %>% 
  group_by(location, marker) %>% dunn_test(percentcd3 ~ donor_type, p.adjust.method = "BH") %>% filter(p.adj <= 0.04)


p3_raw <- tcel %>% filter(marker != "RORγt+/Foxp3+" & marker != "CD3+") %>%
  ggplot(aes(x = donor_type, y = percentcd3)) +
  geom_violin(aes(color = donor_type), position = position_dodge(0.9), scale = "width") + 
  geom_boxplot(position = position_dodge(0.9), width = .1, aes(color = donor_type)) +
  scale_color_manual(values=colmap, name = "Donor Status\n", 
                   labels = c("Healthy - OTU type 1",
                              "Healthy - OTU type 2",
                              "IBD - OTU type 1",
                              "IBD - OTU type 2")) +
  ylim(c(-1, 80)) +
  facet_grid(location~marker, labeller = labeller(location = c(duodenum = "Duodenum", jejunum = "Jejunum", ileum = "Ileum", colon = "Colon"))) +
  stat_compare_means(method = "kruskal.test", label.y = 76, size = 3, aes(group = donor_type)) +
  labs(y = "% of CD3+ cells") +
  theme_pubclean() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  stat_pvalue_manual(tcell.stat, label = "p.adj.signif", 
                     y.position = c(58, 61, 64, 67, 70, 73, 58, 61, 64, 67, 70),
                     tip.length = 0,
                     size = 2)

labelingcd3 <- tcel %>% filter(marker != "RORγt+/Foxp3+" & marker != "CD3+") %>% 
  separate_wider_delim(marker, "|", names = c("CD3","Foxp3", "RORγt"), cols_remove = FALSE) %>%
  mutate_all(~gsub("Foxp3+|RORγt+|CD3+", "", .))

p4 <- labelingcd3 %>%
  pivot_longer(c("CD3","Foxp3", "RORγt")) %>%
  ggplot(aes(x = marker, y = name, label = value)) +
  geom_text(size = 6) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12)) +
  coord_fixed(ratio = .2)

p3_raw/p4

# ggsave("TC_percentofCD3.pdf", width = 6, height = 9.5)


tc_a <- tcel %>% filter(marker != "RORgt+/Foxp3+" & marker != "CD3+") %>% distinct() %>%
  group_by(location, donor_type, marker) %>%
  summarize(meanpercentcd3 = mean(percentcd3, na.rm = TRUE))

tcx <- tc_a %>% group_by(location, donor_type) %>%
  summarize(sumpercentcd3 = sum(meanpercentcd3)) %>%
  mutate("CD3+|Foxp3-|RORgt-" = 100 - sumpercentcd3) %>%
  pivot_longer("CD3+|Foxp3-|RORgt-", names_to = "marker", values_to = "meanpercentcd3") %>%
  select(-sumpercentcd3) %>% bind_rows(tc_a)

colmap4=c("CD3+|Foxp3-|RORgt-"="lightgrey","CD3+|Foxp3-|RORgt+"="#e07a5f","CD3+|Foxp3+|RORgt-"="#fde4cf","CD11c+|CD11b+|CD103+"="#bee9e8")

p4a <- tcx %>%
  ggplot(aes(fill = marker, y = meanpercentcd3, x = donor_type)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(location~., labeller = labeller(location = c(duodenum = "Duodenum", jejunum = "Jejunum", ileum = "Ileum", colon = "Colon"))) +
  scale_x_discrete(name = "", 
                   labels = c("Healthy - OTU type 1",
                              "Healthy - OTU type 2",
                              "IBD - OTU type 1",
                              "IBD - OTU type 2"),
                   guide = guide_axis(angle = 45)) +
  scale_fill_manual(name = "Cell Marker", values = colmap4) +
  labs(y = "Mean % of cell type w.r.t CD3+ cells") +
  theme_pubclean()+
  guides(fill = guide_legend(nrow = 4))

# ggsave("TC_meanpercentcelltype.pdf", width = 3.5, height = 8)

tcell.stat2 <- tcel %>% filter(marker == "RORγt+/Foxp3+") %>% 
  group_by(location) %>%
  mutate(log2ratio = log2(percentcd3)) %>% dunn_test(log2ratio ~ donor_type, p.adjust.method = "BH") %>% filter(p.adj <= 0.05)


p3_ratio <- tcel %>% filter(marker == "RORγt+/Foxp3+" & percentcd3 >0) %>%
  mutate(log2ratio = log2(percentcd3)) %>%
  ggplot(aes(x = donor_type, y = log2ratio)) +
  geom_violin(aes(color = donor_type), position = position_dodge(0.9)) + 
  geom_boxplot(position = position_dodge(0.9), width = .1, aes(color = donor_type)) +
  scale_color_manual(values=colmap, name = "Donor Status\n", 
                     labels = c("Healthy - OTU type 1",
                                "Healthy - OTU type 2",
                                "IBD - OTU type 1",
                                "IBD - OTU type 2")) +
  facet_grid(location~., labeller = labeller(location = c(duodenum = "Duodenum", jejunum = "Jejunum", ileum = "Ileum", colon = "Colon"))) +
  stat_compare_means(method = "kruskal.test", label.y = 9, size = 3, aes(group = donor_type)) +
  labs(y = "Log2 Ratio of RORγt+ and FOXP3+ T-cells") +
  theme_pubclean() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(colour = guide_legend(nrow = 2))
# +
#   stat_pvalue_manual(tcell.stat2, label = "p.adj.signif",
#                        y.position = c(7.2, 7.5, 7.9, 8.2, 8.4),
#                        tip.length = 0,
#                        size = 2)
    

# ggsave("TC_CD3_ratio_nolines.pdf", width = 4.5, height = 9.5)
