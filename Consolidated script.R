#load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(emmeans)
library(rstatix)
library(this.path)
library(DESeq2)
library(Seurat)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(limma)
library(outliers)
library(readr)
library(lubridate)
library(MESS)
library(zoo)

#set working directory to location of this R script
setwd(dirname(this.path::this.path()))

### Figure 1 and S2 ### ----
## Fig 1A ## ----

#clear environment
rm(list = ls())

#import data - filed downloaded from HPAP 13 August 2025
donor_info <- read_excel("data/HPAP/Donor_Summary_192.xlsx",sheet = "donor")
donor_info <- donor_info %>%
  mutate(alpha.endocrine = (alpha_cell_count/(alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)),
         beta.endocrine = (beta_cell_count/(alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)),
         delta.endocrine = (delta_cell_count/(alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)),
         epsilon.endocrine = (epsilon_cell_count/(alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)),
         pp.endocrine = (pp_cell_count/(alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)),
         total_endocrine_count = alpha_cell_count + beta_cell_count + delta_cell_count + epsilon_cell_count + pp_cell_count)

#rename groups
donor_info <- donor_info %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#filter out donors with too few cells or no cells
no_cells <- donor_info %>% filter(total_endocrine_count == 0) %>% dplyr::select(donor_ID)

donor_info %>% #calculate geometric mean
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  summarise(geommean = exp(mean(log(total_endocrine_count)))) %>%
  dplyr::select(geommean) #31512

low_cells <- donor_info %>% #identify those with more than 30-fold below geometric mean
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(total_endocrine_count < 31512/30) %>%
  dplyr::select(donor_ID, total_endocrine_count)

#statistical analysis comparing sexes among donors without diabetes aged 15-39
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(alpha.endocrine ~ age_years + sex) #p=0.028 for sex

donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(beta.endocrine ~ age_years + sex) #p=0.112 for sex

#Graphs
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  ggplot(aes(x=sex, y = alpha.endocrine*100))+
  geom_boxplot(aes(fill = sex), alpha = 0.5)+
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("")+
  ylab("Proportion of endocrine cells\n(%)") +
  ylim(0,120)+
  ggtitle("Alpha-cell proportion") +
  geom_bracket(xmin = 1, xmax = 2, y.position = 100, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Fig1/Fig1A alpha.png", width = 3, height = 3)

donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  ggplot(aes(x=sex, y = beta.endocrine*100))+
  geom_boxplot(aes(fill = sex), alpha = 0.5)+
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("")+
  ylab("Proportion of endocrine cells\n(%)") +
  ggtitle("Beta-cell proportion") +
  ylim(0,110)+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Fig1/Fig1A beta.png", width = 3, height = 3)

## Fig S2A ## ----
donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(alpha.endocrine ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(alpha.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0099 for F-M in Control
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0775 for Control vs T2D in males
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = alpha.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 94, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 99, label = "*", label.size = 7, size = 0.5, tip.length = c(0.04, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,110))+
  ylab("Proportion of endocrine cells\n(%)") +
  labs(fill = "Condition", colour = "Condition")+
  ggtitle("Alpha-cell proportion")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS2/FigS2A alpha.png", width = 5, height = 4)

donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(beta.endocrine ~ age_years + sex*simplified_diagnosis) #p=0.092 for disease, 0.099 for interaction
model <- donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(beta.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0723 for F-M in Control
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0192 for Control vs T2D in males
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = beta.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 91, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.4))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 96, label = "p=0.072", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,100))+
  ylab("Proportion of endocrine cells\n(%)") +
  labs(fill = "Condition", colour = "Condition")+
  ggtitle("Beta-cell proportion")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS2/FigS2A beta.png", width = 5, height = 4)

## Fig 1B ## ----
#clear environment
rm(list = ls())

#import data - filed downloaded from Humanislets.com 12 February 2025
celltypeprop <- read.csv("data/Humanislets.com/cell_pro.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")

#join metadata to data
celltypeprop <- left_join(celltypeprop, donor_data, by = "record_id")

#rename groups
celltypeprop <- celltypeprop %>%
  mutate(simplified_diagnosis = case_when(
    diagnosis == "None" ~ "Control",
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
  ))

#statistical analyses
celltypeprop %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  anova_test(alpha_end ~ donorage + donorsex) #p=0.037 for sex

celltypeprop %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  anova_test(beta_end ~ donorage + donorsex) #p=0.025 for sex

#Graphs
celltypeprop %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y = alpha_end*100))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("")+
  ylab("Proportion of endocrine cells\n(%)") +
  ylim(15,45)+
  ggtitle("Alpha-cell proportion") +
  geom_bracket(xmin = 1, xmax = 2, y.position = 35, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Fig1/Fig1B alpha.png", width = 3, height = 3)

celltypeprop %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y = beta_end*100))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("")+
  ylab("Proportion of endocrine cells\n(%)") +
  ggtitle("Beta-cell proportion") +
  ylim(40,70)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 62, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Fig1/Fig1B beta.png", width = 3, height = 3)

## Fig S2B ## ----
celltypeprop %>%  
  filter(diagnosis != "T1D") %>%
  anova_test(alpha_end ~ donorage + donorsex*diagnosis) #sig for disease and interaction
model <- celltypeprop %>%
  filter(diagnosis != "T1D") %>%
  lm(alpha_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #sig  for F-M in T2D
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for Control vs T2D in females
celltypeprop %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = alpha_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 43, label = "*", label.size = 7, size = 0.5, tip.length = c(0.3, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 42, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(20,50))+
  ylab("Proportion of endocrine cells\n(%)") +
  labs(fill = "Condition", colour = "Condition")+
  ggtitle("Alpha-cell proportion")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/FigS2/FigS2B alpha.png", width = 5, height = 4)

celltypeprop %>%  
  filter(diagnosis != "T1D") %>%
  anova_test(beta_end ~ donorage + donorsex*diagnosis) #sig for disease and interaction
model <- celltypeprop %>%
  filter(diagnosis != "T1D") %>%
  lm(beta_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #sig  for F-M in T2D
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for Control vs T2D in females
celltypeprop %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = beta_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 67, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 66, label = "*", label.size = 7, size = 0.5, tip.length = c(0.53, 0.02))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(40,75))+
  ylab("Proportion of endocrine cells\n(%)") +
  labs(fill = "Condition", colour = "Condition")+
  ggtitle("Beta-cell proportion")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/FigS2/FigS2B beta.png", width = 5, height = 4)



### Figure 2 ### ----
## Fig 2A-B ## ----

#clear environment
rm(list = ls())

#Download and perform DE analysis on pseudobulked beta-cell and alpha-cell scRNAseq data from HPAP
dat <- readRDS("data/HPAP/T1D_T2D_20220428.rds")  #downloaded 7 January 2024
beta <- subset(x = dat, idents = "Beta") #subset by cell type
alpha <- subset(x = dat, idents = "Alpha")
beta_bulk <- AggregateExpression(beta, group.by = c("hpap_id"), return.seurat = TRUE) #pseudobulk
beta_bulk_rawcounts_df <- data.frame(beta_bulk[["RNA"]]$counts) #retrieve raw counts
alpha_bulk <- AggregateExpression(alpha, group.by = c("hpap_id"), return.seurat = TRUE)
alpha_bulk_rawcounts_df <- data.frame(alpha_bulk[["RNA"]]$counts)

donor_all <- read_excel("data/HPAP/Donor_Summary_192.xlsx")
donor_all$donor_ID <- gsub("-","", donor_all$donor_ID) #match donor IDs to those in the Seurat

#rename groups
donor_all <- donor_all %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#DESeq on pseudobulked beta cell data to compare young donors without diabetes
donorlist <- colnames(beta_bulk_rawcounts_df)
betadonor_deseq <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis=="Control") %>%
  filter(age_years > 14, age_years < 40)
donorlist <- betadonor_deseq$donor_ID #including only young donors without diabetes
beta_bulk_rawcounts_df <- beta_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(betadonor_deseq) <- betadonor_deseq$donor_ID
all(rownames(betadonor_deseq) %in% colnames(beta_bulk_rawcounts_df)) #check the samples match

#RunDEseq
betadds <- DESeqDataSetFromMatrix(countData = beta_bulk_rawcounts_df,
                                  colData = betadonor_deseq,
                                  design = ~ sex + age_years) 

#remove low count genes
table(betadonor_deseq$sex, betadonor_deseq$simplified_diagnosis) #7 females, 17 males
smallestGroupSize <- 7
keep <- rowSums(counts(betadds) >= 10) >= smallestGroupSize
betadds <- betadds[keep,]

#make males the control for disease state for comparison
betadds$sex <- relevel(betadds$sex, ref = "Male")
betadds <- DESeq(betadds)
betares <- results(betadds)

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(betadds)
betaresLFCdisease <- lfcShrink(betadds, coef="sex_Female_vs_Male", type="apeglm")
betaresLFCdiseaseOrdered <- betaresLFCdisease[order(betaresLFCdisease$pvalue),] #reorder by p value
betaresLFCdiseaseOrdered

#write csv file
write.csv(betaresLFCdiseaseOrdered, "Output/Fig2/HPAP beta scRNAseq DEseq2 ctrls 15-39.csv")

#GSEA
#GSEA function
gsea_human <- function(df,FC_col, p_col, name){ 
    # remove duplicated entrez
    df <- df[order(df[[p_col]]),]
  df <- df[!duplicated(df$entrez),]
  df <- df[!is.na(df$entrez), ]
  
  #rank by signed -log10 p value
  genelist <- sign(df[[FC_col]]) * (-log10(df[[p_col]]))
  
  # entrez id as names of the gene list
  names(genelist) <- df$entrez
  
  genelist[genelist==Inf] <- max(genelist[is.finite(genelist)])+1
  genelist[genelist==-Inf] <- min(genelist[is.finite(genelist)])-1
  
  genelist = sort(genelist, decreasing = TRUE)
  genelist <- na.omit(genelist)
  
  gseGO.bp <- gseGO(
    geneList=genelist,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <- as.data.frame(gseGO.bp)
  
  gseGO.mf <- gseGO(
    geneList=genelist,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <- as.data.frame(gseGO.mf)
  
  gseGO.cc <- gseGO(
    geneList=genelist,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <- as.data.frame(gseGO.cc)
  
  
  # Add the data frames to separate sheets
  wb <- createWorkbook()
  
  addWorksheet(wb, "GO_BP")
  writeData(wb, "GO_BP", gseGO.bp.df)
  
  addWorksheet(wb, "GO_MF")
  writeData(wb, "GO_MF", gseGO.mf.df)
  
  addWorksheet(wb, "GO_CC")
  writeData(wb, "GO_CC", gseGO.cc.df)
  
  # Save the workbook to a file
  fpath <- paste0("Output/Fig2/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
    cat(paste0("GSEA results saved in - ", fpath))
  
}

betaresLFCdiseaseOrdered <- data.frame(betaresLFCdiseaseOrdered)
betaresLFCdiseaseOrdered$Gene <- rownames(betaresLFCdiseaseOrdered)
betaresLFCdiseaseOrdered$entrez1 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
betaresLFCdiseaseOrdered$entrez2 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

betaresLFCdiseaseOrdered <- betaresLFCdiseaseOrdered %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = betaresLFCdiseaseOrdered,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "beta_scRNAseq_youngcontrols")


#DESeq on pseudobulked alpha cell data to compare young donors without diabetes
donorlist <- colnames(alpha_bulk_rawcounts_df)
alphadonor_deseq <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis=="Control") %>%
  filter(age_years > 14, age_years < 40)
donorlist <- alphadonor_deseq$donor_ID #including only young donors without diabetes
alpha_bulk_rawcounts_df <- alpha_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(alphadonor_deseq) <- alphadonor_deseq$donor_ID
all(rownames(alphadonor_deseq) %in% colnames(alpha_bulk_rawcounts_df)) #check the samples match

#RunDEseq
alphadds <- DESeqDataSetFromMatrix(countData = alpha_bulk_rawcounts_df,
                                  colData = alphadonor_deseq,
                                  design = ~ sex + age_years) 

#remove low count genes
table(alphadonor_deseq$sex, alphadonor_deseq$simplified_diagnosis) #7 females, 17 males
smallestGroupSize <- 7
keep <- rowSums(counts(alphadds) >= 10) >= smallestGroupSize
alphadds <- alphadds[keep,]

#make males the control for disease state for comparison
alphadds$sex <- relevel(alphadds$sex, ref = "Male")
alphadds <- DESeq(alphadds)
alphares <- results(alphadds)

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(alphadds)
alpharesLFCdisease <- lfcShrink(alphadds, coef="sex_Female_vs_Male", type="apeglm")
alpharesLFCdiseaseOrdered <- alpharesLFCdisease[order(alpharesLFCdisease$pvalue),] #reorder by p value
alpharesLFCdiseaseOrdered

#write csv file
write.csv(alpharesLFCdiseaseOrdered, "Output/Fig2/HPAP alpha scRNAseq DEseq2 ctrls 15-39.csv")

#GSEA
alpharesLFCdiseaseOrdered <- data.frame(alpharesLFCdiseaseOrdered)
alpharesLFCdiseaseOrdered$Gene <- rownames(alpharesLFCdiseaseOrdered)
alpharesLFCdiseaseOrdered$entrez1 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
alpharesLFCdiseaseOrdered$entrez2 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

alpharesLFCdiseaseOrdered <- alpharesLFCdiseaseOrdered %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = alpharesLFCdiseaseOrdered,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "alpha_scRNAseq_youngcontrols")

#Comparing young controls for Humanislets.com bulk RNAseq data

#clear environment
rm(list = ls())

#import data downloaded 21 May 2024
proc_rnaseq <- read.csv("data/Humanislets.com/proc_rnaseq.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")

#multiple ANCOVA with age as covariate + BH correction
proc_rnaseq_long <- proc_rnaseq %>%
  pivot_longer(!gene_id, values_to = "logCPM", names_to = "record_id")
proc_rnaseq_long <- inner_join(proc_rnaseq_long, donor_data, by = "record_id") #join metadata
proc_rnaseq_long_15to39 <- proc_rnaseq_long %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None")
table(proc_rnaseq_long_15to39$donorsex, proc_rnaseq_long_15to39$diagnosis)/length(unique(proc_rnaseq_long_15to39$gene_id))
#8 females, 18 males

#filter out genes with NA in more than n.donors in smallest group
smallest_group <- 8
to_remove <- proc_rnaseq_long_15to39 %>%
  group_by(gene_id) %>%
  summarise(na_count = sum(is.na(logCPM))) %>%
  filter(na_count > smallest_group)
proc_rnaseq_long_15to39 <- proc_rnaseq_long_15to39 %>%
  filter(!gene_id %in% to_remove$gene_id)
ngenes <- length(unique(proc_rnaseq_long_15to39$gene_id))

proc_rnaseq_ancovas_15to39 <- list()
for (i in 1:ngenes){ 
  geneid <- unique(proc_rnaseq_long_15to39$gene_id)[i]
  data <- proc_rnaseq_long_15to39 %>%
    filter(gene_id == geneid)
  test.f <- data %>% #to avoid errors in t-test for genes with not enough values
    filter(donorsex == "Female") %>% 
    drop_na(logCPM)
  test.f.val <- dim(test.f)[1]
  test.m <- data %>%
    filter(donorsex == "Male") %>% 
    drop_na(logCPM)
  test.m.val <- dim(test.m)[1]
  if (test.f.val > 1 & test.m.val >1){
    ancova <- data %>%
      anova_test(logCPM ~ donorage + donorsex)
    pval <- ancova[2,5]
    proc_rnaseq_ancovas_15to39[[i]] <- c(geneid, pval)
  } else{
    proc_rnaseq_ancovas_15to39[[i]] <- c(geneid, NA)
  }
}
proc_rnaseq_ancovas_15to39 <- data.frame(do.call(rbind, proc_rnaseq_ancovas_15to39))
colnames(proc_rnaseq_ancovas_15to39) <- c("gene_id","pval")
proc_rnaseq_ancovas_15to39$gene_id <- as.character(proc_rnaseq_ancovas_15to39$gene_id)
proc_rnaseq_ancovas_15to39$pval <- as.numeric(proc_rnaseq_ancovas_15to39$pval)
proc_rnaseq_ancovas_15to39 <- proc_rnaseq_ancovas_15to39 %>%
  filter(is.na(pval) == FALSE) #remove all with no pval

#Benjamini-Hochberg correction for multiple t tests
proc_rnaseq_ancovas_15to39$padj <- p.adjust(proc_rnaseq_ancovas_15to39$pval, method = "BH")

proc_rnaseq_ancovas_15to39.means <- proc_rnaseq_long_15to39 %>%
  group_by(gene_id, donorsex) %>%
  summarise(mean.logCPM = mean(logCPM, na.rm=TRUE))
proc_rnaseq_ancovas_15to39.means <- proc_rnaseq_ancovas_15to39.means %>%
  pivot_wider(names_from = donorsex, names_prefix = "meanlogCPM_", values_from = mean.logCPM)
proc_rnaseq_ancovas_15to39.means$gene_id <- as.character(proc_rnaseq_ancovas_15to39.means$gene_id)
proc_rnaseq_ancovas_15to39 <- inner_join(proc_rnaseq_ancovas_15to39, proc_rnaseq_ancovas_15to39.means, by = "gene_id")
proc_rnaseq_ancovas_15to39 <- proc_rnaseq_ancovas_15to39 %>%
  mutate(logFC = meanlogCPM_Male-meanlogCPM_Female)
head(proc_rnaseq_ancovas_15to39)
proc_rnaseq_ancovas_15to39 <- proc_rnaseq_ancovas_15to39 %>%
  arrange(padj)

#get common names
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
common.names <- as.list(x[mapped_genes])
common.names <- unlist(common.names)
ngenes <- length(unique(proc_rnaseq_ancovas_15to39$gene_id))

proc_rnaseq_ancovas_15to39$common.name <- rep(NA,ngenes)
for (i in 1:ngenes){
  entrez <- proc_rnaseq_ancovas_15to39$gene_id[i]
  proc_rnaseq_ancovas_15to39$common.name[i] <- common.names[entrez]
}
head(proc_rnaseq_ancovas_15to39)

write.csv(proc_rnaseq_ancovas_15to39, "Output/Fig2/Humanisletscom bulk RNAseq DEseq2 ctrls 15-39.csv")

#GSEA
gsea_human <- function(df,FC_col, p_col, name){ 
  # remove duplicated entrez
  df <- df[order(df[[p_col]]),]
  df <- df[!duplicated(df$entrez),]
  df <- df[!is.na(df$entrez), ]
  
  #rank by signed -log10 p value
  genelist <- sign(df[[FC_col]]) * (-log10(df[[p_col]]))
  
  # entrez id as names of the gene list
  names(genelist) <- df$entrez
  
  genelist[genelist==Inf] <- max(genelist[is.finite(genelist)])+1
  genelist[genelist==-Inf] <- min(genelist[is.finite(genelist)])-1
  
  genelist = sort(genelist, decreasing = TRUE)
  genelist <- na.omit(genelist)
  
  gseGO.bp <- gseGO(
    geneList=genelist,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <- as.data.frame(gseGO.bp)
  
  gseGO.mf <- gseGO(
    geneList=genelist,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <- as.data.frame(gseGO.mf)
  
  gseGO.cc <- gseGO(
    geneList=genelist,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <- as.data.frame(gseGO.cc)
  
  
  # Add the data frames to separate sheets
  wb <- createWorkbook()
  
  addWorksheet(wb, "GO_BP")
  writeData(wb, "GO_BP", gseGO.bp.df)
  
  addWorksheet(wb, "GO_MF")
  writeData(wb, "GO_MF", gseGO.mf.df)
  
  addWorksheet(wb, "GO_CC")
  writeData(wb, "GO_CC", gseGO.cc.df)
  
  # Save the workbook to a file
  fpath <- paste0("Output/Fig2/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

proc_rnaseq_ancovas_15to39$entrez1 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_15to39$common.name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
proc_rnaseq_ancovas_15to39$entrez2 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_15to39$common.name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

proc_rnaseq_ancovas_15to39 <- proc_rnaseq_ancovas_15to39 %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = proc_rnaseq_ancovas_15to39,
           FC_col = "logFC",
           p_col = "pval",
           name = "bulk_RNAseq_youngcontrols")

#Combining GSEA results to graph together
HIGO_CC <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_CC")
HIGO_BP <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_BP")
HIGO_MF <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_MF")
HIGO <- bind_rows(HIGO_CC, HIGO_BP, HIGO_MF)
HIGO$NES <- -HIGO$NES #reverse sign to change direction so positive is up in females
HIGO <- HIGO %>% filter(p.adjust < 0.05) #significant only

HPAPbetascGO_CC <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_CC")
HPAPbetascGO_BP <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_BP")
HPAPbetascGO_MF <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_MF")
HPAPbetascGO <- bind_rows(HPAPbetascGO_CC, HPAPbetascGO_BP, HPAPbetascGO_MF)
HPAPbetascGO <- HPAPbetascGO %>% filter(p.adjust < 0.05) #significant only

HPAPalphascGO_CC <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_CC")
HPAPalphascGO_BP <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_BP")
HPAPalphascGO_MF <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_MF")
HPAPalphascGO <- bind_rows(HPAPalphascGO_CC, HPAPalphascGO_BP, HPAPalphascGO_MF)
HPAPalphascGO <- HPAPalphascGO %>% filter(p.adjust < 0.05) #significant only

#Check that column names are consistent
all(colnames(HIGO) == colnames(HPAPbetascGO)) #TRUE
all(colnames(HIGO) == colnames(HPAPalphascGO)) #TRUE

#combine into single dataframe
GO_all <- bind_rows(HIGO, HPAPbetascGO, HPAPalphascGO)
GO_all$dataset <- c(rep("Humanislets.com\nbulk RNAseq",nrow(HIGO)), rep("HPAP beta-cell\nscRNAseq", nrow(HPAPbetascGO)), rep("HPAP alpha-cell\nscRNAseq", nrow(HPAPalphascGO)))
head(GO_all)

#split into female-biased and male-biased pathways
GO_Fbiased <- GO_all %>% filter(NES > 0)
GO_Mbiased <- GO_all %>% filter(NES < 0)

#identify top 40 pathways
GO_Fbiased <- GO_Fbiased %>% arrange(desc(NES))
GO_Fbiased_pathways <- unique(GO_Fbiased$Description)
GO_Fbiased_pathways <- GO_Fbiased_pathways[1:40]

GO_Mbiased <- GO_Mbiased %>% arrange(NES)
GO_Mbiased_pathways <- unique(GO_Mbiased$Description) #only 9 pathways

#Make really long pathway names two lines
GO_Fbiased <- GO_Fbiased %>%
  mutate(Description = case_when(
    Description == "oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor" ~ "oxidoreductase activity, acting on NAD(P)H,\nquinone or similar compound as acceptor",
    .default = Description))
GO_Fbiased_pathways[30] <- "oxidoreductase activity, acting on NAD(P)H,\nquinone or similar compound as acceptor"

#Graph
summary(GO_Fbiased$NES)
summary(GO_Fbiased$p.adjust)
GO_Fbiased %>%
  filter(Description %in% GO_Fbiased_pathways) %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "Female-biased GO pathways", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "darkred", mid = "red", high = "white", midpoint = 0.025, limits = c(0,0.05)) +
  scale_size_continuous(range = c(1,7), limits = c(1,3), breaks = c(1,3)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig2/Fig2A.png", width = 8.5, height = 9)

GO_Mbiased$NES <- -GO_Mbiased$NES
summary(GO_Mbiased$NES)
summary(GO_Mbiased$p.adjust)
GO_Mbiased %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "Male-biased GO pathways", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "blue", high = "white", midpoint = 0.025, limits = c(0,0.05)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.75), breaks = c(1.5,2.75)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig2/Fig2B.png", width = 7.25, height = 3)


## Fig 2C-D ##----
#clear environment
rm(list = ls())

#import proteomics data from Humanislets.com (downloaded 21 May 2024)
prot <- read.csv("data/Humanislets.com/proc_prot.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")
isolation_data <- read.csv("data/Humanislets.com/isolation.csv")

#join donor info with proteomics data
donor_isolation <- inner_join(donor_data, isolation_data, by = "record_id")
prot_long <- pivot_longer(prot, !gene_id, values_to = "Abundance", names_to = "record_id")
prot_data <- inner_join(prot_long, donor_isolation, by = "record_id")

#get protein names from gene id
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
common.names <- as.list(x[mapped_genes])
common.names <- unlist(common.names)
map <- data.frame(names(common.names), common.names)
colnames(map) <- c("gene_id","common_name")
map$gene_id <- as.numeric(map$gene_id)
prot_data <- left_join(prot_data, map, by = "gene_id")

#filter for donors without diabetes aged 15-39
prot_ND15to39 <- prot_data %>%
  filter(diagnosis == "None") %>%
  filter(donorage > 14, donorage < 40)
table(prot_ND15to39$donorsex)/length(unique(prot_ND15to39$gene_id)) #12 females, 21 males

#make separate tables for data and metadata
prot_ND15to39_feature <- prot_ND15to39 %>%
  dplyr::select(gene_id, record_id, Abundance) %>%
  pivot_wider(names_from = "record_id", values_from = "Abundance")
genes <- prot_ND15to39_feature$gene_id
prot_ND15to39_feature <- data.frame(prot_ND15to39_feature)
rownames(prot_ND15to39_feature) <- genes
prot_ND15to39_feature <- prot_ND15to39_feature[,-1] #remove the gene column

prot_ND15to39_metadata <- prot_ND15to39 %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(-Abundance) %>%
  dplyr::select(-common_name)
prot_ND15to39_metadata <- distinct(prot_ND15to39_metadata)
prot_ND15to39_metadata <- data.frame(prot_ND15to39_metadata)

#check order
all(colnames(prot_ND15to39_feature) == prot_ND15to39_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep <- apply(prot_ND15to39_feature, 1, function(x){sum(!is.na(x)) >= 9})
prot_ND15to39_feature <- prot_ND15to39_feature[feature.keep, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(prot_ND15to39_metadata[,"donorsex"]))
fixedEffects <- c("donorage") #correcting for age
all.vars <- c("donorsex", fixedEffects)
design <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = prot_ND15to39_metadata)
colnames(design)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs <- list()
ref <- "Male" #using males as reference
contrasts <- grp.nms[grp.nms != ref]
myargs <- as.list(paste("Female", "-", ref, sep = ""))
myargs[["levels"]] <- design
contrast.matrix <- do.call(makeContrasts, myargs)

# get results
fit <- lmFit(prot_ND15to39_feature, design, trend = TRUE, robust = TRUE)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
res.table <- topTable(fit, number = Inf)

# Remove results rows with NAs
res.table <- res.table[!is.na(res.table$P.Value), ] #none removed

# Process results for output
res.table <- merge(res.table, map, by.x = "row.names", by.y = "gene_id")
res.table <- res.table[,-7]

res.table <- res.table[,c(7,2:6,1,8)]
colnames(res.table)[1] <- "Feature"
colnames(res.table)[3] <- "Average level"
colnames(res.table)[4] <- "T_statistic"
colnames(res.table)[5] <- "P_value"
colnames(res.table)[6] <- "Adjusted p_value"
colnames(res.table)[7] <- "GeneID"

write.csv(res.table, "Output/Fig2/prot_dea_results_correctforage_controls_15to39.csv", row.names = FALSE)

gsea_human <- function(df,FC_col, p_col, name){ 
  # remove duplicated entrez
  df <- df[order(df[[p_col]]),]
  df <- df[!duplicated(df$entrez),]
  df <- df[!is.na(df$entrez), ]
  
  #rank by signed -log10 p value
  genelist <- sign(df[[FC_col]]) * (-log10(df[[p_col]]))
  
  # entrez id as names of the gene list
  names(genelist) <- df$entrez
  
  genelist[genelist==Inf] <- max(genelist[is.finite(genelist)])+1
  genelist[genelist==-Inf] <- min(genelist[is.finite(genelist)])-1
  
  genelist = sort(genelist, decreasing = TRUE)
  genelist <- na.omit(genelist)
  
  gseGO.bp <- gseGO(
    geneList=genelist,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <- as.data.frame(gseGO.bp)
  
  gseGO.mf <- gseGO(
    geneList=genelist,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <- as.data.frame(gseGO.mf)
  
  gseGO.cc <- gseGO(
    geneList=genelist,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <- as.data.frame(gseGO.cc)
  
  
  # Add the data frames to separate sheets
  wb <- createWorkbook()
  
  addWorksheet(wb, "GO_BP")
  writeData(wb, "GO_BP", gseGO.bp.df)
  
  addWorksheet(wb, "GO_MF")
  writeData(wb, "GO_MF", gseGO.mf.df)
  
  addWorksheet(wb, "GO_CC")
  writeData(wb, "GO_CC", gseGO.cc.df)
  
  # Save the workbook to a file
  fpath <- paste0("Output/Fig2/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

res.table$entrez1 <- mapIds(org.Hs.eg.db, keys=c(res.table$Feature), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
res.table$entrez2 <- mapIds(org.Hs.eg.db, keys=c(res.table$Feature), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

res.table <- res.table %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 
gsea_human(df = res.table,
           FC_col = "logFC",
           p_col = "P_value",
           name = "proteomics_youngcontrols")

#Separate into female-biased and male-biased
GO_CC <- read_excel("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)

GO_combined_Fbiased <- GO_combined %>% filter(NES < 0)
GO_combined_Mbiased <- GO_combined %>% filter(NES > 0)

#Plot top 30 pathways
GO_combined_Fbiased$NES <- -GO_combined_Fbiased$NES
GO_combined_Fbiased <- GO_combined_Fbiased[1:30,]
summary(GO_combined_Fbiased$NES)
summary(GO_combined_Fbiased$p.adjust)
#make long pathway name short
GO_combined_Fbiased <- GO_combined_Fbiased %>%
  mutate(Description = case_when(
    Description == "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor" ~ "oxidoreductase activity, acting on the CH-OH\ngroup of donors, NAD or NADP as acceptor",
    .default = Description
  ))

GO_combined_Fbiased %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "Female-biased GO pathways", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "darkred", mid = "red", high = "white", midpoint = 0.005, limits = c(0,0.01)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_continuous(limits = c(1.5,2.5))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig2/Fig2C.png", width = 8, height = 6.5)

GO_combined_Mbiased <- GO_combined_Mbiased[1:30,]
summary(GO_combined_Mbiased$NES)
summary(GO_combined_Mbiased$p.adjust)

GO_combined_Mbiased %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "Male-biased GO pathways", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "blue", high = "white", midpoint = 0.002, limits = c(0,0.004)) +
  scale_size_continuous(range = c(1,7), limits = c(1.4,2.1), breaks = c(1.4,2.1)) +
  scale_x_continuous(limits = c(1.4,2.1))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig2/Fig2D.png", width = 9.5, height = 6.25)

###Figure 3 ### ----
#clear environment
rm(list = ls())

#import data to compare control and T2D donors of all ages in HPAP scRNAseq datasets
dat <- readRDS("data/HPAP/T1D_T2D_20220428.rds")  #downloaded 7 January 2024
beta <- subset(x = dat, idents = "Beta") #subset by cell type
alpha <- subset(x = dat, idents = "Alpha")
beta_bulk <- AggregateExpression(beta, group.by = c("hpap_id"), return.seurat = TRUE) #pseudobulk
beta_bulk_rawcounts_df <- data.frame(beta_bulk[["RNA"]]$counts) #retrieve raw counts
alpha_bulk <- AggregateExpression(alpha, group.by = c("hpap_id"), return.seurat = TRUE)
alpha_bulk_rawcounts_df <- data.frame(alpha_bulk[["RNA"]]$counts)

donor_all <- read_excel("data/HPAP/Donor_Summary_192.xlsx")
donor_all$donor_ID <- gsub("-","", donor_all$donor_ID) #match donor IDs to those in the Seurat

#rename groups
donor_all <- donor_all %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#DESeq on pseudobulked beta cell data to compare donors without diabetes to donors with T2D
#just females
donorlist <- colnames(beta_bulk_rawcounts_df)
betadonor_deseq_F <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis!="T1D") %>%
  filter(sex == "Female")
donorlist <- betadonor_deseq_F$donor_ID #make the donorlist just the females, excluding donors with T1D
betacounts_deseq_F <- beta_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(betadonor_deseq_F) <- betadonor_deseq_F$donor_ID
all(rownames(betadonor_deseq_F) %in% colnames(betacounts_deseq_F)) #check the samples match

#RunDEseq
betadds_F <- DESeqDataSetFromMatrix(countData = betacounts_deseq_F,
                                  colData = betadonor_deseq_F,
                                  design = ~ simplified_diagnosis + age_years) 

#remove low count genes
table(betadonor_deseq_F$simplified_diagnosis) #16 Control, 11 T2D
smallestGroupSize <- 11
keep_F <- rowSums(counts(betadds_F) >= 10) >= smallestGroupSize
betadds_F <- betadds_F[keep_F,]

#make non-diabetic the control for disease state for comparison
betadds_F$simplified_diagnosis <- relevel(betadds_F$simplified_diagnosis, ref = "Control")
betadds_F <- DESeq(betadds_F)
betares_F <- results(betadds_F)
betaresOrdered_F <- betares_F[order(betares_F$pvalue),] #reorder by p value

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(betadds_F)
betaresLFCdisease_F <- lfcShrink(betadds_F, coef="simplified_diagnosis_T2D_vs_Control", type="apeglm")
betaresLFCdiseaseOrdered_F <- betaresLFCdisease_F[order(betaresLFCdisease_F$pvalue),] #reorder by p value

betaresLFCdiseaseOrdered_F <- data.frame(betaresLFCdiseaseOrdered_F)
betaresLFCdiseaseOrdered_F$Gene <- rownames(betaresLFCdiseaseOrdered_F)

#write csv file
write.csv(betaresLFCdiseaseOrdered_F, "Output/Fig3/Female ctrl vs T2D beta scRNAseq.csv")

#GSEA
gsea_human <- function(df,FC_col, p_col, name){ 
  # remove duplicated entrez
  df <- df[order(df[[p_col]]),]
  df <- df[!duplicated(df$entrez),]
  df <- df[!is.na(df$entrez), ]
  
  #rank by signed -log10 p value
  genelist <- sign(df[[FC_col]]) * (-log10(df[[p_col]]))
  
  # entrez id as names of the gene list
  names(genelist) <- df$entrez
  
  genelist[genelist==Inf] <- max(genelist[is.finite(genelist)])+1
  genelist[genelist==-Inf] <- min(genelist[is.finite(genelist)])-1
  
  genelist = sort(genelist, decreasing = TRUE)
  genelist <- na.omit(genelist)
  
  gseGO.bp <- gseGO(
    geneList=genelist,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <- as.data.frame(gseGO.bp)
  
  gseGO.mf <- gseGO(
    geneList=genelist,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <- as.data.frame(gseGO.mf)
  
  gseGO.cc <- gseGO(
    geneList=genelist,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <- as.data.frame(gseGO.cc)
  
  
  # Add the data frames to separate sheets
  wb <- createWorkbook()
  
  addWorksheet(wb, "GO_BP")
  writeData(wb, "GO_BP", gseGO.bp.df)
  
  addWorksheet(wb, "GO_MF")
  writeData(wb, "GO_MF", gseGO.mf.df)
  
  addWorksheet(wb, "GO_CC")
  writeData(wb, "GO_CC", gseGO.cc.df)
  
  # Save the workbook to a file
  fpath <- paste0("Output/Fig3/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

betaresLFCdiseaseOrdered_F$entrez1 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered_F$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
betaresLFCdiseaseOrdered_F$entrez2 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered_F$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

betaresLFCdiseaseOrdered_F <- betaresLFCdiseaseOrdered_F %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = betaresLFCdiseaseOrdered_F,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "beta_scRNAseq_F_ctrlvsT2D")

#same comparison for males
donorlist <- colnames(beta_bulk_rawcounts_df)
betadonor_deseq_M <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis!="T1D") %>%
  filter(sex == "Male")
donorlist <- betadonor_deseq_M$donor_ID #make the donorlist just the females, excluding donors with T1D
betacounts_deseq_M <- beta_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(betadonor_deseq_M) <- betadonor_deseq_M$donor_ID
all(rownames(betadonor_deseq_M) %in% colnames(betacounts_deseq_M)) #check the samples match

#RunDEseq
betadds_M <- DESeqDataSetFromMatrix(countData = betacounts_deseq_M,
                                    colData = betadonor_deseq_M,
                                    design = ~ simplified_diagnosis + age_years) 

#remove low count genes
table(betadonor_deseq_M$simplified_diagnosis) #24 control, 7 T2D
smallestGroupSize <- 7
keep_M <- rowSums(counts(betadds_M) >= 10) >= smallestGroupSize
betadds_M <- betadds_M[keep_M,]

#make non-diabetic the control for disease state for comparison
betadds_M$simplified_diagnosis <- relevel(betadds_M$simplified_diagnosis, ref = "Control")
betadds_M <- DESeq(betadds_M)
betares_M <- results(betadds_M)
betaresOrdered_M <- betares_M[order(betares_M$pvalue),] #reorder by p value

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(betadds_M)
betaresLFCdisease_M <- lfcShrink(betadds_M, coef="simplified_diagnosis_T2D_vs_Control", type="apeglm")
betaresLFCdiseaseOrdered_M <- betaresLFCdisease_M[order(betaresLFCdisease_M$pvalue),] #reorder by p value

betaresLFCdiseaseOrdered_M <- data.frame(betaresLFCdiseaseOrdered_M)
betaresLFCdiseaseOrdered_M$Gene <- rownames(betaresLFCdiseaseOrdered_M)

#write csv file
write.csv(betaresLFCdiseaseOrdered_M, "Output/Fig3/Male ctrl vs T2D beta scRNAseq.csv")

#GSEA
betaresLFCdiseaseOrdered_M$entrez1 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered_M$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
betaresLFCdiseaseOrdered_M$entrez2 <- mapIds(org.Hs.eg.db, keys=c(betaresLFCdiseaseOrdered_M$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

betaresLFCdiseaseOrdered_M <- betaresLFCdiseaseOrdered_M %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = betaresLFCdiseaseOrdered_M,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "beta_scRNAseq_M_ctrlvsT2D")

#alpha-cell scRNAseq females ctrl vs T2D
donorlist <- colnames(alpha_bulk_rawcounts_df)
alphadonor_deseq_F <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis!="T1D") %>%
  filter(sex == "Female")
donorlist <- alphadonor_deseq_F$donor_ID #make the donorlist just the females, excluding donors with T1D
alphacounts_deseq_F <- alpha_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(alphadonor_deseq_F) <- alphadonor_deseq_F$donor_ID
all(rownames(alphadonor_deseq_F) %in% colnames(alphacounts_deseq_F)) #check the samples match

#RunDEseq
alphadds_F <- DESeqDataSetFromMatrix(countData = alphacounts_deseq_F,
                                    colData = alphadonor_deseq_F,
                                    design = ~ simplified_diagnosis + age_years) 

#remove low count genes
table(alphadonor_deseq_F$simplified_diagnosis)
smallestGroupSize <- 11
keep_F <- rowSums(counts(alphadds_F) >= 10) >= smallestGroupSize
alphadds_F <- alphadds_F[keep_F,]

#make non-diabetic the control for disease state for comparison
alphadds_F$simplified_diagnosis <- relevel(alphadds_F$simplified_diagnosis, ref = "Control")
alphadds_F <- DESeq(alphadds_F)
alphares_F <- results(alphadds_F)
alpharesOrdered_F <- alphares_F[order(alphares_F$pvalue),] #reorder by p value

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(alphadds_F)
alpharesLFCdisease_F <- lfcShrink(alphadds_F, coef="simplified_diagnosis_T2D_vs_Control", type="apeglm")
alpharesLFCdiseaseOrdered_F <- alpharesLFCdisease_F[order(alpharesLFCdisease_F$pvalue),] #reorder by p value

alpharesLFCdiseaseOrdered_F <- data.frame(alpharesLFCdiseaseOrdered_F)
alpharesLFCdiseaseOrdered_F$Gene <- rownames(alpharesLFCdiseaseOrdered_F)

#write csv file
write.csv(alpharesLFCdiseaseOrdered_F, "Output/Fig3/Female ctrl vs T2D alpha scRNAseq.csv")

#GSEA
alpharesLFCdiseaseOrdered_F$entrez1 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered_F$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
alpharesLFCdiseaseOrdered_F$entrez2 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered_F$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

alpharesLFCdiseaseOrdered_F <- alpharesLFCdiseaseOrdered_F %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = alpharesLFCdiseaseOrdered_F,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "alpha_scRNAseq_F_ctrlvsT2D")

#alpha-cell scRNAseq males ctrl vs T2D
donorlist <- colnames(alpha_bulk_rawcounts_df)
alphadonor_deseq_M <- donor_all %>%
  filter(donor_ID %in% donorlist) %>% #including only donors that we have scRNAseq data for
  filter(simplified_diagnosis!="T1D") %>%
  filter(sex == "Male")
donorlist <- alphadonor_deseq_M$donor_ID #make the donorlist just the females, excluding donors with T1D
alphacounts_deseq_M <- alpha_bulk_rawcounts_df %>%
  dplyr::select(all_of(donorlist))
rownames(alphadonor_deseq_M) <- alphadonor_deseq_M$donor_ID
all(rownames(alphadonor_deseq_M) %in% colnames(alphacounts_deseq_M)) #check the samples match

#RunDEseq
alphadds_M <- DESeqDataSetFromMatrix(countData = alphacounts_deseq_M,
                                     colData = alphadonor_deseq_M,
                                     design = ~ simplified_diagnosis + age_years) 

#remove low count genes
table(alphadonor_deseq_M$simplified_diagnosis)
smallestGroupSize <- 7
keep_M <- rowSums(counts(alphadds_M) >= 10) >= smallestGroupSize
alphadds_M <- alphadds_M[keep_M,]

#make non-diabetic the control for disease state for comparison
alphadds_M$simplified_diagnosis <- relevel(alphadds_M$simplified_diagnosis, ref = "Control")
alphadds_M <- DESeq(alphadds_M)
alphares_M <- results(alphadds_M)
alpharesOrdered_M <- alphares_M[order(alphares_M$pvalue),] #reorder by p value

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resultsNames(alphadds_M)
alpharesLFCdisease_M <- lfcShrink(alphadds_M, coef="simplified_diagnosis_T2D_vs_Control", type="apeglm")
alpharesLFCdiseaseOrdered_M <- alpharesLFCdisease_M[order(alpharesLFCdisease_M$pvalue),] #reorder by p value

alpharesLFCdiseaseOrdered_M <- data.frame(alpharesLFCdiseaseOrdered_M)
alpharesLFCdiseaseOrdered_M$Gene <- rownames(alpharesLFCdiseaseOrdered_M)

#write csv file
write.csv(alpharesLFCdiseaseOrdered_M, "Output/Fig3/Male ctrl vs T2D alpha scRNAseq.csv")

#GSEA
alpharesLFCdiseaseOrdered_M$entrez1 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered_M$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
alpharesLFCdiseaseOrdered_M$entrez2 <- mapIds(org.Hs.eg.db, keys=c(alpharesLFCdiseaseOrdered_M$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

alpharesLFCdiseaseOrdered_M <- alpharesLFCdiseaseOrdered_M %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = alpharesLFCdiseaseOrdered_M,
           FC_col = "log2FoldChange",
           p_col = "pvalue",
           name = "alpha_scRNAseq_M_ctrlvsT2D")

#Humanislets.com bulk RNAseq data to compare control vs T2D
#import data downloaded 21 May 2024
proc_rnaseq <- read.csv("data/Humanislets.com/proc_rnaseq.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")
proc_rnaseq_long <- proc_rnaseq %>%
  pivot_longer(!gene_id, values_to = "logCPM", names_to = "record_id")
proc_rnaseq_long <- inner_join(proc_rnaseq_long, donor_data, by = "record_id")

#Comparing non-diabetic with T2D for females
proc_rnaseq_long_FctrlvsT2D <- proc_rnaseq_long %>%
  filter(donorsex == "Female")
table(proc_rnaseq_long_FctrlvsT2D$diagnosis)/length(unique(proc_rnaseq_long_FctrlvsT2D$gene_id))
#36 Control, 5 T2D

#remove low count genes
smallest_group <- 5
proc_rnaseq_wide_FctrlvsT2D <- proc_rnaseq_long_FctrlvsT2D %>%
  pivot_wider(names_from = "gene_id", values_from = "logCPM")
proc_rnaseq_wide_FctrlvsT2D_na_counts <- proc_rnaseq_wide_FctrlvsT2D %>%
  summarise_all(~ sum(is.na(.)))
proc_rnaseq_wide_FctrlvsT2D_na_counts <- t(proc_rnaseq_wide_FctrlvsT2D_na_counts)
colnames(proc_rnaseq_wide_FctrlvsT2D_na_counts) <- "NAcount"
proc_rnaseq_wide_FctrlvsT2D_na_counts <- data.frame(proc_rnaseq_wide_FctrlvsT2D_na_counts)
proc_rnaseq_wide_FctrlvsT2D_na_counts <- proc_rnaseq_wide_FctrlvsT2D_na_counts %>%
  filter(NAcount<smallest_group)
proc_rnaseq_long_FctrlvsT2D <- proc_rnaseq_long_FctrlvsT2D %>%
  filter(gene_id %in% rownames(proc_rnaseq_wide_FctrlvsT2D_na_counts))

#multiple ANCOVAs
proc_rnaseq_ancovas_FctrlvsT2D <- list()
for (i in 1:length(unique(proc_rnaseq_long_FctrlvsT2D$gene_id))){ 
  geneid <- unique(proc_rnaseq_long_FctrlvsT2D$gene_id)[i]
  data <- proc_rnaseq_long_FctrlvsT2D %>%
    filter(gene_id == geneid)
  test.c <- data %>% #to avoid errors for genes with not enough values
    filter(diagnosis == "None") %>% 
    drop_na(logCPM)
  test.c.val <- dim(test.c)[1]
  test.d <- data %>%
    filter(diagnosis == "Type2") %>% 
    drop_na(logCPM)
  test.d.val <- dim(test.d)[1]
  if (test.c.val > 1 & test.d.val >1){
    ancova <- data %>%
      anova_test(logCPM ~ donorage + diagnosis)
    pval <- ancova[2,5]
    proc_rnaseq_ancovas_FctrlvsT2D[[i]] <- c(geneid, pval)
  } else{
    proc_rnaseq_ancovas_FctrlvsT2D[[i]] <- c(geneid, NA)
  }
}
proc_rnaseq_ancovas_FctrlvsT2D <- data.frame(do.call(rbind, proc_rnaseq_ancovas_FctrlvsT2D))
colnames(proc_rnaseq_ancovas_FctrlvsT2D) <- c("gene_id","pval")
proc_rnaseq_ancovas_FctrlvsT2D$gene_id <- as.character(proc_rnaseq_ancovas_FctrlvsT2D$gene_id)
proc_rnaseq_ancovas_FctrlvsT2D$pval <- as.numeric(proc_rnaseq_ancovas_FctrlvsT2D$pval)
proc_rnaseq_ancovas_FctrlvsT2D <- proc_rnaseq_ancovas_FctrlvsT2D %>%
  filter(is.na(pval) == FALSE) #remove all with no pval
head(proc_rnaseq_ancovas_FctrlvsT2D)

#Benjamini-Hochberg correction for multiple t tests
proc_rnaseq_ancovas_FctrlvsT2D$padj <- p.adjust(proc_rnaseq_ancovas_FctrlvsT2D$pval, method = "BH")

#get direction of change
proc_rnaseq_ancovas_FctrlvsT2D.means <- proc_rnaseq_long_FctrlvsT2D %>%
  group_by(gene_id, diagnosis) %>%
  summarise(mean.logCPM = mean(logCPM, na.rm=TRUE))
proc_rnaseq_ancovas_FctrlvsT2D.means <- proc_rnaseq_ancovas_FctrlvsT2D.means %>%
  pivot_wider(names_from = diagnosis, names_prefix = "meanlogCPM_", values_from = mean.logCPM)
proc_rnaseq_ancovas_FctrlvsT2D.means$gene_id <- as.character(proc_rnaseq_ancovas_FctrlvsT2D.means$gene_id)
proc_rnaseq_ancovas_FctrlvsT2D <- inner_join(proc_rnaseq_ancovas_FctrlvsT2D, proc_rnaseq_ancovas_FctrlvsT2D.means, by = "gene_id")
proc_rnaseq_ancovas_FctrlvsT2D <- proc_rnaseq_ancovas_FctrlvsT2D %>%
  mutate(logFC = meanlogCPM_None-meanlogCPM_Type2)
head(proc_rnaseq_ancovas_FctrlvsT2D)
proc_rnaseq_ancovas_FctrlvsT2D <- proc_rnaseq_ancovas_FctrlvsT2D %>%
  arrange(padj)

#get common names
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
common.names <- as.list(x[mapped_genes])
common.names <- unlist(common.names)

proc_rnaseq_ancovas_FctrlvsT2D$common.name <- rep(NA,nrow(proc_rnaseq_ancovas_FctrlvsT2D))
for (i in 1:length(proc_rnaseq_ancovas_FctrlvsT2D$gene_id)){
  entrez <- proc_rnaseq_ancovas_FctrlvsT2D$gene_id[i]
  proc_rnaseq_ancovas_FctrlvsT2D$common.name[i] <- common.names[entrez]
}

write.csv(proc_rnaseq_ancovas_FctrlvsT2D, "Humanisletscom bulkRNAseq F Ctrl vs F T2D.csv")

#GSEA
proc_rnaseq_ancovas_FctrlvsT2D$entrez1 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_FctrlvsT2D$common.name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
proc_rnaseq_ancovas_FctrlvsT2D$entrez2 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_FctrlvsT2D$common.name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

proc_rnaseq_ancovas_FctrlvsT2D <- proc_rnaseq_ancovas_FctrlvsT2D %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused")

gsea_human(df = proc_rnaseq_ancovas_FctrlvsT2D,
           FC_col = "logFC",
           p_col = "pval",
           name = "bulkRNAseq_FCtrlvsT2D")

#same comparison in males
proc_rnaseq_long_MctrlvsT2D <- proc_rnaseq_long %>%
  filter(donorsex == "Male")
table(proc_rnaseq_long_MctrlvsT2D$diagnosis)/length(unique(proc_rnaseq_long_MctrlvsT2D$gene_id))
#66 Control, 10 T2D

#remove low count genes
smallest_group <- 10
proc_rnaseq_wide_MctrlvsT2D <- proc_rnaseq_long_MctrlvsT2D %>%
  pivot_wider(names_from = "gene_id", values_from = "logCPM")
proc_rnaseq_wide_MctrlvsT2D_na_counts <- proc_rnaseq_wide_MctrlvsT2D %>%
  summarise_all(~ sum(is.na(.)))
proc_rnaseq_wide_MctrlvsT2D_na_counts <- t(proc_rnaseq_wide_MctrlvsT2D_na_counts)
colnames(proc_rnaseq_wide_MctrlvsT2D_na_counts) <- "NAcount"
proc_rnaseq_wide_MctrlvsT2D_na_counts <- data.frame(proc_rnaseq_wide_MctrlvsT2D_na_counts)
proc_rnaseq_wide_MctrlvsT2D_na_counts <- proc_rnaseq_wide_MctrlvsT2D_na_counts %>%
  filter(NAcount<smallest_group)
proc_rnaseq_long_MctrlvsT2D <- proc_rnaseq_long_MctrlvsT2D %>%
  filter(gene_id %in% rownames(proc_rnaseq_wide_MctrlvsT2D_na_counts))

#multiple ANCOVAs
proc_rnaseq_ancovas_MctrlvsT2D <- list()
for (i in 1:length(unique(proc_rnaseq_long_MctrlvsT2D$gene_id))){ 
  geneid <- unique(proc_rnaseq_long_MctrlvsT2D$gene_id)[i]
  data <- proc_rnaseq_long_MctrlvsT2D %>%
    filter(gene_id == geneid)
  test.c <- data %>% #to avoid errors for genes with not enough values
    filter(diagnosis == "None") %>% 
    drop_na(logCPM)
  test.c.val <- dim(test.c)[1]
  test.d <- data %>%
    filter(diagnosis == "Type2") %>% 
    drop_na(logCPM)
  test.d.val <- dim(test.d)[1]
  if (test.c.val > 1 & test.d.val >1){
    ancova <- data %>%
      anova_test(logCPM ~ donorage + diagnosis)
    pval <- ancova[2,5]
    proc_rnaseq_ancovas_MctrlvsT2D[[i]] <- c(geneid, pval)
  } else{
    proc_rnaseq_ancovas_MctrlvsT2D[[i]] <- c(geneid, NA)
  }
}
proc_rnaseq_ancovas_MctrlvsT2D <- data.frame(do.call(rbind, proc_rnaseq_ancovas_MctrlvsT2D))
colnames(proc_rnaseq_ancovas_MctrlvsT2D) <- c("gene_id","pval")
proc_rnaseq_ancovas_MctrlvsT2D$gene_id <- as.character(proc_rnaseq_ancovas_MctrlvsT2D$gene_id)
proc_rnaseq_ancovas_MctrlvsT2D$pval <- as.numeric(proc_rnaseq_ancovas_MctrlvsT2D$pval)
proc_rnaseq_ancovas_MctrlvsT2D <- proc_rnaseq_ancovas_MctrlvsT2D %>%
  filter(is.na(pval) == FALSE) #remove all with no pval
head(proc_rnaseq_ancovas_MctrlvsT2D)

#Benjamini-Hochberg correction for multiple t tests
proc_rnaseq_ancovas_MctrlvsT2D$padj <- p.adjust(proc_rnaseq_ancovas_MctrlvsT2D$pval, method = "BH")

#get direction of change
proc_rnaseq_ancovas_MctrlvsT2D.means <- proc_rnaseq_long_MctrlvsT2D %>%
  group_by(gene_id, diagnosis) %>%
  summarise(mean.logCPM = mean(logCPM, na.rm=TRUE))
proc_rnaseq_ancovas_MctrlvsT2D.means <- proc_rnaseq_ancovas_MctrlvsT2D.means %>%
  pivot_wider(names_from = diagnosis, names_prefix = "meanlogCPM_", values_from = mean.logCPM)
proc_rnaseq_ancovas_MctrlvsT2D.means$gene_id <- as.character(proc_rnaseq_ancovas_MctrlvsT2D.means$gene_id)
proc_rnaseq_ancovas_MctrlvsT2D <- inner_join(proc_rnaseq_ancovas_MctrlvsT2D, proc_rnaseq_ancovas_MctrlvsT2D.means, by = "gene_id")
proc_rnaseq_ancovas_MctrlvsT2D <- proc_rnaseq_ancovas_MctrlvsT2D %>%
  mutate(logFC = meanlogCPM_None-meanlogCPM_Type2)
head(proc_rnaseq_ancovas_MctrlvsT2D)
proc_rnaseq_ancovas_MctrlvsT2D <- proc_rnaseq_ancovas_MctrlvsT2D %>%
  arrange(padj)

#get common names
proc_rnaseq_ancovas_MctrlvsT2D$common.name <- rep(NA,nrow(proc_rnaseq_ancovas_MctrlvsT2D))
for (i in 1:length(proc_rnaseq_ancovas_MctrlvsT2D$gene_id)){
  entrez <- proc_rnaseq_ancovas_MctrlvsT2D$gene_id[i]
  proc_rnaseq_ancovas_MctrlvsT2D$common.name[i] <- common.names[entrez]
}

write.csv(proc_rnaseq_ancovas_MctrlvsT2D, "Humanisletscom bulkRNAseq M Ctrl vs F T2D.csv")

#GSEA
proc_rnaseq_ancovas_MctrlvsT2D$entrez1 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_MctrlvsT2D$common.name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
proc_rnaseq_ancovas_MctrlvsT2D$entrez2 <- mapIds(org.Hs.eg.db, keys=c(proc_rnaseq_ancovas_MctrlvsT2D$common.name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID

proc_rnaseq_ancovas_MctrlvsT2D <- proc_rnaseq_ancovas_MctrlvsT2D %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused")

gsea_human(df = proc_rnaseq_ancovas_MctrlvsT2D,
           FC_col = "logFC",
           p_col = "pval",
           name = "bulkRNAseq_MCtrlvsT2D")

#Make figures
#females
#Combining GSEA results to graph together
HIGO_F_CC <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_CC")
HIGO_F_BP <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_BP")
HIGO_F_MF <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_MF")
HIGO_F <- bind_rows(HIGO_F_CC, HIGO_F_BP, HIGO_F_MF)
HIGO_F$NES <- -HIGO_F$NES #make it so that negative means down with T2D
HIGO_F <- HIGO_F %>% filter(p.adjust < 0.05) #significant only

HPAPbetascGO_F_CC <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPbetascGO_F_BP <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPbetascGO_F_MF <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPbetascGO_F <- bind_rows(HPAPbetascGO_F_CC, HPAPbetascGO_F_BP, HPAPbetascGO_F_MF)
HPAPbetascGO_F <- HPAPbetascGO_F %>% filter(p.adjust < 0.05) #significant only

HPAPalphascGO_F_CC <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPalphascGO_F_BP <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPalphascGO_F_MF <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPalphascGO_F <- bind_rows(HPAPalphascGO_F_CC, HPAPalphascGO_F_BP, HPAPalphascGO_F_MF)
HPAPalphascGO_F <- HPAPalphascGO_F %>% filter(p.adjust < 0.05) #significant only

#Check that column names are consistent
all(colnames(HIGO_F) == colnames(HPAPbetascGO_F)) #TRUE
all(colnames(HIGO_F) == colnames(HPAPalphascGO_F)) #TRUE

#combine into single dataframe
GO_all_F <- bind_rows(HIGO_F, HPAPbetascGO_F, HPAPalphascGO_F)
GO_all_F$dataset <- c(rep("Humanislets.com\nbulk RNAseq",nrow(HIGO_F)), rep("HPAP beta-cell\nscRNAseq", nrow(HPAPbetascGO_F)), rep("HPAP alpha-cell\nscRNAseq", nrow(HPAPalphascGO_F)))
head(GO_all_F)

#split into control-biased and T2D-biased pathways
GO_F_Ctrlbiased <- GO_all_F %>% filter(NES < 0)
GO_F_T2Dbiased <- GO_all_F %>% filter(NES > 0)

#identify top 30 pathways
GO_F_Ctrlbiased <- GO_F_Ctrlbiased %>% arrange(p.adjust)
GO_F_Ctrlbiased_pathways <- unique(GO_F_Ctrlbiased$Description)
GO_F_Ctrlbiased_pathways <- GO_F_Ctrlbiased_pathways[1:30]

GO_F_T2Dbiased <- GO_F_T2Dbiased %>% arrange(p.adjust)
GO_F_T2Dbiased_pathways <- unique(GO_F_T2Dbiased$Description)
GO_F_T2Dbiased_pathways <- GO_F_T2Dbiased_pathways[1:30]

#Graph
#make long pathway name two lines
GO_F_Ctrlbiased <- GO_F_Ctrlbiased %>% 
  mutate(Description = case_when(
    Description == "oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor" ~ "oxidoreductase activity, acting on NAD(P)H,\nquinone or similar compound as acceptor",
    .default = Description
  ))
GO_F_Ctrlbiased_pathways[24] <- "oxidoreductase activity, acting on NAD(P)H,\nquinone or similar compound as acceptor"

GO_F_Ctrlbiased$NES <- -GO_F_Ctrlbiased$NES #flip sign
GO_F_Ctrlbiased %>%
  filter(Description %in% GO_F_Ctrlbiased_pathways) %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways decreased\nwith T2D in females", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#55066b", mid = "#C337EB", high = "white", midpoint = 2.5e-3, limits = c(0,5e-3)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,3), breaks = c(1.5,3)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig3/Fig3C.png", width = 9, height = 7)

#make long pathway name two lines
GO_F_T2Dbiased %>%
  filter(Description %in% GO_F_T2Dbiased_pathways) %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways increased\nwith T2D in females", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#992002", mid = "#E05A38", high = "white", midpoint = 0.025, limits = c(0,0.05)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig3/Fig3A.png", width = 9.5, height = 7)

#males
#Combining GSEA results to graph together
HIGO_M_CC <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_CC")
HIGO_M_BP <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_BP")
HIGO_M_MF <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_MF")
HIGO_M <- bind_rows(HIGO_M_CC, HIGO_M_BP, HIGO_M_MF)
HIGO_M$NES <- -HIGO_M$NES #make it so that negative means down with T2D
HIGO_M <- HIGO_M %>% filter(p.adjust < 0.05) #significant only

HPAPbetascGO_M_CC <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPbetascGO_M_BP <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPbetascGO_M_MF <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPbetascGO_M <- bind_rows(HPAPbetascGO_M_CC, HPAPbetascGO_M_BP, HPAPbetascGO_M_MF)
HPAPbetascGO_M <- HPAPbetascGO_M %>% filter(p.adjust < 0.05) #significant only

HPAPalphascGO_M_CC <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPalphascGO_M_BP <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPalphascGO_M_MF <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPalphascGO_M <- bind_rows(HPAPalphascGO_M_CC, HPAPalphascGO_M_BP, HPAPalphascGO_M_MF)
HPAPalphascGO_M <- HPAPalphascGO_M %>% filter(p.adjust < 0.05) #significant only

#Check that column names are consistent
all(colnames(HIGO_M) == colnames(HPAPbetascGO_M)) #TRUE
all(colnames(HIGO_M) == colnames(HPAPalphascGO_M)) #TRUE

#combine into single dataframe
GO_all_M <- bind_rows(HIGO_M, HPAPbetascGO_M, HPAPalphascGO_M)
GO_all_M$dataset <- c(rep("Humanislets.com\nbulk RNAseq",nrow(HIGO_M)), rep("HPAP beta-cell\nscRNAseq", nrow(HPAPbetascGO_M)), rep("HPAP alpha-cell\nscRNAseq", nrow(HPAPalphascGO_M)))
head(GO_all_M)

#split into control-biased and T2D-biased pathways
GO_M_Ctrlbiased <- GO_all_M %>% filter(NES < 0)
GO_M_T2Dbiased <- GO_all_M %>% filter(NES > 0)

#identify top 30 pathways
GO_M_Ctrlbiased <- GO_M_Ctrlbiased %>% arrange(p.adjust)
GO_M_Ctrlbiased_pathways <- unique(GO_M_Ctrlbiased$Description)
GO_M_Ctrlbiased_pathways <- GO_M_Ctrlbiased_pathways[1:30]

GO_M_T2Dbiased <- GO_M_T2Dbiased %>% arrange(p.adjust)
GO_M_T2Dbiased_pathways <- unique(GO_M_T2Dbiased$Description)
GO_M_T2Dbiased_pathways <- GO_M_T2Dbiased_pathways[1:30]

#Graph
GO_M_Ctrlbiased$NES <- -GO_M_Ctrlbiased$NES #flip sign
GO_M_Ctrlbiased %>%
  filter(Description %in% GO_M_Ctrlbiased_pathways) %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways decreased\nwith T2D in males", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#0a735b", mid = "#3EDBB8", high = "white", midpoint = 0.025, limits = c(0,0.05)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig3/Fig3D.png", width = 8.35, height = 7)

GO_M_T2Dbiased %>%
  filter(Description %in% GO_M_T2Dbiased_pathways) %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways increased\nwith T2D in males", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#9c6b05", mid = "#DBA126", high = "white", midpoint = 0.025, limits = c(0,0.05)) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.5), breaks = c(1,2.5)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig3/Fig3B.png", width = 8.15, height = 7)

### Figure 4 ### ----

#clear environment
rm(list = ls())

#import proteomics data from Humanislets.com (downloaded 21 May 2024)
prot <- read.csv("data/Humanislets.com/proc_prot.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")
isolation_data <- read.csv("data/Humanislets.com/isolation.csv")

#join donor info with proteomics data
donor_isolation <- inner_join(donor_data, isolation_data, by = "record_id")
prot_long <- pivot_longer(prot, !gene_id, values_to = "Abundance", names_to = "record_id")
prot_data <- inner_join(prot_long, donor_isolation, by = "record_id")

#get protein names from gene id
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
common.names <- as.list(x[mapped_genes])
common.names <- unlist(common.names)
map <- data.frame(names(common.names), common.names)
colnames(map) <- c("gene_id","common_name")
map$gene_id <- as.numeric(map$gene_id)
prot_data <- left_join(prot_data, map, by = "gene_id")

#Comparing control to T2D among females
prot_F <- prot_data %>%
  filter(donorsex == "Female")
table(prot_F$diagnosis)/length(unique(prot_F$gene_id)) #40 control, 7 T2D

#make separate tables for data and metadata
prot_F_feature <- prot_F %>%
  dplyr::select(gene_id, record_id, Abundance) %>%
  pivot_wider(names_from = "record_id", values_from = "Abundance")
genes_F <- prot_F_feature$gene_id
prot_F_feature <- data.frame(prot_F_feature)
rownames(prot_F_feature) <- genes_F
prot_F_feature <- prot_F_feature[,-1] #remove the gene column

prot_F_metadata <- prot_F %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(-Abundance) %>%
  dplyr::select(-common_name)
prot_F_metadata <- distinct(prot_F_metadata)
prot_F_metadata <- data.frame(prot_F_metadata)

#check order
all(colnames(prot_F_feature) == prot_F_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_F <- apply(prot_F_feature, 1, function(x){sum(!is.na(x)) >= 9})
prot_F_feature <- prot_F_feature[feature.keep_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(prot_F_metadata[,"diagnosis"]))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_F <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = prot_F_metadata)
colnames(design_F)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_F <- list()
ref <- "None" #using non-diabetic as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_F <- as.list(paste("Type2", "-", ref, sep = ""))
myargs_F[["levels"]] <- design_F
contrast.matrix_F <- do.call(makeContrasts, myargs_F)

# get results
fit_F <- lmFit(prot_F_feature, design_F, trend = TRUE, robust = TRUE)
fit_F <- contrasts.fit(fit_F, contrast.matrix_F)
fit_F <- eBayes(fit_F)
res.table_F <- topTable(fit_F, number = Inf)

# Remove results rows with NAs
res.table_F <- res.table_F[!is.na(res.table_F$P.Value), ]

# Process results for output
res.table_F <- merge(res.table_F, map, by.x = "row.names", by.y = "gene_id")
res.table_F <- res.table_F[,-7]

write.csv(res.table_F, "Output/Fig4/prot_dea_results_correctforage_F.csv", row.names = FALSE)

#GSEA
res.table_F$entrez1 <- mapIds(org.Hs.eg.db, keys=c(res.table_F$common_name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
res.table_F$entrez2 <- mapIds(org.Hs.eg.db, keys=c(res.table_F$common_name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID
res.table_F <- res.table_F %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human <- function(df,FC_col, p_col, name){ 
  # remove duplicated entrez
  df <- df[order(df[[p_col]]),]
  df <- df[!duplicated(df$entrez),]
  df <- df[!is.na(df$entrez), ]
  
  #rank by signed -log10 p value
  genelist <- sign(df[[FC_col]]) * (-log10(df[[p_col]]))
  
  # entrez id as names of the gene list
  names(genelist) <- df$entrez
  
  genelist[genelist==Inf] <- max(genelist[is.finite(genelist)])+1
  genelist[genelist==-Inf] <- min(genelist[is.finite(genelist)])-1
  
  genelist = sort(genelist, decreasing = TRUE)
  genelist <- na.omit(genelist)
  
  gseGO.bp <- gseGO(
    geneList=genelist,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <- as.data.frame(gseGO.bp)
  
  gseGO.mf <- gseGO(
    geneList=genelist,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <- as.data.frame(gseGO.mf)
  
  gseGO.cc <- gseGO(
    geneList=genelist,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <- as.data.frame(gseGO.cc)
  
  
  # Add the data frames to separate sheets
  wb <- createWorkbook()
  
  addWorksheet(wb, "GO_BP")
  writeData(wb, "GO_BP", gseGO.bp.df)
  
  addWorksheet(wb, "GO_MF")
  writeData(wb, "GO_MF", gseGO.mf.df)
  
  addWorksheet(wb, "GO_CC")
  writeData(wb, "GO_CC", gseGO.cc.df)
  
  # Save the workbook to a file
  fpath <- paste0("Output/Fig4/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

gsea_human(df = res.table_F,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "proteomics_F_controlvsT2D")

#Same comparison for males
prot_M <- prot_data %>%
  filter(donorsex == "Male")
table(prot_M$diagnosis)/length(unique(prot_M$gene_id)) #77 control, 10 T2D

#make separate tables for data and metadata
prot_M_feature <- prot_M %>%
  dplyr::select(gene_id, record_id, Abundance) %>%
  pivot_wider(names_from = "record_id", values_from = "Abundance")
genes_M <- prot_M_feature$gene_id
prot_M_feature <- data.frame(prot_M_feature)
rownames(prot_M_feature) <- genes_M
prot_M_feature <- prot_M_feature[,-1] #remove the gene column

prot_M_metadata <- prot_M %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(-Abundance) %>%
  dplyr::select(-common_name)
prot_M_metadata <- distinct(prot_M_metadata)
prot_M_metadata <- data.frame(prot_M_metadata)

#check order
all(colnames(prot_M_feature) == prot_M_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_M <- apply(prot_M_feature, 1, function(x){sum(!is.na(x)) >= 9})
prot_M_feature <- prot_M_feature[feature.keep_M, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(prot_M_metadata[,"diagnosis"]))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_M <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = prot_M_metadata)
colnames(design_M)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_M <- list()
ref <- "None" #using non-diabetic as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_M <- as.list(paste("Type2", "-", ref, sep = ""))
myargs_M[["levels"]] <- design_M
contrast.matrix_M <- do.call(makeContrasts, myargs_M)

# get results
fit_M <- lmFit(prot_M_feature, design_M, trend = TRUE, robust = TRUE)
fit_M <- contrasts.fit(fit_M, contrast.matrix_M)
fit_M <- eBayes(fit_M)
res.table_M <- topTable(fit_M, number = Inf)

# Remove results rows with NAs
res.table_M <- res.table_M[!is.na(res.table_M$P.Value), ]

# Process results for output
res.table_M <- merge(res.table_M, map, by.x = "row.names", by.y = "gene_id")
res.table_M <- res.table_M[,-7]

write.csv(res.table_M, "Output/Fig4/prot_dea_results_correctforage_M.csv", row.names = FALSE)

#GSEA
res.table_M$entrez1 <- mapIds(org.Hs.eg.db, keys=c(res.table_M$common_name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
res.table_M$entrez2 <- mapIds(org.Hs.eg.db, keys=c(res.table_M$common_name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID
res.table_M <- res.table_M %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

gsea_human(df = res.table_M,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "proteomics_M_controlvsT2D")

#Make graphs
GO_BP_F <- read_excel("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx", sheet = "GO_BP")
GO_BP_M <- read_excel("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx", sheet = "GO_BP")
GO_MF_F <- read_excel("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx", sheet = "GO_MF")
GO_MF_M <- read_excel("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx", sheet = "GO_MF")
GO_CC_F <- read_excel("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx", sheet = "GO_CC")
GO_CC_M <- read_excel("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx", sheet = "GO_CC")

GO_combined_F <- bind_rows(GO_BP_F, GO_MF_F, GO_CC_F)
GO_combined_M <- bind_rows(GO_BP_M, GO_MF_M, GO_CC_M)

GO_combined_F_upinT2D <- GO_combined_F %>% filter(NES > 0)
GO_combined_F_downinT2D <- GO_combined_F %>% filter(NES < 0)
GO_combined_M_upinT2D <- GO_combined_M %>% filter(NES > 0)
GO_combined_M_downinT2D <- GO_combined_M %>% filter(NES < 0)

#Top 30 pathways in females
GO_combined_F_upinT2D <- GO_combined_F_upinT2D %>%
  arrange(p.adjust)
GO_combined_F_upinT2D_select <- GO_combined_F_upinT2D[c(1:30),]

summary(GO_combined_F_upinT2D_select$NES)
summary(GO_combined_F_upinT2D_select$p.adjust)
GO_combined_F_upinT2D_select %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "GO pathways increased\n with T2D in females", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#992002", mid = "#E05A38", high = "white", midpoint = 0.00001, limits = c(0,0.00002)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_continuous(limits = c(1.5,2.5))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig4/Fig4A.png", width = 6, height = 6)

GO_combined_F_downinT2D <- GO_combined_F_downinT2D %>%
  arrange(p.adjust)
GO_combined_F_downinT2D_select <- GO_combined_F_downinT2D[c(1:30),]

GO_combined_F_downinT2D_select$NES <- -GO_combined_F_downinT2D_select$NES
summary(GO_combined_F_downinT2D_select$NES)
summary(GO_combined_F_downinT2D_select$p.adjust)
GO_combined_F_downinT2D_select %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "GO pathways decreased\n with T2D in females", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#55066b", mid = "#C337EB", high = "white", midpoint = 0.00002, limits = c(0,0.00004)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.6), breaks = c(1.5,2.6)) +
  scale_x_continuous(limits = c(1.5,2.6))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig4/Fig4C.png", width = 6.25, height = 6)

#Top 30 pathways in males
GO_combined_M_upinT2D <- GO_combined_M_upinT2D %>%
  arrange(p.adjust)
GO_combined_M_upinT2D_select <- GO_combined_M_upinT2D[c(1:30),]

summary(GO_combined_M_upinT2D_select$NES)
summary(GO_combined_M_upinT2D_select$p.adjust)
GO_combined_M_upinT2D_select %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "GO pathways increased\n with T2D in males", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#9c6b05", mid = "#DBA126", high = "white", midpoint = 0.000015, limits = c(0,0.00003)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,3.5), breaks = c(1.5,3.5)) +
  scale_x_continuous(limits = c(1.5,3.5))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig4/Fig4B.png", width = 6, height = 6)

GO_combined_M_downinT2D <- GO_combined_M_downinT2D %>%
  arrange(p.adjust)
GO_combined_M_downinT2D_select <- GO_combined_M_downinT2D[c(1:30),]

GO_combined_M_downinT2D_select$NES <- -GO_combined_M_downinT2D_select$NES
summary(GO_combined_M_downinT2D_select$NES)
summary(GO_combined_M_downinT2D_select$p.adjust)
GO_combined_M_downinT2D_select %>%
  ggplot(aes(x=NES, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  theme_bw() +
  labs(x = "Normalised Enrichment Score (NES)", y = "", 
       title = "GO pathways decreased\n with T2D in males", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#0a735b", mid = "#3EDBB8", high = "white", midpoint = 0.000002, limits = c(0,0.000004)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_continuous(limits = c(1.5,2.5))+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Fig4/Fig4D.png", width = 6.25, height = 6)

### Figure 5 and S3 ### ----

#clear environment
rm(list = ls())

#UPenn data
#Read in metadata
donors <- read_excel("data/HPAP/Donor_Summary_192.xlsx")

#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/UPenn Perifusion/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}

donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #160

#Read in the data
perifusion_data <- list()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  data <- read.csv(paste0("data/HPAP/UPenn Perifusion/HPAP-",donor,"/Islet Studies/Islet physiology studies/Perifusions/HPAP-",donor,"_Perifusion_data.csv"))
  perifusion_data[[i]] <- data
}

#Do they all have the same colnames?
colnames <- list()
for (i in 1:length(donor_numbers)){
  colnames[[i]] <- colnames(perifusion_data[[i]])
}
colnamesdf <- data.frame(do.call(rbind, colnames))
head(colnamesdf)
colnames(colnamesdf) <- c("C1","C2","C3","C4","C5","C6","C7")
unique(colnamesdf$C1)
unique(colnamesdf$C2)
unique(colnamesdf$C3)
unique(colnamesdf$C4)
unique(colnamesdf$C5)
unique(colnamesdf$C6)
unique(colnamesdf$C7) #Some different versions but at least they all indicate the same thing

#Combine data from different donors into one df
#fix all the colnames so they are the same
for (i in 1:length(donor_numbers)){
  colnames(perifusion_data[[i]]) <- colnames(perifusion_data[[1]])
}

perifusion_all <- do.call(rbind, perifusion_data)
perifusion_all <- data.frame(perifusion_all)
head(perifusion_all)
#add donor_ids
donor_id_col <- c()
for (i in 1:length(donor_numbers)){
  to_add <- rep(paste0("HPAP-",donor_numbers[i]),nrow(perifusion_data[[i]]))
  donor_id_col <- c(donor_id_col, to_add)
}
perifusion_all$DonorID <- donor_id_col

#change to numeric (units written in the column causing variables to be read as character)
perifusion_all$Insulin.Content...ng.islet. <- as.numeric(perifusion_all$Insulin.Content...ng.islet.)

#add metadata
donors <- donors %>%
  mutate(
    simplified_diagnosis = case_when(
      grepl("T2DM", clinical_diagnosis) ~ "T2D",
      grepl("T1DM", clinical_diagnosis) ~ "T1D",
      grepl("control", clinical_diagnosis) ~ "Control",
      .default = clinical_diagnosis
    ))
perifusion_all_with_metadata <- perifusion_all %>%
  left_join(donors, by= c("DonorID" = "donor_ID"))

#refine to focus on insulin content
hormone.content <- perifusion_all_with_metadata %>%
  dplyr::select(-Time..min.) %>%
  dplyr::select(-Insulin.Release..ng.100.islets.min.) %>%
  dplyr::select(-Glucagon.Release..pg.100.islets.min.) %>%
  dplyr::select(-Glucagon.Content..pg.islet.) %>%
  dplyr::select(-Stimulus) %>%
  filter(is.na(Insulin.Content...ng.islet.) == FALSE) %>%
  distinct()

hist(hormone.content$Insulin.Content...ng.islet.) #left-skewed
hist(log(hormone.content$Insulin.Content...ng.islet.)) #right-skewed

## Fig 5A ## ----
hormone.content %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(Insulin.Content...ng.islet. ~ age_years + sex) #p=0.068 for sex
hormone.content %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=Insulin.Content...ng.islet., fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 80, label = "p=0.068", label.size = 4, size = 0.5, tip.length = c(0.02, 0.45))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("Insulin content (ng/islet)") +
  xlab("") +
  labs(colour = "Condition", fill = "Condition") +
  ylim(0,90)+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()     
  )
ggsave("Output/Fig5/Fig5C.png", width = 3, height = 3)

## Fig 5D ## ----
hormone.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(Insulin.Content...ng.islet. ~ age_years + sex*simplified_diagnosis) #sig for sex but not disease
model <- hormone.content %>% filter(simplified_diagnosis != "T1D") %>%
  lm(Insulin.Content...ng.islet. ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0362 for F-M controls
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0475 for ctrl vs T2D in females
hormone.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=Insulin.Content...ng.islet.))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 78, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 90, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Insulin content (ng/islet)") +
  xlab("") +
  labs(colour = "Condition", fill = "Condition") +
  theme_bw() +
  theme(panel.grid = element_blank())+
  ylim(0,100)+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()      
  )
ggsave("Output/Fig5/Fig5F.png", width = 5, height = 4)

#save the hormone content data frame for later
UPenn.insulin.content <- write.csv(hormone.content, "Output/Fig5/UPenn insulin content.csv")

#Vanderbilt data
#clear environment
rm(list = ls())

#Read in metadata
donors <- read_excel("data/HPAP/Donor_Summary_192.xlsx")

#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/Vanderbilt Perifusion/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}

donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #97

#Read in insulin content data for all donors
insulin.content <- list()
for (i in 1:length(donor_numbers)){
  donor <- paste0("HPAP-",donor_numbers[i])
  if (file.exists(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
    if ("Insulin Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
      postperi <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx"), sheet = "Insulin Results", range = "D54:F55")
      colnames(postperi) <- c("Insulin.ng.mL","Insulin.content.ng","Insulin.content.ng.postperiIEQ")
      postperi$DonorID <- rep(donor,nrow(postperi))} else{
        postperi <- NA
      }
  } else{
    if ("Insulin Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx")) == TRUE){
      postperi <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx"), sheet = "Insulin Results", range = "D54:F55")
      colnames(postperi) <- c("Insulin.ng.mL","Insulin.content.ng","Insulin.content.ng.postperiIEQ")
        postperi$DonorID <- rep(donor,nrow(postperi))} else{
        postperi <- NA
      }
  }
  insulin.content[[i]] <- postperi
}

which(is.na(insulin.content)) #18 and 30
donor_numbers[which(is.na(insulin.content))] #HPAP-032 and HPAP-055 don't have insulin data
#remove from the list
insulin.content <- insulin.content [-18]
insulin.content <- insulin.content [-29]

#make dataframe
insulin.content_df <- data.frame(do.call(rbind, insulin.content))
head(insulin.content_df)

#add metadata
donors <- donors %>%
  mutate(
    simplified_diagnosis = case_when(
      grepl("T2DM", clinical_diagnosis) ~ "T2D",
      grepl("T1DM", clinical_diagnosis) ~ "T1D",
      grepl("control", clinical_diagnosis) ~ "Control",
      .default = clinical_diagnosis
    ))

insulin.content_df <- left_join(insulin.content_df, donors, by = c("DonorID" = "donor_ID"))
head(insulin.content_df)

hist(insulin.content_df$Insulin.content.ng.postperiIEQ)

## Fig 5B ## ----
insulin.content_df %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(Insulin.content.ng.postperiIEQ ~ age_years + sex) #ns
insulin.content_df %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=Insulin.content.ng.postperiIEQ, fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("Insulin content (ng/IEQ)") +
  xlab("") +
  labs(colour = "Condition", fill = "Condition") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()     
  )
ggsave("Output/Fig5/Fig5B.png", width = 3, height = 3)

## Fig 5E ## ----
insulin.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(Insulin.content.ng.postperiIEQ ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- insulin.content_df %>% filter(simplified_diagnosis != "T1D") %>%
  lm(Insulin.content.ng.postperiIEQ ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0073 for F-M T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0502 for Ctrl-T2D in females
insulin.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=Insulin.content.ng.postperiIEQ))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 20, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 21, label = "p=0.05", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Insulin content (ng/IEQ)") +
  xlab("") +
  labs(colour = "Condition", fill = "Condition") +
  theme_bw() +
  theme(panel.grid = element_blank())+
  ylim(0,25)+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),  
        axis.line = element_line()       
  )
ggsave("Output/Fig5/Fig5E.png", width = 5, height = 4)

## Fig S3A ## ----
UPenn.insulin.content <- read.csv("Output/Fig5/UPenn insulin content.csv")
UPenn_Vand_inscontent <- insulin.content_df %>%
  inner_join(UPenn.insulin.content, by = "DonorID")
dim(UPenn_Vand_inscontent) #39 donors

#make joint sex + diabetes status variable
UPenn_Vand_inscontent$sex_disease <- paste0(UPenn_Vand_inscontent$sex, " ", UPenn_Vand_inscontent$simplified_diagnosis)

UPenn_Vand_inscontent %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=Insulin.Content...ng.islet., y = Insulin.content.ng.postperiIEQ))+
  geom_point(aes(colour = sex_disease))+
  xlab("Insulin content (ng/islet) measured at UPenn")+
  ylab("Insulin content (ng/IEQ) measured at Vanderbilt")+
  labs(colour = "Group")+
  scale_colour_manual(values = c("#50b5ad","#0f615a","#dbac5a","#b57404"))+
  stat_smooth(method = "lm", colour = "grey50") +
  stat_cor(label.y = 19)+
  theme_classic()
ggsave("Output/FigS3/FigS3A.png", width = 6, height = 4)


#Insulin content from Humanislets.com

#clear environment
rm(list = ls())

#load data
donor_data <- read.csv("data/Humanislets.com/donor.csv")
isolation_data <- read.csv("data/Humanislets.com/isolation.csv")
donor_isolation <- inner_join(donor_data, isolation_data, by = "record_id")

hist(donor_isolation$insulinperieq)
hist(log(donor_isolation$insulinperieq))

## Fig 5C ## ----
donor_isolation %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(insulinperieq ~ donorage + donorsex) #ns
donor_isolation %>%
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "None" ~ "Control",
    .default = diagnosis
  )) %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(insulinperieq), fill= donorsex))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()    
  )
ggsave("Output/Fig5/Fig5A.png", width = 3, height = 3)

## Fig S3B ## ----
donor_isolation %>%
  filter(diagnosis != "Type1") %>%
  anova_test(pdisletparticleindex ~ donorage + donorsex*diagnosis) #ns
model <- donor_isolation %>%   filter(diagnosis != "Type1") %>%
  lm(pdisletparticleindex ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

donor_isolation %>%
  filter(diagnosis != "Type1") %>%
  mutate(diagnosis = case_when(
    diagnosis=="None" ~ "Control",
    diagnosis == "Type2" ~ "T2D",
    .default = diagnosis
  )) %>%
  ggplot(aes(x=donorsex, y=pdisletparticleindex, fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("Islet Particle Index") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS3/FigS3B.png", width = 5, height = 4)


## Fig 5F ## ----
donor_isolation %>%
  filter(diagnosis != "Type1") %>%
  mutate(output = log(insulinperieq)) %>%
  mutate(output = case_when(output == -Inf ~ NA, #remove undetectable results that become -Inf when log-transformed
                            .default = output)) %>%
  filter(is.na(output) == FALSE) %>%
  anova_test(output ~ donorage + donorsex*diagnosis) #p=0.055 for diagnosis
model <- donor_isolation %>% filter(diagnosis != "Type1") %>%
  mutate(output = log(insulinperieq)) %>%
  mutate(output = case_when(output == -Inf ~ NA,
                            .default = output)) %>%
  filter(is.na(output) == FALSE) %>%
  lm(output ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0411 for Ctrl-T2D in females

donor_isolation %>%
  filter(diagnosis != "Type1") %>%
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "None" ~ "Control",
    .default = diagnosis
  )) %>%
  ggplot(aes(x=donorsex, y=log(insulinperieq), fill= diagnosis))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  scale_y_continuous(limits = c(-2.5,7))+ 
  labs(colour = "Condition", fill = "Condition") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Fig5/Fig5D.png", width = 5, height = 4)

### Figure 6 ### ----

#clear environment
rm(list = ls())

#HPAP OCR data
#Read in metadata
donors <- read_excel("data/HPAP/Donor_Summary_192.xlsx")

#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/OCR/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}
donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #96

#Read in the data
all_data <- list()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  data <- read_excel(paste0("data/HPAP/OCR/HPAP-",donor,"/Islet Studies/Islet physiology studies/Oxygen consumption/HPAP-",donor,"_Oxygen-consumption_data.xlsx"))
  all_data[[i]] <- data
}

#Do they all have the same column names?
columnnames <- list()
for (i in 1:length(donor_numbers)){
  columnnames[[i]] <- colnames(all_data[[i]])
}
columnnames <- data.frame(do.call(rbind, columnnames))
head(columnnames)
colnames(columnnames) <- c("C1","C2","C3","C4")
unique(columnnames$C1)
unique(columnnames$C2)
unique(columnnames$C3)
unique(columnnames$C4)
#Yes, all donors have the same column names except for Outflow, where one or more has a slightly different name but indicates same thing

#Fix outflow colnames and non-numeric values (e.g., units written in columns)
for (i in 1:length(donor_numbers)){
  data <- all_data[[i]]
  colnames(data) <- c("Time.min","Inflow","Outflow","OCR.nmol.min.100islets")
  data <- data %>% mutate_all(as.numeric)
  df_filtered <- data %>%
    filter(rowSums(is.na(.)) != ncol(.))
  all_data[[i]] <- df_filtered
}

#Add donorIDs
for (i in 1:length(donor_numbers)){
  data <- all_data[[i]]
  donor <- paste0("HPAP-",donor_numbers[i])
  data$DonorID <- rep(donor,nrow(data))
  all_data[[i]] <- data
}

OCR_all <- do.call(bind_rows, all_data)
OCR_all <- data.frame(OCR_all)
summary(OCR_all)

#Summarise per stimulus
lastbaseline <- OCR_all %>%
  filter(Time.min < 50, Time.min > 49) %>%
  group_by(DonorID) %>%
  summarise(lastbasal = mean(OCR.nmol.min.100islets))
OCR_all <- inner_join(OCR_all, lastbaseline, by = "DonorID")
AAM <- OCR_all %>%
  filter(Time.min > 50, Time.min < 80) %>%
  group_by(DonorID) %>%
  summarise(AAM = mean(OCR.nmol.min.100islets))
OCR_all <- inner_join(OCR_all, AAM, by = "DonorID")
AAM.3mMGlc <- OCR_all %>%
  filter(Time.min > 80, Time.min < 105) %>%
  group_by(DonorID) %>%
  summarise(AAM.3mMGlc = mean(OCR.nmol.min.100islets))
OCR_all <- inner_join(OCR_all, AAM.3mMGlc, by = "DonorID")
lastLG <- OCR_all %>%
  filter(Time.min > 104, Time.min < 105) %>%
  group_by(DonorID) %>%
  summarise(lastLG = mean(OCR.nmol.min.100islets)) #better physiological "baseline"
OCR_all <- inner_join(OCR_all, lastLG, by = "DonorID")
HighGlc <- OCR_all %>%
  filter(Time.min > 105, Time.min < 120) %>% #prevents capturing rise from HG to FCCP
  group_by(DonorID) %>%
  summarise(HighGlc = mean(OCR.nmol.min.100islets), maxHG = max(OCR.nmol.min.100islets))
OCR_all <- inner_join(OCR_all, HighGlc, by = "DonorID")
FCCP <- OCR_all %>%
  filter(Time.min > 140, Time.min < 165) %>% #140 instead of 135 avoids capturing the rise from HG to FCCP
  group_by(DonorID) %>%
  summarise(FCCP = mean(OCR.nmol.min.100islets), maxFCCP = max(OCR.nmol.min.100islets))
OCR_all <- inner_join(OCR_all, FCCP, by = "DonorID")
NaN3 <- OCR_all %>%
  filter(Time.min > 170, Time.min < 200) %>% #170 instead of 165 avoids capturing the drop between FCCP and NaN3
  group_by(DonorID) %>%
  summarise(NaN3 = mean(OCR.nmol.min.100islets),min.NaN3 = min(OCR.nmol.min.100islets, na.rm = TRUE))
OCR_all <- inner_join(OCR_all, NaN3, by = "DonorID")
View(OCR_all)

OCR_all <- OCR_all %>%
  left_join(donors, by = c("DonorID" = "donor_ID")) %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

OCR_summary <- OCR_all %>%
  dplyr::select(-Inflow) %>%
  dplyr::select(-Outflow) %>%
  dplyr::select(-OCR.nmol.min.100islets) %>%
  dplyr::select(-Time.min) %>%
  distinct()

#Make version of summary with calculations for respiratory parameters
OCR_summary_v2 <- OCR_summary %>%
  mutate(spare.resp = FCCP - HighGlc, response.highglc = HighGlc - AAM.3mMGlc, lastbasal = lastbasal, max.resp = FCCP - min.NaN3, response.lowglc = AAM.3mMGlc - AAM, response.AAM = AAM - lastbasal) %>%
  mutate(spare.resp.max = maxFCCP - lastbasal, spare.resp.max.fcbaseline = (maxFCCP - lastbasal)/lastbasal)

#Add metadata
OCR_summary_v2 <- OCR_summary_v2 %>%
  left_join(donors, by = c("DonorID" = "donor_ID"))
OCR_summary_v2 <- OCR_summary_v2 %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

## Fig 6A ## ----
#summarise across time
OCR_all_over_time_summary <- OCR_all %>%
  filter(simplified_diagnosis != "T1D") %>%
  mutate(Time.min = round(Time.min)) %>%
  group_by(DonorID, age_years, Time.min, sex, simplified_diagnosis) %>%
  summarise(
    n = n(),
    donormeanraw = mean(OCR.nmol.min.100islets, na.rm = TRUE),
    .groups = "drop"
  )

#make graph
stimulus_data <- data.frame(
  xmin = c(0,50,80,105,135,165),
  xmax = c(50,200,105,200,200,200),
  label = c("Baseline","AAM","G 3
            ","G 16.7","FCCP","NaN3"),
  ymin = c(0.55,0.55,0.60,0.60,0.65,0.70), 
  ymax = c(0.60,0.60,0.65,0.65,0.70,0.75),
  ylab = c(0.575, 0.575, 0.59, 0.625, 0.675, 0.725)
)
OCR_all_over_time_summary %>%
  filter(Time.min < 200) %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  mutate(Time.min = round(Time.min)) %>%
  ungroup() %>%
  group_by(Time.min, sex) %>%
  summarise(
    n = n(),
    mean = mean(donormeanraw, na.rm = TRUE),
    sd = sd(donormeanraw, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>% 
  ggplot(aes(x=Time.min, y=mean, colour = sex, group = sex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("OCR (nmol/min/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Sex", fill = "Sex") +
  coord_cartesian(xlim = c(0, 210), ylim = c(0,0.75)) + 
  scale_x_continuous(breaks = c(50, 80, 105, 135, 165))+
  geom_vline(xintercept = c(50, 80, 105, 135, 165), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "black",linewidth = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = ylab, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig6/Fig6A.png", width = 5.5, height = 4)

## Fig 6I ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,50,80,105,135,165),   
  xmax = c(50,200,105,200,200,200),
  label = c("Baseline","AAM","G 3","G 16.7","FCCP","NaN3"),
  ymin = c(70,70,77,77,84,91), 
  ymax = c(77,77,84,84,91,98)
)

OCR_all_over_time_summary %>%
  filter(Time.min < 193) %>%
  mutate(Time.min = round(Time.min)) %>%
  ungroup() %>%
  group_by(Time.min, simplified_diagnosis, sex) %>%
  summarise(
    n = n(),
    mean = mean(donormeanraw, na.rm = TRUE),
    sd = sd(donormeanraw, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab ("OCR (nmol/min/100 islets))") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  coord_cartesian(xlim = c(0, 200), ylim = c(0,1)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(50, 80, 105, 135, 165), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin/100, ymax = ymax/100),
    inherit.aes = FALSE,
    fill = NA, colour = "black",size = 0.25
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin/100+ymax/100)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig6/Fig6I.png", width = 5.5, height = 4)

## Fig 6B ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(lastbasal ~ age_years + sex) #ns
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=lastbasal, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("") +
  ylab("Baseline\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6B.png", width = 3, height = 3)

## Fig 6J ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(lastbasal ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(lastbasal ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0231 for F-M T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0068 for Ctrl-T2D in males

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=lastbasal)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.75, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 0.86, label = "*", label.size = 7, size = 0.5, tip.length = c(0.3, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Baseline\n(nmol/min/100 islets)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6J.png", width = 4.5, height = 4)

## Fig 6C ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(maxFCCP ~ age_years + sex) #ns
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=maxFCCP, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("") +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6C.png", width = 3, height = 3)

## Fig 6K ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(spare.resp.max.fcbaseline ~ age_years + sex*simplified_diagnosis) #ns
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(spare.resp.max.fcbaseline ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=ns

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fcbaseline)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6K.png", width = 4.5, height = 4)

## Fig 6D ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(spare.resp.max.fcbaseline ~ age_years + sex) #p=0.055 for sex
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fcbaseline*100, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("") +
  ylab("Spare respiratory capacity\n(% of baseline)") +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 370, label = "p=0.055", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  ylim(0,400) +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6D.png", width = 3, height = 3)

## Fig 6L ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(spare.resp.max.fcbaseline ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(spare.resp.max.fcbaseline ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0466 for F-M in Control, p=0.0403 for F-M in T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0183 for Ctrl-T2D in males

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fcbaseline*100)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 375, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 430, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 360, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.6))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Spare respiratory capacity\n(% of baseline)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,450))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6L.png", width = 4.5, height = 4)

## Fig 6M ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(response.highglc ~ age_years + sex*simplified_diagnosis) #sig for disease and sex*disease interaction
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(response.highglc ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0551 for F-M in Control, p=0.0639 for F-M in T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p<0.0001 for Ctrl-T2D in males

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=response.highglc)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.12, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.4))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("High glucose-stimulated\nOCR (nmol/min/100 islets)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-0.1, 0.15)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6M.png", width = 4.5, height = 4)
  
#Humanislets OCR data
#clear environment
rm(list = ls())

#read metadata and data
donor_data <- read.csv("data/Humanislets.com/donor.csv")
seahorse_norm_dna <- read.csv("data/Humanislets.com/seahorse_norm_dna.csv")
seahorse_norm_dna_baseline <- read.csv("data/Humanislets.com/seahorse_norm_dna_baselineoc.csv")

#first working with baselined data
#remove outlier runs (three runs per donor, check if one of those three should be excluded)
all_time <- seahorse_norm_dna %>%
  pivot_longer(colnames(seahorse_norm_dna)[grepl("^time",colnames(seahorse_norm_dna))], names_to = "time", values_to = "OCR" ) %>%
  group_by(record_id, replicate) %>%
  summarise(mean_all_time = mean(OCR, na.rm = TRUE))

all_time_outliers <- list()
for (i in 1:length(unique(all_time$record_id))){
  donor <- unique(all_time$record_id)[i]
  data <- all_time %>% filter(record_id == donor)
  
  if (nrow(data) >= 3 && all(!is.na(data$mean_all_time))) {
    test_result <- grubbs.test(data$mean_all_time)
    pval <- test_result$p.value
    outlier <- parse_number(test_result$alternative)
  } else {
    pval <- NA
    outlier <- NA
  }
  
  all_time_outliers[[i]] <- c(donor, pval, outlier)
}
all_time_outliers <- do.call(rbind, all_time_outliers)
all_time_outliers <- data.frame(all_time_outliers)
colnames(all_time_outliers) <- c("record_id","pval","outlier")
head(all_time_outliers)
summary(all_time_outliers)
all_time_outliers$outlier <- as.numeric(all_time_outliers$outlier)
all_time_outliers <- all_time_outliers %>%
  filter(pval < 0.05) %>%
  inner_join(all_time, by = c("record_id","outlier" = "mean_all_time")) 
all_time_outliers #6 outliers

seahorse_norm_dna_postgrubbs <- seahorse_norm_dna %>%
  anti_join(all_time_outliers, by = c("record_id", "replicate"))

#add metadata
seahorse_norm_dna_postgrubbs <- left_join(seahorse_norm_dna_postgrubbs, donor_data, by = "record_id")

#summarise per donor
seahorse_norm_dna_summary <- seahorse_norm_dna_postgrubbs %>%
  group_by(record_id, donorage, donorsex, diagnosis) %>%
  summarise(across(where(is.numeric), mean))

## Fig 6E ## ----
#Import raw OCR data (not)
#Make long format
colnames(seahorse_norm_dna_summary)
seahorse_norm_dna_summary_long <- seahorse_norm_dna_summary %>% 
  pivot_longer(cols = 17:48, values_to = "OCR", names_to = "Time.min")
seahorse_norm_dna_summary_long$Time.min <- gsub("time_","",seahorse_norm_dna_summary_long$Time.min)
seahorse_norm_dna_summary_long$Time.min <- as.numeric(seahorse_norm_dna_summary_long$Time.min)

stimulus_data3 <- data.frame(
  xmin = c(0,35,85,145,205), 
  xmax = c(35,85,145,205,260), 
  label = c("G 2.8","G 16.7","Oligomycin","FCCP","Rotenone/\nAntA"),
  ymin = c(rep(300,4),285),  
  ymax = c(rep(330,4),335)
)
seahorse_norm_dna_summary_long %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ungroup() %>%
  group_by(Time.min, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(OCR, na.rm = TRUE),
    sd = sd(OCR, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("OCR (pmol/min/\U03BCg DNA)") +
  xlab ("Time (min)") +
  labs(colour = "Sex", fill = "Sex") +
  coord_cartesian(xlim = c(0, 270), ylim = c(0,330)) +  
  scale_x_continuous(breaks = c(35,85,145,205))+
  geom_vline(xintercept = c(35,85,145,205), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "black",size = 0.5
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "top") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig6/Fig6E.png", width = 5.5, height = 4)

## Fig 6F ## ----
seahorse_norm_dna_summary %>%
  ungroup() %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(calc_basal_resp ~ donorage + donorsex) #ns
seahorse_norm_dna_summary %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ggplot(aes(x=donorsex, y=calc_basal_resp, fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("") +
  ylab("Baseline\n(pmol/min/\U03BCg DNA)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()    
  )
ggsave("Output/Fig6/Fig6F.png", width = 3, height = 3)

## Fig 6G ## ----
seahorse_norm_dna_summary %>%
  ungroup() %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(calc_max_resp ~ donorage + donorsex) #ns
seahorse_norm_dna_summary %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ggplot(aes(x=donorsex, y=calc_max_resp, fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("") +
  ylab("Maximal respiratory capacity\n(pmol/min/\U03BCg DNA)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line()  
  )
ggsave("Output/Fig6/Fig6G.png", width = 3, height = 3)

## Fig 6H ## ----
#Use baselined data and process in the  same way
all_time_baselined <- seahorse_norm_dna_baseline %>%
  pivot_longer(colnames(seahorse_norm_dna_baseline)[grepl("^time",colnames(seahorse_norm_dna_baseline))], names_to = "time", values_to = "baselined_OCR" ) %>%
  group_by(record_id, replicate) %>%
  summarise(mean_all_time = mean(baselined_OCR, na.rm = TRUE))

all_time_outliers_baselined <- list()
for (i in 1:length(unique(all_time_baselined$record_id))){
  donor <- unique(all_time_baselined$record_id)[i]
  data <- all_time_baselined %>%
    filter(record_id == donor)
  pval <- (grubbs.test(data$mean_all_time)[3])$p.value
  outlier <- grubbs.test(data$mean_all_time)[2]
  outlier <- parse_number(outlier$alternative)
  all_time_outliers_baselined[[i]] <- c(donor, pval, outlier)
}
all_time_outliers_baselined <- do.call(rbind, all_time_outliers_baselined)
all_time_outliers_baselined <- data.frame(all_time_outliers_baselined)
colnames(all_time_outliers_baselined) <- c("record_id","pval","outlier")
head(all_time_outliers_baselined)

all_time_outliers_baselined$outlier <- as.numeric(all_time_outliers_baselined$outlier)
all_time_outliers_baselined <- all_time_outliers_baselined %>%
  filter(pval < 0.05) %>%
  inner_join(all_time_baselined, by = c("record_id","outlier" = "mean_all_time")) 
all_time_outliers_baselined #14 outliers

seahorse_norm_dna_baseline_postgrubbs <- seahorse_norm_dna_baseline %>%
  anti_join(all_time_outliers_baselined, by = c("record_id", "replicate"))

seahorse_norm_dna_baseline_postgrubbs <- left_join(seahorse_norm_dna_baseline_postgrubbs, donor_data, by = "record_id")

seahorse_norm_dna_baseline_summary <- seahorse_norm_dna_baseline_postgrubbs %>%
  group_by(record_id, donorage, donorsex, diagnosis, hla_a2) %>%
  summarise(across(where(is.numeric), mean))

seahorse_norm_dna_baseline_summary %>%
  ungroup() %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(calc_spare_cap ~ donorage + donorsex) #p=0.08 for sex
seahorse_norm_dna_baseline_summary %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ggplot(aes(x=donorsex, y=calc_spare_cap*100, fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 72, label = "p=0.08", label.size = 4, size = 0.5,tip.length = c(0.4, 0.02))+
  xlab("") +
  ylab("Spare respiratory capacity\n(% of baseline)") +
  scale_y_continuous(limits = c(0, 80))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line() 
  )
ggsave("Output/Fig6/Fig6H.png", width = 3, height = 3)


## Figure 7 and S4 ## ----

#clear environment
rm(list = ls())

#Read in metadata
donor_info <- read_excel("data/HPAP/Donor_Summary_192.xlsx")
donor_info <- donor_info %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM", clinical_diagnosis) ~ "T2D",
    grepl("T1DM", clinical_diagnosis) ~ "T1D",
    grepl("control", clinical_diagnosis) ~ "Control"
  ))

#Load data
#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/Calcium/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}
donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #82 donors

#Work out which runs are present for which donor (they go up to run 4). File names can have "-" or "_"
not_present1 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run1_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present1 <- c(not_present1,donor)
  }
}

not_present_1 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run1_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present_1 <- c(not_present_1,donor)
  }
}

norun1 <- intersect(not_present1, not_present_1) #just HPAP-069
donor_numbers1 <- setdiff(donor_numbers,norun1)

donor_data_run1 <- list()
for (i in 1:length(donor_numbers1)){
  donor <- donor_numbers1[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run1_data.xlsx")
  if (file.exists(path)) {
    donor_data_run1 [[i]] <- read_excel(paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run1_data.xlsx"), skip=2)
  } else {
    donor_data_run1 [[i]] <- read_excel(paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run1_data.xlsx"), skip=2)
  }
}

not_present2 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run2_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present2 <- c(not_present2,donor)
  }
}
not_present_2 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run2_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present_2 <- c(not_present_2,donor)
  }
}
norun2 <- intersect(not_present2, not_present_2) #9 donors
donor_numbers2 <- setdiff(donor_numbers,norun2)

donor_data_run2 <- list()
for (i in 1:length(donor_numbers2)){
  donor <- donor_numbers2[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run2_data.xlsx")
  if (file.exists(path)) {
    donor_data_run2 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run2_data.xlsx"), skip=2)
  } else {
    donor_data_run2 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run2_data.xlsx"), skip=2)
  }
}

not_present3 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium/HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run3_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present3 <- c(not_present3,donor)
  }
}
not_present_3 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run3_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present_3 <- c(not_present_3,donor)
  }
}
norun3 <- intersect(not_present3, not_present_3) #52 donors
donor_numbers3 <- setdiff(donor_numbers,norun3)

donor_data_run3 <- list()
for (i in 1:length(donor_numbers3)){
  donor <- donor_numbers3[i]
  path <- paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run3_data.xlsx")
  if (file.exists(path)) {
    donor_data_run3 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run3_data.xlsx"), skip=2)
  } else {
    donor_data_run3 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run3_data.xlsx"), skip=2)
  }
}

not_present4 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run4_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present4 <- c(not_present4,donor)
  }
}
not_present_4 <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  path <- paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run4_data.xlsx")
  if (file.exists(path)) {
  } else {
    not_present_4 <- c(not_present_4,donor)
  }
}
norun4 <- intersect(not_present4, not_present_4) #74 donors
donor_numbers4 <- setdiff(donor_numbers,norun4)

donor_data_run4 <- list()
for (i in 1:length(donor_numbers4)){
  donor <- donor_numbers4[i]
  path <- paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run4_data.xlsx")
  if (file.exists(path)) {
    donor_data_run4 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets_Run4_data.xlsx"), skip=2)
  } else {
    donor_data_run4 [[i]] <- read_excel(paste0("data/HPAP/Calcium//HPAP-",donor,"/Islet Studies/Islet physiology studies/Calcium imaging/HPAP-",donor,"_Calcium-imaging-Whole-islets-Run4_data.xlsx"), skip=2)
  }
}

#Join all runs together
all_runs <- c(donor_data_run1, donor_data_run2, donor_data_run3, donor_data_run4)

#Make a vector of the donors in all_runs_X in order.
donor_index <- c(donor_numbers1, donor_numbers2, donor_numbers3, donor_numbers4)
for (i in 1:length(donor_index)){
  ID <- paste0("HPAP_",donor_index[i])
  donor_index[i] <- ID
}

#Fix times to show duration from start rather than Excel times
all_runs_time_fixed <- list()
for (i in 1:length(all_runs)) {
  data <- all_runs[[i]]
  time_diff <- c(0)
  for (j in 2:nrow(data)){
    start <- data$`Time Stamp`[1]
    time <- data$`Time Stamp`[j]
    time_diff[j] <- difftime(time, start, units = "secs")
  }
  data$Time_from_start <- time_diff
  all_runs_time_fixed[[i]] <- data
}


#Exclude regions that had no calcium response to high glucose
all_runs_exclnoglcresp <- list()
for (i in 1:length(all_runs_time_fixed)){
  dataset <- all_runs_time_fixed[[i]]
  y <- ncol(dataset)
  dataset_new <- dataset[,c(1:2,y)]
  colnames_new <- c("Time Stamp","Condition")
  for (j in 3:(y-1)){
    if (is.numeric(dataset[[j]]) == TRUE){
      dataset_new <- cbind(dataset_new, dataset[[j]])
      colnames_new <- c(colnames_new, colnames(dataset)[j])
    }
  } #This ensures any columns that are marked "region omitted" or similar are not included
  y <- ncol(dataset_new)
  nregion <- y-3
  colnames(dataset_new)[4:y] <- paste0("Region",1:nregion)
  #check if HG is even available as some runs didn't have it. Also need a time when AAM is added to get last basal
  if (any(grepl("16.7", unique(dataset_new$Condition))) & any(grepl("4", unique(dataset_new$Condition)))){
    AAMadded <- dataset_new %>%
      filter(grepl("4m",Condition))
    last2minbasal <- dataset_new %>%
      filter(Time_from_start > AAMadded$Time_from_start[1] - 120, Time_from_start < AAMadded$Time_from_start[1]) %>%
      summarise(across(3:(3+nregion), \(x) mean(x, na.rm = TRUE)))
    #baseline the readings
    dataset_baselined <- dataset_new
    for (k in 4:(3+nregion)){
      baseline <- last2minbasal[1,k-2]
      temp <- dataset_new %>%
        mutate(baselined_region = .[[k]] - baseline)
      dataset_baselined <- cbind(dataset_baselined, temp$baselined_region)
    }
    colnames(dataset_baselined)[(4+nregion):(3+nregion+nregion)] <- paste0("Baselined_region",1:nregion)  
    HGadded <- dataset_baselined %>%
      filter(grepl("16.7",Condition))
    HGtime <- HGadded$Time_from_start[1]
    #Delay of ~7 min to peak from where HG occurs. Filter for at least 50% response comparing 1 min before HG added and at 6-8 min.
    HGtimeplus6 <- HGtime + 6*60
    HGtimeplus8 <- HGtime + 8*60
    HGtimeminus1 <- HGtime - 60
    HGaddedplus6to8min <- dataset_baselined %>%
      filter(Time_from_start >= HGtimeplus6, Time_from_start <= HGtimeplus8) %>%
      summarise(across((4+nregion):(3+nregion+nregion), mean, na.rm = TRUE))
    HGaddedminus1min <- dataset_baselined %>%
      filter(Time_from_start >= HGtimeminus1, Time_from_start <= HGtime) %>%
      summarise(across((4+nregion):(3+nregion+nregion), \(x) mean(x, na.rm = TRUE)))
    good_regions <- c()
    for (l in 1:nregion){
      test <- HGaddedplus6to8min[1,l] > (HGaddedminus1min[1,l] + (HGaddedminus1min[1,l]/2))
      if (test == TRUE){
        good_regions <- c(good_regions,l)
      }
    }
    dataset_good <- dataset_baselined[,c(1:3,(good_regions+3),(good_regions+3+nregion))]
    all_runs_exclnoglcresp[[i]] <- dataset_good
  } else{ #if no HG or AAM then exclude that run altogether
    all_runs_exclnoglcresp[[i]] <- NA
  }
  
} 

#Mean of regions per run
all_runs_means <- list()
for (i in 1:length(all_runs_exclnoglcresp)){
  if(all(is.na(all_runs_exclnoglcresp[[i]]))==FALSE){
    dataset <- all_runs_exclnoglcresp[[i]]
    nregion <- (ncol(dataset)-3)/2
    if (nregion >0){
      colnames(dataset)[4:(3+nregion)] <- paste0("Region",1:nregion)
      colnames(dataset)[(4+nregion):(3+nregion*2)] <- paste0("Baselined_region",1:nregion)
      cols_select <- c("Time Stamp","Condition","Time_from_start",paste0("Baselined_region",1:nregion))
      means <- dataset %>%
        dplyr::select(all_of(cols_select)) %>%
        pivot_longer(!c("Time Stamp","Condition","Time_from_start"), names_to = "Region",values_to = "Signal") %>%
        group_by(`Time Stamp`) %>%
        summarise(mean_baselined_signal = mean(Signal))
      dataset_means <- inner_join(dataset, means, by = "Time Stamp")
      dataset_means$Donor <- rep(donor_index[i], nrow(dataset_means))
      all_runs_means[[i]] <- dataset_means}
  } 
}

#Typical solution addition is No substrate, AAM, 3 G, 16.7 G, No substrate, KCl
#find ones that break this pattern
Solution_timings <- list()
for (i in 1:length(all_runs_means)){
  if (is.null(all_runs_means[[i]]) == FALSE){
  data <- all_runs_means[[i]] %>%
    filter(is.na(Condition) == FALSE) %>%
    dplyr::select(Condition, Time_from_start)
  Solution_timings[[i]] <- data
  } else {Solution_timings[[i]] <- NULL}
}

weird <- c()
for (i in 1:length(Solution_timings)){
  if (is.null(Solution_timings[[i]]) == FALSE){
  test1 <- grepl("No",Solution_timings[[i]]$Condition[1])
  test2 <- grepl("AAM",Solution_timings[[i]]$Condition[2])
  test3 <- grepl("3",Solution_timings[[i]]$Condition[3])
  test4 <- grepl("16.7",Solution_timings[[i]]$Condition[4])
  test5 <- grepl("No",Solution_timings[[i]]$Condition[5])
  test6 <- grepl("K",Solution_timings[[i]]$Condition[6])
    if (test1==TRUE & test2 == TRUE & test3 == TRUE & test4 == TRUE & test5==TRUE & test6 == TRUE){
  } else{  
    weird <- c(weird, i)
  }
}
}
weird #indices for unusual runs are 21  22  32  62  97 108 110 126 129 138 192

#find indices of runs with standard timings
standard_timings <- setdiff (1:192, weird)
standard_timings <- setdiff(standard_timings, which(sapply(all_runs_means, is.null)))

#figure out the shortest interval that sample is in each solution so that we can compare across different runs
Comparing_timings <- Solution_timings [[1]]
for (i in 2:length(standard_timings)){
  j <- standard_timings[i]
  Comparing_timings[i+1] <- Solution_timings[[j]][1:6,2]
}
head(Comparing_timings)

donors_standard_timings <- donor_index[standard_timings]
colnames(Comparing_timings) <- c("Condition", donors_standard_timings)
#rename second time no substrate is added to "washout"
Comparing_timings$Condition[5] <- "Washout"
#add run indices to donorIDs as column names
colnames(Comparing_timings)[2:ncol(Comparing_timings)] <- paste0(colnames(Comparing_timings)[2:ncol(Comparing_timings)],"_",standard_timings)

#change to long format
Comparing_timings_long <- Comparing_timings %>%
  pivot_longer(!Condition, names_to = "Donor", values_to = "Time")
ggplot(Comparing_timings_long, aes(x = Condition, y = Time)) +
  geom_point () +
  theme(legend.position = "none") #Timings of each addition roughly clusters together but with some discrepancy

#figure out the shortest interval for each substrate
shortest_interval <- Comparing_timings_long %>%
  pivot_wider(names_from = Condition, values_from = Time) %>%
  mutate(AAMint = `Amino Acid Mix (AAM) 4mM` - `No Substrate`,
         LGint = `3mM Glucose + AAM` - `Amino Acid Mix (AAM) 4mM`,
         HGint = `16.7mM Glucose + AAM` - `3mM Glucose + AAM`,
         Washoutint = Washout - `16.7mM Glucose + AAM`,
         KClint = `KCl, 30mM` - Washout) %>%
  summarise(minAAM = min(AAMint),
            minLG = min(LGint),
            minHG = min(HGint),
            minwash = min(Washoutint),
            minKCl = min(KClint),
            meanAAM = mean(AAMint),
            meanLG = mean(LGint),
            meanHG = mean(HGint),
            meanwash = mean(Washoutint),
            meanKCl = mean(KClint),
            meanminus3SDAAM = mean(AAMint) - 3*sd(AAMint),
            meanminus3SDLG = mean(LGint) - 3*sd(LGint),
            meanminus3SDHG = mean(HGint) - 3*sd(HGint),
            meanminus3SDwash = mean(Washoutint) + 3*sd(Washoutint),
            meanminus3SDKCl = mean(KClint) - 3*sd(KClint))
shortest_interval %>% View #for most substrates, the minimum interval time is not shorter than 3 SD lower than mean. LG is slightly shorted (10 sec) and washout is a fair bit shorter but we are not quantifying that one so it should be okay to just use the min values.

#Take the first min-interval minutes for each of the stimuli for each donor
#First do this for the runs with standard timings
all_runs_min_intervals <- list()
for (j in 1:length(standard_timings)){
  i <- standard_timings[j]
  stimuli_time <- all_runs_means[[i]] %>%
    filter(is.na(Condition) == FALSE)
  extra_basal <- all_runs_means[[i]] %>%
    filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
  extra_AAM <- all_runs_means[[i]] %>%
    filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
  extra_LG <- all_runs_means[[i]] %>%
    filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
  extra_HG <- all_runs_means[[i]] %>%
    filter(Time_from_start > (260 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
  extra_washout <- all_runs_means[[i]] %>%
    filter(Time_from_start > (710 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
  cleaned <- all_runs_means[[i]] %>%
    anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
    anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
    anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
    anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
    anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
  cleaned <- cleaned %>%
    mutate(Time_from_previous_stimulus = case_when(
      Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
      Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
      Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
      Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
      Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
      Time_from_start > stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[6]))
  cleaned <- cleaned %>%
    mutate(adj.time = case_when(
      Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
      Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
      Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
      Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
      Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6]~ Time_from_previous_stimulus + 1500,
      Time_from_start > stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 2210))
  all_runs_min_intervals[[i]] <- cleaned
}

#Fix weird ones
weird

#index 21
Solution_timings[[21]] #There are 2x 10 second pauses between washout and KCl. Ignore these.
i <- 21
all_runs_means[[i]] <- all_runs_means[[i]] %>%
  filter(Time_from_start <2261 | Time_from_start > 2486) #Remove pause period
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[10])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[10] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[10] ~ Time_from_start - stimuli_time$Time_from_start[10]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[10]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[10] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 22
Solution_timings[[22]] #No washout - skip straight to KCl
i <- 22
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[5]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 32
Solution_timings[[32]] #Pause 10s between low and high glucose
i <- 32
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[6])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[7]), Time_from_start < stimuli_time$Time_from_start[8])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[7],
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[8]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 62
Solution_timings[[62]] #No substrate listed twice at start
i <- 62
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[3])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[7]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 97
Solution_timings[[97]] #no KCl, can use other data until this point
i <- 97
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[6]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1500
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 108
Solution_timings[[108]] #Pause 10s between AAM and LG
i <- 108
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[5])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[7]), Time_from_start < stimuli_time$Time_from_start[8])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[7],
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[8]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 110
Solution_timings[[110]] #Extra 10s at start
i <- 110
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[3])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[7]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 126
Solution_timings[[126]] #Extra 10s at start
i <- 126
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[3])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[7]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[6]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 129
Solution_timings[[129]] #No KCl
i <- 129
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[3])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[6]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 1500
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 138
Solution_timings[[138]] #Pause 10s between AAM and LG
i <- 138
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[2])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[2]), Time_from_start < stimuli_time$Time_from_start[5])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[7]), Time_from_start < stimuli_time$Time_from_start[8])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[2],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[7],
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_start - stimuli_time$Time_from_start[8]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[2] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[2] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[7] & Time_from_start < stimuli_time$Time_from_start[8]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[8] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#index 192
Solution_timings[[192]] #Extra 10s at start
i <- 192
stimuli_time <- all_runs_means[[i]] %>%
  filter(is.na(Condition) == FALSE)
extra_basal <- all_runs_means[[i]] %>%
  filter(Time_from_start > 420, Time_from_start < stimuli_time$Time_from_start[3])
extra_AAM <- all_runs_means[[i]] %>%
  filter(Time_from_start > (440 + stimuli_time$Time_from_start[3]), Time_from_start < stimuli_time$Time_from_start[4])
extra_LG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (380 + stimuli_time$Time_from_start[4]), Time_from_start < stimuli_time$Time_from_start[5])
extra_HG <- all_runs_means[[i]] %>%
  filter(Time_from_start > (260 + stimuli_time$Time_from_start[5]), Time_from_start < stimuli_time$Time_from_start[6])
extra_washout <- all_runs_means[[i]] %>%
  filter(Time_from_start > (710 + stimuli_time$Time_from_start[6]), Time_from_start < stimuli_time$Time_from_start[7])
cleaned <- all_runs_means[[i]] %>%
  anti_join(extra_basal, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_AAM, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_LG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_HG, by = c("Time Stamp","Time_from_start","mean_baselined_signal")) %>%
  anti_join(extra_washout, by = c("Time Stamp","Time_from_start","mean_baselined_signal"))
cleaned <- cleaned %>%
  mutate(Time_from_previous_stimulus = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_start,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_start - stimuli_time$Time_from_start[3],
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_start - stimuli_time$Time_from_start[4],
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_start - stimuli_time$Time_from_start[5],
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[6],
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_start - stimuli_time$Time_from_start[7]))
cleaned <- cleaned %>%
  mutate(adj.time = case_when(
    Time_from_start < stimuli_time$Time_from_start[3] ~ Time_from_previous_stimulus,
    Time_from_start > stimuli_time$Time_from_start[3] & Time_from_start < stimuli_time$Time_from_start[4] ~ Time_from_previous_stimulus + 420,
    Time_from_start > stimuli_time$Time_from_start[4] & Time_from_start < stimuli_time$Time_from_start[5] ~ Time_from_previous_stimulus + 860,
    Time_from_start > stimuli_time$Time_from_start[5] & Time_from_start < stimuli_time$Time_from_start[6] ~ Time_from_previous_stimulus + 1240,
    Time_from_start > stimuli_time$Time_from_start[6] & Time_from_start < stimuli_time$Time_from_start[7]~ Time_from_previous_stimulus + 1500,
    Time_from_start > stimuli_time$Time_from_start[7] ~ Time_from_previous_stimulus + 2210
  ))
all_runs_min_intervals[[i]] <- cleaned

#combined into single dataframe
all_runs_min_intervals_real <- Filter(Negate(is.null), all_runs_min_intervals)

all_runs_min_intervalsdf <- list()
for (i in 1:length(all_runs_min_intervals_real)){
  data <- all_runs_min_intervals_real[[i]]
  data <- data %>%
    dplyr::select(Donor, mean_baselined_signal,adj.time, Time_from_start, Time_from_previous_stimulus)
  all_runs_min_intervalsdf[[i]] <- data
}
all_runs_min_intervalsdf <- data.frame(do.call(rbind, all_runs_min_intervalsdf))
summary(all_runs_min_intervalsdf)

#summarise per donor 
per_donor_mininterval <- all_runs_min_intervalsdf %>%
  mutate(adj.time = round(adj.time/60)) %>% #change seconds to minutes and round to nearest minute
  group_by(Donor, adj.time) %>%
  summarise(mean_baselined_signal = mean(mean_baselined_signal, na.rm = TRUE))

#add metadata
per_donor_mininterval$Donor <- gsub("_","-",per_donor_mininterval$Donor)
per_donor_mininterval <- left_join(per_donor_mininterval, donor_info, by = c("Donor" = "donor_ID"))

#Calculate AUCs
#AAM
calculate.auc.AAM <- function(x) {
    n <- length(x)
  if (n < 2) return(NA)  # Not enough points
  time_points <- seq_len(n) + 6 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

AAMdf <- per_donor_mininterval %>%
  filter(adj.time >= 7, adj.time <= 14) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
AAMdf$AUC.AAM <- apply(AAMdf[2:9], 1, calculate.auc.AAM)
AUC_summary <- AAMdf %>% dplyr::select(Donor, AUC.AAM)

#LG
calculate.auc.LG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)  # Not enough points
  time_points <- seq_len(n) + 16 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

LGdf <- per_donor_mininterval %>%
  filter(adj.time >= 17, adj.time <= 23) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
LGdf$AUC.LG <- apply(LGdf[2:8], 1, calculate.auc.LG)
AUC_summary <- LGdf %>% dplyr::select(Donor, AUC.LG) %>% full_join(AUC_summary, by = "Donor")

#HG
calculate.auc.HG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)  # Not enough points
  time_points <- seq_len(n) + 22 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

HGdf <- per_donor_mininterval %>%
  filter(adj.time >= 23, adj.time <= 35) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
HGdf$AUC.HG <- apply(HGdf[2:14], 1, calculate.auc.HG)
AUC_summary <- HGdf %>% dplyr::select(Donor, AUC.HG) %>% full_join(AUC_summary, by = "Donor")

#KCl
calculate.auc.KCl <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)  # Not enough points
  time_points <- seq_len(n) + 39 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
KCldf <- per_donor_mininterval %>%
  filter(adj.time >= 40, adj.time <= 45) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
KCldf$AUC.KCl <- apply(KCldf[2:7], 1, calculate.auc.KCl)
AUC_summary <- KCldf %>% dplyr::select(Donor, AUC.KCl) %>% full_join(AUC_summary, by = "Donor")

#Add metadata
AUC_summary <- left_join(AUC_summary, donor_info, by = c("Donor"="donor_ID"))

## Fig 7A ## ----
stimulus_data <- data.frame(
  xmin = c(0,7,14,21,25,37),
  xmax = c(7,14,21,25,37,45),
  label = c("Base\nline","AAM","G 3","G\n16.7","Washout","KCl"),
  ymin = c(0.9, 0.9 ,1.15, 1.15 ,1.15,1.15),
  ymax = c(1.15, 1.15, 1.4, 1.4, 1.4, 1.4)
)
per_donor_mininterval %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  mutate(adj.time = round(adj.time)) %>%
  filter(adj.time <= 45) %>%
  group_by(adj.time, sex, simplified_diagnosis) %>%
  summarise(mean = mean(mean_baselined_signal, na.rm = TRUE), up = mean(mean_baselined_signal, na.rm = TRUE) + sd(mean_baselined_signal, na.rm = TRUE), down = mean(mean_baselined_signal, na.rm = TRUE) - sd(mean_baselined_signal, na.rm = TRUE)) %>%
  ggplot(aes(x=adj.time, y = mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(colour = NA, alpha = 0.2,(aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  ylab("Baseline-corrected calcium signal") + 
  facet_wrap(~sex)+
  geom_vline(xintercept = stimulus_data$xmin, colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "black",size = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig7/Fig7A.png", width = 7, height = 4)

## Fig 7B ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.AAM ~ age_years + sex*simplified_diagnosis) #ns
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.AAM ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.AAM,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC AAM") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Fig7/Fig7B.png",width = 5, height = 4)

## Fig 7C ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.LG ~ age_years + sex*simplified_diagnosis) #ns
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.LG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.LG,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC 3 mM Glucose") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Fig7/Fig7C.png",width = 5, height = 4)

## Fig 7D ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.HG ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.HG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0763 for Ctrl vs T2D in females, p=0.16 in males

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.HG,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 8.4, label = "p=0.076", label.size = 4, size = 0.5, tip.length = c(0.02, 0.35))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.1, label = "p=0.16", label.size = 4, size = 0.5, tip.length = c(0.02, 0.3))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC 16.7 mM Glucose") +
  scale_y_continuous(limits= c(-1,9))+
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Fig7/Fig7D.png",width = 5, height = 4)

## Fig 7E ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.KCl ~ age_years + sex*simplified_diagnosis) #p=0.057 for disease
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.KCl ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.12 for females, p=0.18 for males

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.KCl,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 7.2, label = "p=0.12", label.size = 4, size = 0.5, tip.length = c(0.02, 0.4))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.2, label = "p=0.18", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC KCl") +
  xlab("")+
  scale_y_continuous(limits= c(0,8))+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Fig7/Fig7E.png",width = 5, height = 4)

## Fig S4A ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,7,14,21,25,37),
  xmax = c(7,14,21,25,37,45),
  label = c("Base\nline","AAM","G 3","G\n16.7","Washout","KCl"),
  ymin = c(0.95, 1.05 ,1.3, 1.2 ,1.3,1.3), 
  ymax = c(1.25, 1.15, 1.4, 1.5, 1.4, 1.4)
)
per_donor_mininterval %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  mutate(adj.time = round(adj.time)) %>%
  filter(adj.time <= 45) %>%
  group_by(adj.time, sex, simplified_diagnosis) %>%
  summarise(mean = mean(mean_baselined_signal, na.rm = TRUE), up = mean(mean_baselined_signal, na.rm = TRUE) + sd(mean_baselined_signal, na.rm = TRUE), down = mean(mean_baselined_signal, na.rm = TRUE) - sd(mean_baselined_signal, na.rm = TRUE)) %>%
  ggplot(aes(x=adj.time, y = mean, colour = sex, group = sex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(colour = NA, alpha = 0.2,(aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  ylab("Baseline-corrected calcium signal") +
  geom_vline(xintercept = stimulus_data2$xmin, colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "black",size = 0.5
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS4/FigS4A.png", width = 5.5, height = 4)

## Fig S4B ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(AUC.AAM ~ age_years + sex) #ns
AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=AUC.AAM,  fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC AAM") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/FigS4/FigS4B.png",width = 3, height = 3)

## Fig S4C ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(AUC.LG ~ age_years + sex) #ns
AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=AUC.LG,  fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC 3 mM Glucose") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/FigS4/FigS4C.png",width = 3, height = 3)

## Fig S4B ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(AUC.HG ~ age_years + sex) #ns
AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=AUC.HG,  fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC 16.7 mM Glucose") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/FigS4/FigS4D.png",width = 3, height = 3)

## Fig S4E ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(AUC.KCl ~ age_years + sex) #ns
AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=AUC.KCl,  fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC KCl") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/FigS4/FigS4E.png",width = 3, height = 3)

### Figure 8 and S5-6,9-10 ### ----

#UPenn perifusion data
#clear environment
rm(list = ls())

#Read in metadata
donors <- read_excel("data/HPAP/Donor_Summary_192.xlsx")

#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/UPenn Perifusion/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}

donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #160

#Read in the new data
perifusion_data <- list()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  data <- read.csv(paste0("data/HPAP/UPenn Perifusion/HPAP-",donor,"/Islet Studies/Islet physiology studies/Perifusions/HPAP-",donor,"_Perifusion_data.csv"))
  perifusion_data[[i]] <- data
}

#fix all the colnames so they are the same
for (i in 1:length(donor_numbers)){
  colnames(perifusion_data[[i]]) <- colnames(perifusion_data[[1]])
}

perifusion_all <- do.call(rbind, perifusion_data)
perifusion_all <- data.frame(perifusion_all)
head(perifusion_all)

#add donor_ids
donor_id_col <- c()
for (i in 1:length(donor_numbers)){
  to_add <- rep(paste0("HPAP-",donor_numbers[i]),nrow(perifusion_data[[i]]))
  donor_id_col <- c(donor_id_col, to_add)
}
perifusion_all$DonorID <- donor_id_col

summary(perifusion_all) #everything except time is character

#Find what is not numeric
insulin.release.char <- which(is.na(as.numeric(perifusion_all$Insulin.Release..ng.100.islets.min.)))
perifusion_all_insulin.release.char <- perifusion_all[insulin.release.char,]
unique(perifusion_all_insulin.release.char$Insulin.Release..ng.100.islets.min.) #seems the cause of the character elements was "NA", the units being written in that column, or "undetectable"
glucagon.release.char <- which(is.na(as.numeric(perifusion_all$Glucagon.Release..pg.100.islets.min.)))
perifusion_all_glucagon.release.char <- perifusion_all[glucagon.release.char,]
unique(perifusion_all_glucagon.release.char$Glucagon.Release..pg.100.islets.min.) #same for glucagon

#Change to numeric
perifusion_all$Insulin.Release..ng.100.islets.min. <- as.numeric(perifusion_all$Insulin.Release..ng.100.islets.min.)
perifusion_all$Glucagon.Release..pg.100.islets.min. <- as.numeric(perifusion_all$Glucagon.Release..pg.100.islets.min.)
perifusion_all$Insulin.Content...ng.islet. <- as.numeric(perifusion_all$Insulin.Content...ng.islet.)
perifusion_all$Glucagon.Content..pg.islet. <- as.numeric(perifusion_all$Glucagon.Content..pg.islet.)

#add metadata
donors <- donors %>%
  mutate(
    simplified_diagnosis = case_when(
      grepl("T2DM", clinical_diagnosis) ~ "T2D",
      grepl("T1DM", clinical_diagnosis) ~ "T1D",
      grepl("control", clinical_diagnosis) ~ "Control",
      .default = clinical_diagnosis
    ))
perifusion_all_with_metadata <- perifusion_all %>%
  left_join(donors, by= c("DonorID" = "donor_ID"))

#calculate AUCs
#AAM
calculate.auc.AAM <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 29 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

AAMins <- perifusion_all_with_metadata %>%
  filter(Stimulus == "Amino acid Mixure, 4 mM") %>%
  dplyr::select(DonorID, Time..min., Insulin.Release..ng.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Insulin.Release..ng.100.islets.min.)
AAMins$AUC.ins.AAM <- apply(AAMins[,c(2:32)], 1, calculate.auc.AAM)

insulinAUCs <- AAMins %>% dplyr::select(DonorID, AUC.ins.AAM)

AAMgcg <- perifusion_all_with_metadata %>%
  filter(Stimulus == "Amino acid Mixure, 4 mM") %>%
  dplyr::select(DonorID, Time..min., Glucagon.Release..pg.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Glucagon.Release..pg.100.islets.min.)
AAMgcg$AUC.gcg.AAM <- apply(AAMgcg[,c(2:32)], 1, calculate.auc.AAM)

glucagonAUCs <- AAMgcg %>% dplyr::select(DonorID, AUC.gcg.AAM)

#LG
calculate.auc.LG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 60 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

LGins <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 3 mM", "Glucose, 3 mM")) %>%
  dplyr::select(DonorID, Time..min., Insulin.Release..ng.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Insulin.Release..ng.100.islets.min.)
LGins$AUC.ins.LG <- apply(LGins[,c(2:21)], 1, calculate.auc.LG)

insulinAUCs <- LGins %>% dplyr::select(DonorID, AUC.ins.LG) %>% full_join(insulinAUCs, by = "DonorID")

LGgcg <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 3 mM", "Glucose, 3 mM")) %>%
  dplyr::select(DonorID, Time..min., Glucagon.Release..pg.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Glucagon.Release..pg.100.islets.min.)
LGgcg$AUC.gcg.LG <- apply(LGgcg[,c(2:21)], 1, calculate.auc.LG)

glucagonAUCs <- LGgcg %>% dplyr::select(DonorID, AUC.gcg.LG) %>% full_join(glucagonAUCs, by = "DonorID")

#HG
calculate.auc.HG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 80 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

HGins <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 16.7 mM", "Glucose, 16.7 mM")) %>%
  dplyr::select(DonorID, Time..min., Insulin.Release..ng.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Insulin.Release..ng.100.islets.min.)
HGins$AUC.ins.HG <- apply(HGins[,c(2:21)], 1, calculate.auc.HG)

insulinAUCs <- HGins %>% dplyr::select(DonorID, AUC.ins.HG) %>% full_join(insulinAUCs, by = "DonorID")

HGgcg <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 16.7 mM", "Glucose, 16.7 mM")) %>%
  dplyr::select(DonorID, Time..min., Glucagon.Release..pg.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Glucagon.Release..pg.100.islets.min.)
HGgcg$AUC.gcg.HG <- apply(HGgcg[,c(2:21)], 1, calculate.auc.HG)

glucagonAUCs <- HGgcg %>% dplyr::select(DonorID, AUC.gcg.HG) %>% full_join(glucagonAUCs, by = "DonorID")

#HG + IBMX
calculate.auc.IBMX <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 100 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

IBMXins <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 16.7 mM + IBMX, 0.1 mM", "Glucose, 16.7 mM + IBMX, 0.1 mM")) %>%
  dplyr::select(DonorID, Time..min., Insulin.Release..ng.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Insulin.Release..ng.100.islets.min.)
IBMXins$AUC.ins.IBMX <- apply(IBMXins[,c(2:21)], 1, calculate.auc.IBMX)

insulinAUCs <- IBMXins %>% dplyr::select(DonorID, AUC.ins.IBMX) %>% full_join(insulinAUCs, by = "DonorID")

IBMXgcg <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("Amino acid Mixure, 4 mM + Glucose, 16.7 mM + IBMX, 0.1 mM", "Glucose, 16.7 mM + IBMX, 0.1 mM")) %>%
  dplyr::select(DonorID, Time..min., Glucagon.Release..pg.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Glucagon.Release..pg.100.islets.min.)
IBMXgcg$AUC.gcg.IBMX <- apply(IBMXgcg[,c(2:21)], 1, calculate.auc.IBMX)

glucagonAUCs <- IBMXgcg %>% dplyr::select(DonorID, AUC.gcg.IBMX) %>% full_join(glucagonAUCs, by = "DonorID")

#KCl
calculate.auc.KCl <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 140 
  if (sum(!is.na(x)) >= 2) {
    x_interp <- zoo::na.approx(x, x = time_points, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = FALSE, na.rm = FALSE)
    x_interp <- zoo::na.locf(x_interp, fromLast = TRUE, na.rm = FALSE)
    auc_value <- MESS::auc(time_points, x_interp, type = "linear")
  } else {
    auc_value <- NA
  }
  return(auc_value)
}

KClins <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("KCl, 30 mM")) %>%
  dplyr::select(DonorID, Time..min., Insulin.Release..ng.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Insulin.Release..ng.100.islets.min.)
KClins$AUC.ins.KCl <- apply(KClins[,c(2:21)], 1, calculate.auc.KCl)

insulinAUCs <- KClins %>% dplyr::select(DonorID, AUC.ins.KCl) %>% full_join(insulinAUCs, by = "DonorID")

KClgcg <- perifusion_all_with_metadata %>%
  filter(Stimulus %in% c("KCl, 30 mM")) %>%
  dplyr::select(DonorID, Time..min., Glucagon.Release..pg.100.islets.min.) %>%
  pivot_wider(names_from = Time..min., values_from = Glucagon.Release..pg.100.islets.min.)
KClgcg$AUC.gcg.KCl <- apply(KClgcg[,c(2:21)], 1, calculate.auc.KCl)

glucagonAUCs <- KClgcg %>% dplyr::select(DonorID, AUC.gcg.KCl) %>% full_join(glucagonAUCs, by = "DonorID")

#mean baseline secretion
baselineins <- perifusion_all_with_metadata %>%
  filter(Stimulus == "No Stimuli", Time..min. < 30) %>%
  group_by(DonorID) %>%
  summarise(mean.ins.baseline = mean(Insulin.Release..ng.100.islets.min.))
insulinAUCs <- baselineins %>% dplyr::select(DonorID, mean.ins.baseline) %>% full_join(insulinAUCs, by = "DonorID")

baselinegcg <- perifusion_all_with_metadata %>%
  filter(Stimulus == "No Stimuli", Time..min. < 30) %>%
  group_by(DonorID) %>%
  summarise(mean.gcg.baseline = mean(Glucagon.Release..pg.100.islets.min.))
glucagonAUCs <- baselinegcg %>% dplyr::select(DonorID, mean.gcg.baseline) %>% full_join(glucagonAUCs, by = "DonorID")

#Add metadata to AUCs
perifusion_summary <- full_join(insulinAUCs, donors, by = c("DonorID" = "donor_ID"))
perifusion_summary <- full_join(perifusion_summary, glucagonAUCs, by = "DonorID")
head(perifusion_summary)
tail(colnames(perifusion_summary))

## Fig 8A ## ----
stimulus_data <- data.frame(
  xmin = c(20,30,61,81,101,121,141), 
  xmax = c(30,121,81,121,121,141,160),
  label = c("Base\nline","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(8.6,8.9,9.5,9.5,10.1,9.8,10.1), 
  ymax = c(9.8,9.5,10.1,10.1,10.7,11,10.7))

perifusion_all_with_metadata %>%
  filter(Time..min. < 160) %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time..min., simplified_diagnosis, sex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.Release..ng.100.islets.min., na.rm = TRUE),
    sd = sd(Insulin.Release..ng.100.islets.min., na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time..min., y=mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab ("Insulin release (ng/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  #coord_cartesian(xlim = c(0, 220), ylim = c(0,520)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(20,30,61,81,101,121,141), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig8/Fig8A.png", width = 8, height = 4)

## Fig 8B ## ----
hist(perifusion_summary$mean.ins.baseline)
hist(log(perifusion_summary$mean.ins.baseline))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(mean.ins.baseline) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(mean.ins.baseline) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0028 for Control vs T2D in females, p=0.0822 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean.ins.baseline), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion \n(ng/100 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-4.5,2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8B.png", width = 4 , height = 4)

## Fig 8C ## ----
hist(perifusion_summary$AUC.ins.AAM)
hist(log(perifusion_summary$AUC.ins.AAM))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.AAM) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.AAM) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0028 for Control vs T2D in females, p=0.0228 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.AAM), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-1,5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8C.png", width = 4 , height = 4)

## Fig 8D ## ----
hist(perifusion_summary$AUC.ins.LG)
hist(log(perifusion_summary$AUC.ins.LG))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.LG) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.LG) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0186 for Control vs T2D in females, p=0.0944 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.LG), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.12))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n3 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-1.5,4)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8D.png", width = 4 , height = 4)

## Fig 8E ## ----
hist(perifusion_summary$AUC.ins.HG)
hist(log(perifusion_summary$AUC.ins.HG))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.HG) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.HG) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0002 for Control vs T2D in females, p=0.0003 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.HG), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.1, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-0.5,5.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8E.png", width = 4 , height = 4)

## Fig 8F ## ----
hist(perifusion_summary$AUC.ins.IBMX)
hist(log(perifusion_summary$AUC.ins.IBMX))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.IBMX) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0001 for Control vs T2D in females, p<0.0001 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,6.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8F.png", width = 4 , height = 4)

## Fig 8G ## ----
hist(perifusion_summary$AUC.ins.KCl)
hist(log(perifusion_summary$AUC.ins.KCl))
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.KCl) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0012 for Control vs T2D in females, p=0.0011 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,5.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Fig8/Fig8G.png", width = 4 , height = 4)

## Fig S5A ## ----
stimulus_data2 <- data.frame(
  xmin = c(20,30,61,81,101,121,141),
  xmax = c(30,121,81,121,121,141,160),
  label = c("Base\nline","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(8.6,8.9,9.5,9.5,10.1,9.8,10.1), 
  ymax = c(9.8,9.5,10.1,10.1,10.7,11,10.7))

perifusion_all_with_metadata %>%
  filter(Time..min. < 160) %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ungroup() %>%
  group_by(Time..min., sex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.Release..ng.100.islets.min., na.rm = TRUE),
    sd = sd(Insulin.Release..ng.100.islets.min., na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time..min., y=mean, colour = sex, group = sex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("Insulin release (ng/100 islets)") +
  xlab ("Time (min)") + 
  geom_vline(xintercept = c(20,30,61,81,101,121,141), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS5/FigS5A.png", width = 5.5, height = 4)

## Fig S5B ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(mean.ins.baseline) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(mean.ins.baseline), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline insulin secretion \n(ng/100 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5B.png", width = 4 , height = 4)

## Fig S5C ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.AAM) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.AAM), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5C.png", width = 4 , height = 4)

## Fig S5D ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.LG) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.LG), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n3 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5D.png", width = 4 , height = 4)

## Fig S5E ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.HG) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.HG), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5E.png", width = 4 , height = 4)

## Fig S5F ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.IBMX) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.IBMX), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5F.png", width = 4 , height = 4)

## Fig S5G ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.KCl) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.KCl), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS5/FigS5G.png", width = 4 , height = 4)

## Fig S9A ## ----
stimulus_data3 <- data.frame(
  xmin = c(18,30,61,81,101,121,141),    
  xmax = c(30,121,81,121,121,141,160),
  label = c("Base\nline","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(335,350,370,370,400,385,400), 
  ymax = c(385,370,400,400,430,445,430))
perifusion_all_with_metadata %>%
  filter(Time..min. < 160) %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  ungroup() %>%
  group_by(Time..min., sex) %>%
  summarise(
    n = n(),
    mean = mean(Glucagon.Release..pg.100.islets.min., na.rm = TRUE),
    sd = sd(Glucagon.Release..pg.100.islets.min., na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time..min., y=mean, colour = sex, group = sex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("Glucagon release (pg/100 islets)") +
  xlab ("Time (min)") +
  geom_vline(xintercept = c(20,30,61,81,101,121,141), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS9/FigS9A.png", width = 5.5, height = 4)

## Fig S9B ## ----
hist(perifusion_summary$mean.gcg.baseline)
hist(log(perifusion_summary$mean.gcg.baseline))
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(mean.gcg.baseline) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(mean.gcg.baseline), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline glucagon secretion\n(pg/100 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS9/FigS9B.png", width = 3 , height = 3)

## Fig S9C ## ----
hist(perifusion_summary$AUC.gcg.AAM)
hist(log(perifusion_summary$AUC.gcg.AAM))
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.gcg.AAM) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.AAM), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS9/FigS9C.png", width = 3 , height = 3)

## Fig S9D ## ----
hist(perifusion_summary$AUC.gcg.IBMX)
hist(log(perifusion_summary$AUC.gcg.IBMX))
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.gcg.IBMX) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.IBMX), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS9/FigS9D.png", width = 3 , height = 3)

## Fig S9E ## ----
hist(perifusion_summary$AUC.gcg.KCl)
hist(log(perifusion_summary$AUC.gcg.KCl))
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.gcg.KCl) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.KCl), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS9/FigS9E.png", width = 3 , height = 3)

## Fig S9K ## ----
glucagon.content <- perifusion_all_with_metadata %>%
  dplyr:: select(Glucagon.Content..pg.islet., DonorID, sex, age_years, simplified_diagnosis) %>%
  distinct()
hist(glucagon.content$Glucagon.Content..pg.islet.)
hist(log(glucagon.content$Glucagon.Content..pg.islet.))
glucagon.content %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(Glucagon.Content..pg.islet.) ~ age_years + sex) #ns

glucagon.content %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(Glucagon.Content..pg.islet.), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(glucagon content\n(pg/islet))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/FigS9/FigS9K.png", width = 3 , height = 3)

## Fig S10A ## ----
perifusion_all_with_metadata %>%
  filter(Time..min. < 160) %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time..min., simplified_diagnosis, sex) %>%
  summarise(
    n = n(),
    mean = mean(Glucagon.Release..pg.100.islets.min., na.rm = TRUE),
    sd = sd(Glucagon.Release..pg.100.islets.min., na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time..min., y=mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab ("Glucagon release (pg/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  #coord_cartesian(xlim = c(0, 220), ylim = c(0,520)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(20,30,61,81,101,121,141), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS10/FigS10A.png", width = 8, height = 4)

## Fig S10B ## ----
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(mean.gcg.baseline) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(mean.gcg.baseline) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0190 for Control vs T2D in females, p=0.1863 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean.gcg.baseline), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion \n(pg/100 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0.5,5.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/FigS10/FigS10B.png", width = 4 , height = 4)

## Fig S10C ## ----
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.AAM) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.AAM) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0156 for Control vs T2D in females, p=0.0934 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.AAM), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.8, label = "p=0.09", label.size = 4, size = 0.5, tip.length = c(0.02, 0.07))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(4.5,9.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/FigS10/FigS10C.png", width = 4 , height = 4)

## Fig S10D ## ----
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.IBMX) ~ age_years + sex*simplified_diagnosis) #ns
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.1152 for Control vs T2D in females, p=0.5490 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/FigS10/FigS10D.png", width = 4 , height = 4)

## Fig S10E ## ----
perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.KCl) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0185 for Control vs T2D in females

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(4.5,9.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/FigS10/FigS10E.png", width = 4 , height = 4)

## Fig S10K ## ----
glucagon.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(Glucagon.Content..pg.islet.) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- glucagon.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(Glucagon.Content..pg.islet.) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0488 for Control vs T2D in females, p=0.0204 for males

glucagon.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Glucagon.Content..pg.islet.), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 10, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.35))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/islet))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(3,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/FigS10/FigS10K.png", width = 5 , height = 4)


#Vanderbilt insulin perifusion data
#clear environment
rm(list = ls())

#Read in metadata
donors <- read_excel("data/HPAP/Donor_Summary_192.xlsx")
donors <- donors %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#Set up donor IDs that exist in this dataset
donor_numbers <- paste0("00",1:9)
donor_numbers <- c(donor_numbers,paste0("0",10:99))
donor_numbers <- c(donor_numbers,100:193)

not_present <- c()
for (i in 1:length(donor_numbers)){
  donor <- donor_numbers[i]
  folder_path <- paste0("data/HPAP/Vanderbilt Perifusion/HPAP-",donor)
  if (dir.exists(folder_path)) {
  } else {
    not_present <- c(not_present,donor)
  }
}

donor_numbers <- setdiff(donor_numbers,not_present)
length(donor_numbers) #97

donor_list <- paste0("HPAP-",donor_numbers)

#Read in data
all_insulin <- list()
for (i in 1:97){
  donor <- donor_list[i]
  if (file.exists(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
    if ("Insulin Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
      data <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx"), sheet = "Insulin Results", range = "A1:F52")
      data <- data[2:51,] #remove first row which is NA because it's merged with the header row
      data$Stimulus[1] <- data$Stimulus[2] #First stimulus entry has donor info so let's remove that
      colnames(data) <- c("Fraction","Time.min","Stimulus","Insulin.ng.mL","Insulin.ng.100IEQs.min","Insulin.pc.content.min")
      data$DonorID <- rep(donor,nrow(data))} else{
        data <- NA
      }
  } else{
    if ("Insulin Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx")) == TRUE){
      data <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx"), sheet = "Insulin Results", range = "A1:F52")
      data <- data[2:51,] #remove first row which is NA because it's merged with the header row
      data$Stimulus[1] <- data$Stimulus[2] #First stimulus entry has donor info so let's remove that
      colnames(data) <- c("Fraction","Time.min","Stimulus","Insulin.ng.mL","Insulin.ng.100IEQs.min","Insulin.pc.content.min")
      data$DonorID <- rep(donor,nrow(data))} else{
        data <- NA
      }
  }
  all_insulin[[i]] <- data
}

which(is.na(all_insulin)) #18, 30
all_insulin <- all_insulin[-18]
all_insulin <- all_insulin[-29]
all_insulin_df <- data.frame(do.call(rbind, all_insulin))

#add metadata
all_insulin_df <- inner_join(all_insulin_df, donors, by = c("DonorID"="donor_ID"))

#Set up to calculte AUCs
all_insulin_wide <- all_insulin_df %>%
  dplyr::select(-Insulin.ng.mL) %>%
  dplyr::select(-Insulin.pc.content.min) %>%
  group_by(DonorID, Time.min) %>%
  summarise(Insulin.ng.100IEQs.min = mean(Insulin.ng.100IEQs.min)) %>% #to deal with duplicate values for a single time point for a single donor
  pivot_wider(values_from = Insulin.ng.100IEQs.min, names_from = Time.min)
head(all_insulin_wide) 
#15 is in the wrong place - move up.
all_insulin_wide <- all_insulin_wide[, c(1:5, 51, 6:50)]
head(all_insulin_wide)

#Calculate AUCs
#High glucose
calculate.auc.HG <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(12,15,18,21,24,27,30,33,36,39)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
all_insulin_wide$AUC.G16.7 <- apply(all_insulin_wide[,c(5:14)], 1, calculate.auc.HG)

#IBMX
calculate.auc.IBMX <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(63,66,69,72,75,78,81,84,87)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
all_insulin_wide$AUC.IBMX <- apply(all_insulin_wide[,c(22:30)], 1, calculate.auc.IBMX)

#KCl
calculate.auc.KCl <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(123,126,129,132,135,138)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
all_insulin_wide$AUC.KCl <- apply(all_insulin_wide[,c(42:47)], 1, calculate.auc.KCl)

#calculate mean baseline secretion
all_insulin_wide <- all_insulin_wide %>%
  mutate(baselineins = rowMeans(across(c(`3`,`6`,`9`))))

#add metadata
all_insulin_wide <- left_join(all_insulin_wide, donors, by = c("DonorID" = "donor_ID"))

## Fig 8H ## ----
stimulus_data <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123),
  xmax = c(12,42,63,72,72,93,102,102,150,132),
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(7,8,7.5,8,9,7.5,8,9,7.5,8),    
  ymax = c(8,8.5,8,9,9.5,8,9,9.5,8,8.5))

all_insulin_df %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, simplified_diagnosis, sex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.ng.100IEQs.min, na.rm = TRUE),
    sd = sd(Insulin.ng.100IEQs.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab ("Insulin release (ng/100 IEQs)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  scale_x_continuous(limits = c(3, 150)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(3,12,42,63,63,72,93,93,102,123), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig8/Fig8H.png", width = 8, height = 4)

## Fig 8I ## ----
hist(all_insulin_wide$baselineins)
hist(log(all_insulin_wide$baselineins))
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(baselineins) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(baselineins) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(baselineins), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion\n(ng/100 IEQs))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8I.png", width = 4, height = 4)

## Fig 8J ## ----
hist(all_insulin_wide$AUC.G16.7)
hist(log(all_insulin_wide$AUC.G16.7))
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.G16.7) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.G16.7) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0252 for ctrl vs T2D in females, 0.0154 in males

all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.G16.7), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.03))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8J.png", width = 4, height = 4)

## Fig 8K ## ----
hist(all_insulin_wide$AUC.IBMX)
hist(log(all_insulin_wide$AUC.IBMX))
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #0.0532 for females, 0.0021 for males.

all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.9, label = "p=0.053", label.size = 4, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.18))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8K.png", width = 4, height = 4)

## Fig 8L ## ----
hist(all_insulin_wide$AUC.KCl)
hist(log(all_insulin_wide$AUC.KCl))
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis) #sig for disease and interaction
model <- all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #, p=0.0828 for F vs M in T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0664 for ctrl vs T2D in females, p=0.0001 in males

all_insulin_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.22))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.8, label = "p=0.07", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 4.6, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.1, 0.07))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0.5,5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8L.png", width = 4, height = 4)

## Fig S5H ## ----
stimulus_data2 <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123), 
  xmax = c(12,42,63,72,72,93,102,102,150,132),
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(7,8,7.5,8,9,7.5,8,9,7.5,8), 
  ymax = c(8,8.5,8,9,9.5,8,9,9.5,8,8.5))

all_insulin_df %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ungroup() %>%
  group_by(Time.min, sex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.ng.100IEQs.min, na.rm = TRUE),
    sd = sd(Insulin.ng.100IEQs.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = sex, group = sex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("Insulin release (ng/100 IEQs)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  scale_x_continuous(limits = c(3, 150)) +  
  geom_vline(xintercept = c(3,12,42,63,63,72,93,93,102,123), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS5/FigS5H.png", width = 5.5, height = 4)

## Fig S5I ## ----
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(baselineins) ~ age_years + sex) #ns

all_insulin_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(baselineins), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline insulin secretion\n(ng/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS5/FigS5I.png", width = 3, height = 3)

## Fig S5J ## ----
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.G16.7) ~ age_years + sex) #ns

all_insulin_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.G16.7), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS5/FigS5J.png", width = 3, height = 3)

## Fig S5K ## ----
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex) #ns

all_insulin_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS5/FigS5K.png", width = 3, height = 3)

## Fig S5L ## ----
all_insulin_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.KCl) ~ age_years + sex) #ns

all_insulin_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS5/FigS5L.png", width = 3, height = 3)

#Vanderbilt glucagon perifusion data
all_glucagon <- list()
for (i in 1:97){
  donor <- donor_list[i]
  if (file.exists(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
    if ("Glucagon Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx")) == TRUE){
      data <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion-data.xlsx"), sheet = "Glucagon Results", range = "A1:F52")
      data <- data[2:51,] #remove first row which is NA because it's merged with the header row
      data$Stimulus[1] <- data$Stimulus[2] #First stimulus entry has donor info so let's remove that
      colnames(data) <- c("Fraction","Time.min","Stimulus","Glucagon.pg.mL","Glucagon.pg.100IEQs.min","Glucagon.pc.content.min")
      data$DonorID <- rep(donor,nrow(data))} else{
        data <- NA
      }
  } else{
    if ("Glucagon Results" %in% excel_sheets(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx")) == TRUE){
      data <- read_excel(paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor ,"_Perifusion_data.xlsx"), sheet = "Glucagon Results", range = "A1:F52")
      data <- data[2:51,] #remove first row which is NA because it's merged with the header row
      data$Stimulus[1] <- data$Stimulus[2] #First stimulus entry has donor info so let's remove that
      colnames(data) <- c("Fraction","Time.min","Stimulus","Glucagon.pg.mL","Glucagon.pg.100IEQs.min","Glucagon.pc.content.min")
      data$DonorID <- rep(donor,nrow(data))} else{
        data <- NA
      }
  }
  all_glucagon[[i]] <- data
}

which(is.na(all_glucagon)) #none
all_glucagon_df <- data.frame(do.call(rbind, all_glucagon))
head(all_glucagon_df)

#add metadata
all_glucagon_df <- inner_join(all_glucagon_df, donors, by = c("DonorID" = "donor_ID")) 
summary(all_glucagon_df)

#Set up to calculte AUCs
all_glucagon_wide <- all_glucagon_df %>%
  dplyr::select(-Glucagon.pg.mL) %>%
  dplyr::select(-Glucagon.pc.content.min) %>%
  group_by(DonorID, Time.min) %>%
  summarise(Glucagon.pg.100IEQs.min = mean(Glucagon.pg.100IEQs.min)) %>% #to deal with duplicate values for a single time point for a single donor
  pivot_wider(values_from = Glucagon.pg.100IEQs.min, names_from = Time.min)
head(all_glucagon_wide) 

#Calculate AUCs
#LG plus epiniphrine
calculate.auc.epi <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(93,96,99,102,105,108,111)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
all_glucagon_wide$AUC.Epi <- apply(all_glucagon_wide[,c(32:38)], 1, calculate.auc.epi)

#IBMX
all_glucagon_wide$AUC.IBMX <- apply(all_glucagon_wide[,c(22:30)], 1, calculate.auc.IBMX)

#KCl
all_glucagon_wide$AUC.KCl <- apply(all_glucagon_wide[,c(42:47)], 1, calculate.auc.KCl)

#calculate mean baseline secretion
all_glucagon_wide <- all_glucagon_wide %>%
  mutate(baselinegcg = rowMeans(across(c(`3`,`6`,`9`))))

#add metadata
all_glucagon_wide <- left_join(all_glucagon_wide, donors, by = c("DonorID" = "donor_ID"))

## Fig S9F ## ----
stimulus_data3 <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123),   
  xmax = c(12,42,63,72,72,93,102,102,150,132),  
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(280,320,300,320,360,300,320,360,300,320),  
  ymax = c(320,340,320,360,380,320,360,380,320,340))

all_glucagon_df %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  filter(age_years > 14, age_years < 40) %>%
  ungroup() %>%
  group_by(Time.min, sex) %>%
  summarise(
    n = n(),
    mean = mean(Glucagon.pg.100IEQs.min, na.rm = TRUE),
    sd = sd(Glucagon.pg.100IEQs.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = sex, group = sex)) +
  #geom_point(size = 0.05) +
  geom_path(linewidth = 1)+
  #geom_smooth(aes(fill = simplified_diagnosis))+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("Glucagon release (pg/100 IEQs)") +
  xlab ("Time (min)") +
  scale_x_continuous(limits = c(3, 150)) +
  scale_y_continuous(limits = c(-40, 390)) +
  geom_vline(xintercept = c(3,12,42,63,63,72,93,93,102,123), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",linewidth = 0.5,
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS9/FigS9F.png", width = 5.5, height = 4)

## Fig S9G ## ----
hist(all_glucagon_wide$baselinegcg)
hist(log(all_glucagon_wide$baselinegcg))
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(baselinegcg) ~ age_years + sex) #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(baselinegcg), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline glucagon\nsecretion (pg/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS9/FigS9G.png", width = 3, height = 3)

## Fig S9H ## ----
hist(all_glucagon_wide$AUC.Epi)
hist(log(all_glucagon_wide$AUC.Epi))
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.Epi) ~ age_years + sex) #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.Epi), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n1.7 mM glucose + Epi)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS9/FigS9H.png", width = 3, height = 3)

## Fig S9I ## ----
hist(all_glucagon_wide$AUC.IBMX)
hist(log(all_glucagon_wide$AUC.IBMX))
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex) #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS9/FigS9I.png", width = 3, height = 3)

## Fig S9J ## ----
hist(all_glucagon_wide$AUC.KCl)
hist(log(all_glucagon_wide$AUC.KCl))
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.KCl) ~ age_years + sex) #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS9/FigS9J.png", width = 3, height = 3)

## Fig S10F ## ----
stimulus_data4 <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123),         
  xmax = c(12,42,63,72,72,93,102,102,150,132), 
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(280,320,300,320,360,300,320,360,300,320), 
  ymax = c(320,340,320,360,380,320,360,380,320,340))

all_glucagon_df %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, simplified_diagnosis, sex) %>%
  summarise(
    n = n(),
    mean = mean(Glucagon.pg.100IEQs.min, na.rm = TRUE),
    sd = sd(Glucagon.pg.100IEQs.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = simplified_diagnosis, group = simplified_diagnosis)) +
  #geom_point(size = 0.05) +
  geom_path(linewidth = 1)+
  #geom_smooth(aes(fill = simplified_diagnosis))+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = simplified_diagnosis)))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab ("Glucagon release (pg/100 IEQs)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  scale_x_continuous(limits = c(3, 150)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(3,12,42,63,63,72,93,93,102,123), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data4,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS10/FigS10F.png", width = 8, height = 4)

## Fig S10G ## ----
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(baselinegcg) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(baselinegcg) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(baselinegcg), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion\n(pg/100 IEQs))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS10/FigS10G.png", width = 4, height = 4)


## Fig S10I ## ----
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.Epi) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.Epi) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.Epi), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n1.7 mM glucose + Epi)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS10/FigS10I.png", width = 4, height = 4)

## Fig S10H ## ----
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS10/FigS10H.png", width = 4, height = 4)

## Fig S10J ## ----
all_glucagon_wide %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS10/FigS10J.png", width = 4, height = 4)

#Vanderbilt glucagon content data
#Load data
glucagon.content <- list()
for (i in 1:length(donor_list)) {
  donor <- donor_list[i]
  
  file_path1 <- paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor, "_Perifusion-data.xlsx")
  file_path2 <- paste0("data/HPAP/Vanderbilt Perifusion/", donor, "/Islet Studies/Islet physiology studies/Islets to Vanderbilt/Islet perifusion at Vanderbilt/", donor, "_Perifusion_data.xlsx")
  
  if (file.exists(file_path1)) {
    file_path <- file_path1
  } else if (file.exists(file_path2)) {
    file_path <- file_path2
  } else {
    glucagon.content[[i]] <- NA
    next
  }
  
  if ("Glucagon Results" %in% excel_sheets(file_path)) {
    data <- read_excel(file_path, sheet = "Glucagon Results", range = "D54:F55")
    colnames(data) <- c("Glucagon.pg.mL", "Glucagon.content.pg", "Glucagon.content.pg.IEQ")
    data$DonorID <- rep(donor, nrow(data))
  } else {
    data <- NA
  }
  
  glucagon.content[[i]] <- data
}

which(is.na(glucagon.content)) #none

glucagon.content_df <- data.frame(do.call(rbind, glucagon.content))
head(glucagon.content_df)

#Add metadata
glucagon.content_df <- inner_join(glucagon.content_df, donors, by = c("DonorID" = "donor_ID")) 

## Fig S9L ## ----
hist(glucagon.content_df$Glucagon.content.pg.IEQ)
hist(log(glucagon.content_df$Glucagon.content.pg.IEQ))
glucagon.content_df %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(age_years >14, age_years < 40) %>%
  anova_test(log(Glucagon.content.pg.IEQ) ~ age_years + sex) #ns

glucagon.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(age_years >14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(Glucagon.content.pg.IEQ), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(glucagon content\n(pg/IEQ))") +
  xlab("") +
  labs(colour = "Condition", fill = "Condition") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line()  
  )
ggsave("Output/FigS9/FigS9L.png", width = 3, height = 3)

## Fig S10L ## ----
glucagon.content_df %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(Glucagon.content.pg.IEQ) ~ age_years + sex*simplified_diagnosis) #ns
model <- glucagon.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(Glucagon.content.pg.IEQ) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

glucagon.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Glucagon.content.pg.IEQ), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = simplified_diagnosis, colour = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/IEQ))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(3,10)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS10/FigS10L.png", width = 5, height = 4)

#Humanislets.com glucose perifusion data
#clear environment
rm(list = ls())

#read data and metadata
peri_gluc <- read.csv("data/Humanislets.com/peri_gluc.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")
donor_data <- donor_data %>%
  mutate(diagnosis = case_when(
    diagnosis == "None" ~ "Control",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "Type2" ~ "T2D"
  ))

#summarise per donor
peri_gluc_means <- peri_gluc %>%
  group_by(record_id) %>%
  summarise_all(mean)

#turn to long format
peri_gluc_long <- peri_gluc_means %>%
  dplyr::select(-replicate) %>%
  pivot_longer(!record_id, names_to = "Time.min", values_to = "Insulin.release")
peri_gluc_long$Time.min <- gsub("time_","",peri_gluc_long$Time.min)
peri_gluc_long$Time.min <- as.numeric(peri_gluc_long$Time.min)
head(peri_gluc_long)

#Add metadata
peri_gluc_long <- inner_join(peri_gluc_long, donor_data, by = "record_id")

#Calculate AUCs
calculate.auc.G15 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(20,22.5,25,27.5,30,32.5,35,40,45,50,55,60,65)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_gluc_means$AUC.G15 <- apply(peri_gluc_means[,c(7:19)], 1, calculate.auc.G15)

calculate.auc.G6 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(85,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,122.5,125,127.5,130,135,140)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_gluc_means$AUC.G6 <- apply(peri_gluc_means[,c(23:42)], 1, calculate.auc.G6)

calculate.auc.KCl <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(160,162.5,165,167.5,170,175,180,185,190)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_gluc_means$AUC.KCl <- apply(peri_gluc_means[,c(46:54)], 1, calculate.auc.KCl)

peri_gluc_means <- peri_gluc_means %>%
  mutate(baseline.mean = rowMeans(across(c(time_5, time_10, time_15, time_20)), na.rm = TRUE))

#add metadata
peri_gluc_means <- left_join(peri_gluc_means, donor_data, by = "record_id")

## Fig 8M ## ----
stimulus_data <- data.frame(
  xmin = c(0,20,90,160), 
  xmax = c(20,70,140,190), 
  label = c("G 3","G 15","G 6","KCl"),
  ymin = c(630,630,630,630),  
  ymax = c(660,660,660,660))

peri_gluc_long %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = diagnosis, group = diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = diagnosis)))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  #coord_cartesian(xlim = c(0, 220), ylim = c(0,520)) +  
  facet_wrap(~donorsex)+
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/Fig8/Fig8M.png", width = 8, height = 4)

## Fig 8N ## ----
hist(peri_gluc_means$baseline.mean)
hist(log(peri_gluc_means$baseline.mean))
peri_gluc_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex*diagnosis) #sig for sex and diagnosis
model <- peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(baseline.mean) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0419 for F-M in controls
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.75, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8N.png", width = 4, height = 4)

## Fig 8O ## ----
hist(peri_gluc_means$AUC.G15)
hist(log(peri_gluc_means$AUC.G15))
peri_gluc_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.G15) ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.G15) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0597 for Ctrl vs T2D in females, p=0.0002 for males

peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.G15), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 10.3, label = "p=0.06", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 11, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n15 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(4,12)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8O.png", width = 4, height = 4)

## Fig 8P ## ----
hist(peri_gluc_means$AUC.G6)
hist(log(peri_gluc_means$AUC.G6))
peri_gluc_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.G6) ~ donorage + donorsex*diagnosis) #sig for sex
model <- peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.G6) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #o=0.0412 for F vs M in controls, 0.0829 in T2D
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.G6), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 10.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n6 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(6,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8P.png", width = 4, height = 4)

## Fig 8Q ## ----
hist(peri_gluc_means$AUC.KCl)
hist(log(peri_gluc_means$AUC.KCl))
peri_gluc_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex*diagnosis) #ns
model <- peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Fig8/Fig8Q.png", width = 4, height = 4)

## Fig S5M ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190),
  label = c("G 3","G 15","G 6","KCl"),
  ymin = c(rep(1000,4)),  
  ymax = c(rep(1100,4)))

peri_gluc_long %>%
  filter(diagnosis %in% c("Control")) %>%
  filter(donorage > 14, donorage < 40) %>%
  ungroup() %>%
  group_by(Time.min, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS5/FigS5M.png", width = 5.5, height = 4)

## Fig S5N ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(baseline.mean ~ donorage + donorsex) #ns

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets))") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/FigS5/FigS5N.png", width = 3, height = 3)

## Fig S5O ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.G15 ~ donorage + donorsex) #p=0.029 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.G15),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 10.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.25, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n15 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(6.5,11))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/FigS5/FigS5O.png", width = 3, height = 3)

## Fig S5P ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.G6 ~ donorage + donorsex) #p=0.031 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.G6),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 10.1, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(6.5,10.5))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/FigS5/FigS5P.png", width = 3, height = 3)

## Fig S5Q ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.KCl ~ donorage + donorsex) #p=0.043 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 9.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.25, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(6.5,10))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/FigS5/FigS5Q.png", width = 3, height = 3)

## Fig S6A ## ----
peri_gluc_means %>% 
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=bodymassindex, y = log(baseline.mean)))+
  geom_point(aes(colour = diagnosis))+
  xlab("BMI")+
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets)) ")+
  labs(colour = "Condition")+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  stat_smooth(method = "lm", colour = "grey50") +
  stat_cor(label.y = 6)+
  theme_classic()
ggsave("Output/FigS6/FigS6A.png", width = 5, height = 4)

## Fig S6B ## ----
peri_gluc_means %>% 
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=predistributionculturetime/24, y = log(baseline.mean)))+
  geom_point(aes(colour = diagnosis))+
  xlab("Culture duration (days)")+
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets)) ")+
  labs(colour = "Condition")+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  stat_smooth(method = "lm", colour = "grey50") +
  stat_cor(label.y = 6)+
  theme_classic()
ggsave("Output/FigS6/FigS6B.png", width = 5, height = 4)

## Fig S6C ## ----
peri_gluc_means %>% 
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=pdisletparticleindex, y = log(baseline.mean)))+
  geom_point(aes(colour = diagnosis))+
  xlab("Islet Particle Index")+
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets)) ")+
  labs(colour = "Condition")+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  stat_smooth(method = "lm", colour = "grey50") +
  stat_cor(label.y = 6)+
  theme_classic()
ggsave("Output/FigS6/FigS6C.png", width = 5, height = 4)

### Figures S7-8 ### ----
#clear environment
rm(list = ls())

#Humanislets.com metadata
donor_data <- read.csv("data/Humanislets.com/donor.csv")
donor_data <- donor_data %>%
  mutate(diagnosis = case_when(
    diagnosis == "None" ~ "Control",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "Type2" ~ "T2D"
  ))


#Humanislets.com leucine perifusion data
peri_leu <- read.csv("data/Humanislets.com/peri_leu.csv")

#summarise per donor
peri_leu_means <- peri_leu %>%
  group_by(record_id) %>%
  summarise_all(mean)

#turn to long format
peri_leu_long <- peri_leu_means %>%
  dplyr::select(-replicate) %>%
  pivot_longer(!record_id, names_to = "Time.min", values_to = "Insulin.release")
peri_leu_long$Time.min <- gsub("time_","",peri_leu_long$Time.min)
peri_leu_long$Time.min <- as.numeric(peri_leu_long$Time.min)
head(peri_leu_long)

#Add metadata
peri_leu_long <- inner_join(peri_leu_long, donor_data, by = "record_id")

#Calculate AUCs
calculate.auc.leu5 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(20,22.5,25,27.5,30,32.5,35,40,45,50,55,60,65,70)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_leu_means$AUC.leu5 <- apply(peri_leu_means[,c(7:20)], 1, calculate.auc.leu5)

calculate.auc.leu5glc6 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,122.5,125,127.5,130,135,140)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_leu_means$AUC.leu5glc6 <- apply(peri_leu_means[,c(24:42)], 1, calculate.auc.leu5glc6)

calculate.auc.KCl <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(160,162.5,165,167.5,170,175,180,185,190)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_leu_means$AUC.KCl <- apply(peri_leu_means[,c(46:54)], 1, calculate.auc.KCl)

peri_leu_means <- peri_leu_means %>%
  mutate(baseline.mean = rowMeans(across(c(time_5, time_10, time_15, time_20)), na.rm = TRUE))

#add metadata
peri_leu_means <- left_join(peri_leu_means, donor_data, by = "record_id")

## Fig S7A ## ----
stimulus_data <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190),
  label = c("G 3","Leu 5","Leu 5 + G 6","KCl"),
  ymin = c(740,740,740,740),
  ymax = c(820,820,820,820))

peri_leu_long %>%
  filter(diagnosis %in% c("Control")) %>%
  filter(donorage > 14, donorage < 40) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS7/FigS7A.png", width = 5.5, height = 4)

## Fig S7B ## ----
hist(peri_leu_means$baseline.mean)
hist(log(peri_leu_means$baseline.mean))
peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex) #p=0.004 for sex

peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(mean baseline insulin\nsecretion (\U03BCU/mL/65 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(-0,6.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7B.png", width = 3, height = 3)

## Fig S7C ## ----
hist(peri_leu_means$AUC.leu5)
hist(log(peri_leu_means$AUC.leu5))
peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.leu5) ~ donorage + donorsex) #p=0.004 for sex

peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.leu5), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 10.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.35, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n5 mM leucine)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7C.png", width = 3, height = 3)

## Fig S7D ## ----
hist(peri_leu_means$AUC.leu5glc6)
hist(log(peri_leu_means$AUC.leu5glc6))
peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.leu5glc6) ~ donorage + donorsex) #p=0.002 for sex

peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.leu5glc6), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 10.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.35, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 5 mM\nleucine + 6 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7D.png", width = 3, height = 3)

## Fig S7E ## ----
hist(peri_leu_means$AUC.KCl)
hist(log(peri_leu_means$AUC.KCl))
peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex) #p=0.045 for sex

peri_leu_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.75, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl) ") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(6,10.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7E.png", width = 3, height = 3)

## Fig S8A ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190), 
  label = c("G 3","Leu 5","Leu 5 + G 6","KCl"),
  ymin = c(600,600,600,600),
  ymax = c(630,630,630,630))

peri_leu_long %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = diagnosis, group = diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = diagnosis)))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  facet_wrap(~donorsex) +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS8/FigS8A.png", width = 8, height = 4)

## Fig S8B ## ----
peri_leu_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex*diagnosis) #sig for sex
model <- peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(baseline.mean) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0187 for F-M in controls
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.17, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,6.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8B.png", width = 5, height = 4)

## Fig S8C ## ----
peri_leu_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.leu5) ~ donorage + donorsex*diagnosis) #ns
model <- peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.leu5) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0355 for F-M in controls
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.leu5), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 10.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 5 mM leucine)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8C.png", width = 5, height = 4)

## Fig S8D ## ----
peri_leu_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.leu5glc6) ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.leu5glc6) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0045 for Ctrl vs T2D in males

peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.leu5glc6), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 10.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.17))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 5 mM leucine\n+ 6 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8D.png", width = 5, height = 4)

## Fig S8E ## ----
peri_leu_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0233 for Ctrl vs T2D in males

peri_leu_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 9.9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.07))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 30 mM KCl) ") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8E.png", width = 5, height = 4)

#Humanislets.com oleate/palmitate perifusion
#metadata
donor_data <- read.csv("data/Humanislets.com/donor.csv")
donor_data <- donor_data %>%
  mutate(diagnosis = case_when(
    diagnosis == "None" ~ "Control",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "Type2" ~ "T2D"
  ))

#Humanislets.com olp perifusion data
peri_olp <- read.csv("data/Humanislets.com/peri_olp.csv")

#summarise per donor
peri_olp_means <- peri_olp %>%
  group_by(record_id) %>%
  summarise_all(mean)

#turn to long format
peri_olp_long <- peri_olp_means %>%
  dplyr::select(-replicate) %>%
  pivot_longer(!record_id, names_to = "Time.min", values_to = "Insulin.release")
peri_olp_long$Time.min <- gsub("time_","",peri_olp_long$Time.min)
peri_olp_long$Time.min <- as.numeric(peri_olp_long$Time.min)
head(peri_olp_long)

#Add metadata
peri_olp_long <- inner_join(peri_olp_long, donor_data, by = "record_id")

#Calculate AUCs
calculate.auc.olp1.5 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(20,22.5,25,27.5,30,32.5,35,40,45,50,55,60,65,70)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_olp_means$AUC.olp1.5 <- apply(peri_olp_means[,c(7:20)], 1, calculate.auc.olp1.5)

calculate.auc.olp1.5glc6 <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,122.5,125,127.5,130,135,140)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_olp_means$AUC.olp1.5glc6 <- apply(peri_olp_means[,c(24:42)], 1, calculate.auc.olp1.5glc6)

calculate.auc.KCl <- function(x){
  if(sum(!is.na(x)) >= 2) {
    x_interpolated <- na.approx(x, na.rm = FALSE)
    time_points <- c(160,162.5,165,167.5,170,175,180,185,190)
    auc_value <- auc(time_points, x_interpolated, type = 'linear')
  } else {
    auc_value <- NA
  }
  return(auc_value)
}
peri_olp_means$AUC.KCl <- apply(peri_olp_means[,c(46:54)], 1, calculate.auc.KCl)

peri_olp_means <- peri_olp_means %>%
  mutate(baseline.mean = rowMeans(across(c(time_5, time_10, time_15, time_20)), na.rm = TRUE))

#add metadata
peri_olp_means <- left_join(peri_olp_means, donor_data, by = "record_id")

## Fig S7F ## ----
stimulus_data <- data.frame(
  xmin = c(0,20,90,160), 
  xmax = c(20,70,140,190), 
  label = c("G 3","OLP 1.5","OLP 1.5 + G 6","KCl"),
  ymin = c(rep(650,4)),     
  ymax = c(rep(700,4)))

peri_olp_long %>%
  filter(diagnosis %in% c("Control")) %>%
  filter(donorage > 14, donorage < 40) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS7/FigS7F.png", width = 5.5, height = 4)

## Fig S7G ## ----
hist(peri_olp_means$baseline.mean)
hist(log(peri_olp_means$baseline.mean))
peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex) #p=0.002 for sex

peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(mean baseline insulin\nsecretion (\U03BCU/mL/65 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(1,6.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7G.png", width = 3, height = 3)

## Fig S7H ## ----
hist(peri_olp_means$AUC.olp1.5)
hist(log(peri_olp_means$AUC.olp1.5))
peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.olp1.5) ~ donorage + donorsex) #p=0.002 for sex

peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.olp1.5), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 10, label = "*", label.size = 7, size = 0.5, tip.length = c(0.4, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7H.png", width = 3, height = 3)

## Fig S7I ## ----
hist(peri_olp_means$AUC.olp1.5glc6)
hist(log(peri_olp_means$AUC.olp1.5glc6))
peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.olp1.5glc6) ~ donorage + donorsex) #p=0.001 for sex

peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.olp1.5glc6), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 10.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate\n+ 6 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7I.png", width = 3, height = 3)

## Fig S7J ## ----
hist(peri_olp_means$AUC.KCl)
hist(log(peri_olp_means$AUC.KCl))
peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex) #p=0.012 for sex

peri_olp_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, aes(colour = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.27, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl) ") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(6,10.2)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/FigS7/FigS7J.png", width = 3, height = 3)

## Fig S8F ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190), 
  label = c("G 3","OLP 1.5","OLP 1.5 + G 6","KCl"),
  ymin = rep(500,4),
  ymax = rep(530,4))

peri_olp_long %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release, na.rm = TRUE),
    sd = sd(Insulin.release, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = diagnosis, group = diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = diagnosis)))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab ("Insulin Release (\U03BCU/mL/65 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  facet_wrap(~donorsex) +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data2,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Output/FigS8/FigS8F.png", width = 8, height = 4)

## Fig S8G ## ----
peri_olp_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex*diagnosis) #sig for sex
model <- peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(baseline.mean) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0117 for F-M in controls
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.18, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(\U03BCU/mL/65 islets))") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(0,6.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8G.png", width = 5, height = 4)

## Fig S8H ## ----
peri_olp_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.olp1.5) ~ donorage + donorsex*diagnosis) #sig for sex and diagnosis
model <- peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.olp1.5) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0275 for F-M in controls
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0784 for Ctrl vs T2D in females

peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.olp1.5), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 10, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5,11)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8H.png", width = 5, height = 4)

## Fig S8I ## ----
peri_olp_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.olp1.5glc6) ~ donorage + donorsex*diagnosis) #sig for sex and diagnosis
model <- peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.olp1.5glc6) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0069 for F vs M in Control
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.018 for Ctrl vs T2D in females, p=0.0483 in males

peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.olp1.5glc6), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 10, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 10.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.4))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 11, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate\n+ 6 mM glucose)") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  scale_y_continuous(limits = c(5.5,12)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8I.png", width = 5, height = 4)

## Fig S8J ## ----
peri_olp_means %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex*diagnosis) #ns
model <- peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #p=0.0768 for Control vs T2D in females

peri_olp_means %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 1, aes(group = diagnosis, colour = diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 30 mM KCl) ") +
  xlab("")+
  labs(colour = "Condition", fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/FigS8/FigS8J.png", width = 5, height = 4)

### Figure S1 ### ---
#clear environment
rm(list = ls())

#Humanislets.com metadata
HI_metadata <- read.csv("data/Humanislets.com/donor.csv")
HI_metadata <- HI_metadata %>%
  mutate(diagnosis = case_when(
    diagnosis == "None" ~ "Control",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "Type2" ~ "T2D"
  ))

#HPAP metadata
HPAP_metadata <- read_excel("data/HPAP/Donor_Summary_192.xlsx",sheet = "donor")
HPAP_metadata <- HPAP_metadata %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#HPAP culture time metadata
HPAP_culturetime <- read_excel("data/HPAP/PancDB_Donors.xlsx")
#convert excel dates to R dates
HPAP_culturetime <- HPAP_culturetime %>%
  mutate(CultureHarvestDate = case_when(
    grepl("^est", CultureHarvestDate) == TRUE ~ mdy(gsub("est ","", CultureHarvestDate)),
    grepl("^est", CultureHarvestDate) == FALSE ~ as.Date(as.numeric(CultureHarvestDate), origin = "1899-12-30")
  ))
HPAP_culturetime <- HPAP_culturetime %>%
  mutate(CultureStartDate = as.numeric(CultureStartDate)) %>%
  mutate(CultureStartDate = as.Date(CultureStartDate, origin = "1899-12-30"))

#check harvest is after start date
all(HPAP_culturetime$CultureHarvestDate > HPAP_culturetime$CultureStartDate) #FALSE

#which ones are problematic?
HPAP_culturetime[which(HPAP_culturetime$CultureHarvestDate < HPAP_culturetime$CultureStartDate),9:10]
#there is a start date of 27 may but harvest date 20 may. Probably they flipped it. Use absolute values to get culture time

#get culture duration
HPAP_culturetime$culturetime <- abs(HPAP_culturetime$CultureHarvestDate - HPAP_culturetime$CultureStartDate)
hist(as.numeric(HPAP_culturetime$culturetime)) #one that's 3000+, doesn't make sense

#which ones are problemtic?
HPAP_culturetime[which(HPAP_culturetime$culturetime>30),9:10]
#It seems there was a typo in one of the months as start is August and harvest is September
#And in another one, 2010 is a typo, should be 2019.

#go back and fix the problematic ones
HPAP_culturetime <- read_excel("data/HPAP/PancDB_Donors.xlsx")
HPAP_culturetime$CultureHarvestDate[16] <- "est 8/17/17" #fixing typo to change harvest date to August in line with start date
HPAP_culturetime$CultureStartDate[45] <- "43679" #fixing typo to change 2010 to 2019

HPAP_culturetime <- HPAP_culturetime %>%
  mutate(CultureHarvestDate = case_when(
    grepl("^est", CultureHarvestDate) == TRUE ~ mdy(gsub("est ","", CultureHarvestDate)),
    grepl("^est", CultureHarvestDate) == FALSE ~ as.Date(as.numeric(CultureHarvestDate), origin = "1899-12-30")
  ))
HPAP_culturetime <- HPAP_culturetime %>%
  mutate(CultureStartDate = as.numeric(CultureStartDate)) %>%
  mutate(CultureStartDate = as.Date(CultureStartDate, origin = "1899-12-30"))

#get culture duration
HPAP_culturetime$culturetime <- abs(HPAP_culturetime$CultureHarvestDate - HPAP_culturetime$CultureStartDate)
hist(as.numeric(HPAP_culturetime$culturetime)) #now it makes sense, distribution 1 to 11 days

#change from difftime to numeric
HPAP_culturetime$culturetime <- as.numeric(HPAP_culturetime$culturetime)

#Join HPAP and Humanislets.com metadata into one dataframe
HPAP_metadataselect <- HPAP_metadata %>%
  filter(simplified_diagnosis != "T1D") %>%
  dplyr::select(donor_ID, sex, simplified_diagnosis, age_years, hba1c, bmi)
HPAP_metadataselect <- HPAP_culturetime %>%
  dplyr::select(DonorID, culturetime) %>%
  right_join(HPAP_metadataselect, by=c("DonorID" = "donor_ID"))

HI_metadataselect <- HI_metadata %>%
  filter(diagnosis != "T1D") %>%
  dplyr::select(record_id, predistributionculturetime, donorsex, diagnosis, donorage, hba1c, bodymassindex) %>%
  mutate(predistributionculturetime = predistributionculturetime/24) #convert from hours to days

head(HPAP_metadataselect)
head(HI_metadataselect)

#make colnames consistent
colnames(HI_metadataselect) <- colnames(HPAP_metadataselect)

#add data source
HPAP_metadataselect$dataset <- rep("HPAP", nrow(HPAP_metadataselect))
HI_metadataselect$dataset <- rep("Humanislets.com", nrow(HI_metadataselect))

#combine
metadata_bothdatasets <- rbind(HPAP_metadataselect, HI_metadataselect)

## Fig S1A ## ----
age <- metadata_bothdatasets %>%
  pivot_wider(names_from = dataset, values_from = age_years)
ks.test(age$HPAP, age$Humanislets.com) #p=1.063e-08
age %>%
  ggplot(aes(x=x)) +
  geom_histogram(bins = 50, aes(x = HPAP, y = after_stat(count)), fill="black" ) +
  annotate("label", x=75, y=5, label="HPAP", color="black") +
  geom_histogram(bins = 50, aes(x = Humanislets.com, y = -after_stat(count)), fill= "#3FA147") +
  annotate("label", x=20, y=-15, label="Humanislets.com", color="#3FA147")+
  geom_hline(yintercept = 0, colour = "grey")+
  xlab("Age (years)")+
  ylab("Frequency")+
  facet_grid(sex~simplified_diagnosis)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("Output/FigS1/FigS1A.png", width = 6, height = 4)

## FIg S1B ## ----
bmi <- metadata_bothdatasets %>%
  pivot_wider(names_from = dataset, values_from = bmi)
ks.test(bmi$HPAP[!is.na(bmi$HPAP)], bmi$Humanislets.com[!is.na(bmi$Humanislets.com)]) #p=4.203e-05
bmi %>%
  ggplot(aes(x=x) ) +
  geom_histogram(bins = 50, aes(x = HPAP, y = after_stat(count)), fill="black" ) +
  annotate("label", x=16, y=8, label="HPAP", color="black") +
  geom_histogram(bins = 50, aes(x = Humanislets.com, y = -after_stat(count)), fill= "#3FA147") +
  annotate("label", x=35, y=-24, label="Humanislets.com", color="#3FA147")+
  geom_hline(yintercept = 0, colour = "grey")+
  xlab("BMI")+
  ylab("Frequency")+
  ylim(-27,10)+
  facet_grid(sex~simplified_diagnosis)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("Output/FigS1/FigS1B.png", width = 6, height = 4)

## Fig S1C ## ----
hba1c <- metadata_bothdatasets %>%
  pivot_wider(names_from = dataset, values_from = hba1c)
ks.test(hba1c$HPAP[!is.na(hba1c$HPAP)], hba1c$Humanislets.com[!is.na(hba1c$Humanislets.com)]) #p=0.003983
hba1c %>%
  ggplot(aes(x=x) ) +
  geom_histogram(bins = 50, aes(x = HPAP, y = after_stat(count)), fill="black" ) +
  annotate("label", x=3.75, y= 8, label="HPAP", color="black") +
  geom_histogram(bins = 50, aes(x = Humanislets.com, y = -after_stat(count)), fill= "#3FA147") +
  annotate("label", x=9, y=-15, label="Humanislets.com", color="#3FA147")+
  geom_hline(yintercept = 0, colour = "grey")+
  xlab("HbA1c")+
  ylab("Frequency")+
  facet_grid(sex~simplified_diagnosis)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("Output/FigS1/FigSC.png", width = 6, height = 4)

## Fig S1D ## ----
culturetime <- metadata_bothdatasets %>%
  pivot_wider(names_from = dataset, values_from = culturetime)
ks.test(culturetime$HPAP[!is.na(culturetime$HPAP)], culturetime$Humanislets.com[!is.na(culturetime$Humanislets.com)]) #p<2.2e-16
culturetime %>%
  ggplot(aes(x=x) ) +
  geom_histogram(bins = 50, aes(x = HPAP, y = after_stat(count)), fill="black" ) +
  annotate("label", x=1.5, y= 13, label="HPAP", color="black") +
  geom_histogram(bins = 50, aes(x = Humanislets.com, y = -after_stat(count)), fill= "#3FA147") +
  annotate("label", x=7, y=-15, label="Humanislets.com", color="#3FA147")+
  geom_hline(yintercept = 0, colour = "grey")+
  xlab("Days of culture")+
  ylab("Frequency")+
  ylim(-60,20)+
  facet_grid(sex~simplified_diagnosis)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("Output/FigS1/FigSD.png", width = 6, height = 4)

### Table S1 ### ----
#clear environment
rm(list = ls())

#read in sheets and merge
HIGO_CC <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_CC")
HIGO_BP <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_BP")
HIGO_MF <- read_excel("Output/Fig2/GSEA_p_bulk_RNAseq_youngcontrols.xlsx",sheet="GO_MF")
HIGO <- bind_rows(HIGO_CC, HIGO_BP, HIGO_MF)
HIGO$NES <- -HIGO$NES #reverse sign to change direction so positive is up in females
HIGO$dataset <- rep("Humanislets.com bulk RNAseq", nrow(HIGO))

HPAPbetascGO_CC <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_CC")
HPAPbetascGO_BP <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_BP")
HPAPbetascGO_MF <- read_excel("Output/Fig2/GSEA_p_beta_scRNAseq_youngcontrols.xlsx",sheet="GO_MF")
HPAPbetascGO <- bind_rows(HPAPbetascGO_CC, HPAPbetascGO_BP, HPAPbetascGO_MF)
HPAPbetascGO$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO))

HPAPalphascGO_CC <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_CC")
HPAPalphascGO_BP <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_BP")
HPAPalphascGO_MF <- read_excel("Output/Fig2/GSEA_p_alpha_scRNAseq_youngcontrols.xlsx",sheet="GO_MF")
HPAPalphascGO <- bind_rows(HPAPalphascGO_CC, HPAPalphascGO_BP, HPAPalphascGO_MF)
HPAPalphascGO$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO))

combined <- bind_rows(HIGO, HPAPbetascGO, HPAPalphascGO)
combined <- combined %>% arrange(p.adjust)
head(combined)

write.csv(combined, "Output/Table S1.csv")

### Table S2 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet$NES <- -sheet$NES #reverse sign to change direction so positive is up in females
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Table S2.csv")

### Table S3 ### ----
#clear environment
rm(list = ls())

#read in sheets and combine
HIGO_F_CC <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_CC")
HIGO_F_BP <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_BP")
HIGO_F_MF <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_FCtrlvsT2D.xlsx",sheet="GO_MF")
HIGO_F <- bind_rows(HIGO_F_CC, HIGO_F_BP, HIGO_F_MF)
HIGO_F$NES <- -HIGO_F$NES #make consistent with scRNAseq direction
HIGO_F$dataset <- rep("Humanislets.com bulk RNAseq", nrow(HIGO_F))


HPAPbetascGO_F_CC <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPbetascGO_F_BP <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPbetascGO_F_MF <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPbetascGO_F <- bind_rows(HPAPbetascGO_F_CC, HPAPbetascGO_F_BP, HPAPbetascGO_F_MF)
HPAPbetascGO_F$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO_F))

HPAPalphascGO_F_CC <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPalphascGO_F_BP <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPalphascGO_F_MF <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_F_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPalphascGO_F <- bind_rows(HPAPalphascGO_F_CC, HPAPalphascGO_F_BP, HPAPalphascGO_F_MF)
HPAPalphascGO_F$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO_F))

combined <- bind_rows(HIGO_F, HPAPbetascGO_F, HPAPalphascGO_F)
combined <- combined %>% arrange(p.adjust)
head(combined)

write.csv(combined, "Output/Table S3.csv")

### Table S4 ### ----
#clear environment
rm(list = ls())

#read in sheets and combine
HIGO_M_CC <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_CC")
HIGO_M_BP <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_BP")
HIGO_M_MF <- read_excel("Output/Fig3/GSEA_p_bulkRNAseq_MCtrlvsT2D.xlsx",sheet="GO_MF")
HIGO_M <- bind_rows(HIGO_M_CC, HIGO_M_BP, HIGO_M_MF)
HIGO_M$NES <- -HIGO_M$NES #make consistent with scRNAseq direction
HIGO_M$dataset <- rep("Humanislets.com bulk RNAseq", nrow(HIGO_M))

HPAPbetascGO_M_CC <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPbetascGO_M_BP <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPbetascGO_M_MF <- read_excel("Output/Fig3/GSEA_p_beta_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPbetascGO_M <- bind_rows(HPAPbetascGO_M_CC, HPAPbetascGO_M_BP, HPAPbetascGO_M_MF)
HPAPbetascGO_M$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO_M))

HPAPalphascGO_M_CC <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_CC")
HPAPalphascGO_M_BP <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_BP")
HPAPalphascGO_M_MF <- read_excel("Output/Fig3/GSEA_p_alpha_scRNAseq_M_ctrlvsT2D.xlsx",sheet="GO_MF")
HPAPalphascGO_M <- bind_rows(HPAPalphascGO_M_CC, HPAPalphascGO_M_BP, HPAPalphascGO_M_MF)
HPAPalphascGO_M$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO_M))


combined <- bind_rows(HIGO_M, HPAPbetascGO_M, HPAPalphascGO_M)
combined <- combined %>% arrange(p.adjust)
head(combined)

write.csv(combined, "Output/Table S4.csv")

### Table S5 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Fig4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Table S5.csv")

### Table S6 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Fig4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Table S6.csv")
