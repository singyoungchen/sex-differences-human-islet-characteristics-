#load packages
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(emmeans)
library(rstatix)
library(this.path)
library(edgeR)
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
library(MatchIt)
source('GOrge.R')

#set working directory to location of this R script
setwd(dirname(this.path::this.path()))

### Figure 1 and S2 ### ----
# HPAP cell proportions
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

donor_info <- donor_info %>%
  mutate(nonab.endocrine = delta.endocrine + epsilon.endocrine + pp.endocrine)
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(nonab.endocrine ~ age_years + sex) #sig for sex

## Fig 1A ## ----
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  ggplot(aes(x=sex, y = alpha.endocrine*100))+
  geom_boxplot(aes(fill = sex), alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("")+
  ylab("% alpha-cells") +
  ylim(0,105)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 100, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1A.tiff", width = 3, height = 3)

## Fig 1C ## ----
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  ggplot(aes(x=sex, y = beta.endocrine*100))+
  geom_boxplot(aes(fill = sex), alpha = 0.5)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 95, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("")+
  ylab("% beta-cells") +
  ylim(0,105)+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1C.tiff", width = 3, height = 3)

## Fig 1E ## ----
donor_info %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  ggplot(aes(x=sex, y = nonab.endocrine*100))+
  geom_boxplot(aes(fill = sex), alpha = 0.5)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 53, label = "*", label.size = 7, size = 0.5, tip.length = c(0.6, 0.02))+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("")+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  ylim(0,105) +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1E.tiff", width = 3, height = 3)

## Fig 1G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 93, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 93, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 102, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 117, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% alpha-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/Fig1/Fig1G.tiff", width = 5, height = 4)

## Fig 1H ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 85, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 94, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 103, label = "p=0.072", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 115, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% beta-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/Fig1/Fig1H beta.tiff", width = 5, height = 4)

## Fig 1I ## ----
donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(nonab.endocrine ~ age_years + sex*simplified_diagnosis) #sig for diagnosis
model <- donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(nonab.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #0.0833 for controls
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for females
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = nonab.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 77.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.25, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 55, label = "*", label.size = 7, size = 0.5, tip.length = c(0.25, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 85, label = "p=0.083", label.size = 4, size = 0.5, tip.length = c(0.25, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 100, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/Fig1/Fig1I.tiff", width = 5, height = 4)

## Fig S2A ## ----
donor_info_matched <- donor_info %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID)

donor_info_matched$Group <- as.logical(donor_info_matched$simplified_diagnosis == "T2D")
donor_info_matched <- matchit(Group ~ age_years + sex,
                              data = donor_info_matched,
                              method = 'nearest',
                              ratio = 1, 
                              exact = ~sex
) 
donor_info_matched <- match.data(donor_info_matched)

donor_info_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #ns
donor_info_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 68, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(25,75))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS2/FigS2A.tiff", width = 5, height = 4)

## Fig S2B ## ----
donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(alpha.endocrine ~ age_years + sex*simplified_diagnosis) #ns
model <- donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(alpha.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = alpha.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 93, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 93, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% alpha-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS2/FigS2B.tiff", width = 5, height = 4)

## Fig S2C ## ----
donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(beta.endocrine ~ age_years + sex*simplified_diagnosis) #ns
model <- donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(beta.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0731 for Control vs T2D in males
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = beta.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 86, label = "p=0.073", label.size = 4, size = 0.5, tip.length = c(0.02, 0.35))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 96, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% beta-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS2/FigS2C.tiff", width = 5, height = 4)

## Fig S2D ## ----
donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  anova_test(nonab.endocrine ~ age_years + sex*simplified_diagnosis) #0.055 for diagnosis
model <- donor_info_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  lm(nonab.endocrine ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #0.0538 for females
donor_info %>%
  filter(!donor_ID %in% no_cells$donor_ID) %>%
  filter(!donor_ID %in% low_cells$donor_ID) %>%
  filter(simplified_diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=sex, y = nonab.endocrine*100))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 80, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.35, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 55, label = "p=0.054", label.size = 4, size = 0.5, tip.length = c(0.3, 0.02))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  xlab("")+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS2/FigS2D.tiff", width = 5, height = 4)


# Humanislets.com cell proportions
#clear environment
rm(list = ls())

#import data - filed downloaded from Humanislets.com 12 February 2025
celltypeprop <- read.csv("data/Humanislets.com/cell_pro.csv")
donor_data <- read.csv("data/Humanislets.com/donor.csv")

#join metadata to data
celltypeprop <- left_join(celltypeprop, donor_data, by = "record_id")
unique(celltypeprop$diagnosis) #no Type1 donors in this dataset

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

celltypeprop <- celltypeprop %>%
  mutate(nonab_end = delta_end + gamma_end)
celltypeprop %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(simplified_diagnosis %in% c("Control")) %>%
  anova_test(nonab_end ~ donorage + donorsex) #ns

## Fig 1B ## ----
celltypeprop %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y = alpha_end*100))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("")+
  ylab("% alpha-cells") +
  ylim(0,105)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 40, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.2))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1B.tiff", width = 3, height = 3)

## Fig 1D ## ----
celltypeprop %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y = beta_end*100))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("")+
  ylab("% beta-cells") +
  ylim(0,105)+
  geom_bracket(xmin = 1, xmax = 2, y.position = 65, label = "*", label.size = 7, size = 0.5, tip.length = c(0.3, 0.1))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1D.tiff", width = 3, height = 3)

## Fig 1F ## ----
celltypeprop %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y = nonab_end*100))+
  geom_boxplot(aes(fill = donorsex), alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.3, 0.1))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  xlab("")+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  ylim(0,105)+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()       
  )
ggsave("Output/Final figures/Fig1/Fig1F.tiff", width = 3, height = 3)

## Fig 1J ## ----
celltypeprop %>%  
  anova_test(alpha_end ~ donorage + donorsex*diagnosis) #sig for disease and interaction
model <- celltypeprop %>%
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 50, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 47, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 55, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% alpha-cells") +
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/Fig1/Fig1J.tiff", width = 5, height = 4)

## Fig 1K ## ----
celltypeprop %>%  
  anova_test(beta_end ~ donorage + donorsex*diagnosis) #sig for disease and interaction
model <- celltypeprop %>%
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 100, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 72, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 85, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% beta-cells") +
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/Fig1/Fig1K.tiff", width = 5, height = 4)

## Fig 1L ## ----
celltypeprop %>%  
  anova_test(nonab_end ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- celltypeprop %>%
  lm(nonab_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ donorsex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for males
celltypeprop %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = nonab_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 25, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 40, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
   geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 50, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/Fig1/Fig1L.tiff", width = 5, height = 4)

## Fig S2E ## ----
unique(celltypeprop$diagnosis) #no Type1 donors in this dataset
celltypeprop$Group <- as.logical(celltypeprop$diagnosis == 'Type2')
celltypeprop_matched <- matchit(Group ~ donorage + donorsex,
                                data = celltypeprop,
                                method = 'nearest',  
                                ratio = 1, 
                                exact = ~donorsex
) 
celltypeprop_matched <- match.data(celltypeprop_matched)
head(celltypeprop_matched)
summary(celltypeprop_matched$subclass)

celltypeprop_matched %>%
  group_by(donorsex) %>%
  t_test(donorage ~ diagnosis) #ns
celltypeprop_matched %>%
  ggplot(aes(x=donorsex, y=donorage))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 68, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(35,75))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS2/FigS2E.tiff", width = 5, height = 4)

#celltype proportion data age-matched
celltypeprop_matched %>%  
  anova_test(alpha_end ~ donorage + donorsex*diagnosis) #sig for disease and interaction
model <- celltypeprop_matched %>%
  lm(alpha_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for Control vs T2D in females
celltypeprop_matched %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = alpha_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 45, label = "*", label.size = 7, size = 0.5, tip.length = c(0.4, 0.02))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 45, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.2, 0.22))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% alpha-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS2/FigS2F.tiff", width = 5, height = 4)

celltypeprop_matched %>%  
  anova_test(beta_end ~ donorage + donorsex*diagnosis) #sig for disease
model <- celltypeprop_matched %>%
  lm(beta_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for Control in females
celltypeprop_matched %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = beta_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.35))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.1))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% beta-cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS2/FigS2G.tiff", width = 5, height = 4)

celltypeprop_matched %>%  
  anova_test(nonab_end ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- celltypeprop_matched %>%
  lm(nonab_end ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for males
celltypeprop_matched %>%  
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control"
  )) %>%
  filter(diagnosis %in% c("Control", "T2D")) %>%
  ggplot(aes(x=donorsex, y = nonab_end*100))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 25, label = "*", label.size = 7, size = 0.5, tip.length = c(0.3, 0.1))+ 
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.1))+ 
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100))+
  ylab("% non-alpha, non-beta\nendocrine cells") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS2/FigS2H.tiff", width = 5, height = 4)


### Figure 2 ### ----
## Fig 2A ## ----

#clear environment
rm(list = ls())

#Pseudobulk scRNAseq from HPAP
dat <- readRDS("data/HPAP/T1D_T2D_20220428.rds")  #downloaded 7 January 2024
dat_bulk <- AggregateExpression(dat, group.by = c("hpap_id"), return.seurat = TRUE) #pseudobulk
dat_bulk_rawcounts_df <- data.frame(dat_bulk[["RNA"]]$counts) #retrieve raw counts
donor_HPAP <- read_excel("data/HPAP/Donor_Summary_192.xlsx")
donor_HPAP$donor_ID <- gsub("-","", donor_HPAP$donor_ID) #match donor IDs to those in the Seurat
donor_HPAP <- donor_HPAP %>%
  mutate(simplified_diagnosis = case_when(
    grepl("T2DM",clinical_diagnosis) == TRUE ~ "T2D",
    grepl("T1DM",clinical_diagnosis) == TRUE ~ "T1D",
    grepl("control", clinical_diagnosis) == TRUE ~ "Control"
  ))

#filter, remove anything with > 80% zeros
dat_bulk_rawcounts_df <- dat_bulk_rawcounts_df %>%
  tibble::rownames_to_column("Gene") %>%
  filter(rowSums(dplyr::select(., -Gene) == 0) / (ncol(.) - 1) <= 0.8) %>%
  tibble::column_to_rownames("Gene")

#logCPM transformation
dat_bulk_rawcounts_df <- DGEList(counts=dat_bulk_rawcounts_df)
dat_bulk_rawcounts_df <- calcNormFactors(dat_bulk_rawcounts_df)
dat_bulk_logCPM <- cpm(dat_bulk_rawcounts_df, log=TRUE, prior.count=2)
dat_bulk_logCPM <- data.frame(dat_bulk_logCPM)
dat_bulk_logCPM_long <- dat_bulk_logCPM %>%
  mutate(Gene = rownames(dat_bulk_logCPM)) %>%
  pivot_longer(!Gene, values_to = "logCPM", names_to = "DonorID")
hist(dat_bulk_logCPM_long$logCPM)

#convert to entrez and sum those with same entrez
dat_bulk_logCPM_long$entrez1 <- mapIds(org.Hs.eg.db, keys=c(dat_bulk_logCPM_long$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()
dat_bulk_logCPM_long$entrez2 <- mapIds(org.Hs.eg.db, keys=c(dat_bulk_logCPM_long$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character()
dat_bulk_logCPM_long <- dat_bulk_logCPM_long %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused")
dat_bulk_logCPM_long <- dat_bulk_logCPM_long %>%
  group_by(entrez, DonorID) %>%
  summarise(entrez_summed = sum(logCPM))

#read in bulk RNAseq from Humanislets.com
proc_rnaseq <- read.csv("data/Humanislets.com/proc_rnaseq.csv")
proc_rnaseq_long <- proc_rnaseq %>%
  pivot_longer(!gene_id, names_to = "DonorID", values_to = "logCPM")
donor_HI <- donor_data %>%
  dplyr::select(record_id, donorage, donorsex, diagnosis) %>%
  mutate(dataset = rep("Humanislets",nrow(donor_data)))
donor_HPAP_select <- donor_HPAP %>%
  dplyr::select(donor_ID, age_years, sex, simplified_diagnosis) %>%
  mutate(dataset = rep("HPAP",nrow(donor_HPAP)))
colnames(donor_HPAP_select) <- colnames(donor_HI)

#combine metadata
metadata_combined <- bind_rows(donor_HPAP_select, donor_HI)
metadata_combined <-metadata_combined %>%
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control",
    .default = diagnosis
  ))


#combine datasets
proc_rnaseq_long <- proc_rnaseq_long %>% dplyr:: select(gene_id, DonorID, logCPM)
dat_bulk_logCPM_long <- dat_bulk_logCPM_long %>% dplyr:: select(entrez, DonorID, entrez_summed)
colnames(dat_bulk_logCPM_long) <- colnames(proc_rnaseq_long)
dat_bulk_logCPM_long$gene_id <- as.numeric(dat_bulk_logCPM_long$gene_id)
bulk_logcpm_combined <- bind_rows(proc_rnaseq_long, dat_bulk_logCPM_long)

#limma
bulk_youngcontrols_feature <- bulk_logcpm_combined %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_bulk_youngcontrols <- bulk_youngcontrols_feature$gene_id
bulk_youngcontrols_feature <- data.frame(bulk_youngcontrols_feature)
rownames(bulk_youngcontrols_feature) <- genes_bulk_youngcontrols
bulk_youngcontrols_feature <- bulk_youngcontrols_feature[,-1] #remove the gene column

bulk_youngcontrols_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control"), donorage > 14, donorage < 40) %>%
  filter(record_id %in% colnames(bulk_youngcontrols_feature))

bulk_youngcontrols_feature <- bulk_youngcontrols_feature %>%
  dplyr::select(bulk_youngcontrols_metadata$record_id)

table(bulk_youngcontrols_metadata$donorsex, bulk_youngcontrols_metadata$dataset)
#7 F, 17 M from HPAP, 8 F, 18 M from humanislets

#check order
all(colnames(bulk_youngcontrols_feature) == bulk_youngcontrols_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_bulk_youngcontrols <- apply(bulk_youngcontrols_feature, 1, function(x){sum(!is.na(x)) >= 9})
bulk_youngcontrols_feature <- bulk_youngcontrols_feature[feature.keep_bulk_youngcontrols, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(bulk_youngcontrols_metadata$donorsex))
fixedEffects <- c("donorage","dataset") #correct for age
all.vars <- c("donorsex", fixedEffects)
design_bulkyoungcontrols <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = bulk_youngcontrols_metadata)
colnames(design_bulkyoungcontrols)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_bulkyoungcontrols <- list()
ref <- "Female" #using female as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_bulkyoungcontrols <- as.list(paste("Male", "-", ref, sep = ""))
myargs_bulkyoungcontrols[["levels"]] <- design_bulkyoungcontrols
contrast.matrix_bulkyoungcontrols <- do.call(makeContrasts, myargs_bulkyoungcontrols)

# get results
fit_bulkyoungcontrols <- lmFit(bulk_youngcontrols_feature, design_bulkyoungcontrols, trend = TRUE, robust = TRUE)
fit_bulkyoungcontrols <- contrasts.fit(fit_bulkyoungcontrols, contrast.matrix_bulkyoungcontrols)
fit_bulkyoungcontrols <- eBayes(fit_bulkyoungcontrols)
res.table_bulkyoungcontrols <- topTable(fit_bulkyoungcontrols, number = Inf)

# Remove results rows with NAs
res.table_bulkyoungcontrols <- res.table_bulkyoungcontrols[!is.na(res.table_bulkyoungcontrols$P.Value), ]

# Save output
res.table_bulkyoungcontrols$entrez <- rownames(res.table_bulkyoungcontrols)
res.table_bulkyoungcontrols$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_bulkyoungcontrols$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_bulkyoungcontrols, "Output/Final figures/Fig2/BulkHI and pbHPAP_combined_youngcontrols_dea_results_correctforage_dataset.csv", row.names = FALSE)

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
  fpath <- paste0("Output/Final figures/Fig2/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

gsea_human(df = res.table_bulkyoungcontrols,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "BulkHI and pbHPAP_youngcontrols")

## Fig 2A ## ----
GO_CC <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined$NES <- -GO_combined$NES #make direction such that positive = female-biased
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
jaccard_similarity <- function(vec1, vec2) {
  intersection <- length(intersect(vec1, vec2))
  union <- length(union(vec1, vec2))
  if (union == 0) return(0)
  return(intersection / union)
}

threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #39

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_F <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_M <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_M$Description), unique(GO_combined_refined_F$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,3.5), breaks = c(1,3.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/Fig2/Fig2A.tiff", width = 8.5, height = 8)

## Fig 2B ##----
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
dim(res.table %>% filter(Adjusted.p_value < 0.05, logFC > 0)) #0 female-biased
dim(res.table %>% filter(Adjusted.p_value < 0.05, logFC < 0)) #0 male-biased

write.csv(res.table, "Output/Final figures/Fig2/prot_dea_results_correctforage_controls_15to39.csv", row.names = FALSE)

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
  fpath <- paste0("Output/Final figures/Fig2/GSEA_p_", name, ".xlsx")
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

#Graph
GO_CC <- read_excel("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined$NES <- -GO_combined$NES #make direction consistent with RNAseq, positive = female-biased

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
jaccard_similarity <- function(vec1, vec2) {
  intersection <- length(intersect(vec1, vec2))
  union <- length(union(vec1, vec2))
  if (union == 0) return(0)
  return(intersection / union)
}

threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]

GO_combined_top60 <- GO_combined_refined[1:60,]
GO_combined_top60$direction <- ifelse(GO_combined_top60$NES>0,-1,1)
GO_combined_top60$directionlogpadj <- GO_combined_top60$direction * log10(GO_combined_top60$p.adjust)

#get pathway order
GO_combined_top60_F <- GO_combined_top60 %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_top60_M <- GO_combined_top60 %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_top60_M$Description), unique(GO_combined_top60_F$Description))

GO_combined_top60$Description <- factor(GO_combined_top60$Description, levels=pathway_order)
GO_combined_top60 %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.5), breaks = c(1,2.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/Fig2/Fig2B.tiff", width = 8, height = 11.5)

### Figure S3 ### ----
beta <- subset(x = dat, idents = "Beta") #subset by cell type
alpha <- subset(x = dat, idents = "Alpha")
beta_bulk <- AggregateExpression(beta, group.by = c("hpap_id"), return.seurat = TRUE) #pseudobulk
beta_bulk_rawcounts_df <- data.frame(beta_bulk[["RNA"]]$counts) #retrieve raw counts
alpha_bulk <- AggregateExpression(alpha, group.by = c("hpap_id"), return.seurat = TRUE)
alpha_bulk_rawcounts_df <- data.frame(alpha_bulk[["RNA"]]$counts)

#filter
beta_bulk_rawcounts_df <- beta_bulk_rawcounts_df %>%
  tibble::rownames_to_column("Gene") %>%
  filter(rowSums(dplyr::select(., -Gene) == 0) / (ncol(.) - 1) <= 0.8) %>%
  tibble::column_to_rownames("Gene")
alpha_bulk_rawcounts_df <- alpha_bulk_rawcounts_df %>%
  tibble::rownames_to_column("Gene") %>%
  filter(rowSums(dplyr::select(., -Gene) == 0) / (ncol(.) - 1) <= 0.8) %>%
  tibble::column_to_rownames("Gene")

#transform to logCPM
beta_bulk_rawcounts_df <- DGEList(counts=beta_bulk_rawcounts_df)
beta_bulk_rawcounts_df <- calcNormFactors(beta_bulk_rawcounts_df) # Defaults to TMM normalization
beta_bulk_logCPM <- cpm(beta_bulk_rawcounts_df, log=TRUE, prior.count=2)
beta_bulk_logCPM <- data.frame(beta_bulk_logCPM)
beta_bulk_logCPM_long <- beta_bulk_logCPM %>%
  mutate(Gene = rownames(beta_bulk_logCPM)) %>%
  pivot_longer(!Gene, values_to = "logCPM", names_to = "DonorID")
hist(beta_bulk_logCPM_long$logCPM)

alpha_bulk_rawcounts_df <- DGEList(counts=alpha_bulk_rawcounts_df)
alpha_bulk_rawcounts_df <- calcNormFactors(alpha_bulk_rawcounts_df) # Defaults to TMM normalization
alpha_bulk_logCPM <- cpm(alpha_bulk_rawcounts_df, log=TRUE, prior.count=2)
alpha_bulk_logCPM <- data.frame(alpha_bulk_logCPM)
alpha_bulk_logCPM_long <- alpha_bulk_logCPM %>%
  mutate(Gene = rownames(alpha_bulk_logCPM)) %>%
  pivot_longer(!Gene, values_to = "logCPM", names_to = "DonorID")
hist(alpha_bulk_logCPM_long$logCPM)

#convert to entrez and sum those with same entrez
beta_bulk_logCPM_long$entrez1 <- mapIds(org.Hs.eg.db, keys=c(beta_bulk_logCPM_long$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
beta_bulk_logCPM_long$entrez2 <- mapIds(org.Hs.eg.db, keys=c(beta_bulk_logCPM_long$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID
beta_bulk_logCPM_long <- beta_bulk_logCPM_long %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused")
beta_bulk_logCPM_long <- beta_bulk_logCPM_long %>%
  group_by(entrez, DonorID) %>%
  summarise(entrez_summed = sum(logCPM))

alpha_bulk_logCPM_long$entrez1 <- mapIds(org.Hs.eg.db, keys=c(alpha_bulk_logCPM_long$Gene), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character()  # use gene symbol to get Entrez ID
alpha_bulk_logCPM_long$entrez2 <- mapIds(org.Hs.eg.db, keys=c(alpha_bulk_logCPM_long$Gene), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() # use gene symbol as alias to get Entrez ID
alpha_bulk_logCPM_long <- alpha_bulk_logCPM_long %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 
alpha_bulk_logCPM_long <- alpha_bulk_logCPM_long %>%
  group_by(entrez, DonorID) %>%
  summarise(entrez_summed = sum(logCPM))

#limma
beta_youngcontrols_feature_HPAP <- beta_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_beta_youngcontrols_HPAP <- beta_youngcontrols_feature_HPAP$gene_id
beta_youngcontrols_feature_HPAP <- data.frame(beta_youngcontrols_feature_HPAP)
rownames(beta_youngcontrols_feature_HPAP) <- genes_beta_youngcontrols_HPAP
beta_youngcontrols_feature_HPAP <- beta_youngcontrols_feature_HPAP[,-1] #remove the gene column

beta_youngcontrols_metadata_HPAP <- donor_HPAP_select %>%
  filter(diagnosis %in% c("Control"), donorage > 14, donorage < 40) %>%
  filter(record_id %in% colnames(beta_youngcontrols_feature_HPAP))

beta_youngcontrols_feature_HPAP <- beta_youngcontrols_feature_HPAP %>%
  dplyr::select(beta_youngcontrols_metadata_HPAP$record_id)

table(beta_youngcontrols_metadata_HPAP$donorsex, beta_youngcontrols_metadata_HPAP$dataset)
#7 F, 17 M

#check order
all(colnames(beta_youngcontrols_feature_HPAP) == beta_youngcontrols_metadata_HPAP$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_beta_youngcontrols_HPAP <- apply(beta_youngcontrols_feature_HPAP, 1, function(x){sum(!is.na(x)) >= 9})
beta_youngcontrols_feature_HPAP <- beta_youngcontrols_feature_HPAP[feature.keep_beta_youngcontrols_HPAP, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(beta_youngcontrols_metadata_HPAP$donorsex))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("donorsex", fixedEffects)
design_betayoungcontrols_HPAP <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = beta_youngcontrols_metadata_HPAP)
colnames(design_betayoungcontrols_HPAP)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_betayoungcontrols_HPAP <- list()
ref <- "Female" #using female as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_betayoungcontrols_HPAP <- as.list(paste("Male", "-", ref, sep = ""))
myargs_betayoungcontrols_HPAP[["levels"]] <- design_betayoungcontrols_HPAP
contrast.matrix_betayoungcontrols_HPAP <- do.call(makeContrasts, myargs_betayoungcontrols_HPAP)

# get results
fit_betayoungcontrols_HPAP <- lmFit(beta_youngcontrols_feature_HPAP, design_betayoungcontrols_HPAP, trend = TRUE, robust = TRUE)
fit_betayoungcontrols_HPAP <- contrasts.fit(fit_betayoungcontrols_HPAP, contrast.matrix_betayoungcontrols_HPAP)
fit_betayoungcontrols_HPAP <- eBayes(fit_betayoungcontrols_HPAP)
res.table_betayoungcontrols_HPAP <- topTable(fit_betayoungcontrols_HPAP, number = Inf)

# Remove results rows with NAs
res.table_betayoungcontrols_HPAP <- res.table_betayoungcontrols_HPAP[!is.na(res.table_betayoungcontrols_HPAP$P.Value), ]

# Save output
res.table_betayoungcontrols_HPAP$entrez <- rownames(res.table_betayoungcontrols_HPAP)
res.table_betayoungcontrols_HPAP$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_betayoungcontrols_HPAP$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_betayoungcontrols_HPAP, "Output/Final figures/FigS3/pbBeta_HPAPonly_youngcontrols_dea_results_correctforage_dataset.csv", row.names = FALSE)

# GSEA
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
  fpath <- paste0("Output/Final figures/FigS3/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}
gsea_human(df = res.table_betayoungcontrols_HPAP,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "HPAPonly_pbBeta_youngcontrols")

#same for alpha-cells
alpha_youngcontrols_feature_HPAP <- alpha_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_alpha_youngcontrols_HPAP <- alpha_youngcontrols_feature_HPAP$gene_id
alpha_youngcontrols_feature_HPAP <- data.frame(alpha_youngcontrols_feature_HPAP)
rownames(alpha_youngcontrols_feature_HPAP) <- genes_alpha_youngcontrols_HPAP
alpha_youngcontrols_feature_HPAP <- alpha_youngcontrols_feature_HPAP[,-1] #remove the gene column

alpha_youngcontrols_metadata_HPAP <- donor_HPAP_select %>%
  filter(diagnosis %in% c("Control"), donorage > 14, donorage < 40) %>%
  filter(record_id %in% colnames(alpha_youngcontrols_feature_HPAP))

alpha_youngcontrols_feature_HPAP <- alpha_youngcontrols_feature_HPAP %>%
  dplyr::select(alpha_youngcontrols_metadata_HPAP$record_id)

table(alpha_youngcontrols_metadata_HPAP$donorsex, alpha_youngcontrols_metadata_HPAP$dataset)
#7 F, 17 M

#check order
all(colnames(alpha_youngcontrols_feature_HPAP) == alpha_youngcontrols_metadata_HPAP$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_alpha_youngcontrols_HPAP <- apply(alpha_youngcontrols_feature_HPAP, 1, function(x){sum(!is.na(x)) >= 9})
alpha_youngcontrols_feature_HPAP <- alpha_youngcontrols_feature_HPAP[feature.keep_alpha_youngcontrols_HPAP, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(alpha_youngcontrols_metadata_HPAP$donorsex))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("donorsex", fixedEffects)
design_alphayoungcontrols_HPAP <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = alpha_youngcontrols_metadata_HPAP)
colnames(design_alphayoungcontrols_HPAP)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_alphayoungcontrols_HPAP <- list()
ref <- "Female" #using female as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_alphayoungcontrols_HPAP <- as.list(paste("Male", "-", ref, sep = ""))
myargs_alphayoungcontrols_HPAP[["levels"]] <- design_alphayoungcontrols_HPAP
contrast.matrix_alphayoungcontrols_HPAP <- do.call(makeContrasts, myargs_alphayoungcontrols_HPAP)

# get results
fit_alphayoungcontrols_HPAP <- lmFit(alpha_youngcontrols_feature_HPAP, design_alphayoungcontrols_HPAP, trend = TRUE, robust = TRUE)
fit_alphayoungcontrols_HPAP <- contrasts.fit(fit_alphayoungcontrols_HPAP, contrast.matrix_alphayoungcontrols_HPAP)
fit_alphayoungcontrols_HPAP <- eBayes(fit_alphayoungcontrols_HPAP)
res.table_alphayoungcontrols_HPAP <- topTable(fit_alphayoungcontrols_HPAP, number = Inf)

# Remove results rows with NAs
res.table_alphayoungcontrols_HPAP <- res.table_alphayoungcontrols_HPAP[!is.na(res.table_alphayoungcontrols_HPAP$P.Value), ]

# Save output
res.table_alphayoungcontrols_HPAP$entrez <- rownames(res.table_alphayoungcontrols_HPAP)
res.table_alphayoungcontrols_HPAP$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_alphayoungcontrols_HPAP$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_alphayoungcontrols_HPAP, "Output/Revisions/pbAlpha_HPAPonly_youngcontrols_dea_results_correctforage_dataset.csv", row.names = FALSE)

# GSEA
gsea_human(df = res.table_betayoungcontrols_HPAP,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "HPAPonly_pbAlpha_youngcontrols")

## Fig S3A ## ----
GO_CC <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined$NES <- -GO_combined$NES #make direction such that positive = female-biased
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway

threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #14

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_F <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_M <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_M$Description), unique(GO_combined_refined_F$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
#make long pathway name two lines
GO_combined_refined <- GO_combined_refined %>%
  mutate(Description = case_when(
    Description == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen" ~ "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen",
    .default = Description
  ))
pathway_order <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen", pathway_order)
GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)

GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.5), breaks = c(1,2.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/FigS3/FigS3A.tiff", width = 6.5, height = 4.5)

## Fig S3B ## ----
GO_CC <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Revisions/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined$NES <- -GO_combined$NES #make direction such that positive = female-biased
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #14

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_F <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_M <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_M$Description), unique(GO_combined_refined_F$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
#make long pathway name two lines
GO_combined_refined <- GO_combined_refined %>%
  mutate(Description = case_when(
    Description == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen" ~ "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen",
    .default = Description
  ))
pathway_order <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen", pathway_order)
GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)

GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.5), breaks = c(1,2.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/FigS3/FigS3B.tiff", width = 6.5, height = 4.5)

###Figure 3-4 ### ----
#Combined transcriptomics for Control vs T2D
#Females
bulk_F_feature <- bulk_logcpm_combined %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_bulk_F <- bulk_F_feature$gene_id
bulk_F_feature <- data.frame(bulk_F_feature)
rownames(bulk_F_feature) <- genes_bulk_F
bulk_F_feature <- bulk_F_feature[,-1] #remove the gene column

bulk_F_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Female") %>%
  filter(record_id %in% colnames(bulk_F_feature))

bulk_F_feature <- bulk_F_feature %>%
  dplyr::select(bulk_F_metadata$record_id)

table(bulk_F_metadata$diagnosis, bulk_F_metadata$dataset)
#16 Control, 11 T2D from HPAP, 36 Control, 5 T2D from humanislets.com

#check order
all(colnames(bulk_F_feature) == bulk_F_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_bulk_F <- apply(bulk_F_feature, 1, function(x){sum(!is.na(x)) >= 9})
bulk_F_feature <- bulk_F_feature[feature.keep_bulk_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(bulk_F_metadata$diagnosis))
fixedEffects <- c("donorage","dataset") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_bulkF <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = bulk_F_metadata)
colnames(design_bulkF)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_bulkF <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_bulkF <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_bulkF[["levels"]] <- design_bulkF
contrast.matrix_bulkF <- do.call(makeContrasts, myargs_bulkF)

# get results
fit_bulkF <- lmFit(bulk_F_feature, design_bulkF, trend = TRUE, robust = TRUE)
fit_bulkF <- contrasts.fit(fit_bulkF, contrast.matrix_bulkF)
fit_bulkF <- eBayes(fit_bulkF)
res.table_bulkF <- topTable(fit_bulkF, number = Inf)

# Remove results rows with NAs
res.table_bulkF <- res.table_bulkF[!is.na(res.table_bulkF$P.Value), ]

# Save output
res.table_bulkF$entrez <- rownames(res.table_bulkF)
res.table_bulkF$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_bulkF$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_bulkF, "Output/Final figures/Fig3-4/BulkHI and pbHPAP_combined_F_dea_results_correctforage_dataset.csv", row.names = FALSE)

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
  fpath <- paste0("Output/Final figures/Fig3-4/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}

gsea_human(df = res.table_bulkF,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "BulkHI and pbHPAP_F")

#Males
bulk_M_feature <- bulk_logcpm_combined %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_bulk_F <- bulk_M_feature$gene_id
bulk_M_feature <- data.frame(bulk_M_feature)
rownames(bulk_M_feature) <- genes_bulk_F
bulk_M_feature <- bulk_M_feature[,-1] #remove the gene column

bulk_M_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Male") %>%
  filter(record_id %in% colnames(bulk_M_feature))

bulk_M_feature <- bulk_M_feature %>%
  dplyr::select(bulk_M_metadata$record_id)

table(bulk_M_metadata$diagnosis, bulk_M_metadata$dataset)
#24 Control, 7 T2D from HPAP, 66 Control, 10 T2D from humanislets.com

#check order
all(colnames(bulk_M_feature) == bulk_M_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_bulk_F <- apply(bulk_M_feature, 1, function(x){sum(!is.na(x)) >= 9})
bulk_M_feature <- bulk_M_feature[feature.keep_bulk_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(bulk_M_metadata$diagnosis))
fixedEffects <- c("donorage","dataset") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_bulkM <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = bulk_M_metadata)
colnames(design_bulkM)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_bulkM <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_bulkM <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_bulkM[["levels"]] <- design_bulkM
contrast.matrix_bulkM <- do.call(makeContrasts, myargs_bulkM)

# get results
fit_bulkM <- lmFit(bulk_M_feature, design_bulkM, trend = TRUE, robust = TRUE)
fit_bulkM <- contrasts.fit(fit_bulkM, contrast.matrix_bulkM)
fit_bulkM <- eBayes(fit_bulkM)
res.table_bulkM <- topTable(fit_bulkM, number = Inf)

# Remove results rows with NAs
res.table_bulkM <- res.table_bulkM[!is.na(res.table_bulkM$P.Value), ]

# Save output
res.table_bulkM$entrez <- rownames(res.table_bulkM)
res.table_bulkM$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_bulkM$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_bulkM, "Output/Final figures/Fig3-4/BulkHI and pbHPAP_combined_M_dea_results_correctforage_dataset.csv", row.names = FALSE)

#GSEA
gsea_human(df = res.table_bulkM,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "BulkHI and pbHPAP_M")

#Proteomics
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

write.csv(res.table_F, "Output/Final figures/Fig3-4/prot_dea_results_correctforage_F.csv", row.names = FALSE)

#GSEA
res.table_F$entrez1 <- mapIds(org.Hs.eg.db, keys=c(res.table_F$common_name), column="ENTREZID", keytype="SYMBOL", multiVals="first") %>% 
  as.character() 
res.table_F$entrez2 <- mapIds(org.Hs.eg.db, keys=c(res.table_F$common_name), column="ENTREZID", keytype="ALIAS", multiVals="first") %>% 
  as.character() 
res.table_F <- res.table_F %>% 
  mutate(entrez = coalesce(entrez1, entrez2), .keep = "unused") 

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

write.csv(res.table_M, "Output/Final figures/Fig3-4/prot_dea_results_correctforage_M.csv", row.names = FALSE)

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

#Refine pathways for visualisation
GO_CC_RNA_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_CC")
GO_MF_RNA_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_MF")
GO_BP_RNA_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_BP")
GO_combined_RNA_F <- bind_rows(GO_CC_RNA_F, GO_MF_RNA_F, GO_BP_RNA_F) %>%
  arrange(pvalue)
GO_combined_RNA_F <-GO_combined_RNA_F %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined_RNA_F$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined_RNA_F <- GO_combined_RNA_F[-repeat_indices,]
dim(GO_combined_refined_RNA_F) #424

GO_CC_Protein_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_CC")
GO_MF_Protein_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_MF")
GO_BP_Protein_F <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_BP")
GO_combined_Protein_F <- bind_rows(GO_CC_Protein_F, GO_MF_Protein_F, GO_BP_Protein_F) %>%
  arrange(pvalue)
GO_combined_Protein_F <-GO_combined_Protein_F %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined_Protein_F$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined_Protein_F <- GO_combined_Protein_F[-repeat_indices,]
dim(GO_combined_refined_Protein_F) #419

length(intersect(GO_combined_refined_RNA_F$Description, GO_combined_refined_Protein_F$Description)) #55 common
hist(GO_combined_refined_RNA_F$p.adjust)
hist(GO_combined_refined_Protein_F$p.adjust)

#make more stringent
GO_combined_refined_RNA_F <- GO_combined_refined_RNA_F %>%
  filter(p.adjust < 0.0001)
dim(GO_combined_refined_RNA_F) #59
GO_combined_refined_Protein_F <- GO_combined_refined_Protein_F %>%
  filter(p.adjust < 0.0001) 
dim(GO_combined_refined_Protein_F) #45
length(intersect(GO_combined_refined_RNA_F$Description, GO_combined_refined_Protein_F$Description)) #4 common
intersect(GO_combined_refined_RNA_F$Description, GO_combined_refined_Protein_F$Description) #mito things and RNP complex biogenesis

#combine into one dataframe
GO_combined_refined_RNA_Protein_F <- bind_rows(GO_combined_refined_RNA_F, GO_combined_refined_Protein_F)
GO_combined_refined_RNA_Protein_F$dataset <- c(rep("RNA", nrow(GO_combined_refined_RNA_F)), rep("Protein",nrow(GO_combined_refined_Protein_F)))
length(unique(GO_combined_refined_RNA_Protein_F$Description)) #100

GO_combined_refined_RNA_Protein_F_up <- GO_combined_refined_RNA_Protein_F %>% filter(NES > 0)
GO_combined_refined_RNA_Protein_F_down <- GO_combined_refined_RNA_Protein_F %>% filter(NES < 0)

## Fig 3A ## ----
#reorder pathways
org.Hs.eg.db::org.Hs.egGO2ALLEGS |> as.list() -> GO2ALLEGS
GO_combined_refined_RNA_Protein_F_up$ID -> sigGOs
length(sigGOs)
CreateGOrge(
  gene_sets = GO2ALLEGS[sigGOs], 
  recursive = T, iter_max = 1
) -> GOrged
GOrged@cluster_n #4 clusters

GO_combined_refined_RNA_Protein_F_up_clusters <- list()
cluster_list <- c()
for (i in 1:GOrged@cluster_n){
  cluster <- GO_combined_refined_RNA_Protein_F_up %>% 
    filter(ID %in% GOrged@cluster_list[[i]])
  GO_combined_refined_RNA_Protein_F_up_clusters[[i]] <- cluster
  cluster_list <- c(cluster_list, rep(i,nrow(cluster)))
}
GO_combined_refined_RNA_Protein_F_up_clusters <- do.call(bind_rows, GO_combined_refined_RNA_Protein_F_up_clusters)
GO_combined_refined_RNA_Protein_F_up_clusters <- data.frame(GO_combined_refined_RNA_Protein_F_up_clusters)
GO_combined_refined_RNA_Protein_F_up_clusters$Cluster <- cluster_list

#use clusters to get order of pathways
GO_combined_refined_RNA_Protein_F_up_clusters <- GO_combined_refined_RNA_Protein_F_up_clusters %>% arrange(desc(Cluster))
GO_combined_refined_RNA_Protein_F_up_clusters$Description <- factor(GO_combined_refined_RNA_Protein_F_up_clusters$Description, levels = unique(GO_combined_refined_RNA_Protein_F_up_clusters$Description))

#graph
summary(GO_combined_refined_RNA_Protein_F_up_clusters$NES)
GO_combined_refined_RNA_Protein_F_up_clusters$dataset <- factor(GO_combined_refined_RNA_Protein_F_up_clusters$dataset, levels = c("RNA","Protein"))
GO_combined_refined_RNA_Protein_F_up_clusters %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways up with T2D", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#992002", mid = "#E05A38", high = "white", midpoint = 0.00005, limits = c(0,0.0001)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Revisions/Fig3A.tiff", width = 8, height = 7)

## Fig 4A ## ----
GO_combined_refined_RNA_Protein_F_down$ID -> sigGOs
length(sigGOs)
CreateGOrge(
  gene_sets = GO2ALLEGS[sigGOs], 
  recursive = T, iter_max = 1
) -> GOrged
GOrged@cluster_n #6 clusters

GO_combined_refined_RNA_Protein_F_down_clusters <- list()
cluster_list <- c()
for (i in 1:GOrged@cluster_n){
  cluster <- GO_combined_refined_RNA_Protein_F_down %>% 
    filter(ID %in% GOrged@cluster_list[[i]])
  GO_combined_refined_RNA_Protein_F_down_clusters[[i]] <- cluster
  cluster_list <- c(cluster_list, rep(i,nrow(cluster)))
}
GO_combined_refined_RNA_Protein_F_down_clusters <- do.call(bind_rows, GO_combined_refined_RNA_Protein_F_down_clusters)
GO_combined_refined_RNA_Protein_F_down_clusters <- data.frame(GO_combined_refined_RNA_Protein_F_down_clusters)
GO_combined_refined_RNA_Protein_F_down_clusters$Cluster <- cluster_list

#use clusters to get order of pathways
GO_combined_refined_RNA_Protein_F_down_clusters <- GO_combined_refined_RNA_Protein_F_down_clusters %>% arrange(desc(Cluster))
GO_combined_refined_RNA_Protein_F_down_clusters$Description <- factor(GO_combined_refined_RNA_Protein_F_down_clusters$Description, levels = unique(GO_combined_refined_RNA_Protein_F_down_clusters$Description))

#graph
summary(GO_combined_refined_RNA_Protein_F_down_clusters$NES)
GO_combined_refined_RNA_Protein_F_down_clusters$dataset <- factor(GO_combined_refined_RNA_Protein_F_down_clusters$dataset, levels = c("RNA","Protein"))
GO_combined_refined_RNA_Protein_F_down_clusters %>%
  ggplot(aes(x=dataset, y=Description, size = -NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways down with T2D", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#55066b", mid = "#C337EB", high = "white", midpoint = 0.00005, limits = c(0,0.0001)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,3), breaks = c(1.5,3)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/Fig3-4/Fig4A.tiff", width = 7.6, height = 7.5)

## Fig 3B ## ----
GO_CC_RNA_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_CC")
GO_MF_RNA_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_MF")
GO_BP_RNA_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_BP")
GO_combined_RNA_M <- bind_rows(GO_CC_RNA_M, GO_MF_RNA_M, GO_BP_RNA_M) %>%
  arrange(pvalue)
GO_combined_RNA_M <-GO_combined_RNA_M %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined_RNA_M$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined_RNA_M <- GO_combined_RNA_M[-repeat_indices,]
dim(GO_combined_refined_RNA_M) 925

GO_CC_Protein_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_CC")
GO_MF_Protein_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_MF")
GO_BP_Protein_M <- read_excel("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_BP")
GO_combined_Protein_M <- bind_rows(GO_CC_Protein_M, GO_MF_Protein_M, GO_BP_Protein_M) %>%
  arrange(pvalue)
GO_combined_Protein_M <-GO_combined_Protein_M %>% filter(p.adjust < 0.05)

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined_Protein_M$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined_Protein_M <- GO_combined_Protein_M[-repeat_indices,]
dim(GO_combined_refined_Protein_M) #324

length(intersect(GO_combined_refined_RNA_M$Description, GO_combined_refined_Protein_M$Description)) #109 common
hist(GO_combined_refined_RNA_M$p.adjust)
hist(GO_combined_refined_Protein_M$p.adjust)

#make more stringent
GO_combined_refined_RNA_M <- GO_combined_refined_RNA_M %>%
  filter(p.adjust < 0.0001)
dim(GO_combined_refined_RNA_M) #144
GO_combined_refined_Protein_M <- GO_combined_refined_Protein_M %>%
  filter(p.adjust < 0.0001) 
dim(GO_combined_refined_Protein_M) #42
length(intersect(GO_combined_refined_RNA_M$Description, GO_combined_refined_Protein_M$Description)) #7 common
intersect(GO_combined_refined_RNA_M$Description, GO_combined_refined_Protein_M$Description) #ribosome, virus, junction, vesicle

#combined into one dataframe
GO_combined_refined_RNA_Protein_M <- bind_rows(GO_combined_refined_RNA_M, GO_combined_refined_Protein_M)
GO_combined_refined_RNA_Protein_M$dataset <- c(rep("RNA", nrow(GO_combined_refined_RNA_M)), rep("Protein",nrow(GO_combined_refined_Protein_M)))
length(unique(GO_combined_refined_RNA_Protein_M$Description)) #180

GO_combined_refined_RNA_Protein_M_up <- GO_combined_refined_RNA_Protein_M %>% filter(NES > 0)
dim(GO_combined_refined_RNA_Protein_M_up) #148
#take top 60
GO_combined_refined_RNA_Protein_M_up <- GO_combined_refined_RNA_Protein_M_up %>% arrange(pvalue)
GO_combined_refined_RNA_Protein_M_up <- GO_combined_refined_RNA_Protein_M_up[1:60,]

GO_combined_refined_RNA_Protein_M_down <- GO_combined_refined_RNA_Protein_M %>% filter(NES < 0)
dim(GO_combined_refined_RNA_Protein_M_down)#39

#cluster for pathway order
GO_combined_refined_RNA_Protein_M_up$ID -> sigGOs
length(sigGOs)
CreateGOrge(
  gene_sets = GO2ALLEGS[sigGOs], 
  recursive = T, iter_max = 1
) -> GOrged
GOrged@cluster_n #4 clusters

GO_combined_refined_RNA_Protein_M_up_clusters <- list()
cluster_list <- c()
for (i in 1:GOrged@cluster_n){
  cluster <- GO_combined_refined_RNA_Protein_M_up %>% 
    filter(ID %in% GOrged@cluster_list[[i]])
  GO_combined_refined_RNA_Protein_M_up_clusters[[i]] <- cluster
  cluster_list <- c(cluster_list, rep(i,nrow(cluster)))
}
GO_combined_refined_RNA_Protein_M_up_clusters <- do.call(bind_rows, GO_combined_refined_RNA_Protein_M_up_clusters)
GO_combined_refined_RNA_Protein_M_up_clusters <- data.frame(GO_combined_refined_RNA_Protein_M_up_clusters)
GO_combined_refined_RNA_Protein_M_up_clusters$Cluster <- cluster_list

#use clusters to get order of pathways
GO_combined_refined_RNA_Protein_M_up_clusters <- GO_combined_refined_RNA_Protein_M_up_clusters %>% arrange(desc(Cluster))
GO_combined_refined_RNA_Protein_M_up_clusters$Description <- factor(GO_combined_refined_RNA_Protein_M_up_clusters$Description, levels = unique(GO_combined_refined_RNA_Protein_M_up_clusters$Description))

#graph
summary(GO_combined_refined_RNA_Protein_M_up_clusters$NES)
GO_combined_refined_RNA_Protein_M_up_clusters$dataset <- factor(GO_combined_refined_RNA_Protein_M_up_clusters$dataset, levels = c("RNA","Protein"))
GO_combined_refined_RNA_Protein_M_up_clusters %>%
  ggplot(aes(x=dataset, y=Description, size = NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways up with T2D", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#9c6b05", mid = "#DBA126", high = "white", midpoint = 0.00005, limits = c(0,0.0001)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,3), breaks = c(1.5,3)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/Fig3-4/Fig3B.tiff", width = 8, height = 7.5)


## Fig 4B ## ----
GO_combined_refined_RNA_Protein_M_down$ID -> sigGOs
length(sigGOs)
CreateGOrge(
  gene_sets = GO2ALLEGS[sigGOs], 
  recursive = T, iter_max = 1
) -> GOrged
GOrged@cluster_n #4 clusters

GO_combined_refined_RNA_Protein_M_down_clusters <- list()
cluster_list <- c()
for (i in 1:GOrged@cluster_n){
  cluster <- GO_combined_refined_RNA_Protein_M_down %>% 
    filter(ID %in% GOrged@cluster_list[[i]])
  GO_combined_refined_RNA_Protein_M_down_clusters[[i]] <- cluster
  cluster_list <- c(cluster_list, rep(i,nrow(cluster)))
}
GO_combined_refined_RNA_Protein_M_down_clusters <- do.call(bind_rows, GO_combined_refined_RNA_Protein_M_down_clusters)
GO_combined_refined_RNA_Protein_M_down_clusters <- data.frame(GO_combined_refined_RNA_Protein_M_down_clusters)
GO_combined_refined_RNA_Protein_M_down_clusters$Cluster <- cluster_list

#use clusters to get order of pathways
GO_combined_refined_RNA_Protein_M_down_clusters <- GO_combined_refined_RNA_Protein_M_down_clusters %>% arrange(desc(Cluster))
GO_combined_refined_RNA_Protein_M_down_clusters$Description <- factor(GO_combined_refined_RNA_Protein_M_down_clusters$Description, levels = unique(GO_combined_refined_RNA_Protein_M_down_clusters$Description))

#graph
summary(GO_combined_refined_RNA_Protein_M_down_clusters$NES)
GO_combined_refined_RNA_Protein_M_down_clusters$dataset <- factor(GO_combined_refined_RNA_Protein_M_down_clusters$dataset, levels = c("RNA","Protein"))
GO_combined_refined_RNA_Protein_M_down_clusters %>%
  ggplot(aes(x=dataset, y=Description, size = -NES, colour = p.adjust))+
  geom_point() +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x = "", y = "", 
       title = "GO pathways down with T2D", colour = "Adjusted P-value", size = "NES") +
  scale_color_gradient2(low = "#0a735b", mid = "#3EDBB8", high = "white", midpoint = 0.00005, limits = c(0,0.0001)) +
  scale_size_continuous(range = c(1,7), limits = c(1.5,2.5), breaks = c(1.5,2.5)) +
  scale_x_discrete(position = "top")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/Fig3-4/Fig4B.tiff", width = 6.65, height = 6.75)

#Venn diagrams
GO_combined_RNA_F_up <- GO_combined_RNA_F %>% filter(NES > 0)
GO_combined_RNA_F_down <- GO_combined_RNA_F %>% filter(NES < 0)
GO_combined_RNA_M_up <- GO_combined_RNA_M %>% filter(NES > 0)
GO_combined_RNA_M_down <- GO_combined_RNA_M %>% filter(NES < 0)
GO_combined_Protein_F_up <- GO_combined_Protein_F %>% filter(NES > 0)
GO_combined_Protein_F_down <- GO_combined_Protein_F %>% filter(NES < 0)
GO_combined_Protein_M_up <- GO_combined_Protein_M %>% filter(NES > 0)
GO_combined_Protein_M_down <- GO_combined_Protein_M %>% filter(NES < 0)

venn.diagram(
  x = list(GO_combined_RNA_F_up$Description, GO_combined_Protein_F_up$Description),
  category.names = c("RNA" , "Protein"),
  filename = "Output/Final figures/Fig3-4/Fig3A Venn.tiff",
  output=TRUE,
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  fill = c("tomato","steelblue"),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(GO_combined_RNA_F_down$Description, GO_combined_Protein_F_down$Description),
  category.names = c("RNA" , "Protein"),
  filename = "Output/Final figures/Fig3-4/Fig4A Venn.tiff",
  output=TRUE,
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  fill = c("tomato","steelblue"),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(GO_combined_RNA_M_up$Description, GO_combined_Protein_M_up$Description),
  category.names = c("RNA" , "Protein"),
  filename = "Output/Final figures/Fig3-4/FigBA Venn.tiff",
  output=TRUE,
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  fill = c("tomato","steelblue"),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(GO_combined_RNA_M_down$Description, GO_combined_Protein_M_down$Description),
  category.names = c("RNA" , "Protein"),
  filename = "Output/Final figures/Fig3-4/Fig4B Venn.tiff",
  output=TRUE,
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  fill = c("tomato","steelblue"),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

### Figure S4 ### ----
## Fig S4A ## ----
beta_F_feature <- beta_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_beta_F <- beta_F_feature$gene_id
beta_F_feature <- data.frame(beta_F_feature)
rownames(beta_F_feature) <- genes_beta_F
beta_F_feature <- beta_F_feature[,-1] #remove the gene column

beta_F_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Female") %>%
  filter(record_id %in% colnames(beta_F_feature))

beta_F_feature <- beta_F_feature %>%
  dplyr::select(beta_F_metadata$record_id)

table(beta_F_metadata$diagnosis, beta_F_metadata$dataset)
#16 Control, 11 T2D from HPAP

#check order
all(colnames(beta_F_feature) == beta_F_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_beta_F <- apply(beta_F_feature, 1, function(x){sum(!is.na(x)) >= 9})
beta_F_feature <- beta_F_feature[feature.keep_beta_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(beta_F_metadata$diagnosis))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_betaF <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = beta_F_metadata)
colnames(design_betaF)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_betaF <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_betaF <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_betaF[["levels"]] <- design_betaF
contrast.matrix_betaF <- do.call(makeContrasts, myargs_betaF)

# get results
fit_betaF <- lmFit(beta_F_feature, design_betaF, trend = TRUE, robust = TRUE)
fit_betaF <- contrasts.fit(fit_betaF, contrast.matrix_betaF)
fit_betaF <- eBayes(fit_betaF)
res.table_betaF <- topTable(fit_betaF, number = Inf)

# Remove results rows with NAs
res.table_betaF <- res.table_betaF[!is.na(res.table_betaF$P.Value), ]

# Save output
res.table_betaF$entrez <- rownames(res.table_betaF)
res.table_betaF$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_betaF$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_betaF, "Output/Final figures/FigS4/beta HPAP only_F_dea_results_correctforage_dataset.csv", row.names = FALSE)

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
  fpath <- paste0("Output/Final figures/FigS4/GSEA_p_", name, ".xlsx")
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(paste0("GSEA results saved in - ", fpath))
  
}
gsea_human(df = res.table_betaF,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "beta_HPAPonly_F")

#Visualise
GO_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)
GO_combined %>% filter(NES > 0) %>% summarise(n=n()) #158
GO_combined %>% filter(NES < 0) %>% summarise(n=n()) #193

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/")
jaccard_similarity <- function(vec1, vec2) {
  intersection <- length(intersect(vec1, vec2))
  union <- length(union(vec1, vec2))
  if (union == 0) return(0)
  return(intersection / union)
}

threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #192
hist(GO_combined_refined$p.adjust)
GO_combined_refined <- GO_combined_refined %>% filter(p.adjust < 0.0001)
dim(GO_combined_refined) #35

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_up <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_down <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_down$Description), unique(GO_combined_refined_up$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "#55066b", mid = "white", high = "#992002", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.5), breaks = c(1,2.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/FigS4/FigS4A.tiff", width = 7.5, height = 6)

## Fig S4C ## ----
beta_M_feature <- beta_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_beta_F <- beta_M_feature$gene_id
beta_M_feature <- data.frame(beta_M_feature)
rownames(beta_M_feature) <- genes_beta_F
beta_M_feature <- beta_M_feature[,-1] #remove the gene column

beta_M_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Male") %>%
  filter(record_id %in% colnames(beta_M_feature))

beta_M_feature <- beta_M_feature %>%
  dplyr::select(beta_M_metadata$record_id)

table(beta_M_metadata$diagnosis, beta_M_metadata$dataset)
#24 Control, 7 T2D from HPAP

#check order
all(colnames(beta_M_feature) == beta_M_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_beta_F <- apply(beta_M_feature, 1, function(x){sum(!is.na(x)) >= 9})
beta_M_feature <- beta_M_feature[feature.keep_beta_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(beta_M_metadata$diagnosis))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_betaM <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = beta_M_metadata)
colnames(design_betaM)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_betaM <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_betaM <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_betaM[["levels"]] <- design_betaM
contrast.matrix_betaM <- do.call(makeContrasts, myargs_betaM)

# get results
fit_betaM <- lmFit(beta_M_feature, design_betaM, trend = TRUE, robust = TRUE)
fit_betaM <- contrasts.fit(fit_betaM, contrast.matrix_betaM)
fit_betaM <- eBayes(fit_betaM)
res.table_betaM <- topTable(fit_betaM, number = Inf)

# Remove results rows with NAs
res.table_betaM <- res.table_betaM[!is.na(res.table_betaM$P.Value), ]

# Save output
res.table_betaM$entrez <- rownames(res.table_betaM)
res.table_betaM$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_betaM$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_betaM, "Output/Final figures/FigS4/beta HPAP only_M_dea_results_correctforage_dataset.csv", row.names = FALSE)

#GSEA
gsea_human(df = res.table_betaM,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "beta_HPAPonly_M")

#Visualise
GO_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)
GO_combined %>% filter(NES > 0) %>% summarise(n=n()) #101
GO_combined %>% filter(NES < 0) %>% summarise(n=n()) #333

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #266
hist(GO_combined_refined$p.adjust)
#take p< 0.0001
GO_combined_refined <- GO_combined_refined %>% filter(p.adjust < 0.0001) #18

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_up <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_down <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_down$Description), unique(GO_combined_refined_up$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "#0a735b", mid = "white", high = "#9c6b05", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,2.6), breaks = c(1,2.6)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/FigS4/FigS4C.tiff", width = 6.25, height = 3.75)

## Fig S4B ## ----
alpha_F_feature <- alpha_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_alpha_F <- alpha_F_feature$gene_id
alpha_F_feature <- data.frame(alpha_F_feature)
rownames(alpha_F_feature) <- genes_alpha_F
alpha_F_feature <- alpha_F_feature[,-1] #remove the gene column

alpha_F_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Female") %>%
  filter(record_id %in% colnames(alpha_F_feature))

alpha_F_feature <- alpha_F_feature %>%
  dplyr::select(alpha_F_metadata$record_id)

table(alpha_F_metadata$diagnosis, alpha_F_metadata$dataset)
#16 Control, 11 T2D from HPAP

#check order
all(colnames(alpha_F_feature) == alpha_F_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_alpha_F <- apply(alpha_F_feature, 1, function(x){sum(!is.na(x)) >= 9})
alpha_F_feature <- alpha_F_feature[feature.keep_alpha_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(alpha_F_metadata$diagnosis))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_alphaF <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = alpha_F_metadata)
colnames(design_alphaF)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_alphaF <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_alphaF <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_alphaF[["levels"]] <- design_alphaF
contrast.matrix_alphaF <- do.call(makeContrasts, myargs_alphaF)

# get results
fit_alphaF <- lmFit(alpha_F_feature, design_alphaF, trend = TRUE, robust = TRUE)
fit_alphaF <- contrasts.fit(fit_alphaF, contrast.matrix_alphaF)
fit_alphaF <- eBayes(fit_alphaF)
res.table_alphaF <- topTable(fit_alphaF, number = Inf)

# Remove results rows with NAs
res.table_alphaF <- res.table_alphaF[!is.na(res.table_alphaF$P.Value), ]

# Save output
res.table_alphaF$entrez <- rownames(res.table_alphaF)
res.table_alphaF$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_alphaF$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_alphaF, "Output/Final figures/FigS4/alpha HPAP only_F_dea_results_correctforage_dataset.csv", row.names = FALSE)

#GSEA
gsea_human(df = res.table_alphaF,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "alpha_HPAPonly_F")

#Visualise
GO_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)
GO_combined %>% filter(NES > 0) %>% summarise(n=n()) #49
GO_combined %>% filter(NES < 0) %>% summarise(n=n()) #252

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
jaccard_similarity <- function(vec1, vec2) {
  intersection <- length(intersect(vec1, vec2))
  union <- length(union(vec1, vec2))
  if (union == 0) return(0)
  return(intersection / union)
}

threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #161
hist(GO_combined_refined$p.adjust)
GO_combined_refined <- GO_combined_refined %>% filter(p.adjust < 0.0001)
dim(GO_combined_refined) #42

GO_combined_refined$direction <- ifelse(GO_combined_refined$NES>0,-1,1)
GO_combined_refined$directionlogpadj <- GO_combined_refined$direction * log10(GO_combined_refined$p.adjust)

#get pathway order
GO_combined_refined_up <- GO_combined_refined %>% filter(NES>0) %>%
  arrange(directionlogpadj)
GO_combined_refined_down <- GO_combined_refined %>% filter(NES<0) %>%
  arrange(directionlogpadj)
pathway_order <- c(unique(GO_combined_refined_down$Description), unique(GO_combined_refined_up$Description))

GO_combined_refined$Description <- factor(GO_combined_refined$Description, levels=pathway_order)
summary(GO_combined_refined$NES)
GO_combined_refined %>%
  ggplot(aes(x=directionlogpadj, y=Description, size = abs(NES), colour = NES))+
  geom_point() +
  theme_bw() +
  labs(x = "Direction signed -log10(p.adjust)", y = "", 
       title = "", colour = "NES", size = "NES") +
  scale_color_gradient2(low = "#55066b", mid = "white", high = "#992002", midpoint = 0) +
  scale_size_continuous(range = c(1,7), limits = c(1,3), breaks = c(1,3)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
ggsave("Output/Final figures/FigS4/FigS4B.tiff", width = 6, height = 6.5)

## Fig S4D ## ----
alpha_M_feature <- alpha_bulk_logCPM_long %>%
  filter(is.na(gene_id)==FALSE) %>%
  dplyr::select(gene_id, DonorID, logCPM) %>%
  pivot_wider(names_from = "DonorID", values_from = "logCPM")
genes_alpha_F <- alpha_M_feature$gene_id
alpha_M_feature <- data.frame(alpha_M_feature)
rownames(alpha_M_feature) <- genes_alpha_F
alpha_M_feature <- alpha_M_feature[,-1] #remove the gene column

alpha_M_metadata <- metadata_combined %>%
  filter(diagnosis %in% c("Control","T2D"), donorsex == "Male") %>%
  filter(record_id %in% colnames(alpha_M_feature))

alpha_M_feature <- alpha_M_feature %>%
  dplyr::select(alpha_M_metadata$record_id)

table(alpha_M_metadata$diagnosis, alpha_M_metadata$dataset)
#24 Control, 7 T2D from HPAP

#check order
all(colnames(alpha_M_feature) == alpha_M_metadata$record_id) #TRUE

#eliminate any feature with fewer than 10 observations
feature.keep_alpha_F <- apply(alpha_M_feature, 1, function(x){sum(!is.na(x)) >= 9})
alpha_M_feature <- alpha_M_feature[feature.keep_alpha_F, ]

# perform analysis
# make design matrix
grp.nms <- sort(unique(alpha_M_metadata$diagnosis))
fixedEffects <- c("donorage") #correct for age
all.vars <- c("diagnosis", fixedEffects)
design_alphaM <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = alpha_M_metadata)
colnames(design_alphaM)[1:length(grp.nms)] <- grp.nms

# make contrast matrix
myargs_alphaM <- list()
ref <- "Control" #using control as reference
contrasts <- grp.nms[grp.nms != ref]
myargs_alphaM <- as.list(paste("T2D", "-", ref, sep = ""))
myargs_alphaM[["levels"]] <- design_alphaM
contrast.matrix_alphaM <- do.call(makeContrasts, myargs_alphaM)

# get results
fit_alphaM <- lmFit(alpha_M_feature, design_alphaM, trend = TRUE, robust = TRUE)
fit_alphaM <- contrasts.fit(fit_alphaM, contrast.matrix_alphaM)
fit_alphaM <- eBayes(fit_alphaM)
res.table_alphaM <- topTable(fit_alphaM, number = Inf)

# Remove results rows with NAs
res.table_alphaM <- res.table_alphaM[!is.na(res.table_alphaM$P.Value), ]

# Save output
res.table_alphaM$entrez <- rownames(res.table_alphaM)
res.table_alphaM$Gene <- mapIds(org.Hs.eg.db, keys=c(res.table_alphaM$entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first") %>% 
  as.character() 
write.csv(res.table_alphaM, "Output/Final figures/FigS4/alpha HPAP only_M_dea_results_correctforage_dataset.csv", row.names = FALSE)

#GSEA
gsea_human(df = res.table_alphaM,
           FC_col = "logFC",
           p_col = "P.Value",
           name = "alpha_HPAPonly_M")

#Visualise
GO_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_CC")
GO_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_MF")
GO_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_BP")
GO_combined <- bind_rows(GO_CC, GO_MF, GO_BP) %>%
  arrange(pvalue)
GO_combined <-GO_combined %>% filter(p.adjust < 0.05)
GO_combined %>% filter(NES > 0) %>% summarise(n=n()) #1
GO_combined %>% filter(NES < 0) %>% summarise(n=n()) #23

#remove redundant pathways
genes_in_pathway <- str_split(GO_combined$core_enrichment, pattern = "/") #returns a list where each list contains a vector of the genes in that pathway
threshold <- 0.6
n <- length(genes_in_pathway)

repeat_indices <- c()
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sim <- jaccard_similarity(genes_in_pathway[[i]], genes_in_pathway[[j]])
    if (sim > threshold) {
      if(length(genes_in_pathway[[i]]) < length(genes_in_pathway[[j]])){
        repeat_indices <- c(repeat_indices,i)} else {repeat_indices <- c(repeat_indices,j)} #keep the one with more genes in it
    }
  }
}

GO_combined_refined <- GO_combined[-repeat_indices,]
dim(GO_combined_refined) #14
hist(GO_combined_refined$p.adjust)
#need to take p< 0.0001 for consistency with the rest
GO_combined_refined <- GO_combined_refined %>% filter(p.adjust < 0.0001) #none


### Figure 5, S5 and S6 ### ----
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

## Fig 5B ## ----
hormone.content %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(Insulin.Content...ng.islet. ~ age_years + sex) #p=0.068 for sex
hormone.content %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(Insulin.Content...ng.islet.), fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 5, label = "p=0.068", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(insulin content (ng/islet))") +
  xlab("") +
  scale_y_continuous(limits = c(-0.5,6))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()     
  )
ggsave("Output/Final figures/Fig5/Fig5B.tiff", width = 3, height = 3)

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
  ggplot(aes(x=sex, y=log(Insulin.Content...ng.islet.)))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(insulin content (ng/islet))") +
  xlab("") +
  scale_y_continuous(limits = c(-5,9))+
  theme_bw() +
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()      
  )
ggsave("Output/Final figures/Fig5/Fig5D.tiff", width = 5, height = 4)

## Fig S6G ## ----
hormone.content_matched <- hormone.content %>%
  filter(simplified_diagnosis != "T1D")
hormone.content_matched$Group <- as.logical(hormone.content_matched$simplified_diagnosis == "T2D")
hormone.content_matched <- matchit(Group ~ age_years + sex,
                                   data = hormone.content_matched,
                                   method = 'nearest',  
                                   ratio = 1, 
                                   exact = ~sex
) 
hormone.content_matched <- match.data(hormone.content_matched)
hormone.content_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #significant for both, switching to optimal doesn't help
hormone.content_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 65, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(20,80))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS6/FigS6G.tiff", width = 5, height = 4)

## Fig S6H ## ----
hormone.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(Insulin.Content...ng.islet. ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- hormone.content_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(Insulin.Content...ng.islet. ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0782 for ctrl vs T2D in females
hormone.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Insulin.Content...ng.islet.)))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 5, label = "p=0.078", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(insulin content (ng/islet))") +
  xlab("") +
  scale_y_continuous(limits = c(-5,9))+
  theme_bw() +
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()      
  )
ggsave("Output/Final figures/FigS6/FigS6H.tiff", width = 5, height = 4)

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

## Fig S5B ## ----
insulin.content_df %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(Insulin.content.ng.postperiIEQ ~ age_years + sex) #ns
insulin.content_df %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(Insulin.content.ng.postperiIEQ), fill = sex))+
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 3.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  scale_y_continuous(limits = c(-0.5,6))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line()     
  )
ggsave("Output/Final figures/FigS5/FigS5B.tiff", width = 3, height = 3)

## Fig S5D ## ----
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
  ggplot(aes(x=sex, y=log(Insulin.content.ng.postperiIEQ)))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 4, label = "p=0.05", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(limits = c(-5,9))+ 
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),  
        axis.line = element_line()       
  )
ggsave("Output/Final figures/FigS5/FigS5D.tiff", width = 5, height = 4)

## Fig S6C ## ----
unique(insulin.content_df$simplified_diagnosis)
insulin.content_df_matched <- insulin.content_df %>%
  filter(simplified_diagnosis != "T1D")
insulin.content_df_matched$Group <- as.logical(insulin.content_df_matched$simplified_diagnosis == "T2D")
insulin.content_df_matched <- matchit(Group ~ age_years + sex,
                                      data = insulin.content_df_matched,
                                      method = 'nearest',  
                                      ratio = 1, 
                                      exact = ~sex
) 
insulin.content_df_matched <- match.data(insulin.content_df_matched)
insulin.content_df_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for males
insulin.content_df_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 68, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  scale_y_continuous(limits = c(20,80))+
  xlab("")+
  ylab("Age (years)") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS6/FigS6C.tiff", width = 5, height = 4)

## Fig S6D ## ----
insulin.content_df_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(Insulin.content.ng.postperiIEQ ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- insulin.content_df_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(Insulin.content.ng.postperiIEQ ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0459 for ctrl vs T2D in females
insulin.content_df_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Insulin.content.ng.postperiIEQ)))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9), aes(fill= simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(limits = c(-5,9))+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),  
        axis.line = element_line()       
  )
ggsave("Output/Final figures/FigS6/FigS6D.tiff", width = 5, height = 4)


#Insulin content from Humanislets.com
#load data
donor_data <- read.csv("data/Humanislets.com/donor.csv")
isolation_data <- read.csv("data/Humanislets.com/isolation.csv")
donor_isolation <- inner_join(donor_data, isolation_data, by = "record_id")

hist(donor_isolation$insulinperieq)
hist(log(donor_isolation$insulinperieq))

## Fig S5A ## ----
donor_isolation %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(log(insulinperieq) ~ donorage + donorsex) #ns
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 5.25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  scale_y_continuous(limits = c(-0.5,6))+
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS5/FigS5A.tiff", width = 3, height = 3)

## Fig S5C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  scale_y_continuous(limits = c(-5,9))+ 
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS5/FigS5C.tiff", width = 5, height = 4)

## Fig S6A ## ----
unique(donor_isolation$diagnosis)
donor_isolation_matched <- donor_isolation %>%
  filter(diagnosis != "Type1") %>%
  filter(insulinperieq > 0)
donor_isolation_matched$Group <- as.logical(donor_isolation_matched$diagnosis == "Type2")
donor_isolation_matched <- matchit(Group ~ donorage + donorsex,
                                   data = donor_isolation_matched,
                                   method = 'nearest',  
                                   ratio = 1, 
                                   exact = ~donorsex
) 

donor_isolation_matched <- match.data(donor_isolation_matched)
table(donor_isolation_matched$donorsex, donor_isolation_matched$diagnosis)

donor_isolation_matched %>%
  group_by(donorsex) %>%
  t_test(donorage ~ diagnosis) #ns

donor_isolation_matched %>%
  ggplot(aes(x=donorsex, y=donorage))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 77, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 80, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(20,80))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS6/FigS6A.tiff", width = 5, height = 4)

## Fig S6B ## ----
donor_isolation_matched %>%
  filter(diagnosis != "Type1") %>%
  mutate(output = log(insulinperieq)) %>%
  mutate(output = case_when(output == -Inf ~ NA, #remove undetectable results that become -Inf when log-transformed
                            .default = output)) %>%
  filter(is.na(output) == FALSE) %>%
  anova_test(output ~ donorage + donorsex*diagnosis) #ns
model <- donor_isolation_matched %>% filter(diagnosis != "Type1") %>%
  mutate(output = log(insulinperieq)) %>%
  mutate(output = case_when(output == -Inf ~ NA,
                            .default = output)) %>%
  filter(is.na(output) == FALSE) %>%
  lm(output ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

donor_isolation_matched %>%
  filter(diagnosis != "Type1") %>%
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "None" ~ "Control",
    .default = diagnosis
  )) %>%
  ggplot(aes(x=donorsex, y=log(insulinperieq), fill= diagnosis))+
  geom_boxplot(alpha = 0.5,position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(insulin content (ng/IEQ))") +
  xlab("") +
  scale_y_continuous(limits = c(-5,9))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS6/FigS6B.tiff", width = 5, height = 4)

# Combined Vanderbilt and Humanislets.com insulin content
Vanderbiltinscontent <- insulin.content_df %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years, Insulin.content.ng.postperiIEQ)
colnames(Vanderbiltinscontent) <- c("DonorID", "sex","diagnosis","age","inscontent")
Vanderbiltinscontent$dataset <- rep("Vanderbilt",nrow(Vanderbiltinscontent))
HIinscontent <- donor_isolation %>%
  dplyr::select(record_id, donorsex, diagnosis, donorage, insulinperieq)
colnames(HIinscontent) <- c("DonorID", "sex","diagnosis","age","inscontent")
HIinscontent$dataset <- rep("Humanislets",nrow(HIinscontent))
combined_inscontent <- combined_inscontent %>%
  mutate(diagnosis = case_when(
    diagnosis == "Type2" ~ "T2D",
    diagnosis == "Type1" ~ "T1D",
    diagnosis == "None" ~ "Control",
    .default = diagnosis
  ))
combined_inscontent <- bind_rows(Vanderbiltinscontent, HIinscontent)
hist(combined_inscontent$inscontent)
hist(log(combined_inscontent$inscontent))

## Fig 5A ## ----
combined_inscontent %>%
  filter(diagnosis == "Control") %>%
  filter(age > 14, age < 40) %>%
  anova_test(log(inscontent) ~ age + dataset + sex) #ns
combined_inscontent %>%
  filter(is.na(inscontent) == FALSE, inscontent > 0) %>%
  filter(dataset %in% c("Vanderbilt", "Humanislets")) %>%
  filter(age > 14, age < 40) %>%
  filter(diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=log(inscontent))) +
  geom_boxplot(alpha = 0.5, aes(fill = sex))+
  scale_fill_manual(values = c("#A7D9DB","#D98B64")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, aes(fill = dataset)) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  scale_y_continuous(limits = c(-0.5,6))+
  xlab("") +
  ylab("log(insulin content (ng/IEQ))") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig5/Fig5A.tiff", width = 3, height = 3)

## Fig 5C ## ----
combined_inscontent %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  filter(dataset %in% c("Vanderbilt","Humanislets")) %>%
  filter(inscontent > 0) %>%
  anova_test(log(inscontent) ~ age + dataset + sex*diagnosis) #significant for dataset only
model <- combined_inscontent %>% 
  filter(diagnosis %in% c("Control","T2D")) %>%
  filter(dataset %in% c("Vanderbilt","Humanislets")) %>%
  filter(inscontent > 0) %>%
  lm(log(inscontent) ~ age + dataset + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns, p=0.0932 for F-M in T2D
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for females

combined_inscontent %>%
  filter(is.na(inscontent) == FALSE, inscontent > 0) %>%
  filter(diagnosis != "T1D") %>%
  filter(dataset %in% c("Vanderbilt","Humanislets")) %>%
  ggplot(aes(x=sex, y=log(inscontent))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  ylab("log(insulin content (ng/IEQ))") +
  xlab("")+
  scale_y_continuous(limits = c(-5,9))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig5/Fig5C.tiff", width = 5, height = 4)

## Fig S6E ## ----
combined_inscontent_matched <- combined_inscontent %>%
  filter(diagnosis != "T1D") %>%
  filter(inscontent > 0)
combined_inscontent_matched$Group <- as.logical(combined_inscontent_matched$diagnosis == "T2D")
combined_inscontent_matched <- matchit(Group ~ age + dataset + sex,
                                       data = combined_inscontent_matched,
                                       method = 'nearest',  
                                       ratio = 1, 
                                       exact = ~sex + dataset
) 

combined_inscontent_matched <- match.data(combined_inscontent_matched)
table(combined_inscontent_matched$sex, combined_inscontent_matched$diagnosis)
combined_inscontent_matched %>%
  group_by(sex) %>%
  anova_test(age ~ dataset + diagnosis) #sig for dataset for both but not age

combined_inscontent_matched %>%
  ggplot(aes(x=sex, y=age))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 77, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 80, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(20,80))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("Output/Final figures/FigS6/FigS6E.tiff", width = 5, height = 4)

## Fig S6F ## ----
combined_inscontent_matched %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  filter(inscontent > 0) %>%
  anova_test(log(inscontent) ~ age + dataset + sex*diagnosis) #significant for dataset only
model <- combined_inscontent_matched %>% 
  filter(diagnosis %in% c("Control","T2D")) %>%
  filter(inscontent > 0) %>%
  lm(log(inscontent) ~ age + dataset + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #no longer sig for females

combined_inscontent_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(inscontent))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  ylab("log(insulin content (ng/IEQ))") +
  xlab("")+
  scale_y_continuous(limits = c(-5,9))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS6/FigS6F.tiff", width = 5, height = 4)


### Figure 6, S7 and S8 ### ----
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
  filter(Time.min > 97, Time.min < 98) %>% #capturing period before the uptick in respiration at transition to high glucose
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
#View(OCR_all)

OCR_all <- OCR_all %>%
  left_join(donors, by = c("DonorID" = "donor_ID"))
OCR_all <- OCR_all %>%
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
  mutate(spare.resp.max = maxFCCP - lastbasal, spare.resp.max.fcbaseline = (maxFCCP - lastbasal)/lastbasal, spare.resp.max.fclastLG = (maxFCCP - lastLG)/lastLG, spare.resp.max.fc.meanLG = (maxFCCP - AAM.3mMGlc)/AAM.3mMGlc)

#Add metadata
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
  label = c("G 0","AAM","G 3
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/Fig6/Fig6A.tiff", width = 5.5, height = 4)

## Fig 6F ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,50,80,105,135,165),   
  xmax = c(50,200,105,200,200,200),
  label = c("G 0","AAM","G 3","G 16.7","FCCP","NaN3"),
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
  ylab ("OCR (nmol/min/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  coord_cartesian(xlim = c(0, 200), ylim = c(0,1)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(50, 80, 105, 135, 165), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin/100, ymax = ymax/100),
    inherit.aes = FALSE,
    fill = NA, colour = "black",linewidth = 0.25
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/Fig6/Fig6F.tiff", width = 5.5, height = 4)

## Fig S7A ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(lastLG ~ age_years + sex) #ns
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=lastLG, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 0.42, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_y_continuous(limits = c(0,0.5))+
  xlab("") +
  ylab("Baseline (3 mM glucose)\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS7/FigS7A.tiff", width = 3, height = 3)

## Fig 6G ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(lastLG ~ age_years + sex*simplified_diagnosis) #0.098 for interaction
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(lastLG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0758 for F-M T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0394 for Ctrl-T2D in males

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=lastLG)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.625, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.8, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 0.925, label = "p=0.076", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 1.05, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Baseline (3 mM glucose)\n(nmol/min/100 islets)") +
  xlab("")+
  scale_y_continuous(limits = c(0,1.1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6G.tiff", width = 4.5, height = 4)

## Fig S7C ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(max.resp ~ age_years + sex) #ns
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=max.resp, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 0.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_y_continuous(limits = c(-0.2,0.5))+
  xlab("") +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS7/FigS7C.tiff", width = 3, height = 3)

## Fig 6H ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(max.resp != -Inf) %>%
  anova_test(max.resp ~ age_years + sex*simplified_diagnosis) #ns
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  filter(max.resp != -Inf) %>%
  lm(max.resp ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=max.resp)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.225, xmax = 2.225, y.position = 1.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 0.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.3, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 0.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.775, y.position = 1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6H.tiff", width = 4.5, height = 4)

## Fig S7E ## ----
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(spare.resp.max.fclastLG ~ age_years + sex) #ns for sex
OCR_summary_v2 %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fclastLG*100, fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  xlab("") +
  ylab("Spare respiratory capacity\n(% of 3 mM glucose)") +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 180, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  ylim(0,200) +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS7/FigS7E.tiff", width = 3, height = 3)

## Fig 6I ## ----
OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(spare.resp.max.fclastLG ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- OCR_summary_v2 %>% filter(simplified_diagnosis != "T1D") %>%
  lm(spare.resp.max.fclastLG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #p=0.0251 for F-M in Control, p=0.0413 for F-M in T2D
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0253 for Ctrl-T2D in males

OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fclastLG*100)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 330, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 400, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.08))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 300, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.5))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 250, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Spare respiratory capacity\n(% of 3 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(0,450))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6I.tiff", width = 4.5, height = 4)

## Fig 6J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.12, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.35))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 0.14, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 0.16, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("High glucose-stimulated\nOCR (nmol/min/100 islets)") +
  xlab("")+
  scale_y_continuous(limits = c(-0.1, 0.16)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6J.tiff", width = 4.5, height = 4)

## Fig S8A ## ----
unique(OCR_summary_v2$simplified_diagnosis)
OCR_summary_v2_matched <- OCR_summary_v2 %>%
  filter(simplified_diagnosis != "T1D")
OCR_summary_v2_matched$Group <- as.logical(OCR_summary_v2_matched$simplified_diagnosis == "T2D")
OCR_summary_v2_matched <- matchit(Group ~ age_years + sex,
                                  data = OCR_summary_v2_matched,
                                  method = 'nearest',  
                                  ratio = 1, 
                                  exact = ~sex
) 
OCR_summary_v2_matched <- match.data(OCR_summary_v2_matched)
OCR_summary_v2_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #ns
OCR_summary_v2_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(20,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS8/FigS8A.tiff", width = 5, height = 4)

## Fig S8B ## ----
stimulus_data3 <- data.frame(
  xmin = c(0,50,80,105,135,165),   
  xmax = c(50,200,105,200,200,200),
  label = c("G 0","AAM","G 3","G 16.7","FCCP","NaN3"),
  ymin = c(70,70,77,77,84,91), 
  ymax = c(77,77,84,84,91,98)
)

OCR_all_over_time_summary %>%
  filter(DonorID %in% OCR_summary_v2_matched$DonorID) %>%
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
  ylab ("OCR (nmol/min/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  coord_cartesian(xlim = c(0, 200), ylim = c(0,1)) +  
  facet_wrap(~sex)+
  geom_vline(xintercept = c(50, 80, 105, 135, 165), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin/100, ymax = ymax/100),
    inherit.aes = FALSE,
    fill = NA, colour = "black",size = 0.25
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin/100+ymax/100)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS8/FigS8B.tiff", width = 5.5, height = 4)

## Fig S8C ## ----
OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(lastLG ~ age_years + sex*simplified_diagnosis) #ns for sex
model <- OCR_summary_v2_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(lastLG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=lastbasal)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.55, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.75, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Baseline (3 mM glucose)\n(nmol/min/100 islets)") +
  xlab("")+
  labs(fill = "Condition")+
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS8/FigS8C.tiff", width = 5, height = 4)

## Fig S8D ## ----
OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(max.resp != -Inf) %>%
  anova_test(max.resp ~ age_years + sex*simplified_diagnosis) #ns
model <- OCR_summary_v2_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(maxFCCP ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=maxFCCP)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 1.05, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.5, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  xlab("")+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS8/FigS8D.tiff", width = 5, height = 4)

## Fig S8E ## ----
OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(spare.resp.max.fclastLG ~ age_years + sex*simplified_diagnosis) #sig for interaction
model <- OCR_summary_v2_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(spare.resp.max.fclastLG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0392 for Ctrl-T2D in males, p=0.0594 for females

OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fclastLG*100)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 170, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 225, label = "p=0.059", label.size = 4, size = 0.5, tip.length = c(0.55, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("Spare respiratory capacity\n(% of 3 mM glucose)") +
  xlab("")+
  labs(fill = "Condition")+
  scale_y_continuous(limits = c(0,250))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS8/FigS8E.tiff", width = 5, height = 4)

## Fig S8F ## ----
OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(response.highglc ~ age_years + sex*simplified_diagnosis) #sig for disease and sex*disease interaction
model <- OCR_summary_v2_matched %>% filter(simplified_diagnosis != "T1D") %>%
  lm(response.highglc ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #still sig for Ctrl-T2D in males

OCR_summary_v2_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=response.highglc)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = simplified_diagnosis))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.12, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.4))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("High glucose-stimulated\nOCR (nmol/min/100 islets)") +
  xlab("")+
  labs(fill = "Condition")+
  scale_y_continuous(limits = c(-0.1, 0.16)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS8/FigS8F.tiff", width = 5, height = 4)

#Humanislets OCR data
#read metadata and data
donor_data <- read.csv("data/Humanislets.com/donor.csv")
seahorse_norm_dna <- read.csv("data/Humanislets.com/seahorse_norm_dna.csv")
seahorse_norm_dna_baseline <- read.csv("data/Humanislets.com/seahorse_norm_dna_baselineoc.csv")

#first working with unbaselined data
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

## Fig 6B ## ----
#Make long format
colnames(seahorse_norm_dna_summary)
seahorse_norm_dna_summary_long <- seahorse_norm_dna_summary %>% 
  pivot_longer(cols = 17:48, values_to = "OCR", names_to = "Time.min")
seahorse_norm_dna_summary_long$Time.min <- gsub("time_","",seahorse_norm_dna_summary_long$Time.min)
seahorse_norm_dna_summary_long$Time.min <- as.numeric(seahorse_norm_dna_summary_long$Time.min)

#Convert raw values to nmol/min/100 islets
seahorse_norm_dna_summary_long <- seahorse_norm_dna_summary_long %>%
  mutate(OCR_nmolmin100islets = (OCR/1000)*(meta_dna_cont/70)*100)

stimulus_data3 <- data.frame(
  xmin = c(0,35,85,145,205), 
  xmax = c(35,85,145,205,260), 
  label = c("G 2.8","G 16.7","Oligomycin","FCCP","Rotenone/\nAntA"),
  ymin = c(rep(0.62,4),0.595),  
  ymax = c(rep(0.67,4),0.695)
)
seahorse_norm_dna_summary_long %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ungroup() %>%
  group_by(Time.min, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(OCR_nmolmin100islets, na.rm = TRUE),
    sd = sd(OCR_nmolmin100islets, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("OCR (nmol/min/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Sex", fill = "Sex") +
  coord_cartesian(xlim = c(0, 270), ylim = c(0,0.75)) +  
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/Fig6/Fig6B.tiff", width = 5.5, height = 4)

# Convert units of summary data
seahorse_norm_dna_summary <- seahorse_norm_dna_summary %>%
  mutate(basalnmolper100islets = (calc_basal_resp/1000)*(meta_dna_cont/70)*100,
         maxnmolper100islets = (calc_max_resp/1000)*(meta_dna_cont/70)*100)

## Fig S7B ## ----
seahorse_norm_dna_summary %>%
  ungroup() %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(basalnmolper100islets ~ donorage + donorsex) #ns
seahorse_norm_dna_summary %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ggplot(aes(x=donorsex, y=basalnmolper100islets, fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 0.35, label = "ns", label.size = 4, size = 0.5,tip.length = c(0.02, 0.1))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  scale_y_continuous(limits = c(0,0.5))+
  xlab("") +
  ylab("Baseline (2.8 mM glucose)\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS7/FigS7B.tiff", width = 3, height = 3)

## Fig S7D ## ----
seahorse_norm_dna_summary %>%
  ungroup() %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  anova_test(maxnmolper100islets ~ donorage + donorsex) #ns
seahorse_norm_dna_summary %>%
  filter(donorage > 14, donorage < 40) %>%
  filter(diagnosis == "None") %>%
  ggplot(aes(x=donorsex, y=maxnmolper100islets, fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 0.4, label = "ns", label.size = 4, size = 0.5,tip.length = c(0.02, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  scale_y_continuous(limits = c(-0.2,0.5))+
  xlab("") +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS7/FigS7D.tiff", width = 3, height = 3)

## Fig S7F ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 80, label = "p=0.08", label.size = 4, size = 0.5,tip.length = c(0.4, 0.1))+
  xlab("") +
  ylab("Spare respiratory capacity\n(% of 2.8 mM glucose)") +
  scale_y_continuous(limits = c(0,400))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  
    axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS7/FigS7F.tiff", width = 3, height = 3)

#Combine HPAP and Humanislets.com OCR data for young control donors
#spare capacity
HPAP_sparecap <- OCR_summary_v2 %>%
  dplyr::select(DonorID, sex, age_years, simplified_diagnosis, spare.resp.max.fclastLG)
HPAP_sparecap$dataset <- rep("HPAP", nrow(HPAP_sparecap))
Humanislets_sparecap <- seahorse_norm_dna_baseline_summary %>%
  dplyr::select(record_id, donorsex, donorage, diagnosis, calc_spare_cap)
Humanislets_sparecap$dataset <- rep("Humanislets", nrow(Humanislets_sparecap))
colnames(Humanislets_sparecap) <- colnames(HPAP_sparecap)
combined_sparecap_summary <- bind_rows(HPAP_sparecap, Humanislets_sparecap)
head(combined_sparecap_summary)
combined_sparecap_summary <- combined_sparecap_summary %>%
  mutate(simplified_diagnosis = case_when(
    simplified_diagnosis == "None" ~ "Control",
    simplified_diagnosis == "Type2" ~ "T2D",
    .default = simplified_diagnosis
  ))
#other parameters
HPAP_OCR <- OCR_summary_v2 %>%
  dplyr::select(DonorID, sex, age_years, simplified_diagnosis, lastLG, max.resp)
HPAP_OCR$dataset <- rep("HPAP", nrow(HPAP_OCR))
Humanislets_OCR <- seahorse_norm_dna_summary %>%
  dplyr::select(record_id, donorsex, donorage, diagnosis, basalnmolper100islets, maxnmolper100islets)
Humanislets_OCR$dataset <- rep("Humanislets", nrow(Humanislets_OCR))
colnames(Humanislets_OCR) <- colnames(HPAP_OCR)
combined_OCR_summary <- bind_rows(HPAP_OCR, Humanislets_OCR)
head(combined_OCR_summary)
combined_OCR_summary <- combined_OCR_summary %>%
  mutate(simplified_diagnosis = case_when(
    simplified_diagnosis == "None" ~ "Control",
    simplified_diagnosis == "Type2" ~ "T2D",
    .default = simplified_diagnosis
  ))

## Fig 6C ## ----
combined_OCR_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(lastLG ~ dataset + age_years + sex) #sig for dataset, ns for sex (p=0.343)
combined_OCR_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=lastLG)) +
  geom_boxplot(alpha = 0.5, aes(fill = sex))+
  scale_fill_manual(values = c("#A7D9DB","#D98B64")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, aes(fill = dataset)) +
  scale_fill_manual(values = c("#50B5AD","#113ED1")) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 0.45, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  xlab("") +
  ylab("Baseline (~3 mM glucose)\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6C.tiff", width = 3, height = 3)

## Fig 6D ## ----
combined_OCR_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(max.resp ~ dataset + age_years + sex) #ns
combined_OCR_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=max.resp)) +
  geom_boxplot(alpha = 0.5, aes(fill = sex))+
  scale_fill_manual(values = c("#A7D9DB","#D98B64")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, aes(fill = dataset)) +
  scale_fill_manual(values = c("#50B5AD","#113ED1")) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 0.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  xlab("") +
  ylab("Maximal respiratory capacity\n(nmol/min/100 islets)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6D.tiff", width = 3, height = 3)

## Fig 6E ## ----
combined_sparecap_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  anova_test(spare.resp.max.fclastLG ~ dataset + age_years + sex) #sig for dataset and sex
combined_sparecap_summary %>%
  filter(age_years > 14, age_years < 40) %>%
  filter(simplified_diagnosis == "Control") %>%
  ggplot(aes(x=sex, y=spare.resp.max.fclastLG*100)) +
  geom_boxplot(alpha = 0.5, aes(fill = sex))+
  scale_fill_manual(values = c("#A7D9DB","#D98B64")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, aes(fill = dataset)) +
  scale_fill_manual(values = c("#50B5AD","#113ED1")) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 180, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_y_continuous(limits = c(0,200))+
  xlab("") +
  ylab("Spare respiratory capacity\n(% of ~3 mM glucose)") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig6/Fig6E.tiff", width = 3, height = 3)


## Figure S9-10 ## ----
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
#AAM
calculate.auc.AAM <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)  
  time_points <- seq_len(n) + 8 
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
  filter(adj.time >= 9, adj.time <= 16) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
AAMdf$AUC.AAM <- apply(AAMdf[2:9], 1, calculate.auc.AAM)
AUC_summary <- AAMdf %>% dplyr::select(Donor, AUC.AAM)

#LG
calculate.auc.LG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  time_points <- seq_len(n) + 15 
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
  filter(adj.time >= 16, adj.time <= 23) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
LGdf$AUC.LG <- apply(LGdf[2:9], 1, calculate.auc.LG)
AUC_summary <- LGdf %>% dplyr::select(Donor, AUC.LG) %>% full_join(AUC_summary, by = "Donor")

#HG
calculate.auc.HG <- function(x) {
  n <- length(x)
  if (n < 2) return(NA) 
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
  filter(adj.time >= 23, adj.time <= 39) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
HGdf$AUC.HG <- apply(HGdf[2:18], 1, calculate.auc.HG)
AUC_summary <- HGdf %>% dplyr::select(Donor, AUC.HG) %>% full_join(AUC_summary, by = "Donor")

#KCl
calculate.auc.KCl <- function(x) {
  n <- length(x)
  if (n < 2) return(NA) 
  time_points <- seq_len(n) + 38 
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
  filter(adj.time >= 39, adj.time <= 45) %>%
  dplyr::select(Donor, adj.time, mean_baselined_signal) %>%
  pivot_wider(names_from = adj.time, values_from = mean_baselined_signal)
KCldf$AUC.KCl <- apply(KCldf[2:7], 1, calculate.auc.KCl)
AUC_summary <- KCldf %>% dplyr::select(Donor, AUC.KCl) %>% full_join(AUC_summary, by = "Donor")

#Add metadata
AUC_summary <- left_join(AUC_summary, donor_info, by = c("Donor"="donor_ID"))

## Fig S10A ## ----
stimulus_data <- data.frame(
  xmin = c(0,9,16,23,27,39),
  xmax = c(9,16,23,27,39,45),
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS10/FigS10A.tiff", width = 7, height = 4)

## Fig S10B ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 2.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.075, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 2.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC AAM") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS10/FigS10B.tiff",width = 5, height = 4)

## Fig S10C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 2.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.3, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 2.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC 3 mM Glucose") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS10/FigS10C.tiff",width = 5, height = 4)

## Fig S10D ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.HG ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.HG ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0736 for Ctrl vs T2D in females, p=0.19 in males

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.HG,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "p=0.074", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 10.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC 16.7 mM Glucose") +
  scale_y_continuous(limits= c(-1,10.7))+
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS10/FigS10D.tiff",width = 5, height = 4)

## Fig S10E ## ----
AUC_summary %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(AUC.KCl ~ age_years + sex*simplified_diagnosis) #p=0.058 for disease
model <- AUC_summary %>% filter(simplified_diagnosis != "T1D") %>%
  lm(AUC.KCl ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | simplified_diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.10 for females, p=0.21 for males

AUC_summary %>%
  filter(simplified_diagnosis %in% c("Control","T2D")) %>%
  ggplot(aes(x=sex, y=AUC.KCl,  fill = simplified_diagnosis))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("AUC KCl") +
  xlab("")+
  scale_y_continuous(limits= c(0,8.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS10/FigS10E.tiff",width = 5, height = 4)

## Fig S9A ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,9,16,23,27,39),
  xmax = c(9,16,23,27,39,45),
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS9/FigS9A.tiff", width = 5.5, height = 4)

## Fig S9B ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 2.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC AAM") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS9/FigS9B.tiff",width = 3, height = 3)

## Fig S9C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 2.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC 3 mM Glucose") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS9/FigS9C.tiff",width = 3, height = 3)

## Fig S9D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 8.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC 16.7 mM Glucose") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS9/FigS9D.tiff",width = 3, height = 3)

## Fig S9E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("AUC KCl") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line()    
  )
ggsave("Output/Final figures/FigS9/FigS9E.tiff",width = 3, height = 3)

### Figure 7 and S11-14, S17-19 ### ----

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
  summarise(mean.ins.baseline = mean(Insulin.Release..ng.100.islets.min., na.rm = TRUE))
insulinAUCs <- baselineins %>% dplyr::select(DonorID, mean.ins.baseline) %>% full_join(insulinAUCs, by = "DonorID")

baselinegcg <- perifusion_all_with_metadata %>%
  filter(Stimulus == "No Stimuli", Time..min. < 30) %>%
  group_by(DonorID) %>%
  summarise(mean.gcg.baseline = mean(Glucagon.Release..pg.100.islets.min., na.rm = TRUE))
glucagonAUCs <- baselinegcg %>% dplyr::select(DonorID, mean.gcg.baseline) %>% full_join(glucagonAUCs, by = "DonorID")

#Add metadata to AUCs
perifusion_summary <- full_join(insulinAUCs, donors, by = c("DonorID" = "donor_ID"))
perifusion_summary <- full_join(perifusion_summary, glucagonAUCs, by = "DonorID")
head(perifusion_summary)
tail(colnames(perifusion_summary))

## Fig 7A and S13A ## ----
stimulus_data <- data.frame(
  xmin = c(20,30,61,81,101,121,141), 
  xmax = c(30,121,81,121,121,141,160),
  label = c("G 0","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(8.9,8.9,9.5,9.5,10.1,9.8,10.1), 
  ymax = c(9.5,9.5,10.1,10.1,10.7,11,10.7))

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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/Fig7/Fig7A.tiff", width = 8, height = 4)
ggsave("Output/Final figures/FigS13/FigS13A.tiff", width = 8, height = 4)

## Fig S13B ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 1.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 3, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion \n(ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13B.tiff", width = 4 , height = 4)

## Fig S13C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  scale_y_continuous(limits = c(-1,5.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13C.tiff", width = 4 , height = 4)

## Fig S13D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 3.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 4.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 4.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n3 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13D.tiff", width = 4 , height = 4)

## Fig S13E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.1, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 6.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(-0.5,6.1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13E.tiff", width = 4 , height = 4)

## Fig S13F ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  scale_y_continuous(limits = c(0,7)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13F.tiff", width = 4 , height = 4)

## Fig S13G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(0,5.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS13/FigS13G.tiff", width = 4 , height = 4)

#Age-matching
## Fig S14B ### ----
perifusion_matched <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(mean.ins.baseline > 0) %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct()
perifusion_matched$Group <- as.logical(perifusion_matched$simplified_diagnosis == "T2D")
perifusion_matched <- matchit(Group ~ age_years + sex,
                              data = perifusion_matched,
                              method = 'nearest',
                              ratio = 1, 
                              exact = ~sex
) 
perifusion_matched <- match.data(perifusion_matched)
table(perifusion_matched$sex, perifusion_matched$simplified_diagnosis)

perifusion_all_with_metadata_matched <- perifusion_all_with_metadata %>%
  filter(DonorID %in% perifusion_matched$DonorID)

perifusion_all_with_metadata_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct() %>%
  group_by(sex, simplified_diagnosis) %>%
  summarise(n=n())

perifusion_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for both
perifusion_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(20,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14B.tiff", width = 4, height = 4)

## Fig S14A ## ----
perifusion_all_with_metadata_matched %>%
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS14/FigS14A.tiff", width = 8, height = 4)

## Fig S14C ## ----
perifusion_summary_matched <- perifusion_summary %>%
  filter(DonorID %in% perifusion_matched$DonorID)

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct() %>%
  group_by(sex, simplified_diagnosis) %>%
  summarise(n=n())

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(mean.ins.baseline) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(mean.ins.baseline) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0114 for Control vs T2D in females, ns for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean.ins.baseline), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 1.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion \n(ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-3,2))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14C.tiff", width = 4, height = 4)

## Fig S14D ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.AAM) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.AAM) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0098 for Control vs T2D in females, p=0.0804 for males

perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.AAM), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.075, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  scale_y_continuous(limits = c(-1,5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14D.tiff", width = 4, height = 4)

## Fig S14E ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.LG) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.LG) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0619 for Control vs T2D in females, ns for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.LG), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.6, label = "p=0.06", label.size = 4, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 3.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n3 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14E.tiff", width = 4, height = 4)

## Fig S14F ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.HG) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.HG) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0002 for Control vs T2D in females, p=0.0093 for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.HG), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.1, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(-0.5,5.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14F.tiff", width = 4, height = 4)

## Fig S14G ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.IBMX) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0002 for Control vs T2D in females, p=0.001 for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  scale_y_continuous(limits = c(0,6))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14G.tiff", width = 4, height = 4)

## Fig S14H ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.ins.KCl) ~ age_years + sex*simplified_diagnosis) #sig for disease and age
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.ins.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0008 for Control vs T2D in females, p=0.0296 for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.ins.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(0,5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS14/FigS14H.tiff", width = 4, height = 4)

## Fig S11A ## ----
stimulus_data2 <- data.frame(
  xmin = c(20,30,61,81,101,121,141),
  xmax = c(30,121,81,121,121,141,160),
  label = c("G 0","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(8.9,8.9,9.5,9.5,10.1,9.8,10.1), 
  ymax = c(9.5,9.5,10.1,10.1,10.7,11,10.7))

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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS11/FigS11A.tiff", width = 5.5, height = 4)

## Fig S11B ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(mean.ins.baseline) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(mean.ins.baseline), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline insulin secretion \n(ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11B.tiff", width = 3 , height = 3)

## Fig S11C ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.AAM) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.AAM), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11C.tiff", width = 3 , height = 3)

## Fig S11D ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.LG) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.LG), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 3.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n3 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11D.tiff", width = 3 , height = 3)

## Fig S11E ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.HG) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.HG), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11E.tiff", width = 3, height = 3)

## Fig S11F ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.IBMX) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.IBMX), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 5.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11F.tiff", width = 3, height = 3)

## Fig S11G ## ----
perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  anova_test(log(AUC.ins.KCl) ~ age_years + sex) #ns

perifusion_summary %>%
  filter(simplified_diagnosis == "Control") %>%
  filter(age_years > 14, age_years < 40) %>%
  ggplot(aes(x=sex, y=log(AUC.ins.KCl), fill = sex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 4.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS11/FigS11G.tiff", width = 3, height = 3)

## Fig S17A ## ----
stimulus_data3 <- data.frame(
  xmin = c(20,30,61,81,101,121,141),    
  xmax = c(30,121,81,121,121,141,160),
  label = c("G 0","AAM","G 3","G 16.7","IBMX","wash\nout","KCl"),
  ymin = c(350,350,370,370,400,385,400), 
  ymax = c(370,370,400,400,430,445,430))
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS17/FigS17A.tiff", width = 5.5, height = 4)

## Fig S17C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.125))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline glucagon\nsecretion (pg/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS17/FigS17C.tiff", width = 3 , height = 3)

## Fig S17E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS17/FigS17E.tiff", width = 3 , height = 3)

## Fig S17G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS17/FigS17G.tiff", width = 3 , height = 3)

## Fig S17I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS17/FigS17I.tiff", width = 3 , height = 3)

## Fig S17K ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(glucagon content\n(pg/islet))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()
  )
ggsave("Output/Final figures/FigS17/FigS17K.tiff", width = 3 , height = 3)

## Fig S18A ## ----
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS18/FigS18A.tiff", width = 8, height = 4)

## Fig S18C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 6.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion \n(pg/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6.2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS18/FigS18C.tiff", width = 4 , height = 4)

## Fig S18E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.07))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.9, label = "p=0.09", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 10.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  scale_y_continuous(limits = c(4.5,10.2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS18/FigS18E.tiff", width = 4 , height = 4)

## Fig S18G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 10.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS18/FigS18G.tiff", width = 4 , height = 4)

## Fig S18I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 10, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 10.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(4.5,10.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS18/FigS18I.tiff", width = 4 , height = 4)

## Fig S18K ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 10, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 10.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 11.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/islet))") +
  xlab("")+
  scale_y_continuous(limits = c(3,11.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS18/FigS18K.tiff", width = 5 , height = 4)

#Age-matched glucagon data
## Fig S19 ## ----
perifusion_matched <- perifusion_summary %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(mean.gcg.baseline > 0) %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct()
perifusion_matched$Group <- as.logical(perifusion_matched$simplified_diagnosis == "T2D")
perifusion_matched <- matchit(Group ~ age_years + sex,
                              data = perifusion_matched,
                              method = 'nearest',  
                              ratio = 1, 
                              exact = ~sex
) 
perifusion_matched <- match.data(perifusion_matched)
table(perifusion_matched$sex, perifusion_matched$simplified_diagnosis)

perifusion_all_with_metadata_matched <- perifusion_all_with_metadata %>%
  filter(DonorID %in% perifusion_matched$DonorID)
perifusion_summary_matched <- perifusion_summary %>%
  filter(DonorID %in% perifusion_matched$DonorID)
perifusion_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for both
perifusion_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(20,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19B.tiff", width = 5, height = 4)

## Fig S19A ## ----
perifusion_all_with_metadata_matched %>%
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS19/FigS19A.tiff", width = 8, height = 4)

## Fig S19E ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(mean.gcg.baseline) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(mean.gcg.baseline) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0255 for Control vs T2D in females, ns for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean.gcg.baseline), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion \n(pg/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6.2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19E.tiff", width = 4, height = 4)

## Fig S19G ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.AAM) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.AAM) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0132 for Control vs T2D in females, ns for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.AAM), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.07))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n4 mM AAM)") +
  xlab("")+
  scale_y_continuous(limits = c(4.5,9.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19G.tiff", width = 4, height = 4)

## Fig S19I ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.IBMX) ~ age_years + sex*simplified_diagnosis) #ns
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0775 for Control vs T2D in females, ns for males

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.4, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.15, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  scale_y_continuous(limits=c(4.8,10))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19I.tiff", width = 4, height = 4)

## Fig S19K ## ----
perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.gcg.KCl) ~ age_years + sex*simplified_diagnosis) #sig for age an disease
model <- perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.gcg.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0216 for Control vs T2D in females

perifusion_summary_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.gcg.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(4.5,9.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19K.tiff", width = 4, height = 4)

## Fig S19M ## ----
glucagon.content_matched <- glucagon.content %>%
  filter(DonorID %in% perifusion_matched$DonorID)
table(glucagon.content_matched$sex, glucagon.content_matched$simplified_diagnosis) #numbers don't line up - need to redo matching

glucagon.content_matched <- glucagon.content %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(Glucagon.Content..pg.islet. > 0) %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years, Glucagon.Content..pg.islet.) %>%
  distinct()

glucagon.content_matched$Group <- as.logical(glucagon.content_matched$simplified_diagnosis == "T2D")
glucagon.content_matched <- matchit(Group ~ age_years + sex,
                                    data = glucagon.content_matched,
                                    method = 'nearest',  
                                    ratio = 1, 
                                    exact = ~sex
) 
glucagon.content_matched <- match.data(glucagon.content_matched)
table(glucagon.content_matched$sex, glucagon.content_matched$simplified_diagnosis)

glucagon.content_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for females
glucagon.content_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 65, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(15,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19M.tiff", width = 4, height = 4)

## Fig S19N ## ----
glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(Glucagon.Content..pg.islet.) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(Glucagon.Content..pg.islet.) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for males only, p=0.0637 for females

glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Glucagon.Content..pg.islet.), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 9.2, label = "p=0.06", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 10.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.4))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/islet))") +
  xlab("")+
  scale_y_continuous(limits = c(5,10.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19N.tiff", width = 4, height = 4)

#Humanislets.com glucose perifusion data
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

#convert units
peri_gluc_long <- peri_gluc_long %>%
  mutate(Insulin.release.uU.min = Insulin.release * (2/5))
peri_gluc_long <- peri_gluc_long %>%
  mutate(Insulin.release.uU.islet.min = Insulin.release.uU.min/65,
         Insulin.release.pmol.islet.min = Insulin.release.uU.islet.min*(6/1000),
         Insulin.release.ng.islet.min = Insulin.release.pmol.islet.min*5.808,
         Insulin.release.ng.100islets.min = Insulin.release.ng.islet.min*100)

peri_gluc_means <- peri_gluc_means %>%
  mutate(across(-c(record_id, replicate), ~ .x * (34.848/1625)))

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

## Fig 7B and S13M ## ----
stimulus_data4 <- data.frame(
  xmin = c(0,20,90,160), 
  xmax = c(20,70,140,190), 
  label = c("G 3","G 15","G 6","KCl"),
  ymin = rep(14.5,4),  
  ymax = rep(16,4))

peri_gluc_long %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release.ng.100islets.min, na.rm = TRUE),
    sd = sd(Insulin.release.ng.100islets.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = diagnosis, group = diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = diagnosis)))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab ("Insulin release (ng/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  facet_wrap(~donorsex)+
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data4,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/Fig7/Fig7B.tiff", width = 8, height = 4)
ggsave("Output/Final figures/FigS13/FigS13M.tiff", width = 8, height = 4)

## Fig S13N ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 3.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13N.tiff", width = 4, height = 4)

## Fig S13O ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "p=0.06", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 7.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n15 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13O.tiff", width = 4, height = 4)

## Fig S13P ## ----
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
  geom_quasirandom(size = 2, shape= 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.3, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n6 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13P.tiff", width = 4, height = 4)

## Fig S13Q ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13Q.tiff", width = 4, height = 4)

## Fig S14P ## ----
peri_gluc_matched <- peri_gluc_means %>%
  filter(baseline.mean > 0) %>%
  dplyr::select(record_id, donorsex, diagnosis, donorage) %>%
  distinct()
peri_gluc_matched$Group <- as.logical(peri_gluc_matched$diagnosis == "T2D")
peri_gluc_matched <- matchit(Group ~ donorage + donorsex,
                             data = peri_gluc_matched,
                             method = 'nearest',
                             ratio = 1, 
                             exact = ~donorsex
) 
peri_gluc_matched <- match.data(peri_gluc_matched)
table(peri_gluc_matched$donorsex, peri_gluc_matched$diagnosis)

peri_gluc_long_matched <- peri_gluc_long %>%
  filter(record_id %in% peri_gluc_matched$record_id)

peri_gluc_long_matched %>%
  dplyr::select(record_id, donorsex, diagnosis, donorage) %>%
  distinct() %>%
  group_by(donorsex, diagnosis) %>%
  summarise(n=n())
peri_gluc_matched %>%
  group_by(donorsex) %>%
  t_test(donorage ~ diagnosis) #ns

peri_gluc_matched %>%
  ggplot(aes(x=donorsex, y=donorage))+
  geom_boxplot(aes(fill = diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 69, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  xlab("")+
  ylab("Age (years)") +
  scale_y_continuous(limits = c(40,75))+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14P.tiff", width = 4, height = 4)

## Fig S14O ## ----
peri_gluc_long_matched %>%
  filter(diagnosis %in% c("Control","T2D")) %>%
  ungroup() %>%
  group_by(Time.min, diagnosis, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release.ng.100islets.min, na.rm = TRUE),
    sd = sd(Insulin.release.ng.100islets.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = diagnosis, group = diagnosis)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = diagnosis)))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab ("Insulin release (ng/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  facet_wrap(~donorsex)+
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data4,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS14/FigS14O.tiff", width = 8, height = 4)

## Fig S14Q ## ----
peri_gluc_means_matched <- peri_gluc_means %>%
  filter(record_id %in% peri_gluc_matched$record_id)
peri_gluc_means_matched %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(baseline.mean) ~ donorage + donorsex*diagnosis) #ns
model <- peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  lm(log(baseline.mean) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 1.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-3.5,1.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14Q.tiff", width = 4, height = 4)

## Fig S14R ## ----
peri_gluc_means_matched %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.G15) ~ donorage + donorsex*diagnosis) #sig for diagnosis
model <- peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.G15) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.G15), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n15 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(0,7)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14R.tiff", width = 4, height = 4)

## Fig S14S ## ----
peri_gluc_means_matched %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.G6) ~ donorage + donorsex*diagnosis) #sig for interaction
model <- peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.G6) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #sig for females

peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.G6), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape= 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.8, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(2.5,6.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14S.tiff", width = 4, height = 4)

## Fig S14T ## ----
peri_gluc_means_matched %>%
  ungroup() %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ donorage + donorsex*diagnosis) #ns
model <- peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ donorage + donorsex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | donorsex)
pairs(emmeans_res, adjust = "tukey") #ns

peri_gluc_means_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl), fill = diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(2.5,6))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14T.tiff", width = 4, height = 4)

## Fig S11M ## ----
stimulus_data5 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190),
  label = c("G 3","G 15","G 6","KCl"),
  ymin = c(rep(21,4)),  
  ymax = c(rep(23,4)))

peri_gluc_long %>%
  filter(diagnosis %in% c("Control")) %>%
  filter(donorage > 14, donorage < 40) %>%
  ungroup() %>%
  group_by(Time.min, donorsex) %>%
  summarise(
    n = n(),
    mean = mean(Insulin.release.ng.100islets.min, na.rm = TRUE),
    sd = sd(Insulin.release.ng.100islets.min, na.rm = TRUE),
    up = if (n > 2) mean + sd else NA_real_,
    down = if (n > 2) mean - sd else NA_real_,
    .groups = "drop"
  ) %>%
  ggplot(aes(x=Time.min, y=mean, colour = donorsex, group = donorsex)) +
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = donorsex)))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab ("Insulin Release (ng/100 islets)") +
  xlab ("Time (min)") +
  labs(colour = "Condition", fill = "Condition") +
  geom_vline(xintercept = c(0,20,70,90,140,160,190), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5
  ) +
  geom_text(
    data = stimulus_data5,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS11/FigS11M.tiff", width = 5.5, height = 4)

## Fig S11N ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(baseline.mean ~ donorage + donorsex) #ns

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(baseline.mean),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(baseline insulin secretion\n(ng/100 islets))") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS11/FigS11N.tiff", width = 3, height = 3)

## Fig S11O ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.G15 ~ donorage + donorsex) #p=0.029 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.G15),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.25, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n15 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(3,7))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS11/FigS11O.tiff", width = 3, height = 3)

## Fig S11P ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.G6 ~ donorage + donorsex) #p=0.031 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.G6),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 6.25, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(2,6.5))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS11/FigS11P.tiff", width = 3, height = 3)

## Fig S11Q ## ----
peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  anova_test(AUC.KCl ~ donorage + donorsex) #p=0.043 for sex

peri_gluc_means %>%
  filter(diagnosis == "Control") %>%
  filter(donorage > 14, donorage < 40) %>%
  ggplot(aes(x=donorsex, y=log(AUC.KCl),  fill = donorsex)) +
  geom_boxplot(alpha = 0.5)+
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE, xmin = 1, xmax = 2, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.25, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(2,6.5))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(), 
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS11/FigS11Q.tiff", width = 3, height = 3)

#Combining HPAP UPenn and Humanislets.com perifusions
UPenn_select <- perifusion_all_with_metadata %>%
  dplyr::select(Time..min.,Insulin.Release..ng.100.islets.min.,DonorID, sex,age_years,simplified_diagnosis)
Humanislets_select <- peri_gluc_long %>%
  dplyr::select(Time.min, Insulin.release.ng.100islets.min, record_id, donorsex, donorage, diagnosis)

colnames(UPenn_select) <- c("Time.min","Insulin.orig.units","DonorID","sex","age","diagnosis")
colnames(Humanislets_select) <- colnames(UPenn_select)
Combined <- bind_rows(UPenn_select, Humanislets_select)
Combined$dataset <- c(rep("UPenn",nrow(UPenn_select)),rep("Humanislets",nrow(Humanislets_select)))

## Fig 7C ## ----
mean_3G <- Combined %>%
  filter(dataset != "Vanderbilt") %>%
  filter(dataset == "Humanislets" & Time.min <= 20 | dataset == "UPenn" & Time.min >= 61 & Time.min <= 80) %>%
  group_by(DonorID, sex, age, diagnosis, dataset) %>%
  summarise(mean_3G = mean(Insulin.orig.units, na.rm = TRUE)) %>%
  ungroup()

hist(mean_3G$mean_3G)
hist(log(mean_3G$mean_3G))

mean_3G %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(mean_3G) ~ dataset + age + sex*diagnosis) #sig for dataset, age, and diagnosis
model <- mean_3G %>% filter(diagnosis != "T1D") %>%
  lm(log(mean_3G) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0327 for Ctrl-T2D in males

mean_3G %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean_3G))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 3.25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_y_continuous(limits = c(-4.5,3.25))+
  ylab("log(mean secretion at\n3 mM glucose (ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig7/Fig7C.tiff", width = 5, height = 4)

## Fig 7D ## ----
max_HG <- Combined %>%
  filter(dataset == "Humanislets" & Time.min >= 20 & Time.min <= 70 | dataset == "UPenn" & Time.min >= 81 & Time.min <= 101) %>%
  group_by(DonorID, sex, age, diagnosis, dataset) %>%
  summarise(max_HG = max(Insulin.orig.units, na.rm = TRUE)) %>%
  ungroup()
unique(max_HG$dataset)

hist(max_HG$max_HG)
hist(log(max_HG$max_HG))

max_HG %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(max_HG) ~ dataset + age + sex*diagnosis) #sig for dataset, age, and diagnosis
model <- max_HG %>% filter(diagnosis != "T1D") %>%
  lm(log(max_HG) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

max_HG %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(max_HG))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 4.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_y_continuous(limits = c(-3,5.5))+
  ylab("log(peak secretion at 15 or\n16.7 mM glucose (ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  ) 
ggsave("Output/Final figures/Fig7/Fig7D.tiff", width = 5, height = 4)

## Fig 7E ## ----
max_KCl <- Combined %>%
  filter(dataset != "Vanderbilt") %>%
  filter(dataset == "Humanislets" & Time.min >= 160 & Time.min <= 190 |  dataset == "UPenn" & Time.min >= 141 & Time.min <= 160) %>%
  group_by(DonorID, sex, age, diagnosis, dataset) %>%
  summarise(max_KCl = max(Insulin.orig.units, na.rm = TRUE)) %>%
  ungroup()
unique(max_KCl$dataset)

hist(max_KCl$max_KCl)
hist(log(max_KCl$max_KCl))

max_KCl %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(max_KCl) ~ dataset + age + sex*diagnosis) #sig for dataset, age, and diagnosis
model <- max_KCl %>% filter(diagnosis != "T1D") %>%
  lm(log(max_KCl) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

max_KCl %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(max_KCl))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 3.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 4.75, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_y_continuous(limits = c(-3,4.75))+
  ylab("log(peak secretion at 30 mM KCl\n(ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  ) 
ggsave("Output/Final figures/Fig7/Fig7E.tiff", width = 5, height = 4)

## Fig 7F ## ----
SI <- inner_join(mean_3G, max_HG, by = c("DonorID", "sex", "age", "diagnosis", "dataset")) %>%
  mutate(SI = max_HG/mean_3G)

hist(SI$SI)
hist(log(SI$SI))

SI %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(SI) ~ dataset + age + sex*diagnosis) #sig for dataset and diagnosis
model <- SI %>% filter(diagnosis != "T1D") %>%
  lm(log(SI) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ sex | diagnosis)
pairs(emmeans_res, adjust = "tukey") #ns
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

SI %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(SI))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_y_continuous(limits = c(0,5.5))+
  ylab("log(stimulation index, 15 or\n16.7 mM glucose/3 mM glucose))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/Fig7/Fig7F.tiff", width = 5, height = 4)

#Age-matched
## Fig S12A ## ----
Combined_matched <- SI %>%
  filter(diagnosis != "T1D") %>%
  filter(is.na(SI) == FALSE) %>%
  dplyr::select(DonorID, sex, diagnosis, age, dataset) %>%
  distinct()
Combined_matched$Group <- as.logical(Combined_matched$diagnosis == "T2D")
Combined_matched <- matchit(Group ~ age + sex + dataset,
                            data = Combined_matched,
                            method = 'nearest', # can switch to 'exact' or 'optimal'
                            ratio = 1, 
                            exact = ~sex + dataset
) 
Combined_matched <- match.data(Combined_matched)
table(Combined_matched$sex, Combined_matched$diagnosis, Combined_matched$dataset)

Combined_matched %>%
  group_by(sex) %>%
  anova_test(age ~ dataset + diagnosis) #sig for diagnosis in both sexes, plus dataset effect in males

Combined_matched %>%
  ggplot(aes(x=sex, y=age))+
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 72, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  xlab("")+
  ylab("Age (years)") +
  ylim(20,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS12/FigS12A.tiff", width = 4, height = 4)

## Fig S12B ## ----
mean_3G_matched <- mean_3G %>%
  filter(DonorID %in% Combined_matched$DonorID)

mean_3G_matched %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(mean_3G) ~ dataset + age + sex*diagnosis) #sig for dataset and diagnosis
model <- mean_3G_matched %>% filter(diagnosis != "T1D") %>%
  lm(log(mean_3G) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0498 for Ctrl-T2D in males

mean_3G_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(mean_3G))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 1.4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_y_continuous(limits = c(-4,2.5))+
  ylab("log(mean secretion at\n3 mM glucose (ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS12/FigS12B.tiff", width = 4, height = 4)

## Fig S12C ## ----
max_HG_matched <- max_HG %>%
  filter(DonorID %in% Combined_matched$DonorID)
max_HG_matched %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(max_HG) ~ dataset + age + sex*diagnosis) #sig for diagnosis
model <- max_HG_matched %>% filter(diagnosis != "T1D") %>%
  filter(max_HG > 0) %>%
  lm(log(max_HG) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

max_HG_matched %>%
  filter(diagnosis != "T1D") %>%
  filter(dataset != "Vanderbilt") %>%
  ggplot(aes(x=sex, y=log(max_HG))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2.9, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.05))+
  scale_y_continuous(limits = c(-3,3.5))+
  ylab("log(peak secretion at 15 or\n16.7 mM glucose (ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  ) 
ggsave("Output/Final figures/FigS12/FigS12C.tiff", width = 4, height = 4)

## Fig S12D ## ----
max_KCl_matched <- max_KCl %>%
  filter(DonorID %in% Combined_matched$DonorID)

max_KCl_matched %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(max_KCl) ~ dataset + age + sex*diagnosis) #sig for diagnosis
model <- max_KCl_matched %>% filter(diagnosis != "T1D") %>%
  lm(log(max_KCl) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

max_KCl_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(max_KCl))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.3, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 3, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_y_continuous(limits = c(-3,4.5))+
  ylab("log(peak secretion at 30 mM KCl\n(ng/100 islets))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  ) 
ggsave("Output/Final figures/FigS12/FigS12D.tiff", width = 4, height = 4)

## Fig S12E ## ----
SI_matched <- SI %>%
  filter(DonorID %in% Combined_matched$DonorID)

SI_matched %>%
  filter(diagnosis != "T1D") %>%
  anova_test(log(SI) ~ dataset + age + sex*diagnosis) #sig diagnosis and dataset
model <- SI_matched %>% filter(diagnosis != "T1D") %>%
  lm(log(SI) ~ dataset + age + sex*diagnosis, data = .)
emmeans_res <- emmeans(model, ~ diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both sexes

SI_matched %>%
  filter(diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(SI))) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9), aes(fill = diagnosis))+
  scale_fill_manual(values = c("grey50", "#FFCE2C")) +
  ggnewscale::new_scale_fill() +
  geom_quasirandom(size = 2, shape = 21, colour = "black", alpha = 0.75, aes(group = diagnosis, fill = dataset), dodge.width=0.9) +
  scale_fill_manual(values = c("#113ED1", "#50B5AD")) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 3.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  scale_y_continuous(limits = c(0,4.))+
  ylab("log(stimulation index, 15 or\n16.7 mM glucose/3 mM glucose))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS12/FigS12E.tiff", width = 4, height = 4)

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


## Fig S13H ## ----
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS13/FigS13H.tiff", width = 8, height = 4)

## Fig S13I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.75, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 1.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion\n(ng/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13I.tiff", width = 4, height = 4)

## Fig S13J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.03))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.1, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.14))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13J.tiff", width = 4, height = 4)

## Fig S13K ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.9, label = "p=0.053", label.size = 4, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.2, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.18))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13K.tiff", width = 4, height = 4)

## Fig S13L ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.22))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 3.8, label = "p=0.07", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 4.6, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.1, 0.07))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 5.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,5.1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS13/FigS13L.tiff", width = 4, height = 4)

## Fig S14J ## ----
Vanderbilt_perifusion_matched <- all_insulin_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct()

Vanderbilt_perifusion_matched$Group <- as.logical(Vanderbilt_perifusion_matched$simplified_diagnosis == "T2D")
Vanderbilt_perifusion_matched <- matchit(Group ~ age_years + sex,
                                         data = Vanderbilt_perifusion_matched,
                                         method = 'nearest', # can switch to 'exact' or 'optimal'
                                         ratio = 1, 
                                         exact = ~sex
) 
Vanderbilt_perifusion_matched <- match.data(Vanderbilt_perifusion_matched)
table(Vanderbilt_perifusion_matched$sex, Vanderbilt_perifusion_matched$simplified_diagnosis)
Vanderbilt_perifusion_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for males

Vanderbilt_perifusion_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 70, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(15,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14J.tiff", width = 4, height = 4)

## Fig S14I
stimulus_data3 <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123),
  xmax = c(12,42,63,72,72,93,102,102,150,132),
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(7,8,7.5,8,9,7.5,8,9,7.5,8),    
  ymax = c(8,8.5,8,9,9.5,8,9,9.5,8,8.5))

all_insulin_df_matched <- all_insulin_df %>%
  filter(DonorID %in% Vanderbilt_perifusion_matched$DonorID)
all_insulin_df_matched %>%
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
    data = stimulus_data3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data3,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS14/FigS14I.tiff", width = 8, height = 4)

## Fig S14K ## ----
all_insulin_wide_matched <- all_insulin_wide %>%
  filter(DonorID %in% Vanderbilt_perifusion_matched$DonorID)

all_insulin_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(baselineins) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(baselineins) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(baselineins), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 0.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 0.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline insulin secretion\n(ng/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14K.tiff", width = 4, height = 4)

## Fig S14L ## ----
all_insulin_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.G16.7) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.G16.7) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0224 for ctrl vs T2D in females, 0.0899 in males

all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.G16.7), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.03))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.6, label = "p=0.1", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14//FigS14L.tiff", width = 4, height = 4)

## Fig S14M ## ----
all_insulin_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #sig for both

all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,6.4)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14M.tiff", width = 4, height = 4)

## Fig S14N ## ----
all_insulin_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis) #sig for disease
model <- all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #p=0.0765 for ctrl vs T2D in females, p=0.0022 in males

all_insulin_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 4, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 4, label = "p=0.08", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  scale_y_continuous(limits = c(0.5,5.1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS14/FigS14N.tiff", width = 4, height = 4)


## Fig S11H ## ----
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS11/FigS11H.tiff", width = 5.5, height = 4)

## Fig S11I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 0.75, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline insulin secretion\n(ng/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS11/FigS11I.tiff", width = 3, height = 3)

## Fig S11J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 4.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n16.7 mM glucose)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS11/FigS11J.tiff", width = 3, height = 3)

## Fig S11K ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 4.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.075, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS11/FigS11K.tiff", width = 3, height = 3)

## Fig S11L ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 3.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS11/FigS11L.tiff", width = 3, height = 3)

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

## Fig S17B ## ----
stimulus_data4 <- data.frame(
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
  geom_path(linewidth = 1)+
  geom_ribbon(alpha = 0.2, colour = NA, (aes(ymax = up, ymin = down, fill = sex)))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab ("Glucagon release (pg/100 IEQs)") +
  xlab ("Time (min)") +
  scale_x_continuous(limits = c(3, 150)) +
  scale_y_continuous(limits = c(-40, 390)) +
  geom_vline(xintercept = c(3,12,42,63,63,72,93,93,102,123), colour = "darkgrey", linetype = 5)+
  geom_rect(
    data = stimulus_data4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",linewidth = 0.5,
  ) +
  geom_text(
    data = stimulus_data4,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS17/FigS17B.tiff", width = 5.5, height = 4)

## Fig S17D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(baseline glucagon\nsecretion (pg/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS17/FigS17D.tiff", width = 3, height = 3)

## Fig S17F ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n1.7 mM glucose + Epi)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS17/FigS17F.tiff", width = 3, height = 3)

## Fig S17H ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 9.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.2))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS17/FigS17H.tiff", width = 3, height = 3)

## Fig S17J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 8.3, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS17/FigS17J.tiff", width = 3, height = 3)

## Fig S18B ## ----
stimulus_data5 <- data.frame(
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
    data = stimulus_data5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data5,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS18/FigS18B.tiff", width = 8, height = 4)

## Fig S18D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion\n(pg/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS18/FigS18D.tiff", width = 4, height = 4)

## Fig S18F ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.12, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 9.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n1.7 mM glucose + Epi)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS18/FigS18F.tiff", width = 4, height = 4)

## Fig S18H ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 9.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS18/FigS18H.tiff", width = 4, height = 4)

## Fig S18J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 8.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 9.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS18/FigS18J.tiff", width = 4, height = 4)

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

## Fig S17L ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = sex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  scale_colour_manual(values = c("#50b5ad", "#dbac5a")) +
  scale_fill_manual(values = c("#50b5ad", "#dbac5a")) +
  ylab("log(glucagon content\n(pg/IEQ))") +
  xlab("") +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank())  +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line()  
  )
ggsave("Output/Final figures/FigS17/FigS17L.tiff", width = 3, height = 3)

## Fig S18L ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.05, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 9.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 9.7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/IEQ))") +
  xlab("")+
  scale_y_continuous(limits = c(3,10)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS18/FigS18L.tiff", width = 5, height = 4)

#Age-matched
## Fig S19D ## ----
Vanderbilt_perifusion_matched_glucagon <- all_glucagon_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years) %>%
  distinct()

Vanderbilt_perifusion_matched_glucagon$Group <- as.logical(Vanderbilt_perifusion_matched_glucagon$simplified_diagnosis == "T2D")
Vanderbilt_perifusion_matched_glucagon <- matchit(Group ~ age_years + sex,
                                                  data = Vanderbilt_perifusion_matched_glucagon,
                                                  method = 'nearest', # can switch to 'exact' or 'optimal'
                                                  ratio = 1, 
                                                  exact = ~sex
) 
Vanderbilt_perifusion_matched_glucagon <- match.data(Vanderbilt_perifusion_matched_glucagon)
table(Vanderbilt_perifusion_matched_glucagon$sex, Vanderbilt_perifusion_matched_glucagon$simplified_diagnosis)
Vanderbilt_perifusion_matched_glucagon %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for males
Vanderbilt_perifusion_matched_glucagon %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(15,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19D.tiff", width = 5, height = 4)

## Fig S19C ## ----
all_glucagon_df_matched <- all_glucagon_df %>%
  filter(DonorID %in% Vanderbilt_perifusion_matched_glucagon$DonorID)

stimulus_data6 <- data.frame(
  xmin = c(3,12,42,63,63,72,93,93,102,123),         
  xmax = c(12,42,63,72,72,93,102,102,150,132), 
  label = c("G\n5.6", "G 16.7", "G 5.6","G\n16.7","IBMX","G 5.6", "G\n1.7","Epi","G 5.6","KCl"),
  ymin = c(280,320,300,320,360,300,320,360,300,320), 
  ymax = c(320,340,320,360,380,320,360,380,320,340))

all_glucagon_df_matched %>%
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
    data = stimulus_data6,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "grey30",size = 0.5,
  ) +
  geom_text(
    data = stimulus_data6,
    aes(x = (xmin + xmax)/2, y = (ymin+ymax)/2, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS19/FigS19C.tiff", width = 8, height = 4)

## Fig S19F ## ----
all_glucagon_wide_matched <- all_glucagon_wide %>%
  filter(DonorID %in% Vanderbilt_perifusion_matched_glucagon$DonorID)

all_glucagon_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(baselinegcg) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(baselinegcg) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(baselinegcg), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 5.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(baseline glucagon secretion\n(pg/100 IEQs))") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19F.tiff", width = 4, height = 4)

## Fig S19J ## ----
all_glucagon_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.IBMX) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.IBMX), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion 16.7 mM\nglucose + 0.1 mM IBMX)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19J.tiff", width = 4, height = 4)

## Fig S19H ## ----
all_glucagon_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.Epi) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.Epi) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.Epi), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 9.1, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.12, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n1.7 mM glucose + Epi)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19H.tiff", width = 4, height = 4)

## Fig S19L ## ----
all_glucagon_wide_matched %>%
  ungroup() %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis) #ns
model <- all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(AUC.KCl) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

all_glucagon_wide_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(AUC.KCl), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 8.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(AUC secretion\n20 mM KCl)") +
  xlab("")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19L.tiff", width = 4, height = 4)

## Fig S19O ## ----
glucagon.content_matched <- glucagon.content_df %>%
  filter(simplified_diagnosis != "T1D") %>%
  filter(Glucagon.content.pg.IEQ > 0) %>%
  dplyr::select(DonorID, sex, simplified_diagnosis, age_years, Glucagon.content.pg.IEQ) %>%
  distinct()

glucagon.content_matched$Group <- as.logical(glucagon.content_matched$simplified_diagnosis == "T2D")
glucagon.content_matched <- matchit(Group ~ age_years + sex,
                                    data = glucagon.content_matched,
                                    method = 'nearest', # can switch to 'exact' or 'optimal'
                                    ratio = 1, 
                                    exact = ~sex
) 
glucagon.content_matched <- match.data(glucagon.content_matched)
table(glucagon.content_matched$sex, glucagon.content_matched$simplified_diagnosis)

glucagon.content_matched %>%
  group_by(sex) %>%
  t_test(age_years ~ simplified_diagnosis) #sig for males
glucagon.content_matched %>%
  ggplot(aes(x=sex, y=age_years))+
  geom_boxplot(aes(fill = simplified_diagnosis), alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 68, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 70, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  xlab("")+
  ylab("Age (years)") +
  ylim(15,75)+
  labs(fill = "Condition")+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS19/FigS19O.tiff", width = 4, height = 4)

## Fig S19P ## ----
glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  anova_test(log(Glucagon.content.pg.IEQ) ~ age_years + sex*simplified_diagnosis) #ns
model <- glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  lm(log(Glucagon.content.pg.IEQ) ~ age_years + sex*simplified_diagnosis, data = .)
emmeans_res <- emmeans(model, ~ simplified_diagnosis | sex)
pairs(emmeans_res, adjust = "tukey") #ns

glucagon.content_matched %>%
  filter(simplified_diagnosis != "T1D") %>%
  ggplot(aes(x=sex, y=log(Glucagon.content.pg.IEQ), fill = simplified_diagnosis)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.9))+
  geom_quasirandom(size = 2, shape = 21, aes(group = simplified_diagnosis, fill = simplified_diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE, xmin = 0.775, xmax = 1.225, y.position = 8.9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE, xmin = 1.775, xmax = 2.225, y.position = 8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey50", "#d65ac9")) +
  scale_fill_manual(values = c("grey50", "#d65ac9")) +
  ylab("log(glucagon content (pg/islet))") +
  xlab("")+
  scale_y_continuous(limits = c(3,9.5))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line()   
  )
ggsave("Output/Final figures/FigS19/FigS19P.tiff", width = 4, height = 4)

### Figures S15-16 ### ----
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

#convert to ng/100 islets
peri_leu_means <- peri_leu_means %>%
  mutate(across(-c(record_id, replicate), ~ .x * (34.848/1625)))

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

## Fig S15A ## ----
stimulus_data <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190),
  label = c("G 3","Leu 5","Leu 5 + G 6","KCl"),
  ymin = rep(16,4),
  ymax = rep(18,4))

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
  ylab ("Insulin Release (ng/100 islets)") +
  xlab ("Time (min)") +
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS15/FigS15A.tiff", width = 5.5, height = 4)

## Fig S15B ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 2.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(mean baseline insulin\nsecretion (ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-3,3)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15B.tiff", width = 3, height = 3)

## Fig S15C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.35, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n5 mM leucine)") +
  xlab("")+
  scale_y_continuous(limits = c(1.5,7)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15C.tiff", width = 3, height = 3)

## Fig S15D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.35, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 5 mM\nleucine + 6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(2.5,7)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15D.tiff", width = 3, height = 3)

## Fig S15E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl) ") +
  xlab("")+
  scale_y_continuous(limits = c(2,6.5)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15E.tiff", width = 3, height = 3)

## Fig S16A ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190), 
  label = c("G 3","Leu 5","Leu 5 + G 6","KCl"),
  ymin = rep(12,4),
  ymax = rep(13,4))

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
  ylab ("Insulin Release (ng/100 islets)") +
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS16/FigS16A.tiff", width = 8, height = 4)

## Fig S16B ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.2, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2.4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 2.7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.17, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 3.6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-4,4)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16B.tiff", width = 5, height = 4)

## Fig S16C ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.05))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.15, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 9, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 5 mM leucine)") +
  xlab("")+
  scale_y_continuous(limits = c(0,10)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16C.tiff", width = 5, height = 4)

## Fig S16D ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.075))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 5 mM leucine\n+ 6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(0,10)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16D.tiff", width = 5, height = 4)

## Fig S16E ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.12))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.07))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 30 mM KCl) ") +
  xlab("")+
  scale_y_continuous(limits = c(-0.5,9)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16E.tiff", width = 5, height = 4)

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

#convert to ng/100 islets
peri_olp_means <- peri_olp_means %>%
  mutate(across(-c(record_id, replicate), ~ .x * (34.848/1625)))


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

## Fig S15F ## ----
stimulus_data <- data.frame(
  xmin = c(0,20,90,160), 
  xmax = c(20,70,140,190), 
  label = c("G 3","OLP 1.5","OLP 1.5 + G 6","KCl"),
  ymin = c(rep(15,4)),     
  ymax = c(rep(17,4)))

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
  ylab ("Insulin Release (ng/100 islets)") +
  xlab ("Time (min)") +
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS15/FigS15F.tiff", width = 5.5, height = 4)

## Fig S15G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 2.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.2, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(mean baseline insulin\nsecretion (ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-3,3)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15G.tiff", width = 3, height = 3)

## Fig S15H ## ----
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
  geom_quasirandom(size = 2,shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.45, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate)") +
  xlab("")+
  scale_y_continuous(limits = c(2,7)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15H.tiff", width = 3, height = 3)

## Fig S15I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6.6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate\n+ 6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(2,7)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15I.tiff", width = 3.25, height = 3)

## Fig S15J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(fill = donorsex)) +
  geom_bracket(inherit.aes = FALSE,xmin = 1, xmax = 2, y.position = 6, label = "*", label.size = 7, size = 0.5, tip.length = c(0.27, 0.02))+
  scale_fill_manual(values = c("#113ed1", "#5e4114")) +
  scale_colour_manual(values = c("#113ed1", "#5e4114")) +
  ylab("log(AUC secretion\n30 mM KCl) ") +
  xlab("")+
  scale_y_continuous(limits = c(2,7)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(), 
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS15/FigS15J.tiff", width = 3, height = 3)

## Fig S16F ## ----
stimulus_data2 <- data.frame(
  xmin = c(0,20,90,160),
  xmax = c(20,70,140,190), 
  label = c("G 3","OLP 1.5","OLP 1.5 + G 6","KCl"),
  ymin = rep(11.5,4),
  ymax = rep(13,4))

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
  ylab ("Insulin Release (ng/100 islets)") +
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
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave("Output/Final figures/FigS16/FigS16F.tiff", width = 8, height = 4)

## Fig S16G ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 1.75, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.15))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 2.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 3, label = "*", label.size = 7, size = 0.5, tip.length = c(0.12, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 4, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(baseline insulin secretion\n(ng/100 islets))") +
  xlab("")+
  scale_y_continuous(limits = c(-2.5,4.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16G.tiff", width = 5, height = 4)

## Fig S16H ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 5.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.1))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 6.25, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.175, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate)") +
  xlab("")+
  scale_y_continuous(limits = c(1.5,8.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16H.tiff", width = 5, height = 4)

## Fig S16I ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.2))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 7, label = "*", label.size = 7, size = 0.5, tip.length = c(0.02, 0.3))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.5, label = "*", label.size = 7, size = 0.5, tip.length = c(0.1, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 0.75 mM\noleate + 0.75 mM palmitate\n+ 6 mM glucose)") +
  xlab("")+
  scale_y_continuous(limits = c(2,9)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16I.tiff", width = 5.25, height = 4)

## Fig S16J ## ----
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
  geom_quasirandom(size = 2, shape = 21, aes(group = diagnosis, fill = diagnosis), dodge.width=0.9) +
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.225, y.position = 6.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.775, xmax = 2.225, y.position = 7, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.25))+
  geom_bracket(inherit.aes = FALSE,xmin = 0.775, xmax = 1.775, y.position = 7.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  geom_bracket(inherit.aes = FALSE,xmin = 1.225, xmax = 2.225, y.position = 8.5, label = "ns", label.size = 4, size = 0.5, tip.length = c(0.02, 0.02))+
  scale_colour_manual(values = c("grey40", "#e04e34"))  +
  scale_fill_manual(values = c("grey40", "#e04e34"))  +
  ylab("log(AUC secretion 30 mM KCl) ") +
  xlab("")+
  scale_y_continuous(limits = c(1.5,9))+
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line() 
  )
ggsave("Output/Final figures/FigS16/FigS16J.tiff", width = 5, height = 4)

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
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
ggsave("Output/Final figures/FigS1/FigS1A.tiff", width = 6, height = 4)

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
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
ggsave("Output/Final figures/FigS1/FigS1B.tiff", width = 6, height = 4)

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
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
ggsave("Output/Final figures/FigS1/FigSC.tiff", width = 6, height = 4)

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
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
ggsave("Output/Final figures/FigS1/FigSD.tiff", width = 6, height = 4)

### Table S1 ### ----
#clear environment
rm(list = ls())

#read in sheets
bulkRNAseq <- read.csv("Output/Final figures/Fig2/BulkHI and pbHPAP_combined_youngcontrols_dea_results_correctforage_dataset.csv")
bulkRNAseq <- bulkRNAseq[,c(8,1:7)]
betascRNAseq <- read.csv("Output/Final figures/FigS3/pbBeta_HPAPonly_youngcontrols_dea_results_correctforage_dataset.csv")
betascRNAseq <- betascRNAseq[,c(8,1:7)]
alphascRNAseq <- read.csv("Output/Final figures/FigS3/pbAlpha_HPAPonly_youngcontrols_dea_results_correctforage_dataset.csv")
alphascRNAseq <- alphascRNAseq[,c(8,1:7)]

#write as new sheet
wb <- createWorkbook()
addWorksheet(wb, "DE results")

writeData(wb, sheet = "DE results", x = bulkRNAseq, startRow = 2, startCol = 1)
writeData(wb, sheet = "DE results", x = betascRNAseq, startRow = 2, startCol = 10)
writeData(wb, sheet = "DE results", x = alphascRNAseq, startRow = 2, startCol = 19)

writeData(wb, sheet = "DE results", x = "Combined bulk RNAseq", startRow = 1, startCol = 1)
writeData(wb, sheet = "DE results", x = "HPAP beta-cell scRNAseq", startRow = 1, startCol = 10)
writeData(wb, sheet = "DE results", x = "HPAP alpha-cell scRNAseq", startRow = 1, startCol = 19)

saveWorkbook(wb, "Output/Final figures/Table S1.xlsx", overwrite = TRUE)

### Table S2 ### ----
#clear environment
rm(list = ls())

#read in sheets and merge
CombinedGO_CC <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_CC")
CombinedGO_BP <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_BP")
CombinedGO_MF <- read_excel("Output/Final figures/Fig2/GSEA_p_BulkHI and pbHPAP_youngcontrols.xlsx",sheet="GO_MF")
CombinedGO <- bind_rows(CombinedGO_CC, CombinedGO_BP, CombinedGO_MF)
CombinedGO$dataset <- rep("Combined bulk RNAseq", nrow(CombinedGO))

HPAPbetascGO_CC <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_CC")
HPAPbetascGO_BP <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_BP")
HPAPbetascGO_MF <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbBeta_youngcontrols.xlsx",sheet="GO_MF")
HPAPbetascGO <- bind_rows(HPAPbetascGO_CC, HPAPbetascGO_BP, HPAPbetascGO_MF)
HPAPbetascGO$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO))

HPAPalphascGO_CC <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_CC")
HPAPalphascGO_BP <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_BP")
HPAPalphascGO_MF <- read_excel("Output/Final figures/FigS3/GSEA_p_HPAPonly_pbAlpha_youngcontrols.xlsx",sheet="GO_MF")
HPAPalphascGO <- bind_rows(HPAPalphascGO_CC, HPAPalphascGO_BP, HPAPalphascGO_MF)
HPAPalphascGO$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO))

all <- bind_rows(CombinedGO, HPAPbetascGO, HPAPalphascGO)
all <- all %>% arrange(p.adjust)
head(all)

write.csv(all, "Output/Final figures/Table S2.csv")

### Table S3 ### ----
#clear environment
rm(list = ls())

#read in sheet
proteomics <- read.csv("Output/Final figures/Fig2/prot_dea_results_correctforage_controls_15to39.csv")

#write as new sheet
write.csv(proteomics, "Output/Final figures/Table S3.csv")

### Table S4 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Final figures/Fig2/GSEA_p_proteomics_youngcontrols.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Final figures/Table S4.csv")

### Table S5 ### ----
#clear environment
rm(list = ls())

#read in sheets
bulkRNAseq <- read.csv("Output/Final figures/Fig3-4/BulkHI and pbHPAP_combined_F_dea_results_correctforage_dataset.csv")
bulkRNAseq <- bulkRNAseq[,c(8,1:7)]
betascRNAseq <- read.csv("Output/Final figures/FigS4/beta HPAP only_F_dea_results_correctforage_dataset.csv")
betascRNAseq <- betascRNAseq[,c(8,1:7)]
alphascRNAseq <- read.csv("Output/Final figures/FigS4/alpha HPAP only_F_dea_results_correctforage_dataset.csv")
alphascRNAseq <- alphascRNAseq[,c(8,1:7)]

#write as new sheet
wb <- createWorkbook()
addWorksheet(wb, "DE results")

writeData(wb, sheet = "DE results", x = bulkRNAseq, startRow = 2, startCol = 1)
writeData(wb, sheet = "DE results", x = betascRNAseq, startRow = 2, startCol = 10)
writeData(wb, sheet = "DE results", x = alphascRNAseq, startRow = 2, startCol = 19)

writeData(wb, sheet = "DE results", x = "Combined bulk RNAseq", startRow = 1, startCol = 1)
writeData(wb, sheet = "DE results", x = "HPAP beta-cell scRNAseq", startRow = 1, startCol = 10)
writeData(wb, sheet = "DE results", x = "HPAP alpha-cell scRNAseq", startRow = 1, startCol = 19)

saveWorkbook(wb, "Output/Final figures/Table S5.xlsx", overwrite = TRUE)

### Table S6 ### ----
#clear environment
rm(list = ls())

#read in sheets
bulkRNAseq <- read.csv("Output/Final figures/Fig3-4/BulkHI and pbHPAP_combined_M_dea_results_correctforage_dataset.csv")
bulkRNAseq <- bulkRNAseq[,c(8,1:7)]
betascRNAseq <- read.csv("Output/Final figures/FigS4/beta HPAP only_M_dea_results_correctforage_dataset.csv")
betascRNAseq <- betascRNAseq[,c(8,1:7)]
alphascRNAseq <- read.csv("Output/Final figures/FigS4/alpha HPAP only_M_dea_results_correctforage_dataset.csv")
alphascRNAseq <- alphascRNAseq[,c(8,1:7)]

#write as new sheet
wb <- createWorkbook()
addWorksheet(wb, "DE results")

writeData(wb, sheet = "DE results", x = bulkRNAseq, startRow = 2, startCol = 1)
writeData(wb, sheet = "DE results", x = betascRNAseq, startRow = 2, startCol = 10)
writeData(wb, sheet = "DE results", x = alphascRNAseq, startRow = 2, startCol = 19)

writeData(wb, sheet = "DE results", x = "Combined bulk RNAseq", startRow = 1, startCol = 1)
writeData(wb, sheet = "DE results", x = "HPAP beta-cell scRNAseq", startRow = 1, startCol = 10)
writeData(wb, sheet = "DE results", x = "HPAP alpha-cell scRNAseq", startRow = 1, startCol = 19)

saveWorkbook(wb, "Output/Final figures/Table S6.xlsx", overwrite = TRUE)


### Table S7 ### ----
#clear environment
rm(list = ls())

#read in sheets and combine
CombinedGO_F_CC <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_CC")
CombinedGO_F_BP <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_BP")
CombinedGO_F_MF <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_F.xlsx",sheet="GO_MF")
CombinedGO_F <- bind_rows(CombinedGO_F_CC, CombinedGO_F_BP, CombinedGO_F_MF)
CombinedGO_F$dataset <- rep("Combined bulk RNAseq", nrow(CombinedGO_F))

HPAPbetascGO_F_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_CC")
HPAPbetascGO_F_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_BP")
HPAPbetascGO_F_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_F.xlsx",sheet="GO_MF")
HPAPbetascGO_F <- bind_rows(HPAPbetascGO_F_CC, HPAPbetascGO_F_BP, HPAPbetascGO_F_MF)
HPAPbetascGO_F$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO_F))

HPAPalphascGO_F_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_CC")
HPAPalphascGO_F_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_BP")
HPAPalphascGO_F_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_F.xlsx",sheet="GO_MF")
HPAPalphascGO_F <- bind_rows(HPAPalphascGO_F_CC, HPAPalphascGO_F_BP, HPAPalphascGO_F_MF)
HPAPalphascGO_F$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO_F))

all <- bind_rows(CombinedGO_F, HPAPbetascGO_F, HPAPalphascGO_F)
all <- all %>% arrange(p.adjust)
head(all)

write.csv(all, "Output/Final figures/Table S7.csv")

### Table S18 ### ----
#clear environment
rm(list = ls())

#read in sheets and combine
CombinedGO_M_CC <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_CC")
CombinedGO_M_BP <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_BP")
CombinedGO_M_MF <- read_excel("Output/Final figures/Fig3-4/GSEA_p_BulkHI and pbHPAP_M.xlsx",sheet="GO_MF")
CombinedGO_M <- bind_rows(CombinedGO_M_CC, CombinedGO_M_BP, CombinedGO_M_MF)
CombinedGO_M$dataset <- rep("Combined bulk RNAseq", nrow(CombinedGO_M))

HPAPbetascGO_M_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_CC")
HPAPbetascGO_M_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_BP")
HPAPbetascGO_M_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_beta_HPAPonly_M.xlsx",sheet="GO_MF")
HPAPbetascGO_M <- bind_rows(HPAPbetascGO_M_CC, HPAPbetascGO_M_BP, HPAPbetascGO_M_MF)
HPAPbetascGO_M$dataset <- rep("HPAP beta-cell scRNAseq", nrow(HPAPbetascGO_M))

HPAPalphascGO_M_CC <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_CC")
HPAPalphascGO_M_BP <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_BP")
HPAPalphascGO_M_MF <- read_excel("Output/Final figures/FigS4/GSEA_p_alpha_HPAPonly_M.xlsx",sheet="GO_MF")
HPAPalphascGO_M <- bind_rows(HPAPalphascGO_M_CC, HPAPalphascGO_M_BP, HPAPalphascGO_M_MF)
HPAPalphascGO_M$dataset <- rep("HPAP alpha-cell scRNAseq", nrow(HPAPalphascGO_M))


all <- bind_rows(CombinedGO_M, HPAPbetascGO_M, HPAPalphascGO_M)
all <- all %>% arrange(p.adjust)
head(all)

write.csv(all, "Output/Final figures/Table S8.csv")

### Table S9 ### ----
#clear environment
rm(list = ls())

#read in sheet
proteomics <- read.csv("Output/Final figures/Fig3-4/prot_dea_results_correctforage_F.csv")

#write as new sheet
write.csv(proteomics, "Output/Final figures/Table S9.csv")

### Table S10 ### ----
#clear environment
rm(list = ls())

#read in sheet
proteomics <- read.csv("Output/Final figures/Fig3-4/prot_dea_results_correctforage_M.csv")

#write as new sheet
write.csv(proteomics, "Output/Final figures/Table S10.csv")

### Table S11 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_F_controlvsT2D.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Final figures/Table S11.csv")

### Table S12 ### ----
#clear environment
rm(list = ls())

#read and write sheets
sheet_CC <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_CC")
sheet_BP <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_BP")
sheet_MF <- read.xlsx("Output/Final figures/Fig3-4/GSEA_p_proteomics_M_controlvsT2D.xlsx",sheet="GO_MF")
sheet <- bind_rows(sheet_CC, sheet_BP, sheet_MF)
sheet <- sheet %>% arrange(p.adjust)

write.csv(sheet, "Output/Final figures/Table S12.csv")
