#################### Figure 2A
setwd("~/Desktop/")
library(vegan)
ab<-read.table("Relative_abundance_table.csv",sep=",",header=TRUE,row.names=1) #supplementary table 5A
ab<-t(ab)
colnames(ab)
ab <- as.data.frame(ab)
samples<-rownames(ab)
str(ab)
metadata<-read.table("Final_GW_metadata.csv", sep = ",",header = TRUE,row.names = 1)
colnames(metadata)
metadata <- as.data.frame(metadata)
metadata$filter <- as.factor(metadata$filter)
str(metadata)
samples<-rownames(metadata)
samples

library(vegan)
ab <- ab [, colSums(ab  != 0) > 0]
ab  <- ab [rowSums(ab  != 0) > 0, ]
ab.log <- decostand(ab,method = "log") #Warning message:non-integer data: divided by smallest positive value
ab.log_matrix<-data.matrix(ab.log)
bcIsoGenie<-vegdist(ab.log,method="bray",na.rm=TRUE)
bcIsoGenie_matrix<-data.matrix(bcIsoGenie)
bc.pcoa<-cmdscale(bcIsoGenie_matrix, k=2, eig=TRUE, add = TRUE)
bc.pcoa$eig
eig2 <- eigenvals(bc.pcoa)
eig2
percentage_variance_explained <- eig2 / sum(eig2)
percentage_variance_explained

#get the x- and y-label PCo1 and 2 %
xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
xlabel= sprintf("%.2f %%", xlabel)
xlabel= paste ("PCo1 (", xlabel, ")")
ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
ylabel= sprintf("%.2f %%", ylabel)
ylabel= paste ("PCo2 (", ylabel, ")")
#all: xlab="PCo1 (14.46%)", ylab="PCo2 (12.39%)"

###year_well
col<- c("2019_H14" = "#ffe119", "2019_H32" = "#fd7401","2019_H41" = "#d32626", "2019_H43" = "#751267","2019_H51" = "#025b40", "2019_H52" = "#011c7c", "2019_H53" = "#4323c4",
        "2022_H14" = "#f2efba", "2022_H32" = "#f3a950","2022_H41" = "#efb1b1", "2022_H43" = "#d31cb0","2022_H51" = "#b1ddd2", "2022_H52" = "#d5e1f7", "2022_H53" = "#d6cfff")
shape <- ifelse(metadata$filter == "01um", 19, ifelse(metadata$filter == "02um", 18, 19))
plot(bc.pcoa$points, col = col[metadata$year_well], pch = shape, xlab = "PCo1 (14.46%)", ylab = "PCo2 (12.39%)", cex = 8, asp=1)
# Add convex hulls for each group
groups <- as.factor(metadata$well)
unique_groups <- unique(groups)

for (group in unique_groups) {
  group_points <- bc.pcoa$points[groups == group, ]
  hull_indices <- chull(group_points)
  polygon(group_points[hull_indices, ], col = adjustcolor(col[group], alpha.f = 0.3), border = col[group])
}

#running adonis
test1<-adonis(ab.log~note,data=metadata,permutations = 999,method = "bray")
test1
# Extract R-squared value
r_squared <- test1$aov.tab$R2[1]
# Extract F-value
f_value <- test1$aov.tab$`F value`[1]
# Extract p-value
p_value <- test1$aov.tab$`Pr(>F)`[1]
# Print the extracted values
print(r_squared)
print(f_value)
print(p_value)

#check variances distances using vegdist
distances.370<-vegdist(ab.log)
test2<-anova(betadisper(distances.370,metadata$note))
test2
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#PERMANOVA statistic well test1: R2 0.3201581, p-value 0.001
#PERMANOVA statistic well test2: p-value 0.1159
#PERMANOVA statistic filter test1: R2 0.03330103, p-value 0.011
#PERMANOVA statistic filter test2: p-value 0.8525
#PERMANOVA statistic year test1: R2 0.09537302, p-value 0.001
#PERMANOVA statistic year test2: p-value 0.04697
#PERMANOVA statistic year_well test1: R2 0.5551084, p-value 0.001
#PERMANOVA statistic year_well test2: p-value 0.7427
#PERMANOVA statistic note-oxic/anoxic test1: R2 0.04344405, p-value 0.001
#PERMANOVA statistic well test2: p-value 0.08918

#anosim
set.seed(123)
anosim.ab.log <- anosim(bcIsoGenie,metadata$note, permutations = 999)
summary(anosim.ab.log)
#ANOSIM statistic well R: R: 0.4426, Significance: 0.001 
#ANOSIM statistic filter R: 0.05512, Significance: 0.022 
#ANOSIM statistic year R: 0.3256, Significance: 0.001  
#ANOSIM statistic well_year R: 0.7332, Significance: 0.001
#ANOSIM statistic note R: 0.08251, Significance: 0.033

#################### Figure 2B
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

##########
setwd("~/Desktop/")

#Macrodiversity-shannonsH
file_name<- "./diversity_table.csv" #supplementary table 6
data<-read.csv(file=file_name,head=T,dec=".",sep=",")

#data summary
summary(data)
View(data)

#check na
is.na(data)

#check normality
## check data
set.seed(1234)
## check normality
## A. Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(data$shannons_H, main = "Density plot of value",xlab = "value")

## B. Q-Q plot: Q-Q plot (or quantile-quantile plot)
ggqqplot(data$shannons_H)

## C. Shapiro-Wilk test of normality
## the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
## In other words, we can assume the normality
shapiro.test(data$shannons_H) 
## Analysed on 8 feb 2024
#shannons_H:W = 0.93216, p-value = 0.001513

#summary statistics
data %>%
  group_by(site) %>%
  get_summary_stats(value, type = "mean_sd")

#statistics: https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
#kruskal test
res.kruskal <- data %>% kruskal_test(shannons_H ~ well)
res.kruskal

#pairwise posthoc: dunn
#A significant Kruskal-Wallis test is generally followed up by Dunn???s test to identify which groups are different
#Compared to the Wilcoxon???s test, the Dunn???s test takes into account the rankings used by the Kruskal-Wallis test. It also does ties adjustments.
pwc<-data%>%dunn_test(shannons_H ~ well, p.adjust.method = "bonferroni")
pwc

#pairwise wilcox
wilcox_test(data, shannons_H ~ well, p.adjust.method = "bonferroni")

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "well")

# Show significance levels
# Hide non-significant tests
s <- ggviolin(data, x = "well", y = "shannons_H", fill = "well") +
  geom_boxplot(width = 0.1, color = "black", fill = "transparent") +  # Add box plot
  stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p.adj.signif") +
  labs(title = "Comparison of Shannons_H between Wells",
       subtitle = get_test_label(res.kruskal, detailed = TRUE),
       caption = get_pwc_label(pwc)) +
  theme_classic() +
  scale_fill_manual(values = c("H14" = "#d7f7ff", "H32" = "#a3cdff",
                               "H41" = "#2dbafc", "H43" = "#2d72fe",
                               "H51" = "#b27afc", "H52" = "#4323c4", 
                               "H53" = "#011c7c"))
s

#################### Figure 2C
## Microdiversity-Avg.pi
## C. Shapiro-Wilk test of normality
shapiro.test(data$pi) 
#data$pi; W = 0.94712, p-value = 0.007637

#kruskal test
res.kruskal <- data %>% kruskal_test(pi ~ well)
res.kruskal

#pairwise posthoc: dunn
#A significant Kruskal-Wallis test is generally followed up by Dunn???s test to identify which groups are different
#Compared to the Wilcoxon???s test, the Dunn???s test takes into account the rankings used by the Kruskal-Wallis test. It also does ties adjustments.
pwc<-data%>%dunn_test(pi ~ well, p.adjust.method = "bonferroni")
pwc

#pairwise wilcox
wilcox_test(data, pi ~ well, p.adjust.method = "bonferroni")

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "well")

# reoder
data$well <- factor(data$well, levels = c("H14", "H32", "H41", "H43", "H51", "H52", "H53"))
pi <- ggviolin(data, x = "well", y = "pi", fill = "well") + 
  geom_boxplot(width = 0.1, color = "black", fill = "transparent") +  # Add box plot
  #stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p.adj.signif") +
  labs(title = "Comparison of pi between Wells",
       subtitle = get_test_label(res.kruskal, detailed = TRUE),
       caption = get_pwc_label(pwc)) +
  theme_classic() +
  scale_fill_manual(values = c("H14" = "#d7f7ff", "H32" = "#a3cdff",
                               "H41" = "#2dbafc", "H43" = "#2d72fe",
                               "H51" = "#b27afc", "H52" = "#4323c4", 
                               "H53" = "#011c7c"))
pi

#################### Figure 2D
# Statistical test: Wilcoxon test for two independent groups
res.wilcox <- data %>%
  wilcox_test(shannons_H ~ note) %>%
  add_y_position()  # note is the description for oxic/anoxic per sample

# Show significance levels
# Hide non-significant tests
s <- ggviolin(data, x = "note", y = "shannons_H", fill = "note") +
  geom_boxplot(width = 0.1, color = "black", fill = "transparent") +  # Add box plot
  stat_pvalue_manual(res.wilcox, hide.ns = TRUE, label = "p") +  # Use "p" as label
  labs(title = "Comparison of Shannons_H between Groups",
       subtitle = get_test_label(res.wilcox, detailed = TRUE)) +
  theme_classic() +
  scale_fill_manual(values = c("oxic" = "#639CBF", "anoxic" = "#A3A646"))
s

#################### Figure 2E
# Statistical test: Wilcoxon test for two independent groups
res.wilcox <- data %>%
  wilcox_test(pi ~ note) %>%
  add_y_position()  # note is the description for oxic/anoxic per sample

pi2 <- ggviolin(data, x = "note", y = "pi", fill = "note") +
  geom_boxplot(width = 0.1, color = "black", fill = "transparent") +  # Add box plot
  stat_pvalue_manual(res.wilcox, hide.ns = TRUE, label = "p") +  # Use "p" as label
  labs(title = "Comparison of Shannons_H between Groups",
       subtitle = get_test_label(res.wilcox, detailed = TRUE)) +
  theme_classic() +
  scale_fill_manual(values = c("oxic" = "#639CBF", "anoxic" = "#A3A646"))
pi2


########################################
#################### Figure 3A&B
# The input data used for this plot is the summed of normalized relative abundance of MAGs (metagenome-assembled genomes) that have virus-host link predictions (supplementary table 5-sheet 4)
# Only taxa contributing >1% relative abundance in any sample were retained as individual groups. Taxa below this threshold were grouped into 'Others' for clarity.
# the virus abundance were colored based on the hosts'

library(ggplot2)
library(cowplot)

setwd("~/Desktop/")

#####host_stacked
file_name<- "./VH_host-perwell.csv" #see Data directory
data<-read.csv(file=file_name,head=T,dec=".",sep=",")

#per_well
category_colors <- c("a_Patescibacteria" =	"#C82127",
                     "b_Nitrospirota" =	"#F27421",
                     "c_Nanoarchaeota" =	"#2D2D2D",
                     "d_Planctomycetota" =	"#7556A3",
                     "e_Proteobacteria" =	"#228970",
                     "f_Omnitrophota" =	"#496EB4",
                     "g_Thermoproteota" =	"#424242",
                     "h_Others"	="grey")

#reordered
# Define the order of levels for the 'sample' variable
well_order_filter2 <- c("H14", "H32", "H41", "H43", "H51", "H52", "H53")

# Use the factor function to reorder 'sample' based on the defined order
#data.h$sample <- factor(data.h$sample, levels = sample_order)
data$well <- factor(data$well, levels = well_order_filter2)

#stacked bar plot_host
h<- ggplot(data, aes(fill=host_taxa, y=value, x=well)) + 
  geom_bar(position="fill", stat="identity",width = 1)+
  scale_fill_manual(values = category_colors)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Host Abundance")

#####virus_stacked
file_name<- "./VH_virus-perwell.csv" #see Data directory
data<-read.csv(file=file_name,head=T,dec=".",sep=",")

#per_well
category_colors <- c("a_Patescibacteria" =	"#C82127",
                     "b_Nitrospirota" =	"#F27421",
                     "c_Nanoarchaeota" =	"#2D2D2D",
                     "d_Planctomycetota" =	"#7556A3",
                     "e_Proteobacteria" =	"#228970",
                     "f_Omnitrophota" =	"#496EB4",
                     "g_Thermoproteota" =	"#424242",
                     "h_Others"	="grey")

#reordered
# Define the order of levels for the 'sample' variable
well_order_filter2 <- c("H14", "H32", "H41", "H43", "H51", "H52", "H53")

# Use the factor function to reorder 'sample' based on the defined order
#data.h$sample <- factor(data.h$sample, levels = sample_order)
data$well <- factor(data$well, levels = well_order_filter2)

#stacked bar plot_host
v<- ggplot(data, aes(fill=host_link, y=value, x=well)) + 
  geom_bar(position="fill", stat="identity",width = 1)+
  scale_fill_manual(values = category_colors)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Virus Abundance (colored by host phylum)")

combined_plot <- plot_grid(h, v, 
                           labels = c("a", "b"),  # Optional: adds panel labels like in your figure
                           label_size = 14,       # Adjust size for publication
                           ncol = 1,              # Stack vertically
                           align = 'v')           # Align vertically
# If you want to see it in R
print(combined_plot)

#################### Figure 3C
# Only viruses and host with the linkage used in this analysis
# The input data used for this plot is the summed of normalized relative abundance of MAGs (metagenome-assembled genomes) that have virus-host link predictions (supplementary table 5-sheet 4)
# And summed of normalized virus relative abundance >= 10 K
# Line plor for VHR
library(ggplot2)

setwd("~/Desktop/")

#vh_stacked
file_name<- "./VHR.csv"
data<-read.csv(file=file_name,head=T,dec=".",sep=",")

# Define the desired order for the y-axis
ordered_levels <- c(
  "1_Actinobacteriota","2_Cyanobacteria","3_Bacteroidota","4_Proteobacteria","5_Myxococcota","6_Acidobacteriota","7_Unclassified_bacteria","8_Verrucomicrobiota","9_Nitrospinota","10_Methylomirabilota","11_UBA9089","12_Thermoplasmatota","13_Chloroflexota","14_Nanoarchaeota","15_CG2-30-53-67","16_Campylobacterota","17_Spirochaetota","18_Lindowbacteria","19_Nitrospirota","20_Desulfobacterota","21_Elusimicrobiota","22_Bdellovibrionota","23_KSB1","24_Thermoproteota","25_Eisenbacteria","26_Gemmatimonadota","27_Fibrobacterota","28_Micrarchaeota","29_Iainarchaeota","30_Planctomycetota","31_Patescibacteria","32_Omnitrophota","33_Firestonebacteria","34_Margulisbacteria","35_Aenigmatarchaeota"
)

# Keep the order
data$vh_links2 <- factor(data$vh_links2, levels = ordered_levels)

# Basic line plot with points
R <- ggplot(data = data, aes(x = vh_ratio_log10, y = vh_links2, group = 1)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 1, color = "red", linetype = "dotted", size = 1) +  # Red dotted line at x = 1
  theme_classic()
R

# Rank abundant for host
library(ggplot2)
library(reshape2)

setwd("~/Desktop/")

file_name<- "./VHR_host_rank_abundance.csv"
data<-read.csv(file=file_name,head=T,dec=".",sep=",")

#data summary
summary(data)

# Sort taxa by their abundance in each sample
sorted_taxa_data <- data[order(-data$abundance), ]

# Calculate the relative abundance by dividing the abundance by the total abundance in each sample
data$Relative_Abundance <- data$abundance / ave(data$abundance, data$sample, FUN = sum)

#for vh_rank abundance
# Define the desired order for the y-axis
ordered_levels <- c(
  "1_Actinobacteriota","2_Cyanobacteria","3_Bacteroidota","4_Proteobacteria","5_Myxococcota","6_Acidobacteriota","7_Unclassified_bacteria","8_Verrucomicrobiota","9_Nitrospinota","10_Methylomirabilota","11_UBA9089","12_Thermoplasmatota","13_Chloroflexota","14_Nanoarchaeota","15_CG2-30-53-67","16_Campylobacterota","17_Spirochaetota","18_Lindowbacteria","19_Nitrospirota","20_Desulfobacterota","21_Elusimicrobiota","22_Bdellovibrionota","23_KSB1","24_Thermoproteota","25_Eisenbacteria","26_Gemmatimonadota","27_Fibrobacterota","28_Micrarchaeota","29_Iainarchaeota","30_Planctomycetota","31_Patescibacteria","32_Omnitrophota","33_Firestonebacteria","34_Margulisbacteria","35_Aenigmatarchaeota")
data$taxa <- factor(data$taxa, levels = ordered_levels)

# Create the rank abundance plot
HR <- rank_abundance_plot <- ggplot(data, aes(x=taxa, y=Relative_Abundance)) +
  geom_bar(stat="identity", fill="grey") +
  coord_flip() +
  labs(x="Taxonomy", y="Relative Abundance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
HR

# Merge for plotting vertically
combined_plot <- plot_grid(R, HR, 
                           ncol = 1,      # or ncol = 2 for side-by-side
                           align = 'hv', 
                           labels = c('A', 'B'))  # Optional: Adds panel labels

# Save the combined figure
ggsave("VHR_combined.png", combined_plot, width = 10, height = 8, dpi = 600)
print(p)

#################### Figure 3D
# The numbers were calculated based on MAGs with virus links and virus with association with MAGs (iPhop results)
library(ggplot2)
library(cowplot)

# Data for the first pie (Host Associations)
data1 <- data.frame(
  Category = c("With associations", "Without"),
  Percentage = c(36.63, 63.37)
)

# Data for the second pie (Virus Links)
data2 <- data.frame(
  Category = c("With associations", "Without"),
  Percentage = c(74.73, 25.27)
)

# Define fill colors
fill_colors <- c("With associations" = "grey80", "Without" = "black")

# First pie chart (Host Associations)
p1 <- ggplot(data1, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = fill_colors) +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            position = position_stack(vjust = 0.5), size = 4) +
  theme_void() +
  labs(title = "Host Associations") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Second pie chart (Virus Links)
p2 <- ggplot(data2, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = fill_colors) +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            position = position_stack(vjust = 0.5), size = 4) +
  theme_void() +
  labs(title = "Virus Links") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Combine plots without individual legends
combined_pie <- plot_grid(
  p1,  # Keep its legend
  p2,  # Keep its legend
  ncol = 1,
  rel_heights = c(1, 1)  # Adjust height ratio if needed
)

combined_pie

#################### Figure 3E
# The co-occurence analysis was done per year
# see: https://github.com/tpq/propr
# devtools::install_github("tpq/propr")
library(propr)

setwd("~/Desktop/")
data<-read.table("Final_adj-table_all_MAGs.csv",sep=",",header=TRUE,row.names=1) # the normalized ALL MAGs relative abundance table per year in "Data" directory
data<-t(data)

pr <- propr(
  data,
  metric = "rho",  # or "phi", "phs", "cor", "vlr"
  ivar = "clr",  # or can use "iqlr" instead
  alpha = NA,  # use to handle zeros
  p = 95  # used for updateCutoffs
) 

# Save the results to a CSV file
getResults(pr)
results <- getResults(pr)
write.csv(results, "CPR_pr_results_*.csv", row.names = FALSE)

# The result was then filtered using criteria described in materials and methods. And the final co-occurrence network was done using cytoscape
# The input for cytoscape can be seen in "Data" directory: CPR_DPANN_final19-consv.ntw and CPR_DPANN_final22-consv.ntw


#################### Figure 4A
# library
library(tidyverse)

setwd("~/Desktop/")

# Load data
file_name <- "./AMG_filtered_manual.csv" # see in "Data" directory
data <- read.csv(file = file_name, head = TRUE, sep = ",")

# Add color group for each functional assignment
data$color <- case_when(
  data$subheader2 %in% c("1_galactose_degradation", "2_pyruvate_metabolism", "3_pentose_pathway",
                         "4_glycolysis", "5_mannose_degradation", "6_fucose_degradation",
                         "7_hydrocarbon_degradation", "8_TCA", "9_Polyphenolics_Cleavage") ~ "#F291A3",
  data$subheader2 %in% c("10_Sulfur", "11_C1", "12_Electron_transport_Chain",
                         "13_Oxygen", "14_Metal_Reduction", "15_Nitrogen") ~ "#C4E1E5",
  data$subheader2 %in% c("17_Nucleotide", "18_Antibiotic_Resistance", "19_aerobic_corrin_ring_synthesis",
                         "20_MISC", "21_ADO-CBL_synthesis", "22_Flagella_Structure") ~ "#F79839",
  data$subheader2 %in% c("23_Peptidase", "24_Amino_Acid") ~ "#357368",
  data$subheader2 == "25_Transporters" ~ "#B3A4CF",
  TRUE ~ "#D3D3D3"  # default grey if anything else
)

# Original scaling step
data <- mutate(data, count.10 = count / 10)
data <- mutate(data, count = -count.10)

# Reorder factor levels
desired_order <- c("1_galactose_degradation","2_pyruvate_metabolism","3_pentose_pathway",
                   "4_glycolysis","5_mannose_degradation","6_fucose_degradation",
                   "7_hydrocarbon_degradation","8_TCA","9_Polyphenolics_Cleavage",
                   "10_Sulfur","11_C1","12_Electron_transport_Chain","13_Oxygen",
                   "14_Metal_Reduction","15_Nitrogen","17_Nucleotide",
                   "18_Antibiotic_Resistance","19_aerobic_corrin_ring_synthesis",
                   "20_MISC","21_ADO-CBL_synthesis","22_Flagella_Structure",
                   "23_Peptidase","24_Amino_Acid","25_Transporters")

data$subheader2 <- factor(data$subheader2, levels = desired_order)

# FINAL PLOT: use assigned colors
ggplot(data, aes(x = subheader2)) +
  geom_col(aes(y = count.10, fill = color), color = "black") +   # Outer colored
  geom_col(aes(y = count, fill = color), color = "black") +      # Inner colored
  scale_fill_identity() +  # Use the hex color directly
  coord_polar() +
  theme_void()+
  geom_text(aes(y = 100, label = subheader2))

#################### Figure 4B
# The input data used for this plot is the summed of normalized relative abundance of vOTUs (>=10kb) (supplementary table 5-sheet 2)
library(ggplot2)
library(reshape2)

setwd("~/Desktop/")

# Load your data from a CSV file
file_name <- "./AMG_bubble_perwell.csv" # see "Data" directory
data <- read.csv(file = file_name, head = TRUE, dec = ".", sep = ",")

# Define the custom order for x-axis categories
custom_order <- c(
  "1_galactose_degradation","2_pyruvate_metabolism","3_pentose_pathway","4_glycolysis",
  "5_mannose_degradation","6_fucose_degradation","7_hydrocarbon_degradation","8_TCA",
  "9_Polyphenolics_Cleavage","10_Sulfur","11_C1","12_Electron_transport_Chain","13_Oxygen",
  "14_Metal_Reduction","15_Nitrogen","17_Nucleotide","18_Antibiotic_Resistance",
  "19_aerobic_corrin_ring_synthesis","20_MISC","21_ADO-CBL_synthesis","22_Flagella_Structure",
  "23_Peptidase","24_Amino_Acid","25_Transporters"
)

# Define the custom order for y-axis categories
y_axis_order <- c("H14", "H32", "H41", "H43", "H51", "H52", "H53")

# Reorder the 'well' and 'subheader2' columns
data$well <- factor(data$well, levels = rev(y_axis_order))
data$subheader2 <- factor(data$subheader2, levels = custom_order)

# Map colors based on subheader2 functional groups
data$color <- case_when(
  data$subheader2 %in% c("1_galactose_degradation","2_pyruvate_metabolism","3_pentose_pathway",
                         "4_glycolysis","5_mannose_degradation","6_fucose_degradation",
                         "7_hydrocarbon_degradation","8_TCA","9_Polyphenolics_Cleavage") ~ "#F291A3",
  data$subheader2 %in% c("10_Sulfur","11_C1","12_Electron_transport_Chain",
                         "13_Oxygen","14_Metal_Reduction","15_Nitrogen") ~ "#C4E1E5",
  data$subheader2 %in% c("17_Nucleotide","18_Antibiotic_Resistance",
                         "19_aerobic_corrin_ring_synthesis","20_MISC",
                         "21_ADO-CBL_synthesis","22_Flagella_Structure") ~ "#F79839",
  data$subheader2 %in% c("23_Peptidase","24_Amino_Acid") ~ "#357368",
  data$subheader2 == "25_Transporters" ~ "#B3A4CF",
  TRUE ~ "#D3D3D3"  # Default gray if anything else
)

# Plot with color mapping
p <- ggplot(data, aes(x = subheader2, y = well, size = ra, color = color)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(5, 15)) +
  scale_color_identity() +  # Tells ggplot to use the hex codes directly
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Show the plot
print(p)

#################### Figure 4C
# Heatmap was visualized using https://github.com/mpg-age-bioinformatics/flaski, with long10 scale. The input file (AMG_cat_metaT_final.csv) can be seen in "Data" directory

#################### Figure 4D
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/")
# Read the data from the CSV file (make sure to provide the correct path to your file)
file_path <- "kegg_recon.csv" # can be found in "Data" directory
data <- read.csv(file_path)

# Define the desired order for the y-axis
ordered_levels <- c(
  "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289"
)

#reorder column
data$order <- factor(data$order, levels = ordered_levels)

# Define a color palette for the heatmap
#no_kegg_ID
color_palette <- scale_fill_gradient2(low = "#FEF4CF", high = "#313695", mid = "#D27F92", 
                                      midpoint = 16.5, 
                                      limits = c(1, 33), 
                                      name = "no_kegg_ID")


# Plotting the heatmap with adjusted legend and separated by 'metabolism'
ggplot(data, aes(x=order, y=type, fill= no_kegg_ID)) + 
  geom_tile() + 
  facet_grid(. ~ metabolism2, scales = "free_x", space = "free_x") + # Separate by metabolism
  color_palette +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # Rotate x-axis text
    legend.position = "bottom", # Move legend to bottom
    strip.background = element_blank(), # Remove facet label background
    strip.text.x = element_text(angle = 0) # Keep facet labels horizontal
  ) +
  guides(fill=guide_legend(nrow=1)) # Ensure legend is a single row


#module_comp_code
# Define a discrete color palette
color_palette <- scale_fill_manual(values = c("1" = "#FCE665", "2" = "#7FB4AF", "3" = "#3A88B9", "4" = "#3D5B81"),
                                   name = "module_comp_code")

data$module_comp_code <- as.character(data$module_comp_code)

# Plotting the heatmap with the new color palette and 'module_comp_code' as the fill
ggplot(data, aes(x=order, y=type, fill= module_comp_code)) + 
  geom_tile() + 
  facet_grid(. ~ metabolism2, scales = "free_x", space = "free_x") + # Separate by metabolism
  color_palette +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # Rotate x-axis text
    legend.position = "bottom", # Move legend to bottom
    strip.background = element_blank(), # Remove facet label background
    strip.text.x = element_text(angle = 0) # Keep facet labels horizontal
  ) +
  guides(fill=guide_legend(nrow=1)) # Ensure legend is a single row

#################### Figure 4E
library(ggplot2)
library(dplyr)

# Data
kegg_data <- data.frame(
  Category = c("Non-targeted by virus", "Targeted by virus"),
  Count = c(198, 91)
)

# Calculate percent and label
kegg_data <- kegg_data %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0(round(Percentage, 1), "%\n(n=", Count, ")"),
         Category = factor(Category, levels = c("Non-targeted by virus", "Targeted by virus")))

# Plot as horizontal stacked bar with scale
ggplot(kegg_data, aes(y = Percentage, x = 1, fill = Category)) +
  geom_bar(stat = "identity", width = 0.4, color = "black") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("Non-targeted by virus" = "grey80",
                               "Targeted by virus" = "#2E5CB8")) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(title = "KEGG pathway modules:"))

