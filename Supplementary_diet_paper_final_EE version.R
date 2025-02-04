--
#Title:Benthic feeding and diet partitioning in Red Sea mesopelagic fish resolved through DNA metabarcoding and ROV footage
#author: Kah Kheng LIM
#Operating System: MacOS

##install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"))

# Installing from GitHub requires the remotes package
install.packages("remotes")
# Windows users will also need to have RTools installed! http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/

# To install the latest version:
remotes::install_github("david-barnett/microViz")

##load packages##
library(microViz)
library(phyloseq)
library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggvenn)
library(patchwork)
library(bipartite)

# Define helper function to remove OTUs by name
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

################################################
#### LOAD THE DATA & PREP A PHYLOSEQ OBJECT ####
################################################
##set the working directory to where your OTU table file is located
otu_mat <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "OTU matrix")
tax_mat <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "Samples")
plot_outdir <- "/Users/limk/Desktop/publication"

# Convert row names
otu_matrix <- otu_mat %>% tibble::column_to_rownames("OTUID")
# Remove the column named "Identity"
tax_mat <- tax_mat %>%
  select(-Identity)
tax_matrix <- tax_mat %>% tibble::column_to_rownames("OTUID")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")

# Convert to matrix
otu_matrix <- as.matrix(otu_matrix)
tax_matrix <- as.matrix(tax_matrix)

# Create phyloseq object
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(tax_matrix)
samples <- sample_data(samples_df)

mesogut <- phyloseq(OTU, TAX, samples)

## Remove host reads
otus_to_remove = c("Otu1", "Otu2", "Otu4") #family Myctophidae, Phosichthyidae
mesogut <- pop_taxa(mesogut, otus_to_remove)
mesogut_nmds<-mesogut #save this original dataframe for betadiversity analysis
mesogut_venn<-mesogut # save this dataframe for Venn diagram

## Validate and fix taxonomy
mesogut <- tax_fix(
  mesogut,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
)

################################################
############### FIGURE 2 STARTS ################
################################################

## Plot the barplot of prey composition relative abundance by host species
lanternfish <- ps_filter(mesogut, SpeciesCode == "BP")
lightfish<-ps_filter(mesogut, SpeciesCode == "VM")

combined <- phyloseq::merge_phyloseq(
  lanternfish %>% tax_agg("Phylum") %>% ps_get(),
  lightfish %>% tax_agg("Phylum") %>% ps_get()
)
combined

prey_barplot<-combined %>%
  tax_transform("compositional") %>%   # Apply compositional transformation
  ps_arrange(desc(Arthropoda), .target = "otu_table") %>%   # Arrange taxa by abundance
  comp_barplot("Phylum", facet_by="SpeciesCode", n_taxa = 10, sample_order = "asis")

prey_barplot <- prey_barplot + theme(
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.ticks.x = element_blank()   # Remove x-axis ticks
)

ggsave(file.path(plot_outdir, "prey rel.abundance_barchart_non_flip.pdf"), plot = prey_barplot, width= 7.04, height= 3.5)

################################################
############### FIGURE 2 ENDS ##################
################################################

################################################
############### FIGURE 3 STARTS ################
################################################

##rarefy abundances to even depth
set.seed=123
mesogut = rarefy_even_depth(mesogut)

##summary of all plots
mesogut@sam_data$Region <- factor(mesogut@sam_data$Province, levels = c("Northen", "Southern"))
a_my_comparisons <- list( c("BP", "VM"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", " "))

figur<-plot_richness(mesogut, x="SpeciesCode", measures= c("Observed", "Shannon","Simpson"),color = "Province")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)
ggsave(file.path(plot_outdir, "alpha diversity plot with p-value_rarefied_province.png"), plot = figur)

################################################
############### FIGURE 3 ENDS ##################
################################################

################################################
############### FIGURE 4 STARTS ################
################################################

## Make an ordination plot to the distribution of OTUs in 2D space
nmds<-mesogut_nmds %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "robust.aitchison") %>%
  ord_calc(
    method = "NMDS"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "SpeciesCode", fill = "SpeciesCode",
    shape = "Province", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = SpeciesCode)
  )
nmds
ggsave(file.path(plot_outdir, "nmds_plot.pdf"), plot = nmds)

## Calculate the stress value of NMDS
nmds_str<-mesogut_nmds %>%
  tax_transform(rank = "unique", trans = "identity") %>% ##check the disssimilarity of OTUs between two hosts
  dist_calc(dist = "robust.aitchison") %>% ## rclr transformation is better at handling zero data
  ord_calc(
    method = "NMDS")
nmds_ord<-ord_get(nmds_str)
nmds_ord$stress

# straight to the betadisp
bd1 <- mesogut_nmds %>%
  tax_agg("unique") %>%
  dist_calc("robust.aitchison") %>%
  dist_bdisp(variables = c("SpeciesCode", "Province")) %>%
  bdisp_get()

# checking the results for each variable
bd1$SpeciesCode #host species
bd1$Province #Province in the Red Sea

################################################
############### FIGURE 4 ENDS ##################
################################################

################################################
############### FIGURE 5 STARTS ################
################################################

## Read the data again from the excel file
otu_mat <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "OTU matrix")
tax_mat <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("/Users/limk/Desktop/Supplementary Table 2.xlsx", sheet = "Samples")

##Bipartite network analysis start
guts_tax <- inner_join(otu_mat,tax_mat, by = "OTUID")
guts_nohost <- guts_tax[!(guts_tax$Phylum %in% c("Apicomplexa","Bigyra", "Cercozoa", "Bacillariophyta", "Bolidophyceae", "Chlorophyta", "Chrysophyceae","Choanoflagellata", "Dinophyceae", "Discosea", "Eustigmatophyceae", "Haptophyta", "Malawimonadida", "Pelagophyceae", "Phaeophyceae", "Prasinodermophyta", "Raphidophyceae", "Rhodophyta", "Streptophyta", "Synurophyceae","Tubulinea", "Xanthophyceae", "Ascomycota", "Basidiomycota", "Oomycota")) & !(guts_tax$Family %in% c("Myctophidae", "Phosichthyidae")), ]
guts_isol <- as_tibble(guts_nohost[c(1:61)]) #isolates column that contains OTUID and the samples (BP_ST1 till VM_ST60)
###transpose otu-by-species table
guts_tr <- guts_isol %>% 
  gather(var, value, -OTUID) %>% 
  spread(OTUID, value)

###calculate relative abundance of OTUS
guts_rel.abu <- guts_tr %>% 
  gather(variable, value, -var) %>% 
  group_by(var) %>% 
  mutate(percentage = value/sum(value)) %>% 
  select(-value) %>% 
  spread(variable, percentage)

# Function to rename samples according to species
rename_samples <- function(sample_name) {
  if (str_detect(sample_name, "^LF-ST([1-9]|[1-2][0-9])$")) {
    return(str_replace(sample_name, "LF-ST", "BP-ST"))
  } else if (str_detect(sample_name, "^LF-ST([3-5][0-9]|60)$")) {
    return(str_replace(sample_name, "LF-ST", "VM-ST"))
  } else {
    return(sample_name)
  }
}

guts_rel.abu <- guts_rel.abu %>%
  mutate(new_var = rename_samples(var)) %>%
  mutate(old_var = var)  # Optionally, save the original variable name in another column

###create metadata vector for species
species_metadata <- as_tibble(sapply(strsplit(guts_rel.abu$new_var,"-"), `[`, 1))

#split species character vector
species_metadata$sp <- substr(species_metadata$value, 1,2)
colnames(species_metadata) <- c("uniqueID", "species","ID")
species_metadata$ID <- guts_rel.abu$old_var

###bipartite network
colnames(guts_tr)[1] <- "ID"
species_metadata
gut.bip <- inner_join(species_metadata, guts_tr)
gut.sum <- gut.bip[,-c(1,3)] %>% group_by(species) %>% summarize_all(sum)

###get relative abundance
gut.bip_rel <- gut.sum %>% 
  gather(variable, value, -species) %>% 
  group_by(species) %>% 
  mutate(percentage = value/sum(value)) %>% 
  select(-value) %>% 
  spread(variable, percentage)

gut.bip_rel <- as.data.frame(gut.bip_rel)
rownames(gut.bip_rel) <- gut.bip_rel$species
bip_motus <- as.data.frame(t(gut.bip_rel[-1]))
bip_motus$OTUID <- rownames(bip_motus)
bip_combo <- inner_join(bip_motus, guts_nohost[,c(1:69)])

########################################
######### FIGURE 5 #####################
########################################

bipartite.n <- bip_combo %>%
  mutate(bipartite.n = case_when(
    Identity >= 0.95 ~ paste(Phylum, Class, Order, sep = "_"),
    Identity >= 0.75 ~ paste(Phylum,Class, sep = "_"),
    Identity < 0.75 ~ paste(Phylum),
  ))
bipartite.n$bipartite.n <- gsub("_NA","",as.character(factor(bipartite.n$bipartite.n)), fixed = TRUE) 
bipartite.n.sub <- bipartite.n[,c(1:2,72)] #new dataframe consisted of columns of the species of interest (1:2) and bipartite column

MOTUs.sum <- as.data.frame(bipartite.n.sub %>% group_by (bipartite.n) %>% summarize_all (sum)) #the dataframe is grouped by the bipartite column, and then all other columns are summarized (summed)
MOTUs.sum.sub <- subset(MOTUs.sum, rowSums(MOTUs.sum[-1]) >= 0.005) #filters MOTUs.sum to include only rows where the sum of all but the first column is greater than or equal to 0.005
split_motus <-  as.data.frame(MOTUs.sum.sub %>% separate(col = bipartite.n, into = c("Cat1","Cat2","Cat3"), sep = '[_]', remove = FALSE))#separates the bipartite.n column into four new columns (Cat1, Cat2, Cat3, Cat4) using underscores as separators. The original column is kept (remove = FALSE).
phyla.new <- as.data.frame(unique(split_motus$Cat1))
colnames(phyla.new) <- "Cat1" #Unique values from Cat1 are extracted and assigned colors using a color ramp palette.

# color function for bipartite, 19 phyla
colfunc4 <- colorRampPalette(c("#4169E1","#008080","gold", "firebrick1"))
phyla.new$colors <- colfunc4(7)
color_vector <- as.data.frame(merge(phyla.new, split_motus, by = "Cat1"))
box_colors <- as.character(color_vector$colors)
color_vec1 <- as.character(rep(color_vector$colors, each = 2)) #Colors are merged with the split_motus data frame, and color vectors are prepared for plotting

# Start the PDF device driver to save the plot
pdf(file = file.path(plot_outdir, "bipartite_network.pdf"), width = 8, height = 5)

# Ensure row names are properly set
rownames(MOTUs.sum.sub) <- as.character(MOTUs.sum.sub$bipartite.n)

# Create the plot
plotweb(MOTUs.sum.sub[-1], method="normal", arrow="down", text.rot=90, 
        col.interaction = color_vec1, col.low = box_colors, bor.col.interaction = FALSE,
        y.width.low=0.03, y.width.high=0.02, labsize = 0.5)

# Close the PDF device
dev.off()

################################################
############### FIGURE 5 ENDS ##################
################################################

################################################
############ SUPPLEMENTARY ANALYSES ############
################################################

################################################
############## Sup.Fig 1 STARTS ################
################################################

# Prepare data for Venn diagram
mesogut_meta_df <- psmelt(mesogut_venn)
str(mesogut_meta_df)

# Function to create Venn data for a given taxonomic level
create_venn_data <- function(df, taxonomic_level) {
  list(
    Bp = unique(df %>% filter(SpeciesCode == "BP", Abundance > 0) %>% pull(!!sym(taxonomic_level))),
    Vm = unique(df %>% filter(SpeciesCode == "VM", Abundance > 0) %>% pull(!!sym(taxonomic_level)))
  )
}

# List of taxonomic levels
taxonomic_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Genspe")

# Create and plot Venn diagrams for each taxonomic level
venn_plots <- lapply(taxonomic_levels, function(level) {
  venn_data <- create_venn_data(mesogut_meta_df, level)
  ggvenn(venn_data, fill_color = c("#FF69B4", "#3B7DB5"), text_size = 3) + 
    ggtitle(paste("Venn Diagram at", level, "Level"))+
    theme(plot.title = element_text(size = 10))
})

# Display all Venn diagrams
venn_plots

combined_plot <- wrap_plots(venn_plots, nrow = 2)
ggsave(file.path(plot_outdir, "Venn diagram.pdf"), plot = combined_plot)

################################################
############### Sup.Fig 1 ENDS #################
################################################