#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)

# CONSTANTS 
MUTATIONAL_SIGNATURES = c("sToL1")
TAXANOMIC_RANKS = c("Fungi", "Metazoa", "Viridiplantae")

# setwd
setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-4")

# Load data
txt = readLines("tree.nwk")
metadata <- read.csv("dtol_all_samples.taxonomic_classification.csv")
somatic_signatures <- read.csv("somatic_mutational_signature_attributions.x0_excluded.rtol_filtered.csv", check.names = F)

# Manipulate data
colnames(somatic_signatures) = paste0('sToL',colnames(somatic_signatures)) # Add column names
rs = rowSums(somatic_signatures[ , -1], na.rm = TRUE) # calculate row sum
somatic_signatures[ , -1] <- somatic_signatures[ , -1] / rs # normalise mutational signature attribution
somatic_signatures[,1] <- ifelse(grepl("\\.", somatic_signatures[,1]), # Change sample names
                                 sub(".*\\.", "", somatic_signatures[,1]),
                                 somatic_signatures[,1])
somatic_signatures$label=NA # Add label consistent with that from the newick tree
somatic_signatures[somatic_signatures<0.035] <- 0 # # Set values < 0.035 to 0
somatic_signatures$label=sapply(somatic_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  } 
  label
})
somatic_signatures$label[somatic_signatures$label=="Hemaris fuciformis"]="Hemaris fuciformis (iHemFuc2)"

# Manipulate newick tree
txt <- gsub('_', '^', txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr   <- read.tree(text = txt2)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)

# Figure 3
pal <- c(
  Fungi = "#8C564B",
  Metazoa = "#001F54",
  Viridiplantae = "#004300"
)

# Get subset of somatic mutational signature attributions
hm <- somatic_signatures %>%
  select(label, sToL1) %>%
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% tr$tip.label) %>%   
  column_to_rownames("label")

# max values
# max(hm[,"sToL1"]) # 0.845

# build annotation
# anno <- tibble(label = org_labels) %>%
anno <- tibble(label = tr$tip.label) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})
anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
    TRUE ~ NA_character_   # << no group = no colour
  )) 

# build data frame for shading phylogenetic tree
grp_list <- split(anno$label, anno$ColorGroup)
hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

# plot
p <- ggtree(tr, layout = "circular", size = 0.25, color = "grey70") +
  geom_hilight(
    data = hilight_df,
    aes(node = node, fill = ColorGroup),
    alpha = 0.6
  ) + 
  scale_fill_manual(
    values = pal, 
    name = "Taxonomic rank", 
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    text = element_text(family = "Helvetica", size = 8)
  )

# return
ggsave("Extended_Data_Figure_8_alpha.pdf", plot = p, width = 7.48, height = 7.48, units = "in")

# define gradient of colours
pal_1 <- c("#FFFFFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")

# layout parameters
tile_offset <- 0.02
tile_width  <- 0.02
gap         <- 0.2

# define grid colour
grid_col  <- "grey70" 

## add sToL1 mutational signature
p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1,
               hm[,"sToL1", drop = F],
               offset = tile_offset,
               width  = tile_width,
               colnames = F,
               colnames_angle = 90,
               colnames_offset_y = 0.5,
               font.size = 8,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_1, 
    name = "sToL1", 
    na.value = "white",
    limits  = c(0, 0.9),
    breaks  = c(0, 0.2, 0.4, 0.6, 0.8),
    labels  = c("0", "0.2", "0.4", "0.6", "0.8"),
    guide = guide_colorbar(
      direction = "horizontal", 
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 1)
  )

# add legend
p1 <- p1 + 
  theme(
    legend.position = "right", 
    legend.box = "vertical", 
    legend.text = element_text(size = 6)
  )

# return plot
ggsave("Extended_Data_Figure_8_beta.pdf", plot = p1, width = 7.48, height = 7.48, units = "in")



