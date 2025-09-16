#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)

# CONSTANTS 
MUTATIONAL_SIGNATURES = c("sToL8", "sToL10", "sToL41")
TAXANOMIC_RANKS = c("Gastropoda", "Platycheirus", "Viridiplantae")

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

# Figure 5
pal <- c(
  Gastropoda = "#000097",
  Platycheirus = "#520066",
  Viridiplantae = "#004300"
)

# pal <- c(
#   # Yellows / Oranges
#   "#FFDC00", "#E5AE38", "#FF8002",
#   
#   # Blues
#   "#30A2DA", "#2B55B9", "#1F77B4", "#4D6EFF", "#6126FF",
#   "#000097", "#0100F8", "#00457B", "#6E7CBB", "#888FFF",
#   
#   # Cyans / Teals / Aquas
#   "#04E3C8", "#17BECF", "#01BE8A", "#00BF00", "#00E83B",
#   "#83D371", "#91FF00", "#B4FF92",
#   
#   # Purples / Violets
#   "#9467BD", "#7D00A0", "#9400F5", "#895DFF", "#520066",
#   "#3A0183", "#A737AE", "#7E7CBB", "#AD8BB1", "#CE85FF",
#   
#   # Pinks / Magentas (strong, not pastel)
#   "#FF55FF", "#DD00FF", "#FF1F83", "#F500C6", "#E377C2",
#   "#B80080", "#84206F", "#FF798F", "#A56089", "#D796AB",
#   
#   # Blue-greens / Blue-purples (bridges)
#   "#95D3FF", "#7CB2FF"
# )

# Get subset of somatic mutational signature attributions
hm <- somatic_signatures %>%
  select(label, sToL8, sToL10, sToL41) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% tr$tip.label) %>%   
  column_to_rownames("label")

# max values
max(hm[,"sToL8"]) # 0.59
max(hm[,"sToL10"]) # 0.31
max(hm[,"sToL41"]) # 0.33

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
anno$Genus <- sapply(anno$Species,function(x){
  unique(metadata$Genus[metadata$Species==x])
})
anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
    Genus  %in% names(pal) ~ Genus,
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
ggsave("Extended_Data_Figure_9_alpha.pdf", plot = p, width = 7.48, height = 7.48, units = "in")

# define gradient of colours
pal_1 <- c("#FFFFFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
pal_2 <- c("#FFFFFF", "#FEE6CE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603", "#7F2704")
pal_4 <- c("#FFFFFF", "#FCE4EC", "#F48FB1", "#EC407A", "#C2185B", "#880E4F")

# layout parameters
tile_offset <- 0.02
tile_width  <- 0.02
gap         <- 0.2

# define grid colour
grid_col  <- "grey70" 

## add sToL4 mutational signature
p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1,
               hm[,"sToL8", drop = F],
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
    name = "sToL8", 
    na.value = "white",
    limits  = c(0, 0.6),
    breaks  = c(0, 0.2, 0.4, 0.6),
    labels  = c("0", "0.2", "0.4", "0.6"),
    guide = guide_colorbar(
      direction = "horizontal", 
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 1)
  )

## add sToL5 mutational signature
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2,
               hm[,"sToL10", drop = F],
               offset = tile_offset + tile_width + gap,
               width  = tile_width,
               colnames = F,
               colnames_angle = 90,
               colnames_offset_y = 0.5,
               font.size = 2,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_2,
    name = "sToL10",
    na.value = "white",
    limits  = c(0, 0.4),
    breaks  = c(0, 0.1, 0.2, 0.3),
    labels  = c("0", "0.1", "0.2", "0.3"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 2)
    )
# p2

## add sToL15 mutational signature
p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3,
               hm["sToL41"],
               offset = tile_offset + 2*(tile_width + gap),
               width  = tile_width,
               colnames = F,
               colnames_angle = 90,
               colnames_offset_y = 0.5,
               font.size = 2,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_4,
    name = "sToL41", 
    na.value = "white",
    limits  = c(0, 0.4),
    breaks  = c(0, 0.1, 0.2, 0.3),
    labels  = c("0", "0.1", "0.2", "0.3"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 3)
    )

# add legend
p3 <- p3 + 
  theme(
    legend.position = "right", 
    legend.box = "vertical", 
    legend.text = element_text(size = 6)
  )

# return plot
ggsave("Extended_Data_Figure_9_beta.pdf", plot = p3, width = 7.48, height = 7.48, units = "in")



