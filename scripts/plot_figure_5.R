#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)

# CONSTANTS 
MUTATIONAL_SIGNATURES = c("gToL1", "gToL3", "gToL4", "gToL6")
TAXANOMIC_RANKS = c("Coleoptera", "Chordata", "Viridiplantae", "Vespidae", "Tenthredinidae")

# setwd
setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-3")

# Load data
nwk = readLines("tree.nwk")
metadata <- read.csv("dtol_all_samples.taxonomic_classification.csv")
germline_signatures <- read.csv("germline_mutational_signature_attributions.x0_excluded.csv", check.names = F)

# Manipulate data
colnames(germline_signatures) = paste0('gToL',colnames(germline_signatures))
germline_signatures[,1] <- ifelse(grepl("\\.", germline_signatures[,1]),
                                  sub(".*\\.", "", germline_signatures[,1]),
                                  germline_signatures[,1])
germline_signatures$label=NA
germline_signatures$label=sapply(germline_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  }
  label
})

# Manipulate newick tree
nwk <- gsub('_', '^', nwk, fixed = TRUE)
nwk <- gsub(" +", "_", nwk)   # turn spaces inside labels into underscores
tr   <- read.tree(text = nwk)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)

# Figure 5
pal <- c(
  Coleoptera = "#A3843F",
  Chordata = "#765FA6",
  Viridiplantae = "#408145",
  Fungi = "#925C35"
)

# Get subset of somatic mutational signature attributions
hm <- germline_signatures %>%
  select(label, gToL1, gToL3, gToL4, gToL6) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% tr$tip.label) %>%   
  column_to_rownames("label")

# max values
# max(hm[,"gToL1"]) # 0.96
# max(hm[,"gToL3"]) # 0.81
# max(hm[,"gToL4"]) # 0.62
# max(hm[,"gToL6"]) # 0.39

# build annotation
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
p <- ggtree(tr, layout = "circular", size = 0.25, color = "#090954") +
  geom_hilight(
    data = hilight_df,
    aes(node = node, fill = ColorGroup),
    alpha=1,
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
ggsave("Figure_5_alpha.pdf", plot = p, width = 7.48, height = 7.48, units = "in")

# define gradient of colours
pal_1 <- c("#F7FCF5","#E1F3DC","#BCE4B5","#8ED08B","#56B567","#2C944C","#05712F","#00441B")
pal_2 <- c("#F7FBFF","#DBE9F6","#BAD6EB","#89BEDC","#539ECD","#2B7BBA","#0B559F","#08306B")
pal_3 <- c("#FFF5EB","#FEE3C8","#FDC692","#FDA057","#F67824","#E05206","#AD3803","#7F2704")
pal_4 <- c("#FCFBFD","#ECEBF4","#D1D2E7","#AFAED4","#8D89C0","#705EAA","#572C92","#3F007D")

# layout parameters
tile_offset <- 0.02
tile_width  <- 0.02
gap         <- 0.2

# define grid colour
grid_col  <- "grey70" 

## add gToL1 mutational signature
p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1,
               hm[,"gToL1", drop = F],
               offset = tile_offset,
               width  = tile_width,
               colnames = F,
               colnames_angle = 90,
               colnames_offset_y = 0.5,
               font.size = 8,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_3,
    name = "gToL1", 
    na.value = "white",
    limits  = c(0, 1.0),
    breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels  = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
    guide = guide_colorbar(
      direction = "horizontal", 
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 1)
  )

## add gToL3 mutational signature
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2,
               hm[,"gToL3", drop = F],
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
    name = "gToL3",
    na.value = "white",
    limits  = c(0, 0.9),
    breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 0.9),
    labels  = c("0", "0.2", "0.4", "0.6", "0.8", "0.9"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 2)
    )
# p2

## add gToL4 mutational signature
p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3,
               hm["gToL4"],
               offset = tile_offset + 2*(tile_width + gap),
               width  = tile_width,
               colnames = F,
               colnames_angle = 90,
               colnames_offset_y = 0.5,
               font.size = 2,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_1,
    name = "gToL4", 
    na.value = "white",
    limits  = c(0, 0.7),
    breaks  = c(0, 0.2, 0.4, 0.6, 0.7),
    labels  = c("0", "0.2", "0.4", "0.6", "0.7"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 3)
    )

## add gToL6 mutational signature
p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4,
               hm["gToL6"],
               offset = tile_offset + 3.03*(tile_width + gap),
               width  = tile_offset,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col
               ) +
  scale_fill_gradientn(
    colours = pal_4,
    name = "gToL6",
    na.value = "white",
    limits  = c(0, 0.4),
    breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
    labels  = c("0", "0.1", "0.2", "0.3", "0.4"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = unit(3, "pt"),
      barwidth = unit(60, "pt"),
      order = 4)
    )

# add legend
p4 <- p4 + 
  theme(
    legend.position = "right", 
    legend.box = "vertical", 
    legend.text = element_text(size = 6)
  )

# return plot
ggsave("Figure_5_beta.pdf", plot = p4, width = 7.48, height = 7.48, units = "in")



