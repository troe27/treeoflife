#!/usr/bin/env Rscript

# Load packages
library(ape)
library(argparse)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(stringr)
library(dplyr)
library(Polychrome)

## setwd
setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-3")

# Load data
nwk = readLines("tree.nwk")
metadata = read.csv("dtol_all_samples.taxonomic_classification.csv")

# ---- CLI ----
parser <- ArgumentParser(description = "Plot Supplementary Figure 3")
parser$add_argument(
  "--newick", 
  required = TRUE,
  help = "Path to Newick tree"
)
parser$add_argument(
  "--taxonomic-classification", 
  required = TRUE,
  help = "Path to dtol_all_samples.taxonomic_classification.csv"
)
args <- parser$parse_args()

# get parsed labels
get_expr_labels <- function(labels){
  esc <- function(s) gsub("'", "\\\\'", s, perl = TRUE)  # escape single quotes
  
  has_paren = grepl("\\(", labels)
  
  # get species name
  species   <- sub("\\s*\\(.*$", "", labels) 
  
  # get sample name
  paren     <- sub("^[^\\(]*", "", labels)
  
  species_e <- esc(species)
  paren_e   <- esc(paren)
  
  ifelse(
    has_paren,
    paste0("italic('", species_e, "')~'", paren_e, "'"),  # italic species + plain "(...)"
    paste0("italic('", species_e, "')")                   # all italic when no "(...)"
  )
}

# Load data
nwk = readLines(args$newick)
metadata = read.csv(args$taxonomic_classification)

# Manipulate data
nwk = gsub("_", "^" , nwk, fixed = TRUE)
nwk <- gsub(" +", "_", nwk)   # turn spaces inside labels into underscores
tr = read.tree(text=nwk)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)

# # Taxanomic rank colours
rank_colours <- c(
  # Kingdom: Viridiplantae
  Bryopsida = "#F8AE15",
  Magnoliopsida = "#066E80",
  
  # Kingdom: Fungi
  Agaricomycetes = "#CE4490",
  
  # Kingdom: Metazoa
  Insecta  = "#48B389",
  Arachnida    = "#B1CB17",
  Actinopteri  = "#355199",
  Mammalia     = "#5EA6DD",
  Aves         = "#3B708A",
  # Mollusca
  Bivalvia = "#BFADE0",
  Gastropoda = "#C4C6E5",
  # Bryozoa
  Gymnolaemata = "#E6AED1",
  # Annelida
  Clitellata = "#401476",
  Polychaeta = "#826DB3"
)


# Build data e annotation
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
anno <- anno %>% mutate(label_expr = get_expr_labels(anno$label))
anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(rank_colours) ~ Phylum,                # use phylum if in rank_coloursette
      Order %in% names(rank_colours) ~ Order,                # use phylum if in rank_coloursette
      Class %in% names(rank_colours) ~ Class,                # use phylum if in rank_coloursette
      Kingdom == "Viridiplantae" ~ "Viridiplantae",
      Kingdom == "Fungi" ~ "Fungi",
      TRUE ~ NA_character_ 
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(rank_colours))  # enforce order
  ) %>% 
  select(label, label_expr, ColorGroup)

# build data frame
tr$tip.label = anno$label
grp_list <- split(anno$label, anno$ColorGroup)
hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

# assign expression label
tr$tip.label = anno$label_expr

# plot
p <- ggtree(tr, layout = "circular", size=0.25, color = "#090954") +
  geom_hilight(
    data = hilight_df,
    aes(node = node, fill = ColorGroup),
    alpha=1,
  ) +
  geom_tiplab(
    aes(label = label),
    parse = TRUE,
    family = "Helvetica",
    size = 0.8
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    text = element_text(family = "Helvetica")
  ) +
  scale_fill_manual(values = rank_colours, name = "Taxanomic rank")
ggsave("SF3.pdf", plot = p, width = 7.48, height = 7.48, units = "in")


