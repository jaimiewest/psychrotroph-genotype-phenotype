### Ordination for Biolog growth substrate use by isolates
### Jaimie West

library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(devtools)
library(data.table)
library(vegan)
library(proxy)     # For custom distance function
library(ape) # for ggtree
#BiocManager::install("ggtree")
library(ggtree) # for nicer looking dendrogram


setwd("path/to/dir")

##### Reading in data etc  ####
df = read.csv("GenIII-6-6.csv",header=TRUE,row.names=1)
metadata<-read.csv("metadata_GenIII-6-6.csv",header=TRUE,row.names=1) 

# Just use the growth substrates
growth.df = subset(df, Substrate_category_general == "growth")

# Remove the metadata columns (last 5 columns)
ncols = ncol(growth.df)
ncols
growth.df2 = growth.df[-c((ncols-4):ncols)]
View(growth.df2)
growth.df2 = t(growth.df2)


# ####Create a distance matrix using Simple Matching Coefficient Gower distance, 
# which moderately penalizes shared absences (shared 0's)
gower_dist <- cluster::daisy(growth.df2, metric = "gower", type = list(symm = colnames(growth.df2)))

# ##### Ordination with metaMDS
nmds <- vegan::metaMDS(gower_dist, distance = "none", k = 5) # report k and stress
nmds$stress # report this stress values. <0.05 is a good fit; <0.1 is fair. increase k to improve stress, see help document.
vegan::stressplot(nmds)


# Plot the ordination
plot(nmds, type = "t")  # type = "t" for text labels (isolate names)


### Make a nicer plot
#extract NMDS scores
dat = as.data.frame(vegan::scores(nmds))

# Pull isolate name back from row.names
dat$Isolate = row.names(dat)
dat$Protocol = metadata$Protocol
dat$Phylum = metadata$Phylum
dat$Source = metadata$Location

# Rename and order factors
dat$Source = recode_factor(dat$Source, "Ice wedge + thermokarst ice mix" = "Ice wedge + \nthermokarst ice mix")
dat$Source = ordered(dat$Source, levels = c("Permafrost", "Active layer soil", "Ice wedge", "Thermokarst ice", "Ice wedge + \nthermokarst ice mix"))

View(dat)

#library(ggpubr)
#show_point_shapes() # Option to explore point shape options

soureshapes = c(22, 21, 25, 24, 23) 
dat$Phylum = ordered(dat$Phylum, levels = c("Actinomycetota",
                                            "Bacillota",
                                            "Pseudomonadota",
                                            "Bacteroidota"))
palette <- c(
  "Actinomycetota" = "#99cc33",
  "Bacillota" = "blue3",
  "Pseudomonadota" = "#cc3366",
  "Bacteroidota" = "gray60"
)

p <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, aes(fill = Phylum, shape = Source), color = "black", stroke = 0.2) +
  
  #  ellipses per Phylum
  stat_ellipse(aes(group = Phylum, color = Phylum),
               type = "t", level = 0.75, linewidth = 0.5, alpha = 0.3, linetype = "solid",
               show.legend = FALSE) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_color_manual(values = palette) +
  scale_fill_manual(
    values = palette,
    guide = guide_legend(override.aes = list(shape = 21, fill = palette, color = "black"))
  ) +
  scale_shape_manual(values = soureshapes, name = "Source material") +
  theme_bw() +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"))
p


#ggsave("NMDS_for_geno_pheno_paper_Gower_ellipses.tiff",width=7, height=4, units = "in", device='tiff', dpi=150)




####### Build dendrogram #######
# Basic dendrogram
hc <- hclust(gower_dist, method = "average") 
plot(hc)

# Make a nicer dendrogram--use ggtree to customize
phylo_tree <- as.phylo(hc)

# Merge tree data with metadata
p <- ggtree(phylo_tree, layout = "rectangular", branch.length = "branch.length") + # could try fan or dendrogram for layout
  coord_cartesian(clip = 'off') +  # allow labels to extend outside
  theme_tree2() +
  theme(legend.position = "right")
 #+ geom_tiplab(align = TRUE, size = 3, aes(label = label))

# Join metadata to ggtree data
p$data <- left_join(p$data, metadata, by = c("label" = "Isolate.1"))

# Add colored labels
p = p + 
  geom_tiplab(aes(color = Phylum), 
              hjust = -0.1, size = 3, angle = 0) +
  scale_color_manual(values = c("#99cc33", "blue3",  "gray30", "#cc3366")) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 5))) +
  xlim(0, 0.28) + #manually adjust--this helps not to clip off labels.
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.position="none") #legend.position="bottom"

p # 435 x 700
ggsave("Dendrogram_for_geno_pheno_paper.tiff",width=3.3, height=5.2, units = "in", device='tiff', dpi=150)



####### PERMANOVA

#Remake gower_dist without Bacteroidota (n = 1 so we don't want to include that in statistics)
growth.df2.2 = growth.df2[rownames(growth.df2) != "AP26", ]
gower_dist.2 <- cluster::daisy(growth.df2.2, metric = "gower", type = list(symm = colnames(growth.df2.2)))

# Rename all the terrestrial ice 'locations' as the same, for statistics
metadata.2 = subset(metadata, Phylum != "Bacteroidota")
View(metadata2)
metadata.3 <- metadata.2 %>%
  mutate(Location = recode(Location,
                           "Thermokarst ice" = "Ice",
                           "Ice wedge" = "Ice",
                           "Ice wedge + thermokarst ice mix" = "Ice"))

View(metadata.2)
# Sicne we have overlapping variance with location and phylum, test by margin so order doesn't matter
P = adonis2(gower_dist.2 ~ Location + Phylum, data = metadata.3, by="margin")
P

# adonis2(formula = gower_dist.2 ~ Location + Phylum, data = metadata.3, by = "margin")
#             Df SumOfSqs      R2      F Pr(>F)   
# Location    2   0.5699 0.11885 3.4343  0.009 **
#   Phylum    2   0.5807 0.12110 3.4995  0.004 **
#   Residual  42   3.4847 0.72673                 
# Total     46   4.7951 1.00000              
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


