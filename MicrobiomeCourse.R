





# Pass to phyloseq for downstream analyses
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())


seqtab.nochim<-readRDS("seqs.rds")
taxa<-readRDS("completed_taxa.rds")


samples.out <- rownames(seqtab.nochim)
sdf<-read.csv("~/sdf.csv")

samples<-data.frame(sdf)

barcodesequence <- samples$barcodesequence
biome<-samples$biome
collection.date<-samples$collection.date
feature<-samples$feature
material<-samples$material
country <- samples$country
latitude <- samples$latitude
longitude<-samples$longitude


sampleinfo<-data.frame(barcodesequence=barcodesequence, biome=biome, collection.date=collection.date, feature=feature, material=material, country=country, latitude=latitude, longitude=longitude)
rownames(sampleinfo) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sampleinfo), tax_table(taxa))
ps

### Identify and remove potential contaminant species
library(decontam)
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=material)) + geom_point()


### Filter mito/chloro/etc
install.packages("readr")
library(readr)
library(stringr)
ps2 <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Family  != "Streptophyta" &
      Class   != "Chloroplast"
  )

taxa_names(ps2)<-paste0("Seq", seq(ntaxa(ps2)))


x <- mean(sample_sums(ps2))/mean(sample_sums(ps))
x
mean(sample_sums(ps2))

readsumsdf = data.frame(nreads = sort(taxa_sums(ps2), TRUE), sorted = 1:ntaxa(ps2), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps2), 
                                                        TRUE), sorted = 1:nsamples(ps2), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")






library(vegan)
tab <- otu_table(ps2)
class(tab)<- "matrix"
tab <- t(tab)
bactcurve<-rarecurve(tab, step=20, label=FALSE, xlim=c(0,250))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome", force = TRUE)
library(microbiome)

write_phyloseq(ps2, 'TAXONOMY')

transrareps2 <- transform_sample_counts(ps.21, function(x) x/sum(x)) #get abundance in %

bact <-psmelt(transrareps2)
head(bact)



aggp<-aggregate(bact$Abundance, by=list(bact$Family), FUN=mean)
aggpsd<-aggregate(bact$Abundance, by=list(bact$Family), FUN=sd)
aggpsd
aggp

bact$flwpart<-sample_data(ps2)$material
bactrich <- estimate_richness(ps2)




library(nlme)
mod<-lme(Observed ~ material, random = ~1|longitude, data=bactrich)
anova(mod)
boxplot(Observed~material,  data=bactrich, colo)

mod1<-lme(Shannon~material,random = ~1|longitude, data=bactrich)
anova(mod1)
boxplot(Shannon~material, data=bactrich)




p<-ggplot(bactrich, aes(x=material, y=Observed))+geom_boxplot()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+xlab("Floral part")+ylab("Observed richness")
p
p2<-ggplot(bactrich, aes(x=material, y=Shannon))+geom_boxplot()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+xlab("Floral part")+ylab("Observed richness")
p2


na.omit(transrareps2)
na.omit(ps2)

ps.prop<-transform_sample_counts(ps.21, function(otu) otu/sum(otu))

na.omit(ps.prop)

ord.ps<-ordinate(ps.prop, method="NMDS", distance="bray", na.rm = TRUE, autotransform = TRUE)



plot_ordination(ps.prop, ord.ps, color="material")+stat_ellipse(position="identity", level=0.95)+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+scale_colour_brewer(palette="paired", labels=c("Bark", "Flower", "Lava", "Leaf", "Nectar", "Root", "Stamen", "Style"), name="Management")








library(vegan)
bact_bray = phyloseq::distance(ps.21, method="bray")
samplebact<-data.frame(sample_data(ps.21))
adonis(bact_bray~material, data=samplebact)
























y1 <- tax_glom(transrareps2, taxrank = 'Family') # agglomerate taxa

y2 = merge_samples(y1, "mgmt") # merge samples on sample variable of interest
sample_data(y2)$mgmt <- levels(sample_data(y1)$mgmt)
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %

# create dataframe from phyloseq object 
y4 <- psmelt(y3) 
#convert to character
y4$Family <- as.character(y4$Family) 
y4$Family[y4$Abundance < 0.01] <- "Family < 1% abund." #rename genera with < 1% abundance
library(grDevices)
library(RColorBrewer)
colourCount = length(unique(y4$Family))
colourCount
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
head(y4)

aggp<-aggregate(y4$Abundance, by=list(y4$Family), FUN=mean)







##################################
## Taxa of interest




install.packages("reshape")
install.packages("reshape2")
library(reshape)
library(reshape2)

list_df<-list(bac2melt, pseud2melt)
merged_df<-merge_all(list_df)
head(merged_df)
merged_df

write.csv(merged_df, file="mergedbact.csv")

aggbact<-aggregate(merged_df$Abundance, by=list(merged_df$mgmt, merged_df$Genus), FUN=mean)
aggbactsd<-aggregate(merged_df$Abundance, by=list(merged_df$mgmt, merged_df$Genus), FUN=sd)
aggbact$sd<-aggbactsd$x
aggbact$n<-(10)
aggbact$se<-aggbact$sd/(sqrt(aggbact$n))
names(aggbact)<-c("mgmt", "Genus", "Abundance", "sd", "n", "se")





# proportional bar graph

bact$group <- paste0(bact$Phylum, "-", bact$Family, sep = "")

colours <-ColourPalleteMulti(bact, "Phylum", "Family")

ggplot(df2, x=material) + 
  geom_bar(aes(fill = Phylum), colour = "grey") +
  scale_fill_manual("Subject", values=colours, guide = "none")



GlobalPatterns_prop = transform_sample_counts(ps2, function(x)       100 * x/sum(x)) 
plot_bar(GlobalPatterns_prop, fill = "Phylum", x= "material")

library(forcats)
library(dplyr)

phylums <- c('Proteobacteria','Bacteroidetes','Firmicutes')
df <- bact
df$Phylum[!df$Phylum %in% phylums] <- "Others"
df$Family[!df$Phylum %in% phylums] <- "Others"

df$Family[df$Phylum=="Proteobacteria" & 
            !df$Family %in% c('Alcaligenaceae','Enterobacteriaceae')] <- "Other Protobacteria"

df$Family[df$Phylum=="Bacteroidetes" &
            !df$Family %in% c('Bacteroidaceae','Rikenellaceae','Porphyromonadaceae')] <- "Other Bacteroidetes"

df$Family[df$Phylum=="Firmicutes" & 
            !df$Family %in% c('Lactobacillaceae','Clostridiaceae','Ruminococcaceae','Lachnospiraceae')] <- "Other Firmicutes"

df2 <- select(df, Sample, Phylum, Family, material) %>%
  mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
         Family=fct_reorder(Family, 10*as.integer(Phylum) + grepl("Others", Family)))
 
ColourPalleteMulti <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

df2<- psmelt(ps.21)

colours <- ColourPalleteMulti(df2, "Phylum", "Family")

ggplot(df2, aes(x=material, fill = Phylum)) + 
  geom_bar(position="fill", colour = "grey") +  # Stacked 100% barplot
  
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance")
 


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

install.packages("ggraph") # for taxatree_plots()
install.packages("DT") # for tax_fix_interactive()
install.packages("corncob") # for example datasets and beta binomial models


# To install the very latest version:
devtools::install_github("david-barnett/microViz")
library(microViz)

pseq <- ps2 %>%
  tax_fix() %>%
  phyloseq_validate()
ord_explore(pseq)




ps.21 %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey") +
 
   coord_flip() +
  facet_wrap("material", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
#> Registered S3 method overwritten by 'seriation':
#>   method         from 
#>   reorder.hclust vegan




# perform ordination
unconstrained_aitchison_pca <- ps.21 %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>%
  tax_agg("Family") %>%
  tax_transform("clr") %>%
  ord_calc()
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.
# ord_calc will automatically infer you want a "PCA" here
# specify explicitly with method = "PCA", or you can pick another method

# create plot
pca_plot <- unconstrained_aitchison_pca %>%
  ord_plot(
    plot_taxa = 1:6, colour = "material", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 0.5),
    auto_caption = 8
  )

# customise plot
customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = material, colour = material), linewidth = 0.3) + # linewidth not size, since ggplot 3.4.0
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 0.5, clip = "off") # makes rotated labels align correctly

# show plot
customised_plot


ps.21 <- ps2%>% subset_taxa(!is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Genus) & !Genus %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Family) & !Family %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Order) & !Order %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Kingdom) & !Kingdom %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Class) & !Class %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Species) & !Species %in% c("", "uncharacterized"))
  
devtools::install_github("gmteunisse/fantaxtic", force = TRUE)
1require("fantaxtic")
require("phyloseq")
library(fantaxtic)
# Load the data
data(ps)
?data
# Get the most abundant phyla and the most abundant families within those phyla
top_nested <- nested_top_taxa(ps.21,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 3, 
                              n_nested_taxa = 3)

# Plot the relative abundances at two levels.
plot_nested_bar(ps_obj = top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Family")





htmp <- ps.21 %>%
  ps_mutate(material = as.character(material)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(ps.21, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(palette = viridis::turbo(11)),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      material = anno_sample_cat(
        var = "material", col = c(flower = "grey35", stamen = "grey85", style = "grey25", lava = "grey45", root = "grey55", nectar = "grey65", bark = "grey95", leaves = "grey15"),
        box_col = NA, legend_title = "material", size = grid::unit(4, "mm")
      )
    )
  )

ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)





# calculate distances
aitchison_dists <- ps.21 %>%
  tax_filter(min_prevalence = 0.1) %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("aitchison")
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99, # you should use at least 999!
    variables = "material"
  )
#> 2022-11-24 01:58:13 - Starting PERMANOVA with 99 perms with 1 processes
#> 2022-11-24 01:58:13 - Finished PERMANOVA

# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()
#>            Df  SumOfSqs         R2        F Pr(>F)
#> bmi_group   2  104.0678 0.04177157 4.773379   0.01
#> Residual  219 2387.2862 0.95822843       NA     NA
#> Total     221 2491.3540 1.00000000       NA     NA

# view the info stored about the distance calculation
info_get(aitchison_perm)
#> psExtra info:
#> tax_agg = "Family" tax_trans = "identity" dist_method = "aitchison"

