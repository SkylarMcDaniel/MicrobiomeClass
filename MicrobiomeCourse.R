# Used Dada2 pipeline on a remote computing cluster to assign taxa, and species 
# now passing to phyloseq for downstream analyses
library(phyloseq)
library(ggplot2)

# Reading in RDS files generated on the remote cluster 
seqtab.nochim<-readRDS("seqs.rds")
taxa<-readRDS("completed_taxa.rds")

# Run biosampleparser to download metadata for the study
library(rentrez)
library(xml2)
source("C:/Users/A02342347/Documents/BioSampleParser.R")
sdf <- BioSampleParser(query = "PRJEB7828")
write.csv(sdf, "sdf.csv")


# Creating a dataframe for the metadata 
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


# Use data frame to create phyloseq object 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sampleinfo), tax_table(taxa))



# Filter mito/chloro/etc
library(readr)
library(stringr)
library(dplyr)
ps2 <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )


# convert taxa_names to an arbitrary "Seq" character 
taxa_names(ps2)<-paste0("Seq", seq(ntaxa(ps2)))


# what percent of samples were lost?
x <- mean(sample_sums(ps2))/mean(sample_sums(ps))
x # we retain 38.6% of samples after this step
# This leads to problems when trying to estimate richness, ordinate the data because this trims out all 
# low abundance sequences
# For this reason I opt to use the unfiltered data (ps)  


# Filter for NAs 
# This allows me to perform ordinations on the filtered data 
ps1 <- ps2%>% subset_taxa(!is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Genus) & !Genus %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Family) & !Family %in% c("", "uncharacterized"))%>%
subset_taxa(!is.na(Order) & !Order %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Kingdom) & !Kingdom %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Class) & !Class %in% c("", "uncharacterized"))%>%
  subset_taxa(!is.na(Species) & !Species %in% c("", "uncharacterized"))
  


# now lets generate plots showing us reads per sample ordered high to low, 
# and reads per observed otu ordered high to low 
# these are for unfiltered data 
readsumsdf = data.frame(nreads = sort(taxa_sums(ps), TRUE), sorted = 1:ntaxa(ps), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps1), 
                                                        TRUE), sorted = 1:nsamples(ps), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


# and for the unfiltered data 
readsumsdf1 = data.frame(nreads = sort(taxa_sums(ps1), TRUE), sorted = 1:ntaxa(ps1), 
                        type = "OTUs")
readsumsdf1 = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps1), 
                                                        TRUE), sorted = 1:nsamples(ps1), type = "Samples"))
title = "Total number of reads"
p1 = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p1 + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")




# make OTU table into a matrix, then transpose to generate a rarefaction curve
tab <- otu_table(ps)
class(tab)<- "matrix"
tab <- t(tab)
bactcurve<-rarecurve(tab, step=20, label=FALSE, xlim=c(0,250))
# Seems that we should not rarefy here, with each sample read we gain resolution
# Same holds true when looking at the filtered data 



# put sample counts in terms of relative abundance 
transrareps <- transform_sample_counts(ps, function(x) x/sum(x))  
transrareps1 <- transform_sample_counts(ps1, function(x) x/sum(x)) 
transrareps2 <- transform_sample_counts(ps2, function(x) x/sum(x))  
# turn this new object into a data frame 
bact <-psmelt(transrareps)
head(bact)
bact1 <-psmelt(transrareps1)
head(bact1)



aggp<-aggregate(bact$Abundance, by=list(bact$Family), FUN=mean)
aggpsd<-aggregate(bact$Abundance, by=list(bact$Family), FUN=sd)



bact$flwpart<-sample_data(ps)$material
bactrich <- estimate_richness(ps)



library(nlme)
library(grDevices)
library(RColorBrewer)
# using a linear model to determine differences in observed richness
mod<-lme(Observed ~ material, random= ~1|collection.date, data=bactrich)
anova(mod) # p-value = 0.0265
p<-ggplot(bactrich, aes(x=material, y=Observed, fill=material))+ scale_fill_brewer(palette = "Set3")+geom_boxplot()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+xlab("Floral part")+ylab("Observed richness")
p # when looking through, it seems observed richness is skewed by the lava samples, and there is 
# not much diversity in the microbiome of the plant samples 

# linear model using shannon index instead 
mod1<-lme(Shannon~material,random = ~1|collection.date, data=bactrich)
anova(mod1) # p-value 0.382
boxplot(Shannon~material, data=bactrich)
p2<-ggplot(bactrich, aes(x=material, y=Shannon, fill = material))+scale_fill_brewer(palette = "Dark2") +geom_boxplot()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+xlab("Floral part")+ylab("Observed richness")
p2 # very little difference between samples here 


# using ordination to hopefully bring out some differences between plant tissue 
ord.ps<-ordinate(transrareps, method="NMDS", distance="bray", autotransform = TRUE)
plot_ordination(transrareps, ord.ps, color="material")+stat_ellipse(position="identity", level=0.95)+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+scale_colour_brewer(palette="paired", labels=c("Bark", "Flower", "Lava", "Leaf", "Nectar", "Root", "Stamen", "Style"), name="Management")

ord.ps1<-ordinate(transrareps, method="PCoA", distance="bray", autotransform = TRUE)
plot_ordination(transrareps, ord.ps1, color="material")+stat_ellipse(position="identity", level=0.95)+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+scale_colour_brewer(palette="Set3", labels=c("Bark", "Flower", "Lava", "Leaf", "Nectar", "Root", "Stamen", "Style"), name="Management")
# These show the same thing, not much dissimilarity between plant tissues 
# For that reason I will not use PERMANOVA/adonis



# function to build color pallates for large datasets
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



# proportional bar graph
bact$group <- paste0(bact$Phylum, "-", bact$Family, sep = "")
colours <- ColourPalleteMulti(bact, "Phylum", "Family")
bact%>%
ggplot() + 
  geom_bar(aes(x=material, y = Abundance, fill = Phylum), colour = "grey", stat = 'identity') +
  scale_fill_manual("Subject", values=colours, guide = "none")

# same for filtered group
bact1$group <- paste0(bact1$Phylum, "-", bact1$Family, sep = "")
colours1 <- ColourPalleteMulti(bact1, "Phylum", "Family")
bact1%>%
  ggplot() + 
  geom_bar(aes(x=material, y = Abundance, fill = Phylum), colour = "grey", stat = 'identity') +
  scale_fill_manual("Subject", values=colours1, guide = "none")



# Another way to visualize 
GlobalPatterns_prop = transform_sample_counts(ps, function(x)  100 * x/sum(x)) 
plot_bar(GlobalPatterns_prop, fill = "Phylum", x= "material")




library(ComplexHeatmap)
# Heatmap to show clustering of ASVs 
htmp <- ps1 %>%
  ps_mutate(material = as.character(material)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(ps1, n = 30), grid_col = NA, name = "Log2p",
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








