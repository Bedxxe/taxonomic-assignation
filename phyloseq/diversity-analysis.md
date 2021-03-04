# Diversity analysis of bacterial metagenomes 
## This file will represent contrived script for use in R with the phyloseq package

~~~
library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
~~~
{: .language-r}

setwd("C:/Users/cairo/Documents/curso-metag")
#-----------------------------
OTUS <- read_delim("JC1A.kraken_ranked-wc","\t", escape_double = FALSE, trim_ws = TRUE)
TAXAS <- read_delim("JC1A.lineage_table-wc", "\t", escape_double = FALSE, 
                    col_types = cols(subspecies = col_character(),  
                                     subspecies_2 = col_character()), trim_ws = TRUE)

names1 = OTUS$OTU
names2 = TAXAS$OTU

OTUS$OTU = NULL
TAXAS$OTU = NULL

abundances = as.matrix(OTUS)
lineages = as.matrix(TAXAS)

row.names(abundances) = names1
row.names(lineages) = names2

OTU = otu_table(abundances, taxa_are_rows = TRUE)
TAX = tax_table(lineages)

metagenome_JC1A = phyloseq(OTU, TAX)
Bacteria <- subset_taxa(metagenome_JC1A, superkingdom == "Bacteria")
##Diego: It is important not to remove any reads before the diversity analysis, because this can
##        take out singletons and doubletons, important factors for the diversity analysis.

#C# metagenome_JC1A <- prune_taxa(taxa_sums(metagenome_JC1A)>10,metagenome_JC1A)

#-----------------------------
OTUS <- read_delim("JP4D.kraken_ranked-wc","\t", escape_double = FALSE, trim_ws = TRUE)
TAXAS <- read_delim("JP4D.lineage_table-wc", "\t", escape_double = FALSE, 
                    col_types = cols(subspecies = col_character(),  
                                     subspecies_2 = col_character()), trim_ws = TRUE)

names1 = OTUS$OTU
names2 = TAXAS$OTU

OTUS$OTU = NULL
TAXAS$OTU = NULL

abundances = as.matrix(OTUS)
lineages = as.matrix(TAXAS)

row.names(abundances) = names1
row.names(lineages) = names2

OTU = otu_table(abundances, taxa_are_rows = TRUE)
TAX = tax_table(lineages)

#OTU<-subset_taxa(metagenome_JP4D, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

metagenome_JP4D = phyloseq(OTU, TAX)
Bacteria <- subset_taxa(metagenome_JP4D, superkingdom == "Bacteria")

#C# metagenome_JP4D <- prune_taxa(taxa_sums(metagenome_JP4D)>10,metagenome_JP4D)

#metagenome_JP4D1<- subset_taxa(metagenome_JP4D, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#___________________________________
merged_metagenomes = merge_phyloseq(metagenome_JC1A, metagenome_JP4D)


##Fs: 
summary(merged_metagenomes@tax_table@.Data== "NA")
summary(merged_metagenomes@tax_table@.Data[,"phylum"]== "NA")


merged_metagenomes <- subset_taxa(merged_metagenomes, phylum != "NA")

p = plot_richness(merged_metagenomes, measures = c("Observed", "Chao1", "Shannon")) 
p + geom_point(size=5, alpha=0.7)  

(sample_sums(merged_metagenomes)

percentages  = transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )

absolute_count = plot_bar(merged_metagenomes, fill="phylum")
absolute_count = absolute_count + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") + ggtitle("Absolute abundance")

percentages = plot_bar(percentages, fill="phylum")
percentages = percentages + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack") + ggtitle("Relative abundance")

absolute_count | percentages
