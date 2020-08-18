#!/usr/bin/Rscript
############################################################################################
# This script produces a cladogram ploting the information of BAPS after RhierBAPS analysis 
# using the single copy genes of E. coli core genome
# Author Arturo Vera
# Aug 2020
# quark.vera@gmail.com
#############################################################################################

#Libraries

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require("ggtree")){
	BiocManager::install("ggtree")
	library(ggtree)
}

if (!require("tidyverse")){
        install.packages("tidyverse")
	library(tidyverse)
}

if (!require("rhierbaps")){
        install.packages("rhierbaps")
	library(rhierbaps)
}

if (!require("phytools")){
        install.packages("phytools")
	library(phytools)
}

if (!require("ape")){
        install.packages("ape")
	library(ape)
}

#Load genome Names

GenomeNames <- read_tsv("GenomeTableNames.tsv.txt")


head(GenomeNames)
GenomeNames <- GenomeNames %>% 
  mutate(Isolate=paste0("genome_",Genome)) %>%
  unite(Grp_Str,sep=" ",c(Group,Strain)) %>%
  select(Isolate,Grp_Str)

#Get clusters by HierBAPS

set.seed(1234)


fasta.file.name <- "all.merged.fasta-gb"

snp.matrix <- load_fasta(fasta.file.name)

hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE)

#Save hierBaps result into a RData
save(hb.results,file="Ecoli.hb.results.RData")
load("Ecoli.hb.results.RData")

head(hb.results$partition.df)

#Teransform a DF with the clusters and lineages partitions
HB.resultsDF <- hb.results$partition.df
HB.resultsDF$Isolate <- gsub(" ","", HB.resultsDF$Isolate)

#Reading tree
iqtree <- phytools::read.newick("all.merged.fasta-gb.treefile.root.nwk")

#Use ggtree to plot the tree
gg <- ggtree(iqtree, layout = "rectangular")
gg <- gg %<+% HB.resultsDF

gg <- gg %<+% GenomeNames

#Add Grp and strain
gg <- gg + 
  geom_tippoint(aes(color = factor(`level 1`))) + 
  geom_tiplab(aes(label=Grp_Str))

#Add lineages
gglevel2 <- gg+ 
  geom_tiplab(aes(label = `level 2`), size = 3, offset = 1,align = TRUE,linetype = "dotted")
gglevel2

ggsave(gglevel2,file="E.coli.Tree.Subclades.pdf",width = 35,height = 35)

#Ading Heatmap for lineage information

#Get the lineages into a DF

DF1 <- data.frame(HB.resultsDF$Isolate,HB.resultsDF$`level 2`)

colnames(DF1) <- c("label","linage")

DF1 <- DF1 %>% mutate(Linage=paste0("linage_",linage))

rownames(DF1) <- DF1$label
DF1 <- select(DF1,Linage)

#Select colors for the lineages

A <- RColorBrewer::brewer.pal(12,"Set3")
B <- RColorBrewer::brewer.pal(12,"Paired")

DF2 <- data.frame(Color=c(A,B),Linage=sort(unique(DF1$Linage)))

DF3 <- DF1 %>% full_join(DF2,by="Linage")

#Create a vector with the colors foreach lineage

heatmap.colors <- as.character(DF3$Color)
names(heatmap.colors) <- DF1$Linage

#Join all data into a new ggtree object

gg2 <- ggtree(iqtree, layout = "circular",branch.length = "none")
gg2 <- gg2 %<+% HB.resultsDF

gg2 <- gg2 %<+% GenomeNames

gg3 <- gg2 + 
  geom_tippoint(aes(color = factor(`level 1`)))+
  scale_color_discrete(name='Cluster')+
  geom_tiplab2(aes(label=Grp_Str),size=2,offset=1)+ 
  geom_tiplab2(aes(label=`level 2`),align=T, linetype=NA, 
               size=2, offset=14)

#Add the heatmap with lineage colors

GHetamp <- gheatmap(gg3, DF1, width=.4, 
         offset=16, 
         colnames=F,
         colnames_offset_y = 1,
         color=NULL) +
  scale_fill_manual(breaks =str_sort(unique(DF3$Linage),
                                     numeric = TRUE),
                    values = heatmap.colors,
                    name="Lineages")

#Saving the tree into pdf

ggsave(GHetamp,file="ecoli.tree.rhierBAPS.final.pdf",width = 12,height = 10)
