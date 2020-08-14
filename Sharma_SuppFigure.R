#Figure S1----------------Alpha diversity analysis on Rarefied Data
# Read in OTU table
otu_table_in <- read.csv("feature_table.txt", sep = "\t", row.names = 1)
# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("Metadata.txt", sep="\t", row.names = 1, header=T)
# Read in tree
phy_tree <- read_tree("tree.nwk")
# Import all as phyloseq objects
OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)
# Same sample names
sample_names(OTU)
sample_names(META)
# Finally merge!
ps <- phyloseq(OTU, TAX, META, phy_tree)

#---- Rarefied at 1000 sequencing depth (4 samples have less then 1k depth)
ps.prune = prune_samples(sample_sums(ps) > 1000, ps)
set.seed(1234)
ps.rare = rarefy_even_depth(ps.prune, sample.size = 1000)
ps_rarefied1k = prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)

p_rarefied <- plot_richness(ps_rarefied1k, "Group", measures = c("Observed", "Chao1", "Shannon"))
p_rarefied <- p_rarefied + geom_boxplot(data = p_rarefied$data, aes(x= Group, y = value, color = Group), alpha = 0.1) + scale_color_manual(values=Colors) + labs(x="",y="Alpha Diversity Measure") + theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("phyloseq_analysis-richness_estimates_rarefiedAt1K.pdf", p_rarefied, width = 6, height = 3)
jpeg("Phyloseq_Richness_rarefiedAt1K.jpg", height = 3, width = 6, units = 'in', res = 600)
p_rarefied
dev.off ()

erich_rarefied1k <- estimate_richness(ps_rarefied1k, measures = c("Observed", "Chao1", "Shannon"))
erich_rarefied1k_group <- merge(erich_rarefied1k, metadata, by=0, all=TRUE)
write.table (erich_rarefied1k_group, file = "Richness_rarefiedat1k.txt", sep = "\t")
wilcox.test(erich_rarefied1k_group$Observed ~ erich_rarefied1k_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 183, p-value = 0.04892
wilcox.test(erich_rarefied1k_group$Chao1 ~ erich_rarefied1k_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 164, p-value = 0.01664
wilcox.test(erich_rarefied1k_group$Shannon ~ erich_rarefied1k_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 223, p-value = 0.2666

#Figure S2----------------------------------- Beta diversity analysis
#Create Phyloseq object at rarefied 1k on filtered data

otu_table_in <- read.csv("feature_table.txt", sep = "\t", row.names = 1)
#-- Additional steps to filter OTUs
otu_table_in_t <- data.frame(t(otu_table_in))
otu.sums <- colSums(otu_table_in_t)
otu_table_filtered <- otu_table_in_t[ , which(otu.sums > 5)]
otu_table_filtered_t <- data.frame(t(otu_table_filtered)) #297 OTUs have sume less then 5

otu_table_in <- as.matrix(otu_table_filtered_t)

OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)

# Finally merge to create Phyloseq object! Filtered for atleast more than 5 OTUs
ps <- phyloseq(OTU, TAX, META, phy_tree)
ps.prune = prune_samples(sample_sums(ps) > 1000, ps)
set.seed(1234)
ps.rare = rarefy_even_depth(ps.prune, sample.size = 1000)
ps_rarefied1k = prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)

Bushman2  = transform_sample_counts(ps_rarefied1k, function(x) x / sum(x) )

# A) --------- UniFrac PCoA at rarefied 1k
UniFrac_distances_rare <- UniFrac(Bushman2, weighted=TRUE)
UniFrac_distances_rare <- as.matrix(UniFrac_distances_rare)
UniFrac_dist_rare_column <- melt(UniFrac_distances_rare)
write.table (UniFrac_dist_rare_column, file = "UniFrac_distances_Rare1k.txt", sep = "\t")

Uni_pcoa_rare <- pcoa(UniFrac_distances_rare)
Uni_pcoa_rare$values[1:2,]
mds.var.per = round(Uni_pcoa_rare$values$Eigenvalues/sum(Uni_pcoa_rare$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("UniFrac_PCoA_rare1k.png", height = 5, width = 5, units = 'in', res = 600)
plot(Uni_pcoa_rare$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[Bushman2@sam_data$Group], pch=c(21,22)[Bushman2@sam_data$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=Bushman2@sam_data$Group, cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa_rare$vectors[,1:2], Bushman2@sam_data$Group, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Uni_pcoa_rare$vectors[,1:2], Bushman2@sam_data$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topleft", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomleft", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis(UniFrac_distances_rare ~ Bushman2@sam_data$Group)
#          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Bushman2@sam_data$Group  1   0.14138 0.141377  1.8122 0.04552  0.087 .


# B) --------- Unweighted UniFrac PCoA at rarefied 1k
Unweigh_UniFrac_distances_rare <- UniFrac(Bushman2, weighted=FALSE)
Unweigh_UniFrac_distances_rare <- as.matrix(Unweigh_UniFrac_distances_rare)
Unweigh_UniFrac_dist_rare_column <- melt(Unweigh_UniFrac_distances_rare)
write.table (Unweigh_UniFrac_dist_rare_column, file = "Unweigh_UniFrac_distnaces_rare1k.txt", sep = "\t")

Unweigh_Uni_pcoa_rare <- pcoa(Unweigh_UniFrac_distances_rare)
Unweigh_Uni_pcoa_rare$values[1:2,]
mds.var.per = round(Unweigh_Uni_pcoa_rare$values$Eigenvalues/sum(Unweigh_Uni_pcoa_rare$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("Unweigh_UniFrac_PCoA_rare1k.png", height = 5, width = 5, units = 'in', res = 600)
plot(Unweigh_Uni_pcoa_rare$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[Bushman2@sam_data$Group], pch=c(21,22)[Bushman2@sam_data$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=Bushman2@sam_data$Group, cex=0.3, font=1, pos=1)
ordiellipse(Unweigh_Uni_pcoa_rare$vectors[,1:2], Bushman2@sam_data$Group, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Unweigh_Uni_pcoa_rare$vectors[,1:2], Bushman2@sam_data$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topright", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomright", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis(Unweigh_UniFrac_distances_rare ~ Bushman2@sam_data$Group)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Bushman2@sam_data$Group  1    0.3343 0.33430  2.1466 0.05347  0.003 **

#C)---------- Bray-Curtis PCoA at rarefied data 1k
OTU1 = as(otu_table(ps_rarefied1k), "matrix")
if(taxa_are_rows(Bushman2)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
write.table(OTUdf, file="checkOTU-table.txt", sep="\t")
metadata_rare1k <- read.csv(file="MetadataForRarefied.txt", sep="\t", row.names = 1, header = T)
write
Bray_pcoa_1k <-pcoa(vegdist(OTUdf, "bray"))
Bray_pcoa_1k$values[1:2,]
mds.var.per = round(Bray_pcoa_1k$values$Eigenvalues/sum(Bray_pcoa_1k$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX_1k <- Bray_pcoa_1k$vectors[,1:2]
Bray_PCoA_MATRIX_1k <- data.frame(Bray_PCoA_MATRIX_1k)
Bray_PCoA_MATRIX_1k_New <- cbind(Bray_PCoA_MATRIX_1k, metadata_rare1k)

pc <- c(1,2)
jpeg("Bray_RareFied1k-PCoA.png", height = 5, width = 5, units = 'in', res = 600)
plot(Bray_pcoa_1k$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[metadata_rare1k$Group], pch=c(21,22)[metadata_rare1k$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Bray_pcoa$vectors[,1:2], labels=metadata$Group, cex=0.3, font=1, pos=1)
ordiellipse(Bray_pcoa_1k$vectors[,1:2], metadata_rare1k$Group, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Bray_pcoa_1k$vectors[,1:2], metadata_rare1k$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topright", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomright", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
Bray_distances_rareFied <-vegdist(OTUdf, "bray")
adonis(Bray_distances_rareFied ~ metadata_rare1k$Group)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#metadata$Group  1     0.489 0.48903  1.9432 0.03814   0.03 *


#Figure S3-------- Enrichment of predicted pathways in Cases and Controls
#--- Controls
Selected_path_co <- read.csv (file="Path_proportions_Controls.txt", row.names = 1, header=T, sep="\t")
Selected_path_co <- data.frame(t(Selected_path_co))
Selected_path_co <- sqrt(sqrt(Selected_path_co))
Group <- metadata$Group
Selected_path_co_Groups <- cbind(Selected_path_co, Group)

path_Co_group_melted <- melt(Selected_path_co_Groups, id.vars = "Group")

library(dplyr)
df.summary_Co <- path_Co_group_melted %>%
  group_by(Group, variable) %>%
  summarise(
    sd = sd(value),
    len = mean(value)
  )
df.summary_Co

jpeg("Sel_controls_Path_new1.jpg", height = 4, width = 8, units = 'in', res = 600)
ggplot(df.summary_Co, aes(x=variable, y=len, fill=Group)) + labs(x="",y="sqrt(pathway's relative abundance)") + theme_classic() +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Group),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = Group), position = position_dodge(0.3)) +
  scale_color_manual(values = c("darkred", "darkgray")) + coord_flip()
dev.off ()

#--- Cases
Selected_path_ca <- read.csv (file="Path_proportions_Cases.txt", row.names = 1, header=T, sep="\t")
Selected_path_ca <- data.frame(t(Selected_path_ca))
Selected_path_ca <- sqrt(sqrt(Selected_path_ca))
Group <- metadata$Group
Selected_path_ca_Groups <- cbind(Selected_path_ca, Group)

path_Ca_group_melted <- melt(Selected_path_ca_Groups, id.vars = "Group")

library(dplyr)
df.summary_Ca <- path_Ca_group_melted %>%
  group_by(Group, variable) %>%
  summarise(
    sd = sd(value),
    len = mean(value)
  )
df.summary_Ca

jpeg("Sel_cases_Path_new1.jpg", height = 8, width = 8, units = 'in', res = 600)
ggplot(df.summary_Ca, aes(x=variable, y=len, fill=Group)) + labs(x="",y="sqrt(pathway's relative abundance)") + theme_classic() +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Group),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = Group), position = position_dodge(0.3)) +
  scale_color_manual(values = c("darkred", "darkgray")) + coord_flip()
dev.off ()
