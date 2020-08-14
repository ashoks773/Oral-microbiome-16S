library(phyloseq)
library (ape)
library (ggplot2)

#Figure1a----------------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in OTU table
otu_table_in <- read.csv("feature_table.txt", sep = "\t", row.names = 1)
total_asv <- colSums(otu_table_in)
range (total_asv)
sd (total_asv)
mean (total_asv)

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

# Finally merge to create Phyloseq object!
ps <- phyloseq(OTU, TAX, META, phy_tree)

Colors <- c("darkred", "darkgray")
p <- plot_richness(ps, "Group", measures = c("Observed", "Chao1", "Shannon"))
p <- p + geom_boxplot(data = p$data, aes(x= Group, y = value, color = Group), alpha = 0.1) + scale_color_manual(values=Colors) + labs(x="",y="Alpha Diversity Measure") + theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("phyloseq_analysis-richness_estimates.pdf", p, width = 6, height = 3)
jpeg("Phyloseq_Richness.jpg", height = 3, width = 6, units = 'in', res = 600)
p
dev.off ()
erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
erich_group <- cbind(erich, metadata)
write.table (erich_group, file = "Richness_total.txt", sep = "\t")
wilcox.test(erich_group$Observed ~ erich_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 162, p-value = 0.002299
wilcox.test(erich_group$Chao1 ~ erich_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 162, p-value = 0.002299
wilcox.test(erich_group$Shannon ~ erich_group$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 281, p-value = 0.4258

#----- For reviewers
erich_group <- read.csv(file="Richness_total.txt", sep = "\t", row.names = 1, header = T)
erich_group_withOutNA <- subset(erich_group, Alcohol=="NA")


#Figure 1b and 1c----------------------------------- Beta diversity analysis
#A) Create Phyloseq object on filtered OTU
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

Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )


#A)-------- UniFrac PCoA
UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
UniFrac_distances <- as.matrix(UniFrac_distances)
UniFrac_dist_column <- melt(UniFrac_distances)
write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")

Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("UniFrac_PCoA.png", height = 5, width = 5, units = 'in', res = 600)
plot(Uni_pcoa$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[Bushman2@sam_data$Group], pch=c(21,22)[Bushman2@sam_data$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=Bushman2@sam_data$Group, cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Group, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topleft", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomleft", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis(UniFrac_distances ~ Bushman2@sam_data$Group)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Bushman2@sam_data$Group  1    0.1672 0.167247  2.0487 0.04013   0.05 *
#----- For reviewers
#adonis(UniFrac_distances ~ Group*Alcohol, data=metadata)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Group          1    0.1672 0.167247 2.05629 0.04013  0.057 .
#Alcohol        4    0.4177 0.104417 1.28380 0.10022  0.161  
#Group:Alcohol  3    0.1664 0.055465 0.68193 0.03993  0.845


Uni_PCoA_MATRIX <- Uni_pcoa$vectors[,1:2]
Uni_PCoA_MATRIX <- data.frame(Uni_PCoA_MATRIX)
Uni_PCoA_MATRIX_New <- cbind(Uni_PCoA_MATRIX, metadata)

pcoa1 <- Uni_PCoA_MATRIX_New[,c(1,8)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group")
jpeg("Uni_PCoA1_Distances.jpg", height = 1, width = 3, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa1_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCoA1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + coord_flip() + theme(legend.position='none')
dev.off ()
wilcox.test(pcoa1$Axis.1 ~ pcoa1$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 302, p-value = 0.6875

pcoa2 <- Uni_PCoA_MATRIX_New[,c(2,8)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group")
jpeg("Uni_PCoA2_Distances.jpg", height = 3, width = 1, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.2, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold"))
ggplot(data = pcoa2_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("UniFrac Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
wilcox.test(pcoa2$Axis.2 ~ pcoa2$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 514, p-value = 0.0002166

Colors <- c("darkred", "darkgray")
UniFrac_dist_groups <- read.csv("UniFrac_distances_Filtered.txt.txt", sep = "\t", header=T)
UniFrac_dist_groups_Melted <- melt(UniFrac_dist_groups, id.vars = "Groups")
jpeg("UniFrac_Distances.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(data = UniFrac_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="UniFrac distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
wilcox.test(UniFrac_dist_groups$Distances ~ UniFrac_dist_groups$Groups, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 264020, p-value = 1.262e-15

#B)------------- Unweighted UniFrac PCoA
Unweigh_UniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
Unweigh_UniFrac_distances <- as.matrix(Unweigh_UniFrac_distances)
Unweigh_UniFrac_dist_column <- melt(Unweigh_UniFrac_distances)
write.table (Unweigh_UniFrac_dist_column, file = "Unweigh_UniFrac_distnaces.txt", sep = "\t")

Unweigh_Uni_pcoa <- pcoa(Unweigh_UniFrac_distances)
Unweigh_Uni_pcoa$values[1:2,]
mds.var.per = round(Unweigh_Uni_pcoa$values$Eigenvalues/sum(Unweigh_Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("Unweigh_UniFrac_PCoA.png", height = 5, width = 5, units = 'in', res = 600)
plot(Unweigh_Uni_pcoa$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[Bushman2@sam_data$Group], pch=c(21,22)[Bushman2@sam_data$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=Bushman2@sam_data$Group, cex=0.3, font=1, pos=1)
ordiellipse(Unweigh_Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Group, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Unweigh_Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topright", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomright", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis(Unweigh_UniFrac_distances ~ Bushman2@sam_data$Group)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Bushman2@sam_data$Group  1    0.4972 0.49716  3.0528 0.05865  0.001 ***

#----- For reviewers
adonis(Unweigh_UniFrac_distances ~ Group*Alcohol, data=metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Group          1    0.4972 0.49716 3.12967 0.05865  0.001 ***
#  Alcohol        4    0.8359 0.20897 1.31551 0.09861  0.054 .  
#Group:Alcohol  3    0.4722 0.15740 0.99083 0.05570  0.484

Unweigh_Uni_PCoA_MATRIX <- Unweigh_Uni_pcoa$vectors[,1:2]
Unweigh_Uni_PCoA_MATRIX <- data.frame(Unweigh_Uni_PCoA_MATRIX)
Unweigh_Uni_PCoA_MATRIX_New <- cbind(Unweigh_Uni_PCoA_MATRIX, metadata)

pcoa1 <- Unweigh_Uni_PCoA_MATRIX_New[,c(1,8)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group")
jpeg("Unweigh_Uni_PCoA1_Distances.jpg", height = 1, width = 3, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa1_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("Unweighted UniFrac Distances") + labs(x="",y="PCoA1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + coord_flip() + theme(legend.position='none')
dev.off ()
wilcox.test(pcoa1$Axis.1 ~ pcoa1$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 493, p-value = 0.001129

pcoa2 <- Unweigh_Uni_PCoA_MATRIX_New[,c(2,8)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group")
jpeg("Unweigh_Uni_PCoA2_Distances.jpg", height = 3, width = 1, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.2, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold"))
ggplot(data = pcoa2_Melted, aes(x=Group, y=value)) + geom_boxplot() + ggtitle("Unweighted UniFrac Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
wilcox.test(pcoa2$Axis.2 ~ pcoa2$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 397, p-value = 0.1727

Colors <- c("darkred", "darkgray")
Unweigh_UniFrac_dist_groups <- read.csv("Unweigh_UniFrac_distnaces_Filtered.txt.txt", sep = "\t", header=T)
Unweigh_UniFrac_dist_groups_Melted <- melt(Unweigh_UniFrac_dist_groups, id.vars = "Groups")
jpeg("Unweigh_UniFrac_Distances.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(data = Unweigh_UniFrac_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="Unweighted UniFrac distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
wilcox.test(Unweigh_UniFrac_dist_groups$Distances ~ Unweigh_UniFrac_dist_groups$Groups, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 311040, p-value < 2.2e-16

#C)--------------- For Bray-Curtis-----------
# Read in OTU table
otu_table_in <- read.csv("feature_table.txt", sep = "\t", row.names = 1)
#-- Additional steps to filter OTUs
otu_table_in_t <- data.frame(t(otu_table_in))
otu.sums <- colSums(otu_table_in_t)
otu_table_filtered <- otu_table_in_t[ , which(otu.sums > 5)]
otu_table_filtered_t <- data.frame(t(otu_table_filtered)) #297 OTUs have sume less then 5
Otu_rel_abundances <- otu_table_filtered_t/colSums(otu_table_filtered_t)[col(otu_table_filtered_t)]
write.table(Otu_rel_abundances, file="OTUs_rel_abundances.txt", sep="\t")
Otu_rel_abundances <- data.frame(t(Otu_rel_abundances))

Bray_pcoa <-pcoa(vegdist(Otu_rel_abundances, "bray"))
Bray_pcoa$values[1:2,]
mds.var.per = round(Bray_pcoa$values$Eigenvalues/sum(Bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- Bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, metadata)

pc <- c(1,2)
jpeg("Bray-PCoA.png", height = 5, width = 5, units = 'in', res = 600)
plot(Bray_pcoa$vectors[,1:2], bg=c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255))[metadata$Group], pch=c(21,22)[metadata$Group], cex=1.3, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Bray_pcoa$vectors[,1:2], labels=metadata$Group, cex=0.3, font=1, pos=1)
ordiellipse(Bray_pcoa$vectors[,1:2], metadata$Sex, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)))
ordispider(Bray_pcoa$vectors[,1:2], metadata$Group, lty=3, spider ="centroid", lwd=1, col="black")
legend("topright", legend = c("Cases", "Controls"), col = c(rgb(154, 0, 0, maxColorValue = 255),rgb(169, 169, 169, maxColorValue = 255)),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomright", legend = c("Cases", "Controls"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
Bray_distances <-vegdist(Otu_rel_abundances, "bray")
adonis(Bray_distances ~ metadata$Group)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#metadata$Group  1     0.489 0.48903  1.9432 0.03814   0.03 *

#----- For reviewers
adonis(Bray_distances ~ Group*Alcohol, metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Group          1    0.4890 0.48903 1.97792 0.03814  0.030 *
#  Alcohol        4    1.2872 0.32181 1.30160 0.10041  0.095 .
#Group:Alcohol  3    0.6600 0.22000 0.88983 0.05148  0.655  

#-------
pcoa1 <- Bray_PCoA_MATRIX_New[,c(1,8)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group")
jpeg("Bray_PCoA1_Distances.jpg", height = 1, width = 3, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.1, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + coord_flip()
ggplot(data = pcoa1_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + coord_flip() + theme(legend.position='none')
dev.off ()
wilcox.test(pcoa1$Axis.1 ~ pcoa1$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 323, p-value = 0.9925

pcoa2 <- Bray_PCoA_MATRIX_New[,c(2,8)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group")
jpeg("Bray_PCoA2_Distances.jpg", height = 3, width = 1, units = 'in', res = 600)
#ggplot(Bray_PCoA_MATRIX_New, aes(x=Groups, y=Axis.2, color=Groups), alpha = 0.1) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2")  + theme_classic() + scale_color_manual(values=Colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold"))
ggplot(data = pcoa2_Melted, aes(x=Group, y=value, fill=Group)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
wilcox.test(pcoa2$Axis.2 ~ pcoa2$Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
#W = 514, p-value = 0.0002166

Bray_dist <- as.matrix(Bray_distances)
Bray_dist_column <- melt(Bray_dist)
write.table (Bray_dist_column, file="Bray-distances.txt", sep="\t")

Colors <- c("darkred", "darkgray")
Bray_dist_groups <- read.csv("Bray-distances_filtered.txt.txt", sep = "\t", header=T)
Bray_dist_groups_Melted <- melt(Bray_dist_groups, id.vars = "Groups")
jpeg("Bray_Distances.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(data = Bray_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="Bray-Curtis distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
wilcox.test(Bray_dist_groups$Distances ~ Bray_dist_groups$Groups, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)


#Figure 2a and 2b  ---------- Taxonomic Analysis---------------------
#qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 7 --o-collapsed-table table-l7.qza
# biom convert -i feature-table.biom -o feature-table-L7.txt --to-tsv

taxa <- read.csv (file="feature-table-L7_filtered.txt", row.names = 1, header=T, sep="\t")
taxa_rel_abundances <- taxa/colSums(taxa)[col(taxa)]
write.table (taxa_rel_abundances, file="taxa_rel_abundances.txt", sep="\t")
taxa_rel_abundances <- data.frame(t(taxa_rel_abundances))
Group <- metadata$Group
taxa_rel_abundances_Groups <- cbind(taxa_rel_abundances, Group)

#---- Machine Learning RandomForest
library (randomForest)
RF1 <- randomForest(Group ~., data=taxa_rel_abundances_Groups, ntree = 500, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
#ntree      OOB      1      2
#ntree      OOB      1      2
#100:  31.37% 25.93% 37.50%
#200:  27.45% 25.93% 29.17%
#300:  27.45% 25.93% 29.17%
#400:  29.41% 29.63% 29.17%
#500:  27.45% 29.63% 25.00%
jpeg("RFImp.jpg", height = 10, width = 15, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
imp <- importance(RF1, type=1)
write.table(imp, file="rf_importance.txt", sep="\t")

OOB.votes <- RF1$votes
write.table (OOB.votes, file = "OOB_pred", sep = "\t")
predictions <- read.csv (file = "OOB_pred", sep = "\t")
roc.multi <- multiclass.roc(taxa_rel_abundances_Groups$Group, predictions$Controls)
auc(roc.multi)
rs <- roc.multi[['rocs']]
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=c("green")))

jpeg("AUC.png", height = 4, width = 4, units = 'in', res = 600)
plot.roc(rs[[1]], print.auc = TRUE, grid=c(0.1, 0.2), lwd=3)
dev.off ()

#---- Reviewers Comments to report the performance on seperate test set.
# Split into Train and Validation sets
train <- sample(nrow(taxa_rel_abundances_Groups), 0.7*nrow(taxa_rel_abundances_Groups), replace = FALSE)
TrainSet <- taxa_rel_abundances_Groups[train,]
ValidSet <- taxa_rel_abundances_Groups[-train,]
summary(TrainSet)
summary(ValidSet)
model1 <- randomForest(Group ~ ., data = TrainSet, importance = TRUE)
#Confusion matrix:
#Cases Controls class.error
#Cases       13        4   0.2352941
#Controls     5       13   0.2777778

predTest<- predict(model1, ValidSet, type = "class")
table(predTest, ValidSet$Group) 

predTest<- predict(model1, ValidSet, type = "prob")

library(ROCR)
perf = prediction(predTest[,2], ValidSet$Group)
# 1. Area under curve
auc = performance(perf, "auc")
auc=0.75
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
jpeg("AUC_onTestSet.png", height = 5, width = 4, units = 'in', res = 600)
plot(pred3,main="ROC Curve for Random Forest",print.auc = TRUE, grid=c(0.1, 0.2), lwd=3)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
dev.off ()

#-- Wilcoxon-rank sum test
source('~/Work/Collaborative_Projects/Irina_Repeated/Wilcox.R')
Wilcoxon_pairwise(taxa_rel_abundances_Groups)

#-- Indicator species Labdsv
library (labdsv)
#Genus2 <- Otu_table_Clusters[ , which(!apply(Otu_table_Clusters==0,2,all))]
#Genus2[is.na(Genus2)] <- 0
iva <- indval(taxa_rel_abundances_Groups[,1:396], taxa_rel_abundances_Groups$Group)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(taxa_rel_abundances_Groups[,1:396]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="indvalsummary.txt", sep = "\t")

#-- SelectedTaxa
#-- Sel_taxa.txt file were prepared by selecting all taxa (p<=0.05,log2FC >=0.2)
#grep -f Sel_taxa.txt taxa_rel_abundances.txt

Selected_taxa <- read.csv (file="Sel_taxa_relative_abun_Filtered_Formatted.txt", row.names = 1, header=T, sep="\t")
Selected_taxa <- data.frame(t(Selected_taxa))
Selected_Taxa_Groups <- cbind(Selected_taxa, metadata)
write.table (Selected_Taxa_Groups, file="Selected_Taxa_rel_abund_Groups.txt", sep="\t")

Group <- metadata$Group
Selected_Taxa_group <- cbind(sqrt(Selected_taxa), Group)
taxa_group_melted <- melt(Selected_Taxa_group, id.vars = "Group")
Colors <- c("darkred", "darkgray")
jpeg("Check_RF_wilcox_sel_taxa.jpg", height = 6, width = 8, units = 'in', res = 600)
p <- ggplot(data = taxa_group_melted, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Group))
p + facet_wrap( ~ variable, scales="free") +scale_color_manual(values=Colors) + scale_fill_manual(values=Colors)
dev.off ()
jpeg("RF_wilcox_sel_taxa_Formatted.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(data = taxa_group_melted, aes(x=variable, y=value, fill=Group)) + geom_boxplot() + ggtitle("") + labs(x="",y="sqrt(Relative abundance)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 5, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='right') + theme(axis.text.x = element_text(angle = 50, hjust = 1))
dev.off ()


#----2c-----Spearman correlations
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)

Selected_taxa <- read.csv (file="../Sel_taxa_relative_abun_Filtered.txt", row.names = 1, header=T, sep="\t")
Selected_taxa <- data.frame(t(Selected_taxa))

Correlation_Pairwise_2(Selected_taxa, Selected_taxa)

SpearCorrTT <- read.csv (file="SpearCorr_Taxa_Taxa.txt", header=T, sep="\t")

graph_taxa_path_Spear_corsTT <- SpearCorrTT %>%
  filter(abs(r) > .4) %>%
  filter(abs(p) < 0.05) %>%
  graph_from_data_frame(directed = FALSE)

#--- This figure is replotted using Cytoscape

jpeg("Taxa_Taxa_Corr_Spear.jpg", height = 7, width = 7, units = 'in', res = 600)
ggraph(graph_taxa_path_Spear_corsTT) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("darkred", "forestgreen")) +
  scale_edge_width(range = c(0.4, 1)) +
  geom_node_point(color = "darkgray", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between disriminating taxa")
dev.off ()


#Figure3a and 3b---------- Diversity correlations with HPB data
erich_group <- read.csv(file="Richness_total.txt", sep="\t", header=T, row.names = 1)
jpeg("HPB.jpg", height = 4, width = 3, units = 'in', res = 600)
boxplot(log10(erich_group$HPB..pmol.mg.DNA.) ~ erich_group$Group, col=c("darkred", "darkgray"))
dev.off ()
wilcox.test(log10(erich_group$HPB..pmol.mg.DNA.) ~ erich_group$Group)

jpeg("Observed_otus_HPB.jpg", height = 4, width = 4, units = 'in', res = 600)
plot (log10(erich_group$HPB..pmol.mg.DNA.), log10(erich_group$Observed), bg=c("darkred","darkgray")[ erich_group$Group], pch=c(21,22)[erich_group$Group])
model<-glm (log10(erich_group$Observed) ~ log10(erich_group$HPB..pmol.mg.DNA.))
abline(model)
dev.off ()
cor.test(x=log10(erich_group$Observed), y=log10(erich_group$HPB..pmol.mg.DNA.), method="spearman")
#rho=-0.3793314, p=0.006047
#-- For Revision
new <- corr.test(x=log10(erich_group$Observed), y=log10(erich_group$HPB..pmol.mg.DNA.), method="spearman", adjust = "fdr", alpha=.05)

#jpeg("Shannon_HPB.jpg", height = 4, width = 4, units = 'in', res = 600)
#plot (log10(erich_group$HPB..pmol.mg.DNA.), log10(erich_group$Shannon), bg=c("darkred","darkgray")[ erich_group$Group], pch=c(21,22)[erich_group$Group])
#model1<-glm (log10(erich_group$Shannon) ~ log10(erich_group$HPB..pmol.mg.DNA.))
#abline(model1)
#dev.off ()
#cor.test(x=log10(erich_group$Shannon), y=log10(erich_group$HPB..pmol.mg.DNA.), method="spearman")
#rho=--0.1081069 , p=0.4502

#--Figure3c---------- Taxonmoic correlations with HPB data
Taxa_Cases <- read.csv(file="Cases_taxa_rel_abundances.txt", sep="\t", header=T, row.names = 1)
Taxa_Cases <- data.frame(t(Taxa_Cases))
Taxa_Cases_Group <- cbind(Taxa_Cases, metadata)

jpeg("Cases_Taxa-HPB-New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot (log10(Taxa_Cases_Group$HPB..pmol.mg.DNA.), log10(Taxa_Cases_Group$Cumulative_abundance+0.0001), bg=c("darkred","darkgray")[ Taxa_Cases_Group$Group], pch=c(21,22)[Taxa_Cases_Group$Group])
model <- lm (log10(Taxa_Cases_Group$Cumulative_abundance+.0001) ~ log10(Taxa_Cases_Group$HPB..pmol.mg.DNA.))
abline(model)
dev.off ()
cor.test(x=log10(Taxa_Cases_Group$Cumulative_abundance+0.0001), y=log10(Taxa_Cases_Group$HPB..pmol.mg.DNA.), method="spearman")
#--S = 14997, p-value = 0.02146, rho = 0.3214068 
#--- New one for revision
new1 <- corr.test(x=log10(Taxa_Cases_Group$Cumulative_abundance+0.0001), y=log10(Taxa_Cases_Group$HPB..pmol.mg.DNA.), method="spearman", adjust = "fdr", alpha=.05)


Taxa_Cases_Group_picked <- log10(Taxa_Cases_Group[c(4,13)])
jpeg("Cases_Taxa-HPB.jpg", height = 4, width = 4, units = 'in', res = 600)
ggscatter(Taxa_Cases_Group_picked, x = "HPB..pmol.mg.DNA.", y = "Cumulative_abundance", 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "darkred"),
          xlab = "log10(HPB, pmol/mg DNA))", ylab = "log10(cumulative taxa abundance)")
dev.off()


Taxa_Control <- read.csv(file="Control_taxa_rel_abundances.txt", sep="\t", header=T, row.names = 1)
Taxa_Control <- data.frame(t(Taxa_Control))
Taxa_Control_Group <- cbind(Taxa_Control, metadata)

jpeg("Control_Taxa-HPB-New.jpg", height = 4, width = 4, units = 'in', res = 600)
plot (log10(Taxa_Control_Group$HPB..pmol.mg.DNA.), log10(Taxa_Control_Group$Cumulative_abundance+0.0001), bg=c("darkred","darkgray")[ Taxa_Cases_Group$Group], pch=c(21,22)[Taxa_Cases_Group$Group])
model <- lm (log10(Taxa_Control_Group$Cumulative_abundance+.0001) ~ log10(Taxa_Control_Group$HPB..pmol.mg.DNA.))
abline(model)
dev.off ()
cor.test(x=log10(Taxa_Control_Group$Cumulative_abundance+0.0001), y=log10(Taxa_Control_Group$HPB..pmol.mg.DNA.), method="spearman")
#--S = 27480, p-value = 0.08519, rho = -0.2434322 
#For revision 
new2 <- corr.test(x=log10(Taxa_Control_Group$Cumulative_abundance+0.0001), y=log10(Taxa_Control_Group$HPB..pmol.mg.DNA.), method="spearman", adjust = "fdr", alpha=.05)

Taxa_Control_Group_picked <- log10(Taxa_Control_Group[c(9,18)])
jpeg("Control_Taxa-HPB.jpg", height = 4, width = 4, units = 'in', res = 600)
ggscatter(Taxa_Control_Group_picked, x = "HPB..pmol.mg.DNA.", y = "Cumulative_abundance", 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "darkgray"),
          xlab = "log10(HPB, pmol/mg DNA))", ylab = "log10(cumulative taxa abundance)")
dev.off()


#---- Figure3d ---- Each taxa and HPB correlations
HPB_actual <- metadata$HPB..pmol.mg.DNA.
HPB_actual <- data.frame(HPB_actual)
Correlation_Pairwise_2(Selected_taxa, HPB_actual)

SpearCorrTaxa_HPB <- read.csv (file="Taxa-HPB_Corr.txt", header=T, sep="\t")

graph_taxa_path_Spear_cors <- SpearCorrTaxa_HPB %>%
  filter(abs(r) > 0.2) %>%
  filter(abs(p) < 0.06) %>%
  graph_from_data_frame(directed = FALSE)

#--- This figure is replotted using Cytoscape

jpeg("Taxa_HPB_Corr_Spear.jpg", height = 6, width = 6, units = 'in', res = 600)
ggraph(graph_taxa_path_Spear_cors) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("darkred", "forestgreen")) +
  scale_edge_width(range = c(0.8, 1)) +
  geom_node_point(color = "darkgray", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between taxa and HPB")
dev.off ()


#------ Figure4 ----- For correlations between taxa and predicted pathways
setwd("~/Work/Collaborative_Projects/Irina_Repeated/Picrust_pathway")

path <- read.csv("pathway_abundance.txt", sep = "\t", row.names = 1)
path_proportions <- path/colSums(path)[col(path)]
path_proportions_t <- data.frame(t(path_proportions))
#path_proportions_to_save <- data.frame(t(path_proportions_t))
#write.table (path_proportions_to_save, "path_proportions.txt", sep="\t")

metadata <- read.table("../Metadata.txt", sep="\t", row.names = 1, header=T)

#--------------------------------Pathway selection
Group <- metadata$Group
path_proportions_Groups <- cbind(path_proportions_t, Group)

#1---- Machine Learning RandomForest
library (randomForest)
RF1 <- randomForest(Group ~., data=path_proportions_Groups, ntree = 500, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
jpeg("RFImp.jpg", height = 10, width = 15, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
imp <- importance(RF1, type=1)
write.table(imp, file="rf_importance.txt", sep="\t")

#2-- Wilcoxon-rank sum test
source('~/Work/Collaborative_Projects/Irina_Repeated/Picrust_pathway/Wilcox.R')
Wilcoxon_pairwise(path_proportions_Groups)

#3-- Indicator species Labdsv (Not Used in this Case)
#library (labdsv)
#path_filtered <- path_proportions_Groups[ , which(!apply(path_proportions_Groups==0,2,all))]
#path_filtered[is.na(path_filtered)] <- 0
#iva <- indval(path_filtered[,1:482], path_filtered$Group)
#gr <- iva$maxcls[iva$pval<=0.05]
#iv <- iva$indcls[iva$pval<=0.05]
#pv <- iva$pval[iva$pval<=0.05]
#fr <- apply(path_filtered[,1:482]>0, 2, sum)[iva$pval<=0.05]
#indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
#indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
#write.table (indvalsummary, file="indvalsummary.txt", sep = "\t")

#-- SelectedPath
# grep -w -f Sel_path path_proportions.txt > path_proportions_Filtered.txt (ORNARGDEG.PWY missing then added manually)

Selected_pathways <- read.csv (file="path_proportions_Filtered_Names.txt", row.names = 1, header=T, sep="\t")
Selected_path <- data.frame(t(Selected_pathways))


#----- Spearman correlations
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)

Selected_taxa <- read.csv (file="../Sel_taxa_relative_abun_Filtered.txt", row.names = 1, header=T, sep="\t")
Selected_taxa <- data.frame(t(Selected_taxa))

source('~/Work/Collaborative_Projects/Irina_Repeated/Picrust_pathway/Corr_spearman.R')

Correlation_Pairwise_2(Selected_taxa, Selected_path)

SpearCorr <- read.csv (file="SpearCorr_Microbe_and_HostGenes.txt", header=T, sep="\t")

graph_taxa_path_Spear_cors <- SpearCorr %>%
  filter(abs(r) > .5) %>%
  filter(abs(p) < 0.05) %>%
  graph_from_data_frame(directed = FALSE)

#--- This figure is replotted using Cytoscape

jpeg("Taxa_Path_Corr_Spear.jpg", height = 12, width = 12, units = 'in', res = 600)
ggraph(graph_taxa_path_Spear_cors) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("darkred", "forestgreen")) +
  scale_edge_width(range = c(0.5, 1)) +
  geom_node_point(color = "darkgray", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between taxa and pathways")
dev.off ()

#------------------------------------------------------------ END-----------------------------------------------------------------------#