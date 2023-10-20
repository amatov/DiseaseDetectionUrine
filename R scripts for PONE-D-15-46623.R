
setwd("~/Google Drive/IBD-690/690urine/PONE") # change accordingly

############
# Figure 1 #
############

# read in data and define factors

myfile="PONE-D-15-46623_spikes.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# create dge object

library(edgeR)
dge <- DGEList(counts = raw_reads)
dge <- calcNormFactors(dge, method="TMM")

# transform counts to log frequencies and asign minimum value to -Inf

logfreq <- t(log2(t(dge$counts)/dge$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)

# load heat map color palette

HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)

# Generate column annotations

annotation = data.frame(Void=void,
                        Fraction=fraction,
                        Subject=subject,
                        Sex=sex)
Void = c("orange", "pink")
Fraction = c("yellow", "purple")
Subject = gray.colors(length(levels(subject)))
Sex=c("green", "blue")
    names(Void) = levels(void)
    names(Fraction) = levels(fraction)
    names(Subject) = levels(subject)
    names(Sex) = levels(sex)
ann_colors = list(Fraction=Fraction,
                  Sex=Sex, Void=Void,
                  Subject=Subject)

# peahtmap

library(pheatmap)
pheatmap (as.matrix(logfreq), scale="none",
color=heatmapColors,
clustering_method="ward",
cellwidth=10, cellhight=4, fontsize=8,
clustering_distance_cols = "correlation",
filename="spike heat map distance=correlation method=ward.pdf",
height=3.5, # may need to adjust hight to allow complete legend
width=20,
border_color=NA,
annotation=annotation, annotation_colors=ann_colors, annotation_legend = TRUE)

############
# Figure 2 #
############

# read in data and define factors

myfile="PONE-D-15-46623_stats.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# filter out low coverage small RNA categories

library(matrixStats)

cutoff <- raw_reads
keep <- (rowMedians(as.matrix(cutoff))>1)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_small_RNA" 

# create dge object

library(edgeR)
dge.filtered <- DGEList(counts = filtered_counts)
dge.filtered <- calcNormFactors(dge.filtered, method="TMM")

# transform counts to log frequencies and asign minimum value to -Inf

logfreq <- t(log2(t(dge.filtered$counts)/dge.filtered$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)

# load heat map color palette

HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)

# Generate column annotations

annotation = data.frame(Void=void,
                        Fraction=fraction,
                        Subject=subject,
                        Sex=sex)
Void = c("orange", "pink")
Fraction = c("yellow", "purple")
Subject = gray.colors(length(levels(subject)))
Sex=c("green", "blue")
names(Void) = levels(void)
names(Fraction) = levels(fraction)
names(Subject) = levels(subject)
names(Sex) = levels(sex)
ann_colors = list(Fraction=Fraction,
                  Sex=Sex, Void=Void,
                  Subject=Subject)

# peahtmap

library(pheatmap)
pheatmap (as.matrix(logfreq), scale="none",
          color=heatmapColors,
          cellwidth=10, cellhight=4, fontsize=8,
          filename="stats heat map.pdf",
          height=5, # may need to adjust hight to allow complete legend
          width=14,
          border_color=NA,
          annotation=annotation, annotation_colors=ann_colors, annotation_legend = TRUE)

###########
# Table B #
###########

# read in data and define factors

myfile="PONE-D-15-46623_stats.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# filter out low coverage small RNA categories

library(matrixStats)

cutoff <- raw_reads
keep <- (rowMedians(as.matrix(cutoff))>1)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_small_RNA" 

# DESeq2

library(DESeq2)
library(plyr)

countData <- as.matrix(round(filtered_counts))
colData <- data.frame (sex, void, fraction, subject)

dds <- DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design = ~ fraction*void + subject)

dds <- DESeq(dds)

# left panel - EV vs. cells (average of voids A and B)

resfraction <-  results(dds, contrast=
                          c(0,1,rep(0,20),0.5))

write.csv(resfraction, file="Table B left panel.csv")

# middle panel - void effect in cells

resvoidCells <-  results(dds, contrast=c("void","B","A"))
write.csv(resvoidCells, file="Table B middle panel.csv")

# right panel - void effect in EV

resvoidEV <-  results(dds, contrast=
                        list( c("void_B_vs_A","fractionEV.voidB")))
write.csv(resvoidEV, file="Table B right panel.csv")

############
# Figure 3 #
############

# read in data and define factors

myfile="PONE-D-15-46623_merged_mir.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# filter out low coverage miRNA

cutoff <- raw_reads
keep <- (rowSums(cutoff)>=5)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs" 

# arrage DESeqDataSet

countData <- as.matrix(round(filtered_counts))
colData <- data.frame (fraction,
  sex, subject, void)

dds <- DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design = ~ sex + fraction)

# transform the counts

vsd <- varianceStabilizingTransformation(dds)

# perform PCA and plot PC2 vs PC1

plotPCA(vsd, intgroup=c("fraction", "sex"))

############
# Figure 4 #
############

# read in data and define factors

myfile="PONE-D-15-46623_clust.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# filter out low coverage miRNA (actually keep only top 15 miRNA clusters)

cutoff <- raw_reads
keep <- rank(-rowSums(raw_reads))<=15
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs" 

# create dge object

library(edgeR)
dge.filtered <- DGEList(counts = filtered_counts)
dge.filtered <- calcNormFactors(dge.filtered, method="TMM")

# transform counts to log frequencies and asign minimum value to -Inf

logfreq <- t(log2(t(dge.filtered$counts)/dge.filtered$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)

# load heat map color palette

HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)

# Generate column annotations

annotation = data.frame(Void=void,
                        Fraction=fraction,
                        Subject=subject,
                        Sex=sex)
Void = c("orange", "pink")
Fraction = c("yellow", "purple")
Subject = gray.colors(length(levels(subject)))
Sex=c("green", "blue")
names(Void) = levels(void)
names(Fraction) = levels(fraction)
names(Subject) = levels(subject)
names(Sex) = levels(sex)
ann_colors = list(Fraction=Fraction,
                  Sex=Sex, Void=Void,
                  Subject=Subject)

# peahtmap

library(pheatmap)
pheatmap (as.matrix(logfreq), scale="none",
          color=heatmapColors,
          clustering_method="complete",
          cellwidth=10, cellhight=8, fontsize=8,
          clustering_distance_cols = "canberra",
          filename="Top 15 miRNA heat map distance=canberra method=complete.pdf",
          height=4, # may need to adjust hight to allow complete legend
          width=15,
          border_color=NA,
          annotation=annotation, annotation_colors=ann_colors,
          annotation_legend = TRUE)

#################################################
# Table F, Figure A, Table G, Table 4, Figure B #
#################################################

# read in data and define factors

myfile="PONE-D-15-46623_merged_mir.txt"

raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)

factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))

# filter out low coverage miRNA

cutoff <- raw_reads
keep <- (rowSums(cutoff)>=5)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs" 

# arrage DESeqDataSet

countData <- as.matrix(round(filtered_counts))
#ID <- factor(paste0(sex,subject))
#library(plyr)
#ID <- mapvalues(ID, from = levels(ID), to = c(1:10, 1:10))
colData <- data.frame (fraction, sex, void)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ fraction*sex + void)

dds <- DESeq(dds)

# Table F - EV vs. cells (average of men and women)

resfraction <-  results(dds, contrast=
                          c(0,1,0,0,0.5))

write.csv(resfraction, file="Table F.csv")

# Figure A

DESeq2::plotMA(resfraction, ylim=c(-15,5),
               ylab="log2 fold change",
               xlab="mean normalized count",
               main="miRNA in EVs vs. cells")
# Table G

resSex <-  results(dds, contrast = # average of cells and EVs
                   c(0,0,1,0,0.5))
write.csv(resSex, file="Table G.csv")

# Table 4

resSexCells <- results(dds, contrast = 
                     c(0,0,1,0,0))
  write.csv(resSexCells, file="Table 4 middle.csv")
resSexEV <- results(dds, contrast = 
                         c(0,0,1,0,1))
  write.csv(resSexEV, file="Table 4 right.csv")

# Figure B

rd <- data.frame(FC_cells=resSexCells$log2FoldChange, 
                 FC_EV=resSexEV$log2FoldChange)
rd$CPM <- resSex$baseMean*79000000/sum(filtered_counts)
rdp <- rd[rd$CPM>=100,]
sr <- round(with(rdp, cor.test(FC_cells, FC_EV, method="pearson"))$estimate,3)

library(ggplot2)
ggplot(rdp, aes(x=FC_cells, y=FC_EV, color=log10(CPM))) +
    geom_point() +
    ylab("log2FC in EVs (F vs. M)") +
    xlab("log2FC in cells (F vs. M)") +
    stat_smooth(method="lm", se = FALSE, color="gold") +
    coord_fixed(ratio = 1) +
    scale_colour_gradientn(name="average \nmiRNA \nabundance \n(log10CPM)",
         colours=rev(rainbow(4))) +
    geom_text(data = NULL, x = 5, y = 20, label = paste0("Pearson r = ",sr),
              color="black", size=3)

