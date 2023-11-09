logfreq <- t(log2(t(dge.filtered$counts)/dge.filtered$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)
HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)
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
library(pheatmap)
pheatmap (as.matrix(logfreq), scale="none",
color=heatmapColors,
cellwidth=10, cellhight=4, fontsize=8,
filename="stats heat map.pdf",
height=5, # may need to adjust hight to allow complete legend
width=14,
border_color=NA,
annotation=annotation, annotation_colors=ann_colors, annotation_legend = TRUE)
myfile="PONE-D-15-46623_stats.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
library(matrixStats)
cutoff <- raw_reads
keep <- (rowMedians(as.matrix(cutoff))>1)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_small_RNA"
library(DESeq2)
library(plyr)
countData <- as.matrix(round(filtered_counts))
colData <- data.frame (sex, void, fraction, subject)
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design = ~ fraction*void + subject)
dds <- DESeq(dds)
resfraction <-  results(dds, contrast=
c(0,1,rep(0,20),0.5))
write.csv(resfraction, file="Table B left panel.csv")
resvoidCells <-  results(dds, contrast=c("void","B","A"))
write.csv(resvoidCells, file="Table B middle panel.csv")
resvoidEV <-  results(dds, contrast=
list( c("void_B_vs_A","fractionEV.voidB")))
write.csv(resvoidEV, file="Table B right panel.csv")
myfile="PONE-D-15-46623_merged_mir.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
cutoff <- raw_reads
keep <- (rowSums(cutoff)>=5)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs"
countData <- as.matrix(round(filtered_counts))
colData <- data.frame (fraction,
sex, subject, void)
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design = ~ sex + fraction)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("fraction", "sex"))
myfile="PONE-D-15-46623_clust.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
cutoff <- raw_reads
keep <- rank(-rowSums(raw_reads))<=15
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs"
library(edgeR)
dge.filtered <- DGEList(counts = filtered_counts)
dge.filtered <- calcNormFactors(dge.filtered, method="TMM")
logfreq <- t(log2(t(dge.filtered$counts)/dge.filtered$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)
HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)
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
myfile="PONE-D-15-46623_merged_mir.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
cutoff <- raw_reads
keep <- (rowSums(cutoff)>=5)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs"
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
resfraction <-  results(dds, contrast=
c(0,1,0,0,0.5))
write.csv(resfraction, file="Table F.csv")
DESeq2::plotMA(resfraction, ylim=c(-15,5),
ylab="log2 fold change",
xlab="mean normalized count",
main="miRNA in EVs vs. cells")
resSex <-  results(dds, contrast = # average of cells and EVs
c(0,0,1,0,0.5))
write.csv(resSex, file="Table G.csv")
resSexCells <- results(dds, contrast =
c(0,0,1,0,0))
write.csv(resSexCells, file="Table 4 middle.csv")
resSexEV <- results(dds, contrast =
c(0,0,1,0,1))
write.csv(resSexEV, file="Table 4 right.csv")
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
plotPCA(vsd, intgroup=c("fraction", "sex"))
plotPCA(vsd, intgroup=c("fraction", "sex"))
myfile="PONE-D-15-46623_spikes.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
library(edgeR)
dge <- DGEList(counts = raw_reads)
dge <- calcNormFactors(dge, method="TMM")
logfreq <- t(log2(t(dge$counts)/dge$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)
HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)
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
myfile="PONE-D-15-46623_stats.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
library(matrixStats)
cutoff <- raw_reads
keep <- (rowMedians(as.matrix(cutoff))>1)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_small_RNA"
library(edgeR)
dge.filtered <- DGEList(counts = filtered_counts)
dge.filtered <- calcNormFactors(dge.filtered, method="TMM")
logfreq <- t(log2(t(dge.filtered$counts)/dge.filtered$samples$lib.size))
d <- density(logfreq)
logfreq[logfreq == -Inf] <- min(d$x)
HC <- read.delim(file="./heatmap_colors.txt", header=TRUE)
HC <- HC[2:4]
heatmapColors <- rgb(HC, maxColorValue=1)
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
library(pheatmap)
pheatmap (as.matrix(logfreq), scale="none",
color=heatmapColors,
cellwidth=10, cellhight=4, fontsize=8,
filename="stats heat map.pdf",
height=5, # may need to adjust hight to allow complete legend
width=14,
border_color=NA,
annotation=annotation, annotation_colors=ann_colors, annotation_legend = TRUE)
myfile="PONE-D-15-46623_stats.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
library(matrixStats)
cutoff <- raw_reads
keep <- (rowMedians(as.matrix(cutoff))>1)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_small_RNA"
library(DESeq2)
library(plyr)
countData <- as.matrix(round(filtered_counts))
colData <- data.frame (sex, void, fraction, subject)
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design = ~ fraction*void + subject)
dds <- DESeq(dds)
resfraction <-  results(dds, contrast=
c(0,1,rep(0,20),0.5))
write.csv(resfraction, file="Table B left panel.csv")
resvoidCells <-  results(dds, contrast=c("void","B","A"))
write.csv(resvoidCells, file="Table B middle panel.csv")
resvoidEV <-  results(dds, contrast=
list( c("void_B_vs_A","fractionEV.voidB")))
write.csv(resvoidEV, file="Table B right panel.csv")
myfile="PONE-D-15-46623_merged_mir.txt"
raw_reads <- read.table(myfile, header=TRUE, skip=5, row.names=1)
raw_reads <- round(raw_reads)
factors <- read.table(myfile, header=TRUE, nrows=4, row.names=1, stringsAsFactors=FALSE)
fraction <- factor(factors[1,], levels=c("cells","EV"))
void <- factor(factors[2,], levels=c("A", "B"))
subject <- factor(factors[3,], levels=1:20)
sex <- factor(factors[4,], levels=c("M", "F"))
cutoff <- raw_reads
keep <- (rowSums(cutoff)>=5)>0
filtered_counts <- raw_reads[keep,]
other_mirs <- as.vector(colSums(cutoff)-colSums(filtered_counts))
filtered_counts <- rbind(filtered_counts, other_mirs)
rownames(filtered_counts)[nrow(filtered_counts)] <- "other_mirs"
countData <- as.matrix(round(filtered_counts))
colData <- data.frame (fraction,
sex, subject, void)
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design = ~ sex + fraction)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("fraction", "sex"))
myfile="PONE-D-15-46623_clust.txt"
install.packages("PTak")
install.packages('PTak')
install.packages(PTak)
library(PTAk)
install.packages("PTAk")
library(PTAk)
install.packages("tensor")
library(PTAk)
install.packages("tensor")
install.packages("tensor")
library(PTAk)
install.packages("tensor")
install.packages("tensor")
library("PTAk)
library("PTAk")
library("PTAk")
help("PTAk")
PTAk(X,nbPT=2,nbPT2=1,minpct=0.1,
smoothing=FALSE,
smoo=list(NA),
verbose=getOption("verbose"),file=NULL,
modesnam=NULL,addedcomment="", ...)
PTAk(X,nbPT=2,nbPT2=1,minpct=0.1,
smoothing=FALSE,
smoo=list(NA),
verbose=getOption("verbose"),file=NULL,
modesnam=NULL,addedcomment="")
data(iris)
iris2 <- as.matrix(iris[,1:4])
dimnames(iris2)[[1]] <- as.character(iris[,5])
D2 <- CauRuimet(iris2,ker=1,withingroup=TRUE)
D2 <- Powmat(D2,(-1))
iris2 <- sweep(iris2,2,apply(iris2,2,mean))
res <- SVDgen(iris2,D2=D2,D1=1)
plot(res,nb1=1,nb2=2,cex=0.5)
summary(res,testvar=0)
library(PTAk)
myt
demo.SVDgen(snr=3,openX11s=T)
data(Zone_climTUN)
library(maptools)
library(RColorBrewer)
Yl=brewer.pal(11,"PuOr")
plot(Zone_climTUN,ol=NA,auxvar=Zone_climTUN$att.data$PREC_OCTO)
Zone_climTUN
don <- array(1:360,c(5,4,6,3))
don <- don + rnorm(360,1,2)
dimnames(don) <- list(paste("s",1:5,sep=""),paste("T",1:4,sep=""),
paste("t",1:6,sep=""),c("young","normal","old"))
# hypothetic data on learning curve at different age and period of year
ones <-list(list(v=rep(1,5)),list(v=rep(1,4)),list(v=rep(1,6)),list(v=rep(1,3)))
don <- PROJOT(don,ones)
don.sol <- PTAk(don,nbPT=1,nbPT2=2,minpct=0.01,
verbose=TRUE,
modesnam=c("Subjects","Trimester","Time","Age"),
addedcomment="centered on each mode")
don.sol[[1]] # mode Subjects results and components
don.sol[[2]] # mode Trimester results and components
don.sol[[3]] # mode Time results and components
don.sol[[4]] # mode Age results and  components with additional information on the call
summary(don.sol,testvar=2)
plot(don.sol,mod=c(1,2,3,4),nb1=1,nb2=NULL,
xlab="Subjects/Trimester/Time/Age",main="Best rank-one approx" )
plot(don.sol,mod=c(1,2,3,4),nb1=4,nb2=NULL,
xlab="Subjects/Trimester/Time/Age",main="Associated to Subject vs1111")
don <- array(1:360,c(5,4,6,3))
don
healthy_miRmedNor_full947 <- read.delim("C:/Users/dsa/Desktop/lung_oncomirs_txt/healthy_miRmedNor_full947.txt", header=FALSE)
View(healthy_miRmedNor_full947)
healthy_miRmedNor_full947 <- read.delim("C:/Users/dsa/Desktop/lung_oncomirs_txt/healthy_miRmedNor_full947.txt", header=TRUE)
View(healthy_miRmedNor_full947)
data("C:/Users/dsa/Desktop/lung_oncomirs_txt/healthy_miRmedNor_full947.txt")
healthy_miRmedNor_full947.txt
healthy_miRmedNor_full947
c(healthy_miRmedNor_full947)
don <- array(1:919,c(healthy_miRmedNor_full947))
don <- c(healthy_miRmedNor_full947)
don <- PROJOT(don)
don <- array(1:360,c(5,4,6,3))
don
don <- don + rnorm(360,1,2)
don
don <- c(healthy_miRmedNor_full947)
don
don(CRJ9)
don(CJR10_S3)
don.sol <- PTAk(don,nbPT=1,nbPT2=2,minpct=0.01,
verbose=TRUE,
modesnam=c("Subjects","Trimester","Time","Age"),
addedcomment="centered on each mode")
don.sol <- PTAk(don
don.sol <- PTAk(don)
x<-don
X<-don
APSOLU3(X,solu,pt3=NULL,nbPT2=1,
smoothing=FALSE,smoo=list(NA),
verbose=getOption("verbose"),file=NULL)
X
healthy_miRmedNor_full947
dim(healthy_miRmedNor_full947)
dim(c)
healthy_miRmedNor_full947[1,]
healthy_miRmedNor_full947[,1]
healthy_miRmedNor_full947[1,]
healthy_miRmedNor_full947[,2]
healthy_miRmedNor_full947[,3]
healthy_miRmedNor_full947[,47]
healthy_miRmedNor_full947[,46]
healthy_miRmedNor_full947[,1:16]
t1<-healthy_miRmedNor_full947[,1:16]
t2<-healthy_miRmedNor_full947[,17:31]
t3<-healthy_miRmedNor_full947[,12:46]
t1<-healthy_miRmedNor_full947[,2:16]
mir<-healthy_miRmedNor_full947[,1]
m<-healthy_miRmedNor_full947[,1]
m
t1
dim(t1)
dim(t2)
dim(t4)
dim(t3)
dim(m)
t3<-healthy_miRmedNor_full947[,31:46]
dim(t3)
t3<-healthy_miRmedNor_full947[,32:46]
dim(t3)
m<-healthy_miRmedNor_full947[,1]
m
dim(m)
p<-healthy_miRmedNor_full947[1,]
p
p<-healthy_miRmedNor_full947[0,]
p
data(iris)
iris.data <- iris[,1:4]
iris.data
agnes.mod <- agnes(iris.data) # create cluster tree
install.packages('seqinr')
library('seqinr')
install.packages('seqinr')
library('seqinr')
install.packages("seqinr")
install.packages('seqinr')
library('seqinr')
devtools::install_github("rstudio/keras")
install.packages('devtools')
devtools::install_github("rstudio/keras")
install.packages("devtools")
devtools::install_github("rstudio/keras")
install.packages("devtools")
install.packages("devtools")
install.packages("devtools")
install.packages("rlang")
install.packages("rlang")
install.packages("devtools")
install.packages("devtools")
library(keras)
devtools::install_github("rstudio/keras")
library(rlang)
update.packages()
devtools::install_github("rstudio/keras")
install.packages("rlang")
install.packages("rlang")
library(rlang)
install.packages("rlang")
install.packages("rlang")
install.packages("rlang")
install.packages("rlang")
install.packages("rlang")
install.packages("rlang")
remove.packages("rlang")
install.packages("rlang")
install.packages("rlang")
remove.packages("rlang")

