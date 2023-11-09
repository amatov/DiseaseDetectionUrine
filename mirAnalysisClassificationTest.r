bagsvm <- classify(data = data.trainS4, method = "bagsvm", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
bagsvm
trained(bagsvm)
confusionMat(bagsvm)
library(DESeq2)
library(DESeq2)
source('C:/Users/dsa/Desktop/micro_RNA_analyses_1_12_2017.R', echo=TRUE)
## Load the necessary libraries
library(DESeq2)
library(edgeR)
library(DESeq2)
install.packages("DESeq2")
source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
library(DESeq2)
library(edgeR)
temp <- read.table("C:/Users/dsa/Desktop/Data_1_12_2017_2019.csv", sep = ";", header = TRUE, row.names = 1, nrows = 687)
RNA <- temp[-1, seq(from = 1, to = 76, by = 5)]
colnames(RNA)
group <- factor(c(rep(1,8),rep(2,8)))
contrasts(group) <- c(-1,1)
design <- model.matrix( ~ group)
y <- DGEList(counts = RNA, group = group)
## show on the screen what is in the list object y.
y
et <- exactTest(y)
topTags(et, n = 15, adjust.method = "hochberg")
y
et <- exactTest(y)
y_norm <- calcNormFactors(y, method = "TMM")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(d,design)
result <- glmQLFTest(fit)
topTags(result, n = 30)
## Same, but now with robust = TRUE
## robust = TRUE: whether to estimate the prior QL dispersion distribution robustly.
fit <- glmQLFit(d, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results, n = 30)
## Same, but now with abundance.trend = FALSE
fit <- glmQLFit(d, design, abundance.trend=FALSE)
results <- glmQLFTest(fit)
topTags(results, n = 30)
y_norm <- calcNormFactors(y, method = "RLE")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_RLE <- glmQLFTest(fit)
topTags(result_RLE, n = 30)
## Different normalisation method
## 2. method = "upperquartile"
y_norm <- calcNormFactors(y, method = "upperquartile")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_upperquartile <- glmQLFTest(fit)
topTags(result_upperquartile, n = 30)
## Different normalisation method
## 3. method = "none"
y_norm <- calcNormFactors(y, method = "none")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_none <- glmQLFTest(fit)
topTags(result_none, n = 30)
ratio <- rowMeans(RNA[,9:16]+1)/rowMeans(RNA[,1:8]+1)
ratio_mean <- ratio[order(ratio, decreasing = TRUE)]
## plot of log(mean counts + 1) over groups
## points above the line indicate
plot(log(rowMeans(RNA[,9:16]+1)),rowMeans(log(RNA[,1:8]+1)),
     xlab = "log mean(RNA_counts + 1) control",
     ylab = "log mean(RNA_counts + 1) patients")
abline(a = 0, b = 1)
## Assess whether there are outliers, with comparison mean / median.
ratio <- rowMedians(as.matrix(RNA[,1:3])+1)/rowMedians(as.matrix(RNA[,4:6])+1)
ratio_median <- ratio[order(ratio, decreasing = TRUE)]
temp <- cbind.data.frame(RNA,ratio_mean = ratio_mean, ratio_median = ratio_median)
plot(temp$ratio_mean,temp$ratio_median,
     xlab = "ratio patient/control mean counts",
     ylab = "ratio patient/control median counts")
abline(a = 0, b = 1)
abline(a = 0, b = 2, col = "red")
temp[temp$ratio_median/temp$ratio_mean>2,]
# weggooien lowexpressed microRNAs (is that relevant?)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
class<-data.frame(condition = factor(rep(c("H","P"),c(8,8))))
as.factor(class[,1])
data.train <- RNA
data.train <- as.matrix(round(data.train + 1))
classtr <- data.frame(condition = class)
dim(data.train)
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train, colData = classtr, formula(~condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")
data.trainS4
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
cart
trained(cart)
confusionMat(cart)
## method = 'randomforest'
rf <- classify(data = data.trainS4, method = "randomforest", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
rf
trained(rf)
confusionMat(rf)
## method = 'randomforest'
rf_tmm <- classify(data = data.trainS4, method = "randomforest", normalize = "tmm", deseqTransform = "voom_CPM", cv = 5, rpt = 3, ref = "H")
rf_tmm
trained(rf_tmm)
confusionMat(rf_tmm)
## method = 'svm'
svm <- classify(data = data.trainS4, method = "svm", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
svm
trained(svm)
confusionMat(svm)
## method = 'bagsvm'
bagsvm <- classify(data = data.trainS4, method = "bagsvm", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
bagsvm
trained(bagsvm)
confusionMat(bagsvm)
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
cart
rf_tmm <- classify(data = data.trainS4, method = "randomforest", normalize = "tmm", deseqTransform = "voom_CPM", cv = 5, rpt = 3, ref = "H")
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
data
data.trainS4
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
y
et <- exactTest(y)
install.packages("FactoMineR")
library(FactoMineR)
df <- read_excel("C:\\Users\\dsa\\Desktop\\LungCancerPanel25of947miRs28pts.xlsx", sheet=1)
library(readxl)
df <- read_excel("C:\\Users\\dsa\\Desktop\\LungCancerPanel25of947miRs28pts.xlsx", sheet=1)
df
res <- PCA(df)
df <- read_excel("C:\\Users\\dsa\\Desktop\\LungCancerPanel25of947miRs28ptsClean.xlsx", sheet=1)
res <- PCA(df)
res <- PCA(df')
LungCancerPanel20of743miRs13ptsH2listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH2listT.txt", header=TRUE, row.names=1)
>   View(LungCancerPanel20of743miRs13ptsH2listT)
LungCancerPanel20of743miRs13ptsH2listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH2listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH2listT)
library(FactoMineR)
res <- PCA(LungCancerPanel20of743miRs13ptsH2listT)
df<-scale(LungCancerPanel20of743miRs13ptsH2listT)
heatmap(df,scale="none")
df<-scale(LungCancerPanel20of743miRs13ptsH2listT'')
df<-scale(LungCancerPanel20of743miRs13ptsH2listT')
LungCancerPanel20of743miRs13ptsH2list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH2list.txt", header=TRUE, row.names=1)
df1<-scale(LungCancerPanel20of743miRs13ptsH2list')
df1<-scale(LungCancerPanel20of743miRs13ptsH2list)
heatmap(df1,scale="none")
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df1, scale = "none", col =  col,
RowSideColors = rep(c("blue", "pink"), each = 16),
ColSideColors = c(rep("purple", 5), rep("orange", 6)))
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col,
RowSideColors = rep(c("blue", "pink"), each = 28),
ColSideColors = c(rep("purple", 5), rep("orange", 6)))
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df1, scale = "none", col =  col,
RowSideColors = c(rep("purple", 13), rep("orange", 15)))
library("RColorBrewer")
heatmap(df1,scale="none")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df1, scale = "none", col =  col)
heatmap(df2, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
## Compute genewise EXACT tests for differences in the means between two groups of
## negative-binomially distributed counts.
## logFC = log Fold change
## Display the 15 microRNAs with the lowest p-value
## In the column Pvalue the p-value of the differential test is displayed
## Give 15th most diferential markers and calculate family wise error based on method = 'hochberg'.
## Karin to do: what method to calculate family wise error is best?
## alternatives are: holm, hommel, bonferroni, BH, BY and fdr.
y <- estimateDisp(y, design)
et <- exactTest(y)
topTags(et, n = 15, adjust.method = "hochberg")
y_norm <- calcNormFactors(y, method = "TMM")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(d,design)
result <- glmQLFTest(fit)
topTags(result, n = 30)
## Same, but now with robust = TRUE
## robust = TRUE: whether to estimate the prior QL dispersion distribution robustly.
fit <- glmQLFit(d, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results, n = 30)
## Same, but now with abundance.trend = FALSE
fit <- glmQLFit(d, design, abundance.trend=FALSE)
results <- glmQLFTest(fit)
topTags(results, n = 30)
## Test other normalisation factors
## 1. method = "RLE
y_norm <- calcNormFactors(y, method = "RLE")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_RLE <- glmQLFTest(fit)
topTags(result_RLE, n = 30)
## Different normalisation method
## 2. method = "upperquartile"
y_norm <- calcNormFactors(y, method = "upperquartile")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_upperquartile <- glmQLFTest(fit)
topTags(result_upperquartile, n = 30)
## Different normalisation method
## 3. method = "none"
y_norm <- calcNormFactors(y, method = "none")
y_norm$samples
d <- estimateDisp(y_norm, design)
fit <- glmQLFit(y_norm,design)
result_none <- glmQLFTest(fit)
topTags(result_none, n = 30)
###############################################################################
###############################################################################
## Ruwe data analyse
###############################################################################
###############################################################################
ratio <- rowMeans(RNA[,9:16]+1)/rowMeans(RNA[,1:8]+1)
ratio_mean <- ratio[order(ratio, decreasing = TRUE)]
## plot of log(mean counts + 1) over groups
## points above the line indicate
plot(log(rowMeans(RNA[,9:16]+1)),rowMeans(log(RNA[,1:8]+1)),
xlab = "log mean(RNA_counts + 1) control",
ylab = "log mean(RNA_counts + 1) patients")
abline(a = 0, b = 1)
## Assess whether there are outliers, with comparison mean / median.
ratio <- rowMedians(as.matrix(RNA[,1:3])+1)/rowMedians(as.matrix(RNA[,4:6])+1)
ratio_median <- ratio[order(ratio, decreasing = TRUE)]
temp <- cbind.data.frame(RNA,ratio_mean = ratio_mean, ratio_median = ratio_median)
plot(temp$ratio_mean,temp$ratio_median,
xlab = "ratio patient/control mean counts",
ylab = "ratio patient/control median counts")
abline(a = 0, b = 1)
abline(a = 0, b = 2, col = "red")
temp[temp$ratio_median/temp$ratio_mean>2,]
# weggooien lowexpressed microRNAs (is that relevant?)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
class<-data.frame(condition = factor(rep(c("H","P"),c(8,8))))
as.factor(class[,1])
data.train <- RNA
data.train <- as.matrix(round(data.train + 1))
classtr <- data.frame(condition = class)
dim(data.train)
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train, colData = classtr, formula(~condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")
data.trainS4
## method = 'cart'
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
cart
trained(cart)
library(DESeq2)
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train, colData = classtr, formula(~condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")
data.trainS4
## method = 'cart'
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
cart
trained(cart)
confusionMat(cart)
library(classify)
BiocManager::install("MLSeq")
## method = 'cart'
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
acart
trained(cart)
confusionMat(cart)
## method = 'randomforest'
rf <- classify(data = data.trainS4, method = "randomforest", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
## different possible algorithms: svm (=support vector machines),
## bagsvm (= support vector machines with bagging ensemble),
## randomForest (= random forest algorithm), cart (classification
## and regression trees algorithm).
## Normalization possibilities: 'none', 'tmm' (trimmed mean of M values),
## 'deseq' (deseq normalization).
## other parameters:
## cv = number of cross-validation folds
## rpt = number of complete sets of folds for computation
## ref = reference class (choose healthy samples)
BiocManager::install("MLSeq")
library(MLSeq)
## different possible algorithms: svm (=support vector machines),
## bagsvm (= support vector machines with bagging ensemble),
## randomForest (= random forest algorithm), cart (classification
## and regression trees algorithm).
## Normalization possibilities: 'none', 'tmm' (trimmed mean of M values),
## 'deseq' (deseq normalization).
## other parameters:
## cv = number of cross-validation folds
## rpt = number of complete sets of folds for computation
## ref = reference class (choose healthy samples)
BiocManager::install("MLSeq")
library(MLSeq)
## method = 'cart'
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
acart
trained(cart)
LungCancerPanel20of743miRs13ptsH1list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1list.txt", header=FALSE)
View(LungCancerPanel20of743miRs13ptsH1list)
LungCancerPanel20of743miRs13ptsH1list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1list.txt", header=TRUE)
View(LungCancerPanel20of743miRs13ptsH1list)
data(LungCancerPanel20of743miRs13ptsH1list)
LungCancerPanel20of743miRs13ptsH1list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1list.txt", header=TRUE)
View(LungCancerPanel20of743miRs13ptsH1list)
data("LungCancerPanel20of743miRs13ptsH1list")
res.pca<-(LungCancerPanel20of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1list)
LungCancerPanel20of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1listT.txt", header=FALSE)
View(LungCancerPanel20of743miRs13ptsH1listT)
LungCancerPanel20of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
>   View(LungCancerPanel20of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1listT)
LungCancerPanel20of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1listT)
LungCancerPanel20of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1listT)
LungCancerPanel20of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1listT)
LungCancerPanel20of743miRs13ptsH1list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH1list.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH1list)
d <- dist(as.matrix(LungCancerPanel20of743miRs13ptsH1listT))
hc <- hclust(d)
plot(hc)
heatmap(df, scale = "none", col =  col)
heatmap(df1, scale = "none", col =  col)
heatmap(df2, scale = "none", col =  col)
df_new<-scale(LungCancerPanel20of743miRs13ptsH1listT')
df_new<-scale(LungCancerPanel20of743miRs13ptsH1listT)
heatmap(df_new, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
LungCancerPanel17of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel17of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel17of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel17of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel17of743miRs13ptsH1listT)
df_new2<-scale(LungCancerPanel17of743miRs13ptsH1listT)
heatmap(df_new2, scale = "none", col =  col)
heatmap(df_new, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
LungCancerPanel17of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel17of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel17of743miRs13ptsH3listT)
LungCancerPanel20of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH3listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH3listT)
df_new3<-scale(LungCancerPanel20of743miRs13ptsH3listT)
heatmap(df_new3, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
lung_oncomirs_medNor_722H3.xlsx <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/lung_oncomirs_medNor_722H3.xlsx", header=TRUE, row.names=1)
View(lung_oncomirs_medNor_722H3.xlsx)
lung_oncomirs_medNor_722H3.xlsx <- read.delim("C:\Users\dsa\Documents\MATLAB\lung_oncomirs_medNor_722H3.xlsx", header=TRUE, row.names=1)
View(lung_oncomirs_medNor_722H3.xlsx)
lung_oncomirs_medNor_722H3.xlsx <- read.delim("C:\Users\dsa\Documents\MATLAB\lung_oncomirs_medNor_722H3.xlsx", header=TRUE, row.names=1)
View(lung_oncomirs_medNor_722H3.xlsx)
LungCancerPanel20of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH3listT)
res.pca<-PCA(LungCancerPanel20of743miRs13ptsH3listT)
df_new3<-scale(LungCancerPanel20of743miRs13ptsH3listT)
heatmap(df_new3, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
LungCancerPanel5of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel5of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel5of743miRs13ptsH3listT)
res.pca<-PCA(LungCancerPanel5of743miRs13ptsH3listT)
df_new5<-scale(LungCancerPanel5of743miRs13ptsH3listT)
heatmap(df_new5, scale = "none", col =  col)
LungCancerPanel5of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel5of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel5of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel5of743miRs13ptsH1listT)
df_new6<-scale(LungCancerPanel5of743miRs13ptsH1listT)
heatmap(df_new6, scale = "none", col =  col)
heatmap(df_new6, scale = "none", col =  col)
heatmap(df_new5, scale = "none", col =  col)
LungCancerPanel5of743miRs13ptsH2listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel5of743miRs13ptsH2listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel5of743miRs13ptsH2listT)
df_new7<-scale(LungCancerPanel5of743miRs13ptsH2listT)
heatmap(df_new7, scale = "none", col =  col)
heatmap(df, scale = "none", col =  col)
LungCancerPanel3of743miRs13ptsH2listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH2listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH2listT)
d23<-scale(LungCancerPanel3of743miRs13ptsH2listT)
heatmap(df23, scale = "none", col =  col)
df23<-scale(LungCancerPanel3of743miRs13ptsH2listT)
heatmap(df23, scale = "none", col =  col)
LungCancerPanel3of743miRs13ptsH1listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH1listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH1listT)
df13<-scale(LungCancerPanel3of743miRs13ptsH1listT)
heatmap(df13, scale = "none", col =  col)
LungCancerPanel3of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH3listT)
df33<-scale(LungCancerPanel3of743miRs13ptsH3listT)
heatmap(df33, scale = "none", col =  col)
heatmap(df23, scale = "none", col =  col)
res.pca<-PCA(LungCancerPanel2of743miRs13ptsH1listT)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH1listT)
LungCancerPanel3of743miRs13ptsH3list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH3list.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH3list)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH3list)
LungCancerPanel3of743miRs13ptsH2list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH2list.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH2list)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH2list)
LungCancerPanel3of743miRs13ptsH1list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel3of743miRs13ptsH1list.txt", header=TRUE, row.names=1)
View(LungCancerPanel3of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH1list)
LungCancerPanel1of743miRs13ptsH2list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel1of743miRs13ptsH2list.txt", header=TRUE, row.names=1)
View(LungCancerPanel1of743miRs13ptsH2list)
res.pca<-PCA(LungCancerPanel3of743miRs13ptsH1list)
res.pca<-PCA(LungCancerPanel1of743miRs13ptsH2list)
LungCancerPanel1of743miRs13ptsH2listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel1of743miRs13ptsH2listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel1of743miRs13ptsH2listT)
df12<-scale(LungCancerPanel1of743miRs13ptsH3listT)
df12<-scale(LungCancerPanel1of743miRs13ptsH2listT)
heatmap(df12, scale = "none", col =  col)
library(clv)
install.packages(clv)
library(clusterCrit)
BiocManager::install("clv")
BiocManager::install("clusterCrit")
library(clusterCrit)
library(clv)
plot(hc)
clv.Dunn(hc)
Dunn(hc)
plot(hclust(hc),  main="Dissimilarity = 1 - Correlation", xlab="")
plot(hc,  main="Dissimilarity = 1 - Correlation", xlab="")
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH2listT)
distance <- as.dist(dissimilarity)
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH2list)
distance <- as.dist(dissimilarity)
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH1list)
distance <- as.dist(dissimilarity)
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH3list)
dissimilarity <- 1 - cor(LungCancerPanel20of722miRs13ptsH3list)
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH3list)
LungCancerPanel20of743miRs13ptsH3listT <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH3listT.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH3list)
LungCancerPanel20of743miRs13ptsH3list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH3list.txt", header=FALSE)
View(LungCancerPanel20of743miRs13ptsH3list)
LungCancerPanel20of743miRs13ptsH3list <- read.delim("D:/ProgramFiles/R-3.6.2/library/FactoMineR/data/LungCancerPanel20of743miRs13ptsH3list.txt", header=TRUE, row.names=1)
View(LungCancerPanel20of743miRs13ptsH3list)
dissimilarity <- 1 - cor(LungCancerPanel20of743miRs13ptsH3list)
distance <- as.dist(dissimilarity)
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
devtools::install_github("rstudio/keras")
install.packages('devtools')
install.packages("devtools")
devtools::install_github("rstudio/keras")
exit
devtools::install_github("rstudio/keras")
update.packages('rlang')
library('rlang')
library(rlang)
detach("package:rlang", unload = TRUE)
library(rlang)
update.packages()
install.packages('rland')
install.packages('installr')
install.Rtools()
installr::install.Rtools()
devtools::install_github("rstudio/keras")
library("rlang")
update.packages()
save.image(file='LungCancer_miRs.RData')
load("C:/Users/dsa/Documents/LungCancer_miRs.RData")
devtools::install_github("rstudio/keras")
install.packages('rland')
install.packages('rlang')
devtools::install_github("rstudio/keras")
library('keras')
devtools::install_github("rstudio/keras")
install.packages('ps')
devtools::install_github("rstudio/keras")
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
plot(hc,  main="Dissimilarity = 1 - Correlation", xlab="")
heatmap(df12, scale = "none", col =  col)
df12<-scale(LungCancerPanel1of743miRs13ptsH2listT)
heatmap(df12, scale = "none", col =  col)
plot(hc)
plot(hclust(hc),  main="Dissimilarity = 1 - Correlation", xlab="")
plot(hc,  main="Dissimilarity = 1 - Correlation", xlab="")
df<-scale(LungCancerPanel20of743miRs13ptsH2listT')
df<-scale(LungCancerPanel20of743miRs13ptsH2listT)
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col)
