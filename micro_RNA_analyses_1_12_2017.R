###############################################################################
###############################################################################
## Analysis microRNA data
## December 1st, 2017
## Karin Groothuis-Oudshoorn
###############################################################################
###############################################################################
## Load the necessary libraries
library(DESeq2)
library(edgeR)

## December 1st, 2017: next dataset
## File: Data_1_12_2017.csv
## From this dataset we only need the names of the microRNAs (= first column) and the raw counts.
## The raw counts are in the columnnumbers 2, 7, 12, 17 etc (=seq(from = 2, by = 5, to = 77))
##temp <- read.table("C:/Users/oudshoornc/Dropbox/UriWatch/Data/Data_1_12_2017.csv", sep = ";", header = TRUE, row.names = 1, nrows = 687)

temp <- read.table("C:/Users/dsa/Desktop/Data_1_12_2017_2019.csv", sep = ";", header = TRUE, row.names = 1, nrows = 687)

RNA <- temp[-1, seq(from = 1, to = 76, by = 5)]

## What are the names of the samples?
colnames(RNA)
## First column is the microRNA name.
## Column 2 to 9 are control data, rest is cases.
## The factor 'group' indicates which sample belongs to which group (two groups: 2 = case (patient), 1 = healthy)
## 'design' is a matrix that is needed later on to test whether there are differential expression between the samples due to the factor 'group'.
## y is a list where the data is stored together with some extra information
## 1 december 2017: first 8 are healthy, second 8 are patient.
## To compare patients with controls we define the contrast such that effects are estimated as differences
## from controls (statement: contrast(group) <- c(-1,1))
group <- factor(c(rep(1,8),rep(2,8)))
contrasts(group) <- c(-1,1)
design <- model.matrix( ~ group)
y <- DGEList(counts = RNA, group = group)

## show on the screen what is in the list object y.
y

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

## Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
## Conduct microRNAwise statistical tests for differenes between the two groups.
## Remark: the result of the 15th most differential microRNAs depends on the normalisation factor!
## Calculate normalization factors to scale the raw library sizes
## Default normalisation factor is "TMM". 

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
## Conclusion: differential expression depends on normalisation method.

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


###############################################################################
###############################################################################
##### classify
##### This is how the code is. 
##### Beware: six samples is not enough..
##### Below all samples are used to derive the classifier...
##### cross validation is applied.

class<-data.frame(condition = factor(rep(c("H","P"),c(8,8))))
as.factor(class[,1])

data.train <- RNA
data.train <- as.matrix(round(data.train + 1))
classtr <- data.frame(condition = class)
dim(data.train)

data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train, colData = classtr, formula(~condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")
data.trainS4

## classification part
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
## method = 'cart'
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", deseqTransform = "vst", cv = 5, rpt = 3, ref = "H")
acart
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

### TO DO (Karin)
### 1. find out how to graphically present results of ML algorithm (dendrogram, VarImp)
### 2. how to present results over time (differential expressions change after adding more samples)
### 3. What can be done with PCA analysis? (on samples and on microRNA)
### 4. How to identify what microRNA markers are important in the classifier (e.g. for RandomForest)