help(cor)
x=c
x=c(1,2,3,4)
cor(x, x, method="spearman")
library("gplots")
# Read in data
#HeatMap_Data = read.csv("Heatmap.csv")
HeatMap_Data = read.csv("MSK-PCa1-7_70Genes_v3.txt", sep="\t")
rownames(HeatMap_Data) = HeatMap_Data$X
HeatMap_Data = HeatMap_Data[,-c(1)]
Group = substring(colnames(HeatMap_Data),1,2)
Group.Color = c("yellow","cyan","blue","black","red","green")
ColSideColors = Group.Color[match(Group, unique(Group))];
###################################################
# Create Heap Maps with raw values and no "dendrogram" clustering
win.graph()
heatmap.2(as.matrix(HeatMap_Data), cexRow = 0.8, dendrogram = "none", scale ="none", Rowv = NA, Colv = NA, col=greenred(75), key=TRUE, symkey=FALSE, density.info="none", trace="none")
win.graph()
# Create Heap Maps with z score value on row, and no "dendrogram" clustering, and ColSideColors
heatmap.2(as.matrix(HeatMap_Data), na.rm = T, cexRow = 0.8, dendrogram = "none", scale ="row", Rowv = NA, Colv = NA, col=greenred(75), key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=ColSideColors)
win.graph()
# Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
heatmap.2(as.matrix(HeatMap_Data), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="row", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=ColSideColors)
win.graph()
# Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
getwd(
  getwd()
  getwd()
  library("gplots")
  # Read in data
  #HeatMap_Data = read.csv("Heatmap.csv")
  HeatMap_Data = read.csv("MSK-PCa1-7_70Genes_v3.txt", sep="\t")
  rownames(HeatMap_Data) = HeatMap_Data$X
  HeatMap_Data = HeatMap_Data[,-c(1)]
  Group = substring(colnames(HeatMap_Data),1,2)
  Group.Color = c("yellow","cyan","blue","black","red","green")
  ColSideColors = Group.Color[match(Group, unique(Group))];
  ###################################################
  # Create Heap Maps with raw values and no "dendrogram" clustering
  win.graph()
  heatmap.2(as.matrix(HeatMap_Data), cexRow = 0.8, dendrogram = "none", scale ="none", Rowv = NA, Colv = NA, col=greenred(75), key=TRUE, symkey=FALSE, density.info="none", trace="none")
  win.graph()
  # Create Heap Maps with z score value on row, and no "dendrogram" clustering, and ColSideColors
  heatmap.2(as.matrix(HeatMap_Data), na.rm = T, cexRow = 0.8, dendrogram = "none", scale ="row", Rowv = NA, Colv = NA, col=greenred(75), key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=ColSideColors)
  win.graph()
  # Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
  heatmap.2(as.matrix(HeatMap_Data), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="row", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none", ColSideColors=ColSideColors)
  win.graph()
  # Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
  heatmap.2(as.matrix(t(HeatMap_Data)), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
  library("gplots")
  # Read in data
  #HeatMap_Data = read.csv("Heatmap.csv")
  HeatMap_Data = read.csv("MSK-PCa1-PCa2-PCa3-PCa5_70Genes.txt", sep="\t")
  rownames(HeatMap_Data) = HeatMap_Data$X
  HeatMap_Data = HeatMap_Data[,-c(1)]
  Group = substring(colnames(HeatMap_Data),1,2)
  Group.Color = c("yellow","cyan","blue","black","red","green")
  ColSideColors = Group.Color[match(Group, unique(Group))];
  ###################################################
  win.graph()
  # Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
  heatmap.2(as.matrix(t(HeatMap_Data)), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
  library("gplots")
  # Read in data
  #HeatMap_Data = read.csv("Heatmap.csv")
  HeatMap_Data = read.csv("MSK-PCa1-PCa2-PCa3-PCa5_70Genes.txt", sep="\t")
  rownames(HeatMap_Data) = HeatMap_Data$X
  HeatMap_Data = HeatMap_Data[,-c(1)]
  Group = substring(colnames(HeatMap_Data),1,2)
  Group.Color = c("yellow","cyan","blue","black","red","green")
  ColSideColors = Group.Color[match(Group, unique(Group))];
  ###################################################
  win.graph()
  # Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
  heatmap.2(as.matrix(t(HeatMap_Data)), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
  library("gplots")
  # Read in data
  #HeatMap_Data = read.csv("Heatmap.csv")
  HeatMap_Data = read.csv("MSK-PCa1-PCa2-PCa3-PCa5_70Genes.txt", sep="\t")
  rownames(HeatMap_Data) = HeatMap_Data$X
  HeatMap_Data = HeatMap_Data[,-c(1)]
  Group = substring(colnames(HeatMap_Data),1,2)
  Group.Color = c("yellow","cyan","blue","black","red","green")
  ColSideColors = Group.Color[match(Group, unique(Group))];
  ###################################################
  speed = c(15.5, 19, 19.3, 18.7)
  win.graph()
  # Create Heap Maps with z score value on row, with clustering on rows and columns, and ColSideColors
  heatmap.2(as.matrix(t(HeatMap_Data)), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
  dim(HeatMap_Data)
  HeatMap_Data[1,]
  HeatMap_Data[2]
  HeatMap_Data[2]
  HeatMap_Data[2]
  HeatMap_Data[2,]
  dim(HeatMap_Data)
  x = ""
  for i = 1:66
  x = c(x, cor(speed, HeatMap_Data[i,], method="spearman"))
  end
  x = ""
  for i in 1:66
  x = c(x, cor(speed, HeatMap_Data[i,], method="spearman"))
  end
  x = NULL
  for i in 1:66
  x = c(x, cor(speed, HeatMap_Data[i,], method="spearman"))
  end
  for i in 1 to 66
  help for
  x = NULL
  for (i in 1:66) {
    x = c(x, cor(speed, HeatMap_Data[i,], method="spearman"))
  }
  }i
i
cor(speed, HeatMap_Data[i,], method="spearman")
speed
HeatMap_Data[1,]
as.numeric(HeatMap_Data[1,])
cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
x = NULL
for (i in 1:66) {
  cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
}
x = NULL
for (i in 1:66) {
  x = c(x, cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman"))
}
x
plot(x)
plot(x, xlab=rownames(HeatMap_Data))
plot(x, xlab=rownames(HeatMap_Data))
rownames(HeatMap_Data)
rbind(rownames(HeatMap_Data), x)
plot(speed, HeatMap_Data[49,])
plot(speed, as.numeric(HeatMap_Data[49,]))
plot(speed, as.numeric(HeatMap_Data[48,]))
plot(speed, as.numeric(HeatMap_Data[48,]), ylab="RNASeq FPKM")
write.table(file= "Speed_Correlation.txt", plot(speed, as.numeric(HeatMap_Data[48,]), ylab="RNASeq FPKM"), row.names=F, colnames=F, quote=F)
write.table(file= "Speed_Correlation.txt", plot(speed, as.numeric(HeatMap_Data[48,]), ylab="RNASeq FPKM"), row.names=F, col.names=F, quote=F)
write.table(file= "Speed_Correlation.txt", rbind(speed, HeatMap_Data[48,]), ylab="RNASeq FPKM"), row.names=F, col.names=F, quote=F)
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data), x), ylab="RNASeq FPKM"), row.names=F, col.names=F, quote=F)
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data), ylab="RNASeq FPKM"), row.names=F, col.names=F, quote=F)
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data),x), row.names=F, col.names=F, quote=F)
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data),x), row.names=F, col.names=F, quote=F, sep="\t")
plot(speed, as.numeric(HeatMap_Data[48,]), ylab="RNASeq FPKM")
heatmap.2( rbind(rownames(HeatMap_Data),x))
heatmap.2(x)
x
heatmap.2(as.numeric(x))
heatmap.2(rbind(x,x))
heatmap.2(as.matrix(t(HeatMap_Data)), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
heatmap.2(rbind(x,x), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
heatmap.2(x, na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
heatmap.2(rbind(x,x), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
dim(HeatMap_Data)
length(x)
rbind(HeatMap_Data,x)
cbind(HeatMap_Data,x)
cbind(HeatMap_Data,x)
heatmap.2(as.matrix(t(cbind(HeatMap_Data,x))), na.rm = T, cexRow = 0.8, dendrogram = "both", scale ="col", Rowv = TRUE, Colv = TRUE, col=greenred(75),  key=TRUE, symkey=FALSE, density.info="none", trace="none")
help(cor)
cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
[x1,x2] = cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")$p.value
x = cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
x
x = cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")$rho
x
cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")
x = NULL
p = NULL
for (i in 1:66) {
  x = c(x, cor(speed, as.numeric(HeatMap_Data[i,]), method="spearman"))
  p = c(p, cor.test(speed, as.numeric(HeatMap_Data[i,]), method="spearman")$p.value)
}
plot(speed, as.numeric(HeatMap_Data[48,]), ylab="RNASeq FPKM")
heatmap.2( rbind(rownames(HeatMap_Data),x))
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data),x,p), row.names=F, col.names=F, quote=F, sep="\t")
write.table(file= "Speed_Correlation.txt", rbind(rownames(HeatMap_Data),x,p), row.names=F, col.names=F, quote=F, sep="\t")
plot(speed, as.numeric(HeatMap_Data[41,]), ylab="RNASeq FPKM of NUSAP1")
q()
