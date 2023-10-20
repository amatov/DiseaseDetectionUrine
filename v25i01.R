## Load the package
library("FactoMineR")

## Load the dataset
data("decathlon")					       

## Perform a PCA with the decathlon dataset, 2 quantitative
## supplementary variables and 1 qualitative supplementary variable
res.pca <- PCA(decathlon, quanti.sup = 11:12, quali.sup = 13)     

## Plot the individuals graph; the individuals are colored
## according to the qualitative variable
plot(res.pca, habillage = 13)				       

## Plot the barplot of the eigenvalues
barplot(res.pca$eig[,1], main = "Eigenvalues", 	       
  names.arg = paste("Dim", 1:nrow(res.pca$eig), sep = ""))

## Plot the variables graph with dimensions 3 and 4
plot(res.pca, choix = "var", axes = c(3, 4), lim.cos2.var = 0) 

## Print the results
print(res.pca)  					       

## Perform the description of the dimensions
dimdesc(res.pca, proba = 0.2)				       

## Load the dataset
data("children")					       

## Perform a correspondence analysis with the children
## dataset, and supplementary columns and rows
res.ca <- CA(children, col.sup = 6:8, row.sup = 15:18)         

# Plot the graph with only the active data
plot(res.ca, invisible = c("row.sup", "col.sup"))	       
