remove(mirs)
install.packages("rlang")

install.packages("readxl")
library(readxl)

df <- read_excel("C:\\Users\\dsa\\Desktop\\UriWatchDataSet001.xlsx", sheet=1 , col_types = "numeric")

p24=df[4:613,11]
p25=df[4:613,20]
p26=df[4:613,29]
c7=df[4:613,38]
c8=df[4:613,47]
c9=df[4:613,56]

df <- read_excel("C:\\Users\\dsa\\Desktop\\UriWatchDataSet001.xlsx", sheet=1 )

mirs <- matrix(0,610, 7)

mirs[1:610,1]=t(mirList)

mirs[1:610,2]= as.numeric(
  t(c7))
mirs[1:610,3]= as.numeric(
  t(c8))
mirs[1:610,4]=t(c9)
mirs[1:610,5]=t(p24)
mirs[1:610,6]=t(p25)


typeof(mirs)
#[1] "double"
dim(mirs)
#[1] 610   6
install.packages("FactoMineR") 
library(FactoMineR)
res <- PCA(mirs)

df <- read_excel("C:\\Users\\dsa\\Desktop\\UriWatchDataSet001.xlsx", sheet=1)
mirList=df[4:613,1]
mirList
# A tibble: 610 x 1
X__1
#<chr>
#  1 hsa-miR-10b-5p
#2 hsa-miR-10a-5p
##3 hsa-miR-486-5p
#4  hsa-let-7f-5p
#5  hsa-let-7a-5p
#6 hsa-miR-99b-5p
#7 hsa-miR-191-5p
#8 hsa-miR-100-5p
#9 hsa-miR-30a-5p
#10  hsa-let-7b-5p
# ... with 600 more rows

mirs[1:610,2:7]=as.numeric(mirs[1:610,2:7])


PCA(df[,2:7])

PCA(as.numeric(df3[,1:610]))
