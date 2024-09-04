library("plyr")
library("data.table")


## this script performs the preprocessing of GDSC data

## gene expression and mutation data are kept as two sources, with cancer type as modifying variable




#===============
# The user needs load the datasets from https://www.cancerrxgene.org/downloads/bulk_download archived data folder 'release-5.0'. 
# Downloading the three datasets used for our analysis.
#===============

features <- data.frame(read.csv("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/gdsc data/gdsc_en_input_w5.csv", head=T))
names.fea <- strsplit(rownames(features), "")
features <- t(features)
p <- c(13321, 13747-13321, 13818-13747)
Cell.Line <- rownames(features)
features <- data.frame(Cell.Line, features)

ic50_00 <- data.frame(read.csv("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/gdsc data/gdsc_drug_sensitivity_fitted_data_w5.csv", head=T))
ic50_0 <- ic50_00[,c(1,4,7)]
drug.id <- data.frame(read.csv("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/gdsc data/gdsc_tissue_output_w5.csv", head=T))[,c(1,3)]
drug.id2 <- drug.id[!duplicated(drug.id$drug.id),]
# delete drug.id=1066 since ID1066 and ID156 both correspond drug AZD6482, 
# and no ID1066 in the "suppl.Data1" by Garnett et al. (2012)
drug.id2 <- drug.id2[drug.id2$drug.id!=1066,] 
drug.id2$drug.name <- as.character(drug.id2$drug.name)
drug.id2$drug.name <- substr(drug.id2$drug.name, 1, nchar(drug.id2$drug.name)-6)
drug.id2$drug.name <- gsub(" ", "-", drug.id2$drug.name)

ic50 <- ic50_0
# mapping the drug_id to drug names in drug sensitivity data set
ic50$drug_id <- plyr::mapvalues(ic50$drug_id, from = drug.id2[,2], to = drug.id2[,1])
colnames(ic50) <- c("Cell.Line", "compound", "IC50")

# transform drug sensitivity overall cell lines to a data matrix
y0 <- reshape(ic50, v.names="IC50", timevar="compound", idvar="Cell.Line", direction="wide")
y0$Cell.Line <- gsub("-", ".", y0$Cell.Line)

#===============
# select nonmissing pharmacological data
#===============
y00 <- y0
m0 <- dim(y0)[2]-1
eps <- 0.05
# r1.na is better to be not smaller than r2.na
r1.na <- 0.3
r2.na <- 0.2
k <- 1
while(sum(is.na(y0[,2:(1+m0)]))>0){
  r1.na <- r1.na - eps/k
  r2.na <- r1.na - eps/k
  k <- k + 1
  ## select drugs with <30% (decreasing with k) missing data overall cell lines
  na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y<r1.na)<m0){
    y0 <- y0[,-c(1+which(na.y>=r1.na))]
    m0 <- sum(na.y<r1.na)
    na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  }
  
  ## select cell lines with treatment of at least 80% (increasing with k) drugs
  na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y0<r2.na)<(dim(y0)[1])){
    y0 <- y0[na.y0<r2.na,]
    na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  }
  num.na <- sum(is.na(y0[,2:(1+m0)]))
  message("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
}

#===============
# combine drug sensitivity, tissues and molecular features
#===============
yx <- merge(y0, features, by="Cell.Line")
names.cell.line <- yx$Cell.Line
names.drug <- colnames(yx)[2:(dim(y0)[2])]
names.drug <- substr(names.drug, 6, nchar(names.drug))

# numbers of gene expression features, copy number features and mutation features
p <- c(13321, 13747-13321, 13818-13747) 
num.nonpen <- 13
yx <- data.matrix(yx[,-1])

# we select x1 as gene expression and x2 as mutations
y <- yx[,1:(dim(y0)[2]-1)]
x <- cbind(yx[,dim(y0)[2]-1+sum(p)+1:num.nonpen], yx[,dim(y0)[2]-1+1:sum(p)])

## preselect gene expression explaining 10%/30%/50% variations over all samples
var.x1 <- apply(log(x[,num.nonpen+1:p[1]]), 2, var)
var.sortx1 <- sort(var.x1, decreasing=TRUE)
sum.ax1 <- cumsum(var.sortx1)
half.ax1 <- var.sortx1[which(sum.ax1>(sum.ax1[length(var.x1)]*0.5))[1]] # explaining 50% variations
#half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]*0.3))[1]] # explaining 30% variations
#half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]*0.1))[1]] # explaining 10% variations
x1 <- x[,num.nonpen+1:p[1]][,var.x1>=half.ax1]

## preselect gene expression explaining 10%/30%/50% variations over all samples
x2 <- x[,num.nonpen+1+p[1]+p[2]:(p[2]+p[3]-1)]


z <- cbind(x[,1:num.nonpen]) ##
#p[1] <- dim(x1)[2]

# delete genes with only one mutated cell line
#x <- x[,-c(num.nonpen+p[1]+p[2]+which(colSums(x[,num.nonpen+p[1]+p[2]+1:p[3]])<=1))]
#p[3] <- ncol(x) - num.nonpen - p[1] - p[2]

GDSC <- list(y=y,z=z, x1=x1, x2=x2, p=dim(x1)[2], num.nonpen=num.nonpen, names.cell.line=names.cell.line, 
             names.drug=names.drug)
saveRDS(GDSC, file="GDSC_data_processed.rds")