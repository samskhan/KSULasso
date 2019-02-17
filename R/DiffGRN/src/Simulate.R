setwd("~/KSULasso/R/DiffGRN/")
source("src/RandomGraph.R")
source("src/RandomLasso.R")

library(readxl)
library(igraph)
library(ggraph)
library(ggnetwork)
library(lars)
library(pROC)
library(gplots)
library(RColorBrewer)
library(roxygen2)

# Loads xlsx data into memory. -Matthew
expData <- read_excel("res/AsthmaDNA.xlsx", col_names = TRUE)
# Removes entire row if any NA data is present. -Matthew
expData <- data.frame(na.omit(expData))
# Name of the genes are located in column one and two. -Matthew
gene_names <- as.character(as.vector(expData[, 2]))
# Removing columns one and two. -Matthew
expGene_Asthma <- expData[ ,-c(1,2)]
# Alt code that removes all columns with non-numeric values. -Matthew
# expGene_Asthma <- expData[ , !is.na(as.numeric(as.character(expData[1, ])))]
# Printing the dimensions of the data -Matthew
dim(expGene_Asthma)

# Converting gene names to numeric. -Matthew
tmpGene <- t(as.matrix(as.character(gene_names)))
n_r <- length(tmpGene)
dim(tmpGene)

# Attching gene names to gene numbers as a list. -Matthew
geneList <- list()
for (i in 1:n_r) {
    geneList[[i]] <- unlist(strsplit(tmpGene[i], "\\|"))
}
length(geneList)

tmp1 <- matrix(FALSE, nrow = n_r, ncol = 1 )
tmp_names <- matrix(FALSE, nrow = n_r, ncol = 1 )
for (i in 1:n_r) {
   tmp1[i] <- length(geneList[[i]]) == 1
}
dim(tmp1)

expGene_Asthma <- t(expGene_Asthma[tmp1, ])
gene_names <- gene_names[tmp1]
colnames(expGene_Asthma) <- gene_names
dim(expGene_Asthma)

gene_names  <- colnames(expGene_Asthma)
# Changing "." to "_" in gene names. -Matthew
gene_names <- chartr(".", "_", gene_names)
# Changing "-" to "_" in gene names. -Matthew
gene_names <- chartr("-", "_", gene_names)
length(gene_names) # 3139

# Checking for repeat gene names -Matthew
length(gene_names) # 3139
length(unique(gene_names)) # 3016

unique_gene <- unique(colnames(expGene_Asthma))
length(unique_gene) # 3016 
nc_gene <- length(unique_gene)
nr_gene <- nrow(expGene_Asthma)

AsthmaData <- matrix(0, nrow = nr_gene, ncol = nc_gene)
dim(AsthmaData)
colnames(AsthmaData) <- unique_gene
rownames(AsthmaData) <- rownames(expGene_Asthma)
for (i in 1:nc_gene) {
    same_gene <- which(!is.na(match(gene_names, unique_gene[i])))
    if (length(same_gene) == 1) {
        AsthmaData[, i] <- expGene_Asthma[, same_gene]
    } else if (length(same_gene) > 1) { 
        AsthmaData[, i] <- apply(expGene_Asthma[, same_gene], 1, max) 
    }
}
dim(AsthmaData) # 194 by 3016

# The first three letters of the row names. -Matthew
Group <- as.matrix(substr(rownames(AsthmaData), 1, 3))
# Dividing the data into two groups. -Matthew
# AST are those with Asthma, AST == Asthma. -Matthew
expGene_AST <- AsthmaData[(Group == "AST"), ] # 97 by 3016
# CON are those without Asthma, CON == Control. -Matthew
expGene_CON <- AsthmaData[(Group == "CON"), ] # 97 by 3016

# Number of columns. -Matthew
nc <- ncol(AsthmaData) # 3016
# Creating matrix for Pvalues. -Matthew
Pval_Of_ttest <- matrix(0, ncol = nc)
# Assigning gene names to matrix for Pvalues. -Matthew
names(Pval_Of_ttest) <- colnames(AsthmaData) 
# Running T-test between AST and CON for every column. -Matthew
for (ii in 1:nc) {
    tmp <- t.test(expGene_AST[, ii], expGene_CON[, ii]) 
    Pval_Of_ttest[ii] <- tmp$p.value
}
# Only P-values greater than 0.05 are kept. -Matthew
id_gene <- which(Pval_Of_ttest < 0.05)
Asthma_gene <- AsthmaData[, id_gene]

# Adding row names as a row. -Matthew
write.csv(cbind(rownames(Asthma_gene), Asthma_gene),
          "res/AsthmaGenes.csv", row.names = FALSE)

# Normalization. -Matthew
Norm_Asthma <- scale(Asthma_gene)

# The first three letters of the row names. -Matthew
Group <- as.matrix(substr(rownames(Norm_Asthma), 1, 3))
# Dividing the data into two groups. -Matthew
# AST are those with Asthma, AST == Asthma. -Matthew
Norm_Asthma_TRT <- Norm_Asthma[(Group == "AST"), ]
# CON are those without Asthma, CON == Control. -Matthew
Norm_Asthma_CON <- Norm_Asthma[(Group == "CON"), ]
dim(Norm_Asthma_TRT)  # 97 by 214

# Alt Normalization seperately by groups. -Matthew
# Group <- as.matrix(substr(rownames(Asthma_gene), 1, 3))
# Norm_Asthma_TRT <- Asthma_gene[(Group == "AST"), ]
# Norm_Asthma_CON <- Asthma_gene[(Group == "CON"), ]
# Norm_Asthma_TRT <- scale(Norm_Asthma_TRT)
# Norm_Asthma_CON <- scale(Norm_Asthma_CON) 
# dim(Norm_Asthma_TRT)  # 97 by 214

# Seed for random number generation. -Matthew
set.seed(34567)

# See RandomLasso.R for function. -Matthew
#size <- Determin_Feature_size(Data = Norm_Asthma_TRT, C_I = 0.95, P = 0.5, S_E = 0.05, REP = 20) 
#q1 <- size[[1]] # 120
#q2 <- size[[2]] # 120 # This is not used. -Matthew
#B <- size[[3]] # 36


# Check paper on Random Lasso page 472 for on Step 1. -Matthew
# http://dept.stat.lsa.umich.edu/~jizhu/pubs/Wang-AOAS11.pdf
#library("future") # Temp parrelle programming.
#plan(multiprocess)
#library(parallel)
library(glmnet)
setwd("~/KSULasso/R/DiffGRN/")
source("src/DiffGRN.R") # Temp for reloading functions for testing. -Matthew
Importance_genes_asthma <- DiffGRN(Norm_Asthma_TRT, Norm_Asthma_CON)
#Importance_genes_asthma <- RandomLasso(Norm_Asthma_TRT, Norm_Asthma_CON, suppress = TRUE)

# Check paper on Random Lasso page 472 for on Step 2. -Matthew
# http://dept.stat.lsa.umich.edu/~jizhu/pubs/Wang-AOAS11.pdf
#Importance_TRT_asthma2 <- abs(Importance_genes_asthma[[1]])
#Importance_CON_asthma2 <- abs(Importance_genes_asthma[[3]])
#Adj_Asthma_b_RandomLASSO <- RandomLasso(Norm_Asthma_TRT, Norm_Asthma_CON) 


Adj_gene_net_b_TRT_asthma <- Adj_Asthma_b_RandomLASSO[[1]]
Adj_gene_SE_b_TRT_asthma  <- Adj_Asthma_b_RandomLASSO[[2]]
Adj_gene_net_b_CON_asthma <- Adj_Asthma_b_RandomLASSO[[3]]  
Adj_gene_SE_b_CON_asthma  <- Adj_Asthma_b_RandomLASSO[[4]]

# This is the end. Everthing past here is mostly visualization. -Matthew
# END -Matthew

## save result
# write.csv(Adj_gene_net_b_TRT_asthma, "C:/Research/R/asthma/src/result/Adj_gene_net_b_TRT_asthma.csv", row.names=FALSE)
# write.csv(Adj_gene_SE_b_TRT_asthma, "C:/Research/R/asthma/src/result/Adj_gene_SE_b_TRT_asthma.csv", row.names=FALSE)
# write.csv(Adj_gene_net_b_CON_asthma, "C:/Research/R/asthma/src/result/Adj_gene_net_b_CON_asthma.csv", row.names=FALSE)
# write.csv(Adj_gene_SE_b_CON_asthma, "C:/Research/R/asthma/src/result/Adj_gene_SE_b_CON_asthma.csv", row.names=FALSE)

# # re-load result : sometimes we need this procedure
# Adj_gene_net_b_TRT_asthma <- as.matrix(data.frame(read.csv("C:/Research/R/asthma/src/result/Adj_gene_net_b_TRT_asthma.csv")))
# rownames(Adj_gene_net_b_TRT_asthma) <- colnames(Adj_gene_net_b_TRT_asthma)
# Adj_gene_SE_b_TRT_asthma<- as.matrix(data.frame(read.csv("C:/Research/R/asthma/src/result/Adj_gene_SE_b_TRT_asthma.csv")))
# rownames(Adj_gene_SE_b_TRT_asthma ) <- colnames(Adj_gene_SE_b_TRT_asthma)
# Adj_gene_net_b_CON_asthma <- as.matrix(data.frame(read.csv("C:/Research/R/asthma/src/result/Adj_gene_net_b_CON_asthma.csv")))
# rownames(Adj_gene_net_b_CON_asthma ) <- colnames(Adj_gene_net_b_CON_asthma)
# Adj_gene_SE_b_CON_asthma <- as.matrix(data.frame(read.csv("C:/Research/R/asthma/src/result/Adj_gene_SE_b_CON_asthma.csv")))
# rownames(Adj_gene_SE_b_CON_asthma ) <- colnames(Adj_gene_SE_b_CON_asthma)

## compute differential network our DGRN

##=========================================================================================================
## compute differentail network
diff_beta_asthma <- Adj_gene_net_b_TRT_asthma - Adj_gene_net_b_CON_asthma
typeof(diff_beta_asthma )
z_score_asthma <- diff_beta_asthma / sqrt(Adj_gene_SE_b_TRT_asthma^2 + Adj_gene_SE_b_CON_asthma^2)
z_score_asthma[z_score_asthma == 'NaN'] <- 0
z_score_asthma <- as.matrix(data.frame(z_score_asthma))

p_value_asthma <- pnorm(-abs(z_score_asthma)) * 2   #  p value of z-score 
Prediction_Log10_pval_asthma <- -log10(p_value_asthma)  ##  
Prediction_Log10_pval_asthma[Prediction_Log10_pval_asthma == 'NaN'] <- 0

dim(Prediction_Log10_pval_asthma)
typeof(Prediction_Log10_pval_asthma)
sum(diag(Prediction_Log10_pval_asthma))
sum(p_value_asthma < .05)

sig_diff_asthma <- as.matrix(Prediction_Log10_pval_asthma > -log10(1e-8))
dim(sig_diff_asthma)  
sum(diag(sig_diff_asthma))  # diagonal elements equal zero
sum(sig_diff_asthma == TRUE) #  1826 

##=========================================================================================================
## making dataset for graph in cyptoscape : Adj  treatment group : save estimated coefficients

nr <- ncol(sig_diff_asthma)
sign_trt_con <- matrix('', ncol = nr, nrow = nr)

Graph_asthma <- matrix(0, ncol = 8)
for (i in 1:nr)
{
    for (j in 1:nr)
    {
        if (sig_diff_asthma[j, i] == TRUE)
        {
            end_Gene <- colnames(sig_diff_asthma)[i]
            start_Gene <- rownames(sig_diff_asthma)[j]
            weight <- Prediction_Log10_pval_asthma[j, i]   # -log10
            diff_trt_con <- diff_beta_asthma[j, i]
            sign_trt <- Adj_gene_net_b_TRT_asthma[j, i]   # -log10
            sign_con <- Adj_gene_net_b_CON_asthma[j, i]
            # group : sign direction
            if (Adj_gene_net_b_TRT_asthma[j, i] == 0 & Adj_gene_net_b_CON_asthma[j, i]  == 0) {sign_trt_con <- 'ZZ'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  > 0 & Adj_gene_net_b_CON_asthma[j, i]  > 0) {sign_trt_con <- 'PP'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  > 0 & Adj_gene_net_b_CON_asthma[j, i]  < 0) {sign_trt_con <- 'PN'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  < 0 & Adj_gene_net_b_CON_asthma[j, i]  > 0) {sign_trt_con <- 'NP'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  < 0 & Adj_gene_net_b_CON_asthma[j, i]  < 0) {sign_trt_con <- 'NN'}
            else if (Adj_gene_net_b_TRT_asthma[j, i] == 0 & Adj_gene_net_b_CON_asthma[j, i]  > 0) {sign_trt_con <- 'ZP'}
            else if (Adj_gene_net_b_TRT_asthma[j, i] == 0 & Adj_gene_net_b_CON_asthma[j, i]  < 0) {sign_trt_con <- 'ZN'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  > 0 & Adj_gene_net_b_CON_asthma[j, i] == 0) {sign_trt_con <- 'PZ'}
            else if (Adj_gene_net_b_TRT_asthma[j, i]  < 0 & Adj_gene_net_b_CON_asthma[j, i] == 0) {sign_trt_con <- 'NZ'}
            Graph_asthma <- rbind(Graph_asthma, cbind(start_Gene, end_Gene, weight,  diff_trt_con, diff_trt_con, sign_trt, sign_con, sign_trt_con))
        }
    }
}

Graph_asthma <- Graph_asthma[-1, ]
dim(Graph_asthma)   # 1826 by 6
table(Graph_asthma[,8]) 

# table results 
# NP  NZ  PN  PP  PZ  ZN  ZP 
# 234 328 256   1 363 340 304 

sum(table(Graph_asthma[, 8]) ) # 1826

NN_asthma <- Graph_asthma[Graph_asthma[, 8] == 'NN', ] # 0
NP_asthma <- Graph_asthma[Graph_asthma[, 8] == 'NP', ] # 234
NZ_asthma <- Graph_asthma[Graph_asthma[, 8] == 'NZ', ] # 328
PN_asthma <- Graph_asthma[Graph_asthma[, 8] == 'PN', ] # 256
PP_asthma <- Graph_asthma[Graph_asthma[, 8] == 'PP', ] # 1
PZ_asthma <- Graph_asthma[Graph_asthma[, 8] == 'PZ', ] # 363
ZN_asthma <- Graph_asthma[Graph_asthma[, 8] == 'ZN', ] # 340
ZP_asthma <- Graph_asthma[Graph_asthma[, 8] == 'ZP', ] # 304

write.csv(Graph_asthma, "res/GraphAsthma.csv",    row.names=FALSE)
write.csv(NN_asthma,    "res/GraphAsthmaNN.csv", row.names=FALSE)
write.csv(NP_asthma,    "res/GraphAsthmaNP.csv", row.names=FALSE)
write.csv(NZ_asthma,    "res/GraphAsthmaNZ.csv", row.names=FALSE)
write.csv(PN_asthma,    "res/GraphAsthmaPN.csv", row.names=FALSE)
write.csv(PP_asthma,    "res/GraphAsthmaPP.csv", row.names=FALSE)
write.csv(PZ_asthma,    "res/GraphAsthmaPZ.csv", row.names=FALSE)
write.csv(ZN_asthma,    "res/GraphAsthmaZN.csv", row.names=FALSE)
write.csv(ZP_asthma,    "res/GraphAsthmaZP.csv", row.names=FALSE)

##==============confirm our significant genes with string DB====================================
sig_gene_names1 <- Graph_asthma[, 1]
sig_gene_names2 <- Graph_asthma[, 2]

unique_sig_gene_names <- unique(sig_gene_names1, sig_gene_names2)
length(unique_sig_gene_names)  # 1477

write.csv(unique_sig_gene_names, "res/unique_sig_gene_names.csv", row.names=FALSE)

sig_gene_names <- cbind(sig_gene_names1, sig_gene_names2 )

match_gene <- readCSVdata(pathTodata="res/match_gene.csv") #, sheet=1, verbose=FALSE, perl="c:/Research/R/asthma/perl.exe", header = TRUE)
dim(match_gene)

n_sig_gene <- length(sig_gene_names1)   #244
nr_match_gene <- nrow(match_gene)  # 61
match_string_sig_gene <- matrix("", nrow = n_sig_gene, ncol=2)

dim(match_string_sig_gene)

for (i in 1:n_sig_gene)
{
    for (j in 1: nr_match_gene)
    {
        if (((sig_gene_names[i, 1] == match_gene[ j, 1]) && (sig_gene_names[i, 2] == match_gene[ j, 2])) || 
            ((sig_gene_names[i, 1] == match_gene[ j, 2]) && (sig_gene_names[i, 1] == match_gene[ j, 1])))
        { match_string_sig_gene[i, 1] <- sig_gene_names[i, 1]
          match_string_sig_gene[i, 2] <- sig_gene_names[i, 2] }
    }
}

match_gene <- as.matrix(data.frame(match_gene), nrow= 39, ncol=2)
match_gene[, 1]

##==============================================================================================
## heat map for GRN of two groups

Lg10_pvalue <- Prediction_Log10_pval_asthma
sum(Lg10_pvalue <= 8 )
Lg10_pvalue[Lg10_pvalue <= 8] <- 0
sum(Lg10_pvalue !=0)  # 1826

Adj_TRT_asthma <- Adj_gene_net_b_TRT_asthma
Adj_CON_asthma <- Adj_gene_net_b_CON_asthma

Adj_TRT_asthma[Lg10_pvalue <= 8] <- 0
Adj_CON_asthma[Lg10_pvalue <= 8] <- 0
diff_asthma <- Adj_TRT_asthma - Adj_CON_asthma

max(c(Adj_gene_net_b_TRT_asthma, Adj_gene_net_b_CON_asthma) )  # 1.048338
min(c(Adj_gene_net_b_TRT_asthma, Adj_gene_net_b_CON_asthma) )  # -0.8278301

max(c(Adj_TRT_asthma, Adj_CON_asthma) )  # 0.868562
min(c(Adj_TRT_asthma, Adj_CON_asthma) )  # -0.795268

# color define
myCol <- c("black", colorRampPalette(c("green","darkgreen"))(18),"white", colorRampPalette(c("red","darkred"))(20)) 
length(myCol)  # 40

# defining breaks, the number of breask is one greater than color size
myBreaks <- seq(-1, 1, 0.05) 
length(myBreaks)  # 41

heatmap.2(Adj_gene_net_b_TRT_asthma, Colv=NA,
          col = myCol,        # using my colors
          breaks = myBreaks,  # using my breaks
          dendrogram = "none",  ## non dendograms
          cexRow=1, cexCol=1, key=FALSE,
          margins = c(2, 12),trace="none")
legend("topleft", fill = c("black","green","white","red", "black"),
       legend = c("-1", "-1 < to < 0", "0", "0< to <1", "1"), cex=1, horiz =TRUE)

heatmap.2(Adj_gene_net_b_CON_asthma, Colv=NA,
          col = myCol,        # using my colors
          breaks = myBreaks,  # using my breaks
          dendrogram = "none",  ## non dendograms
          cexRow=1, cexCol=1, key=FALSE,
          margins = c(2, 12),trace="none")
legend("topleft", fill = c("black","green","white","red", "black"),
       legend = c("-1", "-1 < to < 0", "0", "0< to <1", "1"), cex=1, horiz =TRUE)

# Build the palette and heatmap again using leveplot
pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
levelplot(Adj_TRT_asthma, xlab="", ylab="", col.regions=pal(41), cuts=40, at=seq(-1,1,0.05))

pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
levelplot(Adj_CON_asthma, xlab="", ylab="", col.regions=pal(41), cuts=40, at=seq(-1,1,0.05))

pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
levelplot(diff_asthma, xlab="", ylab="", col.regions=pal(41), cuts=40, at=seq(-1,1,0.05))


# again heatmap for subset from upper 150 of pvalue
select_value <- 150 
selected <- Lg10_pvalue <= select_value
Lg10_pvalue_sel <- Lg10_pvalue
Lg10_pvalue_sel[selected] <- 0
Lg10_pvalue_sel[Lg10_pvalue_sel == Inf] <- 500
Adj_TRT_asthma_sel <- Adj_TRT_asthma
Adj_TRT_asthma_sel[selected] <- 0
Adj_CON_asthma_sel <- Adj_CON_asthma
Adj_CON_asthma_sel[selected] <- 0

filteredRowCol <- colSums(Lg10_pvalue_sel) > 0
Lg10_pvalue_sel <- Lg10_pvalue_sel[filteredRowCol, filteredRowCol]
Adj_TRT_asthma_sel <- Adj_TRT_asthma_sel[filteredRowCol, filteredRowCol]
Adj_CON_asthma_sel <- Adj_CON_asthma_sel[filteredRowCol, filteredRowCol]
diff_sel = Adj_TRT_asthma_sel - Adj_CON_asthma_sel

# colnames(Lg10_pvalue_sel)[Lg10_pvalue_sel > 150]
# par(pin = c(4.5, 4.5), cex.lab = 1.2 )   # graph size pin=c(width, height)

Title <- "Heatmap for significant differential coefficients (# of Paire = 1826) "
pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(Lg10_pvalue_sel, xlab="", ylab="", col.regions=pal(51), cuts=50, at=seq(0,500,10), 
          scales=list(x=list(rot=90, cex = .6), y=list(cex=0.6)), 
          main = list(Title, font = 0.8, side=1,line=0.5) )

Title <- "Heatmap for significant Adjacency matrix of Treatment Group"
pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
myCol <- c(colorRampPalette(c("darkgreen","white"))(20), colorRampPalette(c("white","darkred"))(20)) 
levelplot(Adj_TRT_asthma_sel, xlab="", ylab="", col.regions=myCol, cuts=40, at=seq(-1,1,0.05), 
          scales=list(x=list(rot=90, cex = .6), y=list(cex=0.6)), 
          main = list(Title, font = 0.8, side=1,line=0.5) )

Title <- "Heatmap for significant Adjacency matrix of Control Group"
levelplot(Adj_CON_asthma_sel, xlab="", ylab="", col.regions=myCol, cuts=40, at=seq(-1,1,0.05), 
          scales=list(x=list(rot=90, cex = .6), y=list(cex=0.6)), 
          main = list(Title, font = 0.8, side=1,line=0.5) )


# ======================================================================================================
### ==================== Figure 5 : heatmap for upper 150 ===============================
#              NP  NZ  PN  PP  PZ  ZN  ZP 
#             234 328 256   1 363 340 304 

# ======== heatmap higher than 150 ============
 # NN_asthma    (Asthma = Negative, Control = Negative)"
 # NP_asthma    (Asthma = Negative, Control = Positive)"
 # NZ_asthma    (Asthma = Negative, Control = Zero)"
 # PN_asthma    (Asthma = Positive, Control = Negative)"
 # PP_asthma    (Asthma = Positive, Control = Positive)"
 # PZ_asthma    (Asthma = Positive, Control = Zero)"
 # ZN_asthma    (Asthma = Zero, Control = Negative)"
 # ZP_asthma    (Asthma = Zero, Control = Positive)"

select_value <- 150 
setwd("C:/Users/Phoebe/Dropbox/IJDMB_DiffGRN/IJDMB 2018/Figures")

##==== NZ graph : Figure 5
asthma <- NZ_asthma
tmp <- asthma[as.numeric(asthma[, 3]) > select_value ,  ]
# 
# Title <- "Heatmap (Asthma = Zero, Control = Negative)"

rnames <- unique(tmp[, 1])
cnames <- unique(tmp[, 2])

TRT <- Adj_TRT_asthma[rnames, cnames] 
CON <- Adj_CON_asthma[rnames, cnames] 
log10_pavlue <- Prediction_Log10_pval_asthma[rnames, cnames] 

TRT <-                   TRT * (TRT < 0) * (CON == 0)  # NZ
CON <-                   CON * (TRT < 0) * (CON == 0)  # 
log10_pavlue <- log10_pavlue * (TRT < 0) * (CON == 0)  # 

selected <- log10_pavlue <= select_value
log10_pavlue_sel <- log10_pavlue
log10_pavlue_sel[selected] <- 0
log10_pavlue_sel[log10_pavlue_sel == Inf] <- 500

TRT_sel <- TRT
TRT_sel[selected] <- 0
CON_sel <- CON
CON_sel[selected] <- 0

pdf(file = "./heatmap_NZ_GT_150.pdf", width=4, height=4)
par(mar=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(log10_pavlue_sel, xlab="", ylab="", col.regions=pal(51), cuts=50, at=seq(0,500,10), 
          scales=list(x=list(rot=90, cex = .8), y=list(cex=0.8)) )# , aspect = 0.25 )
#main = list(Title, font = 0.8, side=1,line=0.5) )
dev.off()


##==== PZ graph : Figure 5
asthma <- PZ_asthma
tmp <- asthma[as.numeric(asthma[, 3]) > select_value ,  ]
# 
# Title <- "Heatmap (Asthma = Zero, Control = Negative)"

rnames <- unique(tmp[, 1])
cnames <- unique(tmp[, 2])

TRT <- Adj_TRT_asthma[rnames, cnames] 
CON <- Adj_CON_asthma[rnames, cnames] 
log10_pavlue <- Prediction_Log10_pval_asthma[rnames, cnames] 

TRT <-                   TRT * (TRT > 0) * (CON == 0)  # NZ
CON <-                   CON * (TRT > 0) * (CON == 0)  # 
log10_pavlue <- log10_pavlue * (TRT > 0) * (CON == 0)  # 

selected <- log10_pavlue <= select_value
log10_pavlue_sel <- log10_pavlue
log10_pavlue_sel[selected] <- 0
log10_pavlue_sel[log10_pavlue_sel == Inf] <- 500

TRT_sel <- TRT
TRT_sel[selected] <- 0
CON_sel <- CON
CON_sel[selected] <- 0

pdf(file = "./heatmap_PZ_GT_150.pdf", width=4.5, height=4)
par(mar=c(0,0,0,0), mai=c(0,0,0,0))
pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(log10_pavlue_sel, xlab="", ylab="", col.regions=pal(51), cuts=50, at=seq(0,500,10), 
          scales=list(x=list(rot=90, cex = .8), y=list(cex=0.8)) )#, aspect = 0.25 )
#main = list(Title, font = 0.8, side=1,line=0.5) )
dev.off()


##==== ZN graph : Figure 5
asthma <- ZN_asthma
tmp <- asthma[as.numeric(asthma[, 3]) > select_value ,  ]
# 
# Title <- "Heatmap (Asthma = Zero, Control = Negative)"

rnames <- unique(tmp[, 1])
cnames <- unique(tmp[, 2])

TRT <- Adj_TRT_asthma[rnames, cnames] 
CON <- Adj_CON_asthma[rnames, cnames] 
log10_pavlue <- Prediction_Log10_pval_asthma[rnames, cnames] 

TRT <-                   TRT * (TRT == 0) * (CON < 0)  # ZP
CON <-                   CON * (TRT == 0) * (CON < 0)  # 
log10_pavlue <- log10_pavlue * (TRT == 0) * (CON < 0)  # 

selected <- log10_pavlue <= select_value
log10_pavlue_sel <- log10_pavlue
log10_pavlue_sel[selected] <- 0
log10_pavlue_sel[log10_pavlue_sel == Inf] <- 500

TRT_sel <- TRT
TRT_sel[selected] <- 0
CON_sel <- CON
CON_sel[selected] <- 0

pdf(file = "./heatmap_ZN_GT_150.pdf", width=8, height=4)
par(mar=c(0,0,0,0))
pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(log10_pavlue_sel, xlab="", ylab="", col.regions=pal(51), cuts=50, at=seq(0,500,10), 
          scales=list(x=list(rot=90, cex = .8), y=list(cex=0.8)) )#, aspect = 0.25 )
          #main = list(Title, font = 0.8, side=1,line=0.5) )
dev.off()

##==== ZP graph  : Figure 5

asthma <- ZP_asthma
tmp <- asthma[as.numeric(asthma[, 3]) > select_value ,  ]
# 
# Title <- "Heatmap (Asthma = Zero, Control = Negative)"

rnames <- unique(tmp[, 1])
cnames <- unique(tmp[, 2])

TRT <- Adj_TRT_asthma[rnames, cnames] 
CON <- Adj_CON_asthma[rnames, cnames] 
log10_pavlue <- Prediction_Log10_pval_asthma[rnames, cnames] 

TRT <-                   TRT * (TRT == 0) * (CON > 0)  # ZP
CON <-                   CON * (TRT == 0) * (CON > 0)  # 
log10_pavlue <- log10_pavlue * (TRT == 0) * (CON > 0)  # 

selected <- log10_pavlue <= select_value
log10_pavlue_sel <- log10_pavlue
log10_pavlue_sel[selected] <- 0
log10_pavlue_sel[log10_pavlue_sel == Inf] <- 500

TRT_sel <- TRT
TRT_sel[selected] <- 0
CON_sel <- CON
CON_sel[selected] <- 0

pdf(file = "./heatmap_ZP_GT_150.pdf", width=8, height=4)
par(mar=c(0,0,0,0))

pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(log10_pavlue_sel, xlab="", ylab="", col.regions=pal(51), cuts=50, at=seq(0,500,10), 
          scales=list(x=list(rot=90, cex = .8), y=list(cex=0.8)) )#, aspect = 0.25 )
#main = list(Title, font = 0.8, side=1,line=0.5) )
dev.off()


# ======================================================================================================

# pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
# myCol <- c(colorRampPalette(c("darkgreen","white"))(100), colorRampPalette(c("white","darkred"))(100)) 
# levelplot(TRT_sel, xlab="", ylab="", col.regions=myCol, cuts=20, at=seq(-0.5,0.5,0.005))
# levelplot(CON_sel, xlab="", ylab="", col.regions=myCol, cuts=20, at=seq(-0.5,0.5,0.005))

## save file with upper 150 of pvalue=============================================================================
# == upper rank  ============
# NP NZ PZ ZN ZP 
# 1 34 49 89 71 

Upper_all_asthma <- Graph_asthma[as.numeric(Graph_asthma[, 3]) > 150, ]
Upper_NN_asthma <- NN_asthma[as.numeric(NN_asthma[, 3]) > 150, ]
Upper_NP_asthma <- NP_asthma[as.numeric(NP_asthma[, 3]) > 150, ]
Upper_NZ_asthma <- NZ_asthma[as.numeric(NZ_asthma[, 3]) > 150, ]
Upper_PN_asthma <- PN_asthma[as.numeric(PN_asthma[, 3]) > 150, ]
Upper_PP_asthma <- PP_asthma[as.numeric(PP_asthma[, 3]) > 150, ]
Upper_PZ_asthma <- PZ_asthma[as.numeric(PZ_asthma[, 3]) > 150, ]
Upper_ZN_asthma <- ZN_asthma[as.numeric(ZN_asthma[, 3]) > 150, ]
Upper_ZP_asthma <- ZP_asthma[as.numeric(ZP_asthma[, 3]) > 150, ]

write.csv(Upper_all_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_all_asthma.csv", row.names=FALSE)
write.csv(Upper_NN_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_NN_asthma.csv", row.names=FALSE)
write.csv(Upper_NP_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_Np_asthma.csv", row.names=FALSE)
write.csv(Upper_NZ_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_NZ_asthma.csv", row.names=FALSE)
write.csv(Upper_PN_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_PN_asthma.csv", row.names=FALSE)
write.csv(Upper_PP_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_PP_asthma.csv", row.names=FALSE)
write.csv(Upper_PZ_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_PZ_asthma.csv", row.names=FALSE)
write.csv(Upper_ZN_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_ZN_asthma.csv", row.names=FALSE)
write.csv(Upper_ZP_asthma, "C:/Research/R/asthma/src/result/Upper_pairs_GT_150/Upper_ZP_asthma.csv", row.names=FALSE)


##=============Pathway_gene_Adjacency matrix for asthma =====================================================
## pathway for gene
#setwd("C:/Users/Phoebe/Dropbox/Youngsoon Kim/Publications/Scientific Reports/Src/result")
setwd("C:/Research/R/Publications/Scientific Reports/Src/result")

## include 4 kinds of pathway : "KEGG", "REACTOME", "PID", "BIOCARTA"
adj_pathway_gene <- readCSVdata(pathTodata=".\\adj_pathway_gene.csv")
rownames(adj_pathway_gene) <- adj_pathway_gene[, 1]
adj_pathway_gene <- adj_pathway_gene[, -1]  # 1273 by 12042

path_gene_names <- colnames(adj_pathway_gene)
length(path_gene_names)   # 12042

# Graph_asthma <- readCSVdata("C:/Research/R/asthma/src/result/network_graph/Graph_asthma.csv")

sig_gene_names1 <- Graph_asthma[, 1]
sig_gene_names2 <- Graph_asthma[, 2]

unique_sig_gene_names_asthma <- unique(sig_gene_names1, sig_gene_names2)
length(unique_sig_gene_names_asthma)  # 1477

adj_pathway_gene_asthma <- Choose_same_gene(adj_pathway_gene, unique_sig_gene_names_asthma)

dim(adj_pathway_gene_asthma)  # 1273 by 148

length(rowSums(adj_pathway_gene_asthma) > 0)

upper_pathway_N_includedGenes <- tail(sort(rowSums(adj_pathway_gene_asthma)), 10)

upper_pathway_name <- names(upper_pathway_N_includedGenes )   # 10

upper_adj_pathway_gene_asthma <- Choose_same_sample(adj_pathway_gene_asthma, upper_pathway_name)
dim(upper_adj_pathway_gene_asthma)   # 10 by 148
    
n_upper_pathway <- nrow(upper_adj_pathway_gene_asthma)

for (i in 1 : n_upper_pathway)
{
    print(colnames(upper_adj_pathway_gene_asthma[i, which(upper_adj_pathway_gene_asthma[i, ] != 0)]))
}

for (i in 1 : n_upper_pathway)
{
    print(rownames(upper_adj_pathway_gene_asthma[i, which(upper_adj_pathway_gene_asthma[i, ] != 0)]))
}
# until here, we use print result to table2, 




# upper_pathway_gene_asthma <- cbind(upper_pathway_name, upper_pathway_N_includedGenes, geneList_temp)

write.csv(upper_pathway_N_includedGenes, "C:/Research/R/asthma/src/result/pathway_gene/Upper_pathway_4gene.csv", row.names = TRUE)
 
## computer upper pathway size using "upper_pathway_name"
adj_pathway_gene <- readCSVdata(pathTodata=".\\adj_pathway_gene.csv")
rownames(adj_pathway_gene) <- adj_pathway_gene[, 1]
adj_pathway_gene <- adj_pathway_gene[, -1]  # 1273 by 12042

tmp_rowSums <- rowSums(adj_pathway_gene)
tmp_rowSums[names(tmp_rowSums) == upper_pathway_name]
tmp_rowSums[which(!is.na(match(names(tmp_rowSums), upper_pathway_name)))]

pathwaygeneNames <- gene_names[which(colSums(adj_pathway_gene) > 0)] 
length(pathwaygeneNames)   # 6103 genes

InPathGene <- which(!is.na(match(colnames(Adj_gene_net_b_TRT_asthma) , pathwaygeneNames )))   #  91
NonPathGene <- (which(is.na(match(colnames(Adj_gene_net_b_TRT_asthma) , pathwaygeneNames )))) # 123

InPath_Adj_gene_net_b_TRT_asthma <- Adj_gene_net_b_TRT_asthma[InPathGene,  InPathGene]  # 91 by 91
InPath_Adj_gene_net_b_CON_asthma <- Adj_gene_net_b_CON_asthma[InPathGene,  InPathGene]  # 91 by 91

Path_gene_names <- colnames(InPath_Adj_gene_net_b_TRT_asthma)

adj_pathway_gene_c <- adj_pathway_gene[, Path_gene_names] 
dim(adj_pathway_gene_c )
hist(rowSums(adj_pathway_gene_c))

numPathway <- nrow(adj_pathway_gene_c)
order_gene <- vector(length=1)
order_pathway_num <- vector(length=1)
order_gene_name <- vector(length=1)

for (i in 1 : numPathway)
{
    Inpathway <- which(adj_pathway_gene_c[i, ] == 1)
    temp_names <- paste("P", i, "_", colnames(adj_pathway_gene_c)[Inpathway], sep='')
    
    if ((i == 1) & (length(Inpathway) != 0))
    {
        order_gene <- Inpathway
        order_gene_name <- temp_names
        order_pathway_num <- rep(i, length = length(Inpathway))
    }
    else if ((i == 1) & (length(Inpathway) == 0))
    {
        next
    }
    else if ((i > 1) & (length(Inpathway) != 0)) 
    {
        order_gene <- append(order_gene, Inpathway)
        order_gene_name <- append(order_gene_name, temp_names)
        order_pathway_num <- append(order_pathway_num,  rep(i, length = length(Inpathway)))
    }
    else if ((i > 1) & (length(Inpathway) == 0))  
    {
        next  
    }
}
length(order_gene) # 777
length(order_gene_name)  #777
length(order_pathway_num)  #777


# make Adj_network considering pathway
Path_Adj_gene_net_b_TRT_asthma <- Adj_gene_net_b_TRT_asthma[order_gene,  order_gene]  # 777 by 777
Path_Adj_gene_net_b_CON_asthma <- Adj_gene_net_b_CON_asthma[order_gene,  order_gene]  # 777 by 777
colnames(Path_Adj_gene_net_b_TRT_asthma) <- order_gene_name
rownames(Path_Adj_gene_net_b_TRT_asthma) <- order_gene_name
colnames(Path_Adj_gene_net_b_CON_asthma) <- order_gene_name
rownames(Path_Adj_gene_net_b_CON_asthma) <- order_gene_name


## include only "REACTOME" pathway
adj_pathway_gene_REACTOME <- readCSVdata(pathTodata=".\\adj_pathway_gene_REACTOME.csv")
rownames(adj_pathway_gene_REACTOME) <- adj_pathway_gene_REACTOME[, 1]
adj_pathway_gene_REACTOME <- adj_pathway_gene_REACTOME[, -1]  
dim(adj_pathway_gene_REACTOME) # 674 by 12042
gene_names <- colnames(adj_pathway_gene_REACTOME)

# select all pathways gene names
pathwaygeneNames <- gene_names[which(colSums(adj_pathway_gene_REACTOME) > 0)] 
length(pathwaygeneNames)   # 4632 genes

InPathGene <- which(!is.na(match(colnames(Adj_gene_net_b_TRT_asthma) , pathwaygeneNames )))   #  63
NonPathGene <- (which(is.na(match(colnames(Adj_gene_net_b_TRT_asthma) , pathwaygeneNames )))) # 151
length(InPathGene)
length(NonPathGene)

InPath_Adj_gene_net_b_TRT_asthma <- Adj_gene_net_b_TRT_asthma[InPathGene,  InPathGene]  # 63 by 63
dim(InPath_Adj_gene_net_b_TRT_asthma)
InPath_Adj_gene_net_b_CON_asthma <- Adj_gene_net_b_CON_asthma[InPathGene,  InPathGene]  # 63 by 63
dim(InPath_Adj_gene_net_b_CON_asthma)

Path_gene_names <- colnames(InPath_Adj_gene_net_b_TRT_asthma)

adj_pathway_gene_c <- adj_pathway_gene_REACTOME[, Path_gene_names]   # 674 by 63
dim(adj_pathway_gene_c )

numPathway <- nrow(adj_pathway_gene_c)
order_gene <- vector(length=1)
order_pathway_num <- vector(length=1)
order_gene_name <- vector(length=1)

for (i in 1 : numPathway)
{
    Inpathway <- which(adj_pathway_gene_c[i, ] == 1)
    temp_names <- paste("P", i, "_", colnames(adj_pathway_gene_c)[Inpathway], sep='')
    print(i)
    print(temp_names)
    
    if ((i == 1) & (length(Inpathway) != 0))
    {
        order_gene <- Inpathway
        order_gene_name <- temp_names
        order_pathway_num <- rep(i, length = length(Inpathway))
    }
    else if ((i == 1) & (length(Inpathway) == 0))
    {
        next
    }
     else if ((i > 1) & (length(Inpathway) != 0)) 
    {
        order_gene <- append(order_gene, Inpathway)
        order_gene_name <- append(order_gene_name, temp_names)
        order_pathway_num <- append(order_pathway_num,  rep(i, length = length(Inpathway)))
    }
    else if ((i > 1) & (length(Inpathway) == 0))  
    {
        next  
    }
}
order_gene <- order_gene[-1]
order_gene_name <- order_gene_name[-1]
order_pathway_num <- order_pathway_num[-1]

length(order_gene) # 416
length(order_gene_name)  # 416
length(order_pathway_num)  # 416

# make Adj_network considering pathway
Path_Adj_gene_net_b_TRT_asthma <- Adj_gene_net_b_TRT_asthma[order_gene,  order_gene]  # 416 by 416
Path_Adj_gene_net_b_CON_asthma <- Adj_gene_net_b_CON_asthma[order_gene,  order_gene]  # 416 by 416
colnames(Path_Adj_gene_net_b_TRT_asthma) <- order_gene_name
rownames(Path_Adj_gene_net_b_TRT_asthma) <- order_gene_name
colnames(Path_Adj_gene_net_b_CON_asthma) <- order_gene_name
rownames(Path_Adj_gene_net_b_CON_asthma) <- order_gene_name

## heat map for DiffGRN considering pathway 
# Build the palette and plot it
pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
levelplot(Path_Adj_gene_net_b_TRT_asthma, xlab="", ylab="", col.regions=pal(41), cuts=40, at=seq(-1,1,0.05))

# Build the palette and plot it
pal <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
levelplot(Path_Adj_gene_net_b_CON_asthma, xlab="", ylab="", col.regions=pal(41), cuts=40, at=seq(-1,1,0.05))




## ======================================================================================================
## compute differential network our DINGO for asthma
n_r1 <- nrow(Norm_Asthma_TRT)

Group_asthma <- c(rep(1, n_r1), rep(2, n_r1))
dat_dingo_asthma <- rbind(Norm_Asthma_CON, Norm_Asthma_TRT)
fit_dingo_asthma <- dingo(dat_dingo_asthma, Group_asthma, B= 100)

Dingo_Log10_pval_asthma <- fit_dingo_asthma$p
Dingo_Log10_pval_asthma[is.nan(Dingo_Log10_pval_asthma)] <- 1
Prediction_Dingo_Log10_pval_asthma <- -log10(Dingo_Log10_pval_asthma)

sigNum_dingo <- which((Prediction_Dingo_Log10_pval_asthma > 5))

sig_genename_dingo <- fit_dingo_asthma$genepair[sigNum_dingo,]
sign_diff_dingo <- fit_dingo_asthma$diff.score[sigNum_dingo]
sig_pval <- fit_dingo_asthma$p.val[sigNum_dingo]
TRT_CORR <- fit_dingo_asthma$R1[sigNum_dingo]
CON_CORR <- fit_dingo_asthma$R2[sigNum_dingo]

dim(sig_genename_dingo)  # 61 by 2

aa <- sig_genename_dingo[, 1]
bb<- sig_genename_dingo[, 2]
sig_con_mean_a <- apply(Norm_Asthma_CON[, aa], 2, mean)
sig_con_mean_b <- apply(Norm_Asthma_CON[, bb], 2, mean)
sig_con_mean <- sig_con_mean_a / sig_con_mean_b
sig_con_mean
sig_trt_mean_a <- apply(Norm_Asthma_TRT[, aa], 2, mean)
sig_trt_mean_b <- apply(Norm_Asthma_TRT[, bb], 2, mean)
sig_trt_mean <- sig_trt_mean_a / sig_trt_mean_b
sig_trt_mean

dingo_result <- cbind(sig_genename_dingo, sign_diff_dingo, sig_pval, TRT_CORR, CON_CORR,sig_con_mean, sig_trt_mean)
write.csv(dingo_result,    "C:/Research/R/asthma/src/result/dingo_asthma.csv", row.names=FALSE)
# ======================================================================================================