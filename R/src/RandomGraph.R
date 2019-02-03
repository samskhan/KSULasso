
#Generate random graph for gene expression according to the Erdos-Renyi model.
#the arguments are, numbr of samples, number of genes, and the Sparce rate.
Rnetwork_GRN <- function(n, p, x)
{
    flag <- FALSE
    while(flag == FALSE)
    { 
        graph <- erdos.renyi.game(p, x, directed = TRUE, loops = FALSE)
        degree.distribution(graph)
        flag <- is.dag(graph)
    } 
    # Ground truth of gene expression
    Adj_Rnetwork <- as_adjacency_matrix(graph)
    #E has gaussian distribution with mean=0 and sd=0.01
    E <- matrix(rnorm(n * p, mean = 0, sd = 0.01), n, p)
    Z <- as.matrix(E %*% solve(diag(p) - Adj_Rnetwork))  
    # T <- scale(Z)
    return(list(Adj_Rnetwork, Z))
} 

##==============================================================================
#function to find confusion matrix
Confusion <- function(GT, Adj_Predition)
{
    TP <- sum(GT != 0 & Adj_Predition != 0)
    FP <- sum(GT == 0 & Adj_Predition != 0)
    FN <- sum(GT != 0 & Adj_Predition == 0)
    TN <- sum(GT == 0 & Adj_Predition == 0) - length(diag(GT))

    return(c(TP, FP, FN, TN))
}

#==============================================================================
# AUC and FDR
AUC_FDR_From_confusion <- function(threshold_pvalue, GT_diff, Log10_pvalue) 
{
    n_row <- length(threshold_pvalue)
    ROC_data <- matrix(0, nrow = n_row, ncol = 6)
    ROC_data[, 1] <- threshold_pvalue
    colnames(ROC_data) <- c('pvalue', 'TPR', 'PPV', 'FPR', 'AUC', 'FDR')
    n_r <- nrow(GT_diff)
    
    for (i in threshold_pvalue)
    {
        ConFusion <- Confusion(GT_diff, Log10_pvalue > i )
        
        TP <- ConFusion[1]
        FP <- ConFusion[2]
        FN <- ConFusion[3]
        TN <- ConFusion[4]
        
        ROC_data[(ROC_data[, 1] == i), 2] <- (TP / (TP + FN))  # TPR, sensitivity for sd = .01
        ROC_data[(ROC_data[, 1] == i), 3] <- (TP / (TP + FP))  # PPV, precision   for sd = .01
        ROC_data[(ROC_data[, 1] == i), 4] <- (FP / (TN + FP))   #FPR
        
        GT_subset <- as.vector(GT_diff[-seq(1,n_r^2,n_r+1)])  #  remove diagonal elements
        Log10_pvalue_subset <- as.vector(Log10_pvalue[-seq(1,n_r^2,n_r+1)])
        roc_obj <- roc(GT_subset, Log10_pvalue_subset)
        ROC_data[(ROC_data[, 1] == i), 5] <- auc(roc_obj)
        
        # FDR <- FP/(TP + FP) # original equation
        length_diff_FDR <- sum(ConFusion)
        FDR_value <-  FP / length_diff_FDR  #simulation data included population truth. denominator sames # of comparison pairs
        ROC_data[(ROC_data[, 1] == i), 6] <- FDR_value
    }
    return(ROC_data)
}

#==============================================================================
# FDR
FDR_From_confusion <- function(threshold_pvalue, GT_diff, Log10_pvalue) 
{

    n_row <- length(threshold_pvalue)
    FDR <- matrix(0, nrow = n_row, ncol=2)
    colnames(FDR) <- c('pvalue', 'FDR')
    FDR[, 1] <- threshold_pvalue
    
    for (i in threshold_pvalue)
    {
        ConFusion_FDR <- Confusion(GT_diff, Log10_pvalue > i)
        TP <- ConFusion_FDR[1]
        FP <- ConFusion_FDR[2]
        FN <- ConFusion_FDR[3]
        TN <- ConFusion_FDR[4]
        # TPR_DGRN <- FP / (FP + TN)  # sensitivity for sd = .01
        
        # FDR <- FP/(TP + FP) # original equation
        length_diff_FDR <- sum(ConFusion_FDR)
        FDR_value <- FP/length_diff_FDR  #simulation data included population truth. denominator sames # of comparison pairs
        FDR[(FDR[, 1] == i), 2] <-  FDR_value
    }
    return(FDR)
}
 
