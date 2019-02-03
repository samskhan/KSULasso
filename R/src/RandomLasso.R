
####### !!attention : feature name have to not included '-'. So feature names which are in source data are changed by '_'. 
##==============================================================================
# read csv file
readCSVdata <- function (pathTodata )
{
    data <- read.csv(pathTodata, sep = ',', header = TRUE)
}
##==============================================================================
# read XLS file
readXLSdata <- function (pathTodata )
{
    data <- read.xlsx(pathTodata, header = TRUE)
}
##==============================================================================
# determine sample size
Determin_Feature_size <- function(Data, C_I = 0.95, P = 0.5, S_E = 0.05, REP = 10) 
{
    # Determining sample size of finite population for RANDOM LASSO
    numGenes <- ncol(Data) # population size, i.e. the numer of genes
    z_score <- qnorm(C_I, 0, 1) # z score from N(0, 1) 
    p <- P # proportion of population, unknown
    e <- S_E # sampling error, we can determine by ourselves.
    n0 <- (z_score^2 * p * (1 - p)) / e^2
    n  <- floor((n0 * numGenes) / (n0 + (numGenes - 1))) + 1 # sample size of 90% confidence interval
    
    q1 = n   # feature size for step1 in RandomLASSO
    q2 = n   # feature size for step2 in RandomLASSO
    rep_Boostrap = round(numGenes / n * REP) # repeated numer of bootstraping 
                                             # a gene to be included in overal boostraping on average REP times
    size <- list(q1, q2, rep_Boostrap)
    return(size)
}
##==============================================================================
# Random Lasso with selction weight
RandomLasso <- function(TRT_Data, CON_Data, Importance_TRT = Importance_weight_TRT, Importance_CON = Importance_weight_CON,
                        NumOfFeatures = q1, repeat_Boostrapping = B, PVAL = 0.05, STEP2) 
{    
    expGene_names <- colnames(TRT_Data)  # 3139, same between Ast and con
    numGenes  <- ncol(TRT_Data)          # number of genes
    
    Adj_gene_net_b_TRT <- matrix(0, nrow = numGenes, ncol = numGenes) # beta hat matrix 
    rownames(Adj_gene_net_b_TRT) <- expGene_names
    colnames(Adj_gene_net_b_TRT) <- expGene_names
    
    Adj_gene_net_b_CON <- matrix(0, nrow = numGenes, ncol = numGenes) # beta hat matrix 
    rownames(Adj_gene_net_b_CON) <- expGene_names
    colnames(Adj_gene_net_b_CON) <- expGene_names
    
    Adj_gene_SE_b_TRT <- matrix(0, nrow = numGenes, ncol = numGenes) # beta hat matrix 
    rownames(Adj_gene_SE_b_TRT) <- expGene_names
    colnames(Adj_gene_SE_b_TRT) <- expGene_names
    
    Adj_gene_SE_b_CON <- matrix(0, nrow = numGenes, ncol = numGenes) # beta hat matrix 
    rownames(Adj_gene_SE_b_CON) <- expGene_names
    colnames(Adj_gene_SE_b_CON) <- expGene_names

    for (i in 1 : numGenes)  # i is index of dependent variable, 
    {
        print("i= ")
        print(i)
        ## data preparation
        # AST data
        X_TRT <- as.matrix(TRT_Data)
        y_TRT <- X_TRT[, i]
        X_TRT[, i] <- 1 # intercept
        
        # CON data
        X_CON <- as.matrix(CON_Data)
        y_CON <- X_CON[, i]
        X_CON[, i] <- 1 # intercept
        
        # step 1 of Random LASSO < making improtance measures of genes >
        ## call LASSO ftn for step 1 ( ) : how to process W_p and importance 
        Adj_temp_b_TRT <- matrix(0, nrow = numGenes, ncol = repeat_Boostrapping)  # temporary storage for coefficient, AST
        rownames(Adj_temp_b_TRT) <- expGene_names
        
        Adj_temp_b_CON <- matrix(0, nrow = numGenes, ncol = repeat_Boostrapping)  # temporary storage for coefficient, CON
        rownames(Adj_temp_b_CON) <- expGene_names
        
        Adj_temp_SE_TRT <- matrix(0, nrow = numGenes, ncol = repeat_Boostrapping)  # temporary storage for coefficient, AST
        rownames(Adj_temp_SE_TRT) <- expGene_names
        
        Adj_temp_SE_CON <- matrix(0, nrow = numGenes, ncol = repeat_Boostrapping)  # temporary storage for coefficient, CON
        rownames(Adj_temp_SE_CON) <- expGene_names
        
        for (j in 1 : repeat_Boostrapping)  # For ith dependent gene, randomly repeat lasso method of B times.
        {
            print("j= ")
            print(j)

            # select genes by sampling with selection probability
            ## gene selection with gene selection weight
            if (STEP2 == FALSE) 
            { 
                select_prob <- 1/numGenes + Importance_TRT[, i]  ## In setp1, Importance = zero of two groups
                sampleGene <- sample(expGene_names[-i], NumOfFeatures, replace = FALSE, prob = select_prob[-i])

                XX_TRT <- Choose_same_gene(X_TRT, sampleGene)
                XX_TRT <- cbind(1, XX_TRT)  # add intercept first column in xx matrix
                colnames(XX_TRT)[1] <- 'Y'
                
                XX_CON <- Choose_same_gene(X_CON, sampleGene)
                XX_CON <- cbind(1, XX_CON)  # add intercept first column in xx matrix
                colnames(XX_CON)[1] <- 'Y'
                
            } else if (STEP2 == TRUE)
            {
                select_prob_TRT <- 1/numGenes + Importance_TRT[, i]
                sampleGene_TRT <- sample(expGene_names[-i], NumOfFeatures, replace = FALSE, prob = select_prob_TRT[-i])

                XX_TRT <- Choose_same_gene(X_TRT, sampleGene_TRT)
                XX_TRT <- cbind(1, XX_TRT)  # add intercept first column in xx matrix
                colnames(XX_TRT)[1] <- 'Y'

                select_prob_CON <- 1/numGenes + Importance_CON[, i]
                sampleGene_CON <- sample(expGene_names[-i], NumOfFeatures, replace = FALSE, prob = select_prob_CON[-i])

                XX_CON <- Choose_same_gene(X_CON, sampleGene_CON)
                XX_CON <- cbind(1, XX_CON)  # add intercept first column in xx matrix
                colnames(XX_CON)[1] <- 'Y'
            } 
                
            ## Fit in LASSO model : AST  from here
            XX <- XX_TRT
            y  <- y_TRT
            Adj_temp_TRT <- LarsForRandomLASSO(XX, y, expGene_names, j, Adj_temp_b_TRT, Adj_temp_SE_TRT, PVAL = 0.01, STEP2)
            Adj_temp_b_TRT[, j]  <- Adj_temp_TRT[[1]][, j]
            Adj_temp_SE_TRT[, j] <- Adj_temp_TRT[[2]][, j]

            XX <- XX_CON
            y  <- y_CON
            Adj_temp_CON <- LarsForRandomLASSO(XX, y, expGene_names, j, Adj_temp_b_CON, Adj_temp_SE_CON, PVAL = 0.01, STEP2)
            Adj_temp_b_CON[, j]  <- Adj_temp_CON[[1]][, j]
            Adj_temp_SE_CON[, j] <- Adj_temp_CON[[2]][, j]

        } # end of for bootstrapping    to here  AST and CON respectively repeat
        
        # average the coefficient 
        nonzero_repeat_Boostrapping1 <- rowSums(Adj_temp_b_TRT != 0)
        nonzero_repeat_Boostrapping1[nonzero_repeat_Boostrapping1 == 0] = 1
        print(nonzero_repeat_Boostrapping1)
        nonzero_repeat_Boostrapping2 <- rowSums(Adj_temp_b_CON != 0)
        nonzero_repeat_Boostrapping2[nonzero_repeat_Boostrapping2 == 0] = 1
        print(nonzero_repeat_Boostrapping2)
        
        Adj_gene_net_b_TRT[, i] <- rowSums(Adj_temp_b_TRT)  / nonzero_repeat_Boostrapping1
        Adj_gene_SE_b_TRT[, i]  <- sqrt(rowSums(Adj_temp_SE_TRT^2)  / nonzero_repeat_Boostrapping1 )
        Adj_gene_net_b_CON[, i] <- rowSums(Adj_temp_b_CON)  / nonzero_repeat_Boostrapping2
        Adj_gene_SE_b_CON[, i]  <- sqrt(rowSums(Adj_temp_SE_CON^2)  / nonzero_repeat_Boostrapping2)
    }# end of i, every column feature  

    Adj_gene_net_b <- list(Adj_gene_net_b_TRT, Adj_gene_SE_b_TRT, Adj_gene_net_b_CON, Adj_gene_SE_b_CON)    
    return(Adj_gene_net_b) 
} # end of randomLasso ftn  
##==============================================================================
# fit lasso
LarsForRandomLASSO <- function(XX, y, expGene_names, j, Adj_temp_b,  Adj_temp_SE, PVAL, STEP2)
{
    las       <- lars(as.matrix(XX), y, type = "lasso", use.Gram = FALSE)
    cvlas     <- cv.lars(as.matrix(XX), y, type = "lasso", plot.it = FALSE, use.Gram = FALSE)
    frac      <- cvlas$index[which.min(cvlas$cv)]
    las.coef  <- predict.lars(las, type = "coefficients", mode = "fraction", s = frac)
    hat_beta  <- las.coef$coefficients

    idx <- which(hat_beta != 0)

    if (length(idx) == 0) 
    { 
        Adj_temp_b[, j] <- 0
        Adj_temp_SE[, j] <- 0

    } else if (STEP2 == FALSE)
    {
        gene_id <- match(names(hat_beta[-1]), expGene_names)  # 1th is intercept, so remove
        Adj_temp_b[gene_id, j] <- hat_beta[-1]
        Adj_temp_SE[gene_id, j] <- 0 
        print("fighting")
    } else if (STEP2 == TRUE)
    {
        # find lm for nonzero columns of gene expression after applying lasso for
        # each gene i and consider those columns that corresponging p-value is less than 0.01
        
        nonzero_xx <- Choose_same_gene(XX, names(idx)) #nonZeroCoeff 
        print("There is nonzero xx")
        
        # linear regression for p-value
        fit.lm    <- lm(y ~. , data = data.frame(nonzero_xx))

        pval <- summary(fit.lm)$coefficients[, 4]
        pval <- pval[-1]
        sig_gene <- which(pval < PVAL)

        if (length(sig_gene) == 0 )
        {
            Adj_temp_b[, j] <- 0
            Adj_temp_SE[, j] <- 0
        } else
        {
            print("There is significant xx")
            gene_id <- match(names(sig_gene), expGene_names)
            Adj_temp_b[gene_id, j]  <- summary(fit.lm)$coefficients[sig_gene + 1, 1]   # hat beta
            Adj_temp_SE[gene_id, j] <- summary(fit.lm)$coefficients[sig_gene + 1, 2]  # SE of hat beta
        }
        
    } 
    Adj_temp <- list(Adj_temp_b, Adj_temp_SE)
    return(Adj_temp)
}
##==============================================================================
# choose data from same gene names
Choose_same_gene <- function(data1, data2)  # column is gene
{
    data1 <- data1[, which(!is.na(match(colnames(data1), data2)))]
    return(data1)
}
##==============================================================================
# choose data from same sample names
Choose_same_sample <- function(data1, data2)      # row is sample
{
    data1 <- data1[which(!is.na(match(rownames(data1), data2))), ]
    return(data1)
}

