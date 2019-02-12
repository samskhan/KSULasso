RandomLasso <- function(treatment.data, control.data, features, bootstrap,
                        PVAL = 0.05, STEP2, C_I = 0.95, P = 0.5, S_E = 0.05,
                        REP = 10) {

  source("src/EstimateTime.R")
  
  gene.names <- colnames(treatment.data)
  total.genes <- ncol(treatment.data)
  
  if (missing(features) | missing(bootstrap)) {
    un <- ((qnorm(C_I, 0, 1))^2 * P * (1 - P)) / S_E^2
    features <- floor((un * total.genes) / (un + (total.genes - 1))) + 1
    bootstrap <- round(total.genes / features * REP) * 2
  }
  
  unit.test.genes <- ncol(control.data) # 214
  try( if (total.genes != unit.test.genes) {
    stop("Fatal: Number of columns in Treatment and Control must be equal.")
    })
  
  Importance_TRT <- matrix(0, nrow = total.genes, ncol = total.genes)
  Importance_CON <- matrix(0, nrow = total.genes, ncol = total.genes)
    
  adj.genes.net.treatment <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.net.treatment) <- gene.names
  colnames(adj.genes.net.treatment) <- gene.names
  
  adj.genes.net.control <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.net.control) <- gene.names
  colnames(adj.genes.net.control) <- gene.names
  
  adj.genes.se.treatment <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.se.treatment) <- gene.names
  colnames(adj.genes.se.treatment) <- gene.names
  
  adj.genes.se.control <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.se.control) <- gene.names
  colnames(adj.genes.se.control) <- gene.names
  
  tracked.time <- matrix(0, nrow = total.genes, ncol = 6)
  colnames(tracked.time) <- c("Remaining", "Complete",
                              "Percent", "Passed", "Average", "Total")
  start.time = Sys.time()
  
  for (ii in 1:total.genes) {
    cat(paste("Step 1: ", ii, " / ", total.genes, "\n"))
    if (ii > 1) {
      tracked.time = PrintEstimateTime(tracked.time, start.time,
                                       ii, total.genes, TRUE)
    }
    x.treatment <- treatment.data
    y.treatment <- x.treatment[, ii] # Y value for lasso. -Matthew
    x.treatment[, ii] <- 1
    
    x.control <- control.data
    y.control <- x.control[, ii]  # Y value for lasso. -Matthew
    x.control[, ii] <- 1
    
    # step 1 of Random LASSO < making improtance measures of genes >
    ## call LASSO ftn for step 1 ( ) : how to process W_p and importance 
    Adj_temp_b_TRT <- matrix(0, nrow = total.genes, ncol = bootstrap)
    rownames(Adj_temp_b_TRT) <- gene.names
    
    Adj_temp_b_CON <- matrix(0, nrow = total.genes, ncol = bootstrap)
    rownames(Adj_temp_b_CON) <- gene.names
    
    Adj_temp_SE_TRT <- matrix(0, nrow = total.genes, ncol = bootstrap)
    rownames(Adj_temp_SE_TRT) <- gene.names
    
    Adj_temp_SE_CON <- matrix(0, nrow = total.genes, ncol = bootstrap)
    rownames(Adj_temp_SE_CON) <- gene.names
    # For ith dependent gene, randomly repeat lasso method of B times.
    for (jj in 1:bootstrap) {
      if (jj < 10) {
        cat(paste( "  ", jj, ":", sep = ""))
      } else {
        cat(paste( " ", jj, ":", sep = ""))
      }
      if (STEP2 == FALSE) { 
          select_prob <- 1/total.genes + Importance_TRT[, ii]
          # Random sample portion. Probabiltiy is determined above. -Matthew
          gene.sample <- sample(gene.names[-ii], features, replace = FALSE,
                                prob = select_prob[-ii])

          xx.treatment <- x.treatment[, which(!is.na(match(colnames(x.treatment), gene.sample)))]
          xx.control <- x.control[, which(!is.na(match(colnames(x.control), gene.sample)))]
          
          xx.treatment <- cbind(1, xx.treatment)  # add intercept first column in xx matrix
          colnames(xx.treatment)[1] <- 'Y'
          
          xx.control <- cbind(1, xx.control)  # add intercept first column in xx matrix
          colnames(xx.control)[1] <- 'Y'
          
      } else if (STEP2 == TRUE) {
          select_prob_TRT <- 1/total.genes + Importance_TRT[, ii]
          gene.sample.treatment <- sample(gene.names[-ii], features, replace = FALSE, prob = select_prob_TRT[-ii])

          xx.treatment <- x.treatment[, which(!is.na(match(colnames(x.treatment), gene.sample.treatment)))]
          xx.treatment <- cbind(1, xx.treatment)  # add intercept first column in xx matrix
          colnames(xx.treatment)[1] <- 'Y'

          select_prob_CON <- 1/total.genes + Importance_CON[, ii]
          gene.sample.control <- sample(gene.names[-ii], features, replace = FALSE, prob = select_prob_CON[-ii])

          xx.control <- x.control[, which(!is.na(match(colnames(x.control), gene.sample.control)))]
          xx.control <- cbind(1, xx.control)  # add intercept first column in xx matrix
          colnames(xx.control)[1] <- 'Y'
      }
      #xx <- list(xx.treatment, xx.control)
      #yy <- list(y.treatment, y.control)
      #bb <- list(Adj_temp_b_TRT, Adj_temp_b_CON)
      #ss <- list(Adj_temp_SE_TRT, Adj_temp_SE_CON)
      #gg <- list(gene.names)
      
      #Vectorization -Matthew
      #mcmapply(LarsForRandomLASSO, xx, yy, gg, jj, bb, ss, 0.01, STEP2, mc.cores = 2)
      
      Adj_temp_TRT <- LarsForRandomLASSO(xx.treatment, y.treatment, gene.names, jj, Adj_temp_b_TRT, Adj_temp_SE_TRT, PVAL = 0.01, STEP2)
      Adj_temp_CON <- LarsForRandomLASSO(xx.control, y.control, gene.names, jj, Adj_temp_b_CON, Adj_temp_SE_CON, PVAL = 0.01, STEP2)
      
      Adj_temp_b_TRT[, jj]  <- Adj_temp_TRT[[1]][, jj] # 36 by 214
      Adj_temp_SE_TRT[, jj] <- Adj_temp_TRT[[2]][, jj]
      
      Adj_temp_b_CON[, jj]  <- Adj_temp_CON[[1]][, jj]
      Adj_temp_SE_CON[, jj] <- Adj_temp_CON[[2]][, jj]
      if ((jj %% 6) == 0) {
        cat("\n")
      }
    }
    
    bootstrap.1 <- rowSums(Adj_temp_b_TRT != 0)
    bootstrap.1[bootstrap.1 == 0] = 1
    print(bootstrap.1)
    
    bootstrap.2 <- rowSums(Adj_temp_b_CON != 0)
    bootstrap.2[bootstrap.2 == 0] = 1
    print(bootstrap.2)
    
    adj.genes.net.treatment[, ii] <- rowSums(Adj_temp_b_TRT)  / bootstrap.1
    adj.genes.se.treatment[, ii]  <- sqrt(rowSums(Adj_temp_SE_TRT^2)  / bootstrap.1 )
    adj.genes.net.control[, ii] <- rowSums(Adj_temp_b_CON)  / bootstrap.2
    adj.genes.se.control[, ii]  <- sqrt(rowSums(Adj_temp_SE_CON^2)  / bootstrap.2)
  }
  Adj_gene_net_b <- list(adj.genes.net.treatment, adj.genes.se.treatment,
                         adj.genes.net.control, adj.genes.se.control)   
  write.csv(adj.genes.net.treatment, paste("res/adj.genes.net.treatment", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.treatment, paste("res/adj.genes.se.treatment", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.net.control, paste("res/adj.genes.net.control", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.control, paste("res/adj.genes.se.control", Sys.Date(), ".csv", sep = ""))
  write.csv(tracked.time, paste("res/TrackedTime", Sys.Date(), ".csv", sep = ""))
  return(Adj_gene_net_b)
}

LarsForRandomLASSO <- function(XX, y, gene.names, jj, Adj_temp_b,  Adj_temp_SE, PVAL, STEP2) {
    cat("A")  
    las       <- lars(as.matrix(XX), y, type = "lasso", use.Gram = FALSE)
    cat("B")
    cvlas     <- cv.lars(as.matrix(XX), y, type = "lasso", plot.it = FALSE, use.Gram = FALSE)
    cat("C")
    frac      <- cvlas$index[which.min(cvlas$cv)]
    cat("D")
    hat_beta  <- predict.lars(las, type = "coefficients", mode = "fraction", s = frac)$coefficients
    cat("E")
    #hat_beta  <- las.coef$coefficients

    idx <- which(hat_beta != 0)
    cat("F")

    if (length(idx) == 0) { 
        Adj_temp_b[, jj] <- 0
        Adj_temp_SE[, jj] <- 0

    } else if (STEP2 == FALSE) {
        gene_id <- match(names(hat_beta[-1]), gene.names)  # 1th is intercept, so remove
        Adj_temp_b[gene_id, jj] <- hat_beta[-1]
        Adj_temp_SE[gene_id, jj] <- 0
        #cat("fighting")
    } else if (STEP2 == TRUE) {

        nonzero_xx <- XX[, which(!is.na(match(colnames(XX), names(idx))))]
        print("There is nonzero xx")
        # linear regression for p-value
        fit.lm    <- lm(y ~. , data = data.frame(nonzero_xx))
        pval <- summary(fit.lm)$coefficients[, 4]
        pval <- pval[-1]
        sig_gene <- which(pval < PVAL)

        if (length(sig_gene) == 0 ) {
            Adj_temp_b[, jj] <- 0
            Adj_temp_SE[, jj] <- 0
        } else {
            print("There is significant xx")
            gene_id <- match(names(sig_gene), gene.names)
            Adj_temp_b[gene_id, jj]  <- summary(fit.lm)$coefficients[sig_gene + 1, 1]   # hat beta
            Adj_temp_SE[gene_id, jj] <- summary(fit.lm)$coefficients[sig_gene + 1, 2]  # SE of hat beta
        }
    } 

    Adj_temp <- list(Adj_temp_b, Adj_temp_SE)
    cat("G")
    return(Adj_temp)
}