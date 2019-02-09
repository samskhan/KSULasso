Determin_Feature_size <- function(Data, C_I = 0.95, P = 0.5, S_E = 0.05, REP = 10) {
    # P is the proportion of population that is unknown. -Matthew
    # S_E is sampling error. -Matthew
    # Determining sample size of finite population for RANDOM LASSO
    total.genes <- ncol(Data)
    # Returns z score, mean = 0 std = 1. -Matthew
    z_score <- qnorm(C_I, 0, 1)
    n0 <- (z_score^2 * P * (1 - P)) / S_E^2
    # Floor rounds down, e.g. 3.942 -> 3. -Matthew
    n <- floor((n0 * total.genes) / (n0 + (total.genes - 1))) + 1
    # Round rounds like normal. -Matthew
    rep_Boostrap <- round(total.genes / n * REP)
    # Generic vector filled with three values. -Matthew
    size <- list(n, n, rep_Boostrap)
    return(size)
}

RandomLasso <- function(treatment.data, control.data, Importance_TRT, Importance_CON,
                        features, bootstrap, PVAL, STEP2) {
    
  gene.names <- colnames(treatment.data)
  total.genes <- ncol(treatment.data) # 214
  unit.test.genes <- ncol(control.data) #214
  try( if (total.genes != unit.test.genes) {
    stop("Fatal: Number of columns in Treatment and control.data must be equal.")
    })
    
  # Projection matrix, Beta Hat. -Matthew
  adj.genes.net.treatment <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.net.treatment) <- gene.names # Why name rows?
  colnames(adj.genes.net.treatment) <- gene.names
  
  # Projection matrix, Beta Hat. -Matthew
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
    print("Hi!!!")
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
    Adj_temp_b_TRT <- matrix(0, nrow = total.genes, ncol = bootstrap)  # temporary storage for coefficient, AST
    rownames(Adj_temp_b_TRT) <- gene.names
    
    Adj_temp_b_CON <- matrix(0, nrow = total.genes, ncol = bootstrap)  # temporary storage for coefficient, CON
    rownames(Adj_temp_b_CON) <- gene.names
    
    Adj_temp_SE_TRT <- matrix(0, nrow = total.genes, ncol = bootstrap)  # temporary storage for coefficient, AST
    rownames(Adj_temp_SE_TRT) <- gene.names
    
    Adj_temp_SE_CON <- matrix(0, nrow = total.genes, ncol = bootstrap)  # temporary storage for coefficient, CON
    rownames(Adj_temp_SE_CON) <- gene.names
    # For ith dependent gene, randomly repeat lasso method of B times.
    for (jj in 1:bootstrap) {
      if (jj < 10) {
        cat(paste("...", jj, "  "))
      } else {
        if ((jj %% 10) == 0)
          cat(paste("...", jj, " \n"))
        else {
          cat(paste("...", jj, " "))
        }
      }
      if (STEP2 == FALSE) { 
          select_prob <- 1/total.genes + Importance_TRT[, ii]
          # Random sample portion. Probabiltiy is determined above. -Matthew
          gene.sample <- sample(gene.names[-ii], features, replace = FALSE, prob = select_prob[-ii])

          xx.treatment %<-% x.treatment[, which(!is.na(match(colnames(x.treatment), gene.sample)))]
          xx.control %<-% x.control[, which(!is.na(match(colnames(x.control), gene.sample)))]
          
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
      ##Fit in LASSO model : AST  from here
      Adj_temp_TRT %<-% LarsForRandomLASSO(xx.treatment, y.treatment, gene.names, jj, Adj_temp_b_TRT, Adj_temp_SE_TRT, PVAL = 0.01, STEP2)
      Adj_temp_CON %<-% LarsForRandomLASSO(xx.control, y.control, gene.names, jj, Adj_temp_b_CON, Adj_temp_SE_CON, PVAL = 0.01, STEP2)
      
      Adj_temp_b_TRT[, jj]  <- Adj_temp_TRT[[1]][, jj] # 36 by 214
      Adj_temp_SE_TRT[, jj] <- Adj_temp_TRT[[2]][, jj]
      
      Adj_temp_b_CON[, jj]  <- Adj_temp_CON[[1]][, jj]
      Adj_temp_SE_CON[, jj] <- Adj_temp_CON[[2]][, jj]

    } # end of for bootstrapping    to here  AST and CON respectively repeat
    
    # average the coefficient 
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
  write.csv(adj.genes.net.treatment, paste("bin/adj.genes.net.treatment", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.treatment, paste("bin/adj.genes.se.treatment", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.net.control, paste("bin/adj.genes.net.control", Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.control, paste("bin/adj.genes.se.control", Sys.Date(), ".csv", sep = ""))
  write.csv(tracked.time, paste("bin/TrackedTime", Sys.Date(), ".csv", sep = ""))
  return(Adj_gene_net_b)
}

LarsForRandomLASSO <- function(XX, y, gene.names, jj, Adj_temp_b,  Adj_temp_SE, PVAL, STEP2) {
    las       %<-% lars(as.matrix(XX), y, type = "lasso", use.Gram = FALSE)
    cvlas     %<-% cv.lars(as.matrix(XX), y, type = "lasso", plot.it = FALSE, use.Gram = FALSE)
    frac      <- cvlas$index[which.min(cvlas$cv)]
    hat_beta  <- predict.lars(las, type = "coefficients", mode = "fraction", s = frac)$coefficients
    #hat_beta  <- las.coef$coefficients

    idx <- which(hat_beta != 0)

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
    return(Adj_temp)
}

PrintEstimateTime <- function(tracked.time, start.time, current.increment,
                              end.increment, save = FALSE) {
  time.snip <- EstimateTime(start.time, current.increment, end.increment)
  cat(paste("Time Remaining:", Seconds2Time(time.snip[1]), "\n"))
  cat(paste("Estimated Completion:", as.POSIXct(time.snip[2], origin = "1970-01-01"), "\n"))
  cat(paste("Percent Complete:", round((100 * time.snip[3]), 2),"%\n"))
  cat(paste("Time Passed:", floor(time.snip[4]), "Seconds", "\n"))
  cat(paste("Average Time:", round(time.snip[5], 1), "Seconds"), "\n")
  cat(paste("Estimated Total Time:", Seconds2Time(time.snip[6]), "\n"))
  if (save) {
    tracked.time[current.increment - 1, 1] = time.snip[1]
    tracked.time[current.increment - 1, 2] = time.snip[2]
    tracked.time[current.increment - 1, 3] = time.snip[3]
    tracked.time[current.increment - 1, 4] = time.snip[4]
    tracked.time[current.increment - 1, 5] = time.snip[5]
    tracked.time[current.increment - 1, 6] = time.snip[6]
    return(tracked.time)
  }
}

EstimateTime <- function(start.time, current.increment, end.increment) {
  start.time <- as.numeric(start.time)
  current.time <- as.numeric(Sys.time())
  percent <- current.increment / end.increment
  passed <- current.time - start.time
  average <- passed / current.increment
  remaining <- (passed / percent) - passed
  complete <- current.time + remaining
  total <- remaining + passed
  values <-  c(remaining, complete, percent, passed, average, total)
  return(values)
}

Seconds2Time <- function(sec) {
  day <- floor(sec / 86400)
  hr <- floor(sec / 3600) - (day * 24)
  min <- floor(sec / 60) - (day * 1440) - (hr * 60)
  sec <- floor(sec) - (day * 86400) - (hr * 3600) - (min * 60)
  return(paste(day, " Days  ", hr, " Hours  ", min, " Minutes  ",
               sec, " Seconds", sep = ""))
}