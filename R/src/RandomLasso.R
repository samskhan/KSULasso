#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' @param treatment.data Matrix of treatment data.
#' @param control.data Matrix of control data.
#' @param bootstrap
#' @param suppress Supresses all printing and time estimation.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(asthma.treatment, asthma.control)

source("src/EstimateTime.R")

RandomLasso <- function(treatment.data, control.data, bootstrap,
                        PVAL = 0.05, C_I = 0.95, P = 0.5, S_E = 0.05,
                        REP = 10, suppress = FALSE) {

  gene.names <- colnames(treatment.data)
  total.genes <- ncol(treatment.data)
  STEP2 = FALSE
  
  un <- ((qnorm(C_I, 0, 1))^2 * P * (1 - P)) / S_E^2
  features <- floor((un * total.genes) / (un + (total.genes - 1))) + 1
  if (missing(bootstrap)) {
    bootstrap <- round(total.genes / features * REP) * 2
  }
  
  # Error Handling
  unit.test.genes <- ncol(control.data)
  try( if (total.genes != unit.test.genes) {
    stop("Fatal: Argument ", treatment.data, " and ", unit.test.genes, " have
         different columns.")
    })
  
  Importance_TRT <- matrix(0, nrow = total.genes, ncol = total.genes)
  Importance_CON <- matrix(0, nrow = total.genes, ncol = total.genes)
  select.prob <-  matrix(1 / total.genes, nrow = total.genes, ncol = total.genes)
  
  adj.genes.net.treatment <- matrix(0, nrow = total.genes, ncol = total.genes)
  rownames(adj.genes.net.treatment) <- gene.names
  colnames(adj.genes.net.treatment) <- gene.names
  
  adj.genes.net.control <- adj.genes.net.treatment
  adj.genes.se.treatment <- adj.genes.net.treatment
  adj.genes.se.control <- adj.genes.net.treatment
  
  temp.b.treatment <- matrix(0, nrow = total.genes, ncol = bootstrap)
  rownames(temp.b.treatment) <- gene.names

  temp.b.control <- temp.b.treatment
  temp.se.treatment <- temp.b.treatment
  temp.se.control <- temp.b.treatment

  if (!suppress) {
    tracked.time <- matrix(0, nrow = total.genes, ncol = 6)
    colnames(tracked.time) <- c("Remaining", "Complete",
                                "Percent", "Passed", "Average", "Total")
  }
  
  start.time = Sys.time()
  
  for (ii in 1:total.genes) {
    
    if (!suppress) {
      cat(paste("Step 1: ", ii, " / ", total.genes, "\n"))
      if (ii > 1) {
        tracked.time = PrintEstimateTime(tracked.time, start.time,
                                         ii, total.genes, TRUE)
      }
    }
   
    x.treatment <- treatment.data
    x.treatment[, ii] <- 1 # Why make a column of 1's?

    x.control <- control.data
    x.control[, ii] <- 1 # Why make a column of 1's?

    temp.b.treatment[] <- 0
    temp.b.control[] <- 0
    temp.se.treatment[] <- 0
    temp.se.control[] <- 0
   
    for (jj in 1:bootstrap) {
      
      if (!suppress) {
        if (jj < 10) {
          cat(paste( "  ", jj, ":", sep = ""))
        } else {
          cat(paste( " ", jj, ":", sep = ""))
        }
      }
      
      gene.sample <- sample(gene.names[-ii], features, replace = FALSE,
                            prob = select.prob[-ii, ii])
      xx.treatment <- x.treatment[, which(!is.na(match(colnames(x.treatment),
                                                       gene.sample)))]
      xx.control <- x.control[, which(!is.na(match(colnames(x.control),
                                                   gene.sample)))]
      xx.treatment <- cbind(1, xx.treatment)
      colnames(xx.treatment)[1] <- 'Y'

      xx.control <- cbind(1, xx.control)
      colnames(xx.control)[1] <- 'Y'

      #xx <- list(xx.treatment, xx.control)
      #yy <- list(treatment.data[, ii], control.data[, ii])
      #bb <- list(temp.b.treatment, temp.b.control)
      #ss <- list(temp.se.treatment, temp.se.control)
      #gg <- list(gene.names)
      
      # Early Attempts at Vectorization -Matthew
      #a = mapply(LeastAngleSquare, xx, yy, gg, jj, bb, ss, 0.01, STEP2)
      
      bootstrapped.treatment <- LeastAngleSquare(xx.treatment, treatment.data[, ii],
                                         gene.names, jj, temp.b.treatment,
                                         temp.se.treatment, PVAL = 0.01,
                                         STEP2 = STEP2)
      bootstrapped.control <- LeastAngleSquare(xx.control, control.data[, ii],
                                         gene.names, jj, temp.b.control,
                                         temp.se.control, PVAL = 0.01,
                                         STEP2 = STEP2)
      
      temp.b.treatment[, jj]  <- bootstrapped.treatment[[1]][, jj] # 36 by 214
      temp.se.treatment[, jj] <- bootstrapped.treatment[[2]][, jj]
      temp.b.control[, jj]  <- bootstrapped.control[[1]][, jj]
      temp.se.control[, jj] <- bootstrapped.control[[2]][, jj]
      
      if (!suppress) {
        if ((jj %% 12) == 0) {
          cat("\n")
        }
      }
    }
    
    bootstrap.1 <- rowSums(temp.b.treatment != 0)
    bootstrap.1[bootstrap.1 == 0] = 1
    
    bootstrap.2 <- rowSums(temp.b.control != 0)
    bootstrap.2[bootstrap.2 == 0] = 1
    
    if (!suppress) {
      print(bootstrap.1)
      print(bootstrap.2)
    }
    
    adj.genes.net.treatment[, ii] <- rowSums(temp.b.treatment)  / bootstrap.1
    adj.genes.se.treatment[, ii] <- sqrt(rowSums(temp.se.treatment^2) / bootstrap.1)
    adj.genes.net.control[, ii] <- rowSums(temp.b.control)  / bootstrap.2
    adj.genes.se.control[, ii] <- sqrt(rowSums(temp.se.control^2) / bootstrap.2)
  }
  
  # Step II | To Be Continued
  while (FALSE) {
    Importance_TRT <- adj.genes.net.treatment
    Importance_CON <- adj.genes.net.control
    STEP2 = TRUE
    
    for (ii in 1:total.genes) {
      
      if (!suppress) {
        cat(paste("Step 2: ", ii, " / ", total.genes, "\n"))
        if (ii > 1) {
          tracked.time = PrintEstimateTime(tracked.time, start.time,
                                           ii, total.genes, TRUE)
        }
      }
      
      x.treatment <- treatment.data
      y.treatment <- x.treatment[, ii]
      x.treatment[, ii] <- 1
      
      x.control <- control.data
      y.control <- x.control[, ii]
      x.control[, ii] <- 1
  
      temp.b.treatment <- matrix(0, nrow = total.genes, ncol = bootstrap)
      rownames(temp.b.treatment) <- gene.names
      
      temp.b.control <- matrix(0, nrow = total.genes, ncol = bootstrap)
      rownames(temp.b.control) <- gene.names
      
      temp.se.treatment <- matrix(0, nrow = total.genes, ncol = bootstrap)
      rownames(temp.se.treatment) <- gene.names
      
      temp.se.control <- matrix(0, nrow = total.genes, ncol = bootstrap)
      rownames(temp.se.control) <- gene.names
      
      for (jj in 1:bootstrap) {
        
        if (!suppress) {
          if (jj < 10) {
            cat(paste( "  ", jj, ":", sep = ""))
          } else {
            cat(paste( " ", jj, ":", sep = ""))
          }
        }
        
        select.prob_TRT <- 1/total.genes + Importance_TRT[, ii]
        gene.sample.treatment <- sample(gene.names[-ii], features,replace = FALSE,
                                        prob = select.prob_TRT[-ii])
        
        xx.treatment <- x.treatment[, which(!is.na(match(colnames(x.treatment),
                                                         gene.sample.treatment)))]
        xx.treatment <- cbind(1, xx.treatment)
        colnames(xx.treatment)[1] <- 'Y'
        
        select.prob_CON <- 1/total.genes + Importance_CON[, ii]
        gene.sample.control <- sample(gene.names[-ii], features, replace = FALSE,
                                      prob = select.prob_CON[-ii])
        
        xx.control <- x.control[, which(!is.na(match(colnames(x.control),
                                                     gene.sample.control)))]
        xx.control <- cbind(1, xx.control)
        colnames(xx.control)[1] <- 'Y'
      
        bootstrapped.treatment <- LeastAngleSquare(xx.treatment, y.treatment, gene.names,
                                           jj, temp.b.treatment, temp.se.treatment,
                                           PVAL = 0.01, STEP2)
        bootstrapped.control <- LeastAngleSquare(xx.control, y.control, gene.names, jj,
                                           temp.b.control, temp.se.control,
                                           PVAL = 0.01, STEP2)
        
        temp.b.treatment[, jj]  <- bootstrapped.treatment[[1]][, jj]
        temp.se.treatment[, jj] <- bootstrapped.treatment[[2]][, jj]
        
        temp.b.control[, jj]  <- bootstrapped.control[[1]][, jj]
        temp.se.control[, jj] <- bootstrapped.control[[2]][, jj]
        if ((jj %% 6) == 0) {
          cat("\n")
        }
      }
      
      bootstrap.1 <- rowSums(temp.b.treatment != 0)
      bootstrap.1[bootstrap.1 == 0] = 1
      
      bootstrap.2 <- rowSums(temp.b.control != 0)
      bootstrap.2[bootstrap.2 == 0] = 1
      
      if (!suppress) {
        print(bootstrap.1)
        print(bootstrap.2)
      }
      
      adj.genes.net.treatment[, ii] <- rowSums(temp.b.treatment)  / bootstrap.1
      adj.genes.se.treatment[, ii]  <- sqrt(rowSums(temp.se.treatment^2)  / bootstrap.1 )
      adj.genes.net.control[, ii] <- rowSums(temp.b.control)  / bootstrap.2
      adj.genes.se.control[, ii]  <- sqrt(rowSums(temp.se.control^2)  / bootstrap.2)
    }
  }
  
  results <- list(adj.genes.net.treatment, adj.genes.se.treatment,
                         adj.genes.net.control, adj.genes.se.control)   
  write.csv(adj.genes.net.treatment, paste("res/adj.genes.net.treatment",
                                           Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.treatment, paste("res/adj.genes.se.treatment",
                                          Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.net.control, paste("res/adj.genes.net.control",
                                         Sys.Date(), ".csv", sep = ""))
  write.csv(adj.genes.se.control, paste("res/adj.genes.se.control",
                                        Sys.Date(), ".csv", sep = ""))
  write.csv(tracked.time, paste("res/TrackedTime",
                                Sys.Date(), ".csv", sep = ""))
  return(results)
}

LeastAngleSquare <- function(XX, y, gene.names, jj, Adj_temp_b,  Adj_temp_SE, PVAL, STEP2) {
    #las <- lars(as.matrix(XX), y, type = "lasso", use.Gram = FALSE)
    las2 <- glmnet(XX, y, alpha = 1, family = "gaussian")
    #cvlas <- cv.lars(as.matrix(XX), y, type = "lasso", plot.it = FALSE, use.Gram = FALSE)
    cvlas2 <- cv.glmnet(XX, y, alpha = 1)
    #frac <- cvlas$index[which.min(cvlas$cv)]
    #hat_beta <- predict.lars(las, type = "coefficients", mode = "fraction", s = frac)$coefficients
    hat_beta2 <- predict.glmnet(las2, type = "coefficients", s = cvlas2$lambda.1se)
    idx <- which(hat_beta2 != 0)

    if (length(idx) == 0) { 
        Adj_temp_b[, jj] <- 0
        Adj_temp_SE[, jj] <- 0
    } else if (STEP2 == FALSE) {
        gene_id2 <- match(hat_beta2@Dimnames[[1]][-c(1:2)], gene.names)
        Adj_temp_b[gene_id2, jj] <- hat_beta2[-c(1,2)]
        Adj_temp_SE[, jj] <- 0
    } else if (STEP2 == TRUE) {
        nonzero_xx <- XX[, which(!is.na(match(colnames(XX), names(idx))))]
        print("There is nonzero xx")
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
            Adj_temp_b[gene_id, jj]  <- summary(fit.lm)$coefficients[sig_gene + 1, 1]
            Adj_temp_SE[gene_id, jj] <- summary(fit.lm)$coefficients[sig_gene + 1, 2]
        }
    } 

    Adj_temp <- list(Adj_temp_b, Adj_temp_SE)
    return(Adj_temp)
}