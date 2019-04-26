# Follow these instructions https://rpubs.com/wendyu/sparkr

Sys.setenv(SPARK_HOME = "/Users/MatthewHamilton/Dropbox/KSULasso/R/SparkSandbox/spark-2.4.1-bin-hadoop2.7/")
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
library("SparkR", lib.loc = "/Users/MatthewHamilton/Dropbox/KSULasso/R/SparkSandbox/spark-2.4.1-bin-hadoop2.7/")
library(SparkR)
sc <- sparkR.session(master = "local",
                     sparkEnvir = list(spark.driver.memory = "2g",
                                       spark.num.executors = "2",
                                       spark.executor.cores = "2"))
sc <- sparkR.session(sc) #?



setwd("~/Dropbox/KSULasso/R/SparkSandbox/")
x = read.csv("x")[,-1]
y = read.csv("y")[,-1]
all <- cbind(y, x)
features <- ncol(x)
samples <- nrow(x)
bootstraps <- round(features / samples) * 80

# df <- createDataFrame(sqlContext, all) # SPARK DATAFRAME

# >>>>>>>>>> Sandbox <<<<<<<<<<
head(df)
head(select(df, df$y, df$P1))
head(select(df, columns(df)))
# >>>>>>>>>> Sandbox <<<<<<<<<<

#random.features <- c(1, base::sample(features, samples, replace = FALSE))
#df_sample <- sample(df, withReplacement = TRUE, fraction = 1, seed = 1)
#df_sample <- subset(df, select = random.features)
#df_coef <- glm(dependent ~ ., data = df_sample)
#summary(df_coef)

# rdd <- SparkR:::parallelize(sqlContext, 1:2)
verbose = TRUE
alpha = 1
number.of.features <- ncol(x)
number.of.samples <- nrow(x)

Selper.part.a <- function(ii, x, y) {
  library(glmnet)
  # Sample features column numbers equal to the number of samples.

  random.features <- base::sample(number.of.features, number.of.samples,
                                  replace = FALSE)
  # Mix up the rows.
  random.samples <- base::sample(number.of.samples, replace = TRUE)
  # Subset the columns and rows from the independent data.
  random.x <- x[random.samples, random.features]
  # print(random.x)
  # Sys.sleep(1)
  # Subset the rows from the dependent data.
  random.y <- y[random.samples, ]

  # Centering the dependent variable.
  random.y.mean <- mean(random.y)
  random.y.scale <- random.y - random.y.mean

  # Centering the independent variable.
  random.x.mean <- apply(random.x, 2, mean)
  random.x.scale <- scale(random.x, random.x.mean, FALSE)
  standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
  random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)
  # Filling in the empty rows one by one with results from lasso.
  cv.lasso.results <- cv.glmnet(random.x.scale, random.y.scale, type.measure = "mse",
                                nfold = 5, alpha = alpha)
  lasso.results <- glmnet(random.x.scale, random.y.scale, lambda = cv.lasso.results$lambda.min,
                          standardize = FALSE, alpha = alpha, intercept = FALSE)
  # Coeffiecents of random sample.
  hat.beta <- coef(lasso.results)[-1]
  return(hat.beta)
}
start <- Sys.time()
beta.hatbeta.hat <- lapply(1:12, Selper.part.a, x, y)
test1 <- Sys.time() - start

start <- Sys.time()
beta.hat <- spark.lapply(1:12, function(ii) Selper.part.a(ii, x, y))
test2 <- Sys.time() - start



map <- SparkR:::map()
reduce <- SparkR:::reduce()
