# Follow these instructions https://rpubs.com/wendyu/sparkr

Sys.setenv(SPARK_HOME = "/Users/MatthewHamilton/Dropbox/KSULasso/R/SparkSandbox/spark-2.4.1-bin-hadoop2.7/")
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
library(SparkR)
sc <- sparkR.session(master = "local", sparkEnvir = list(spark.driver.memory = "2g"))
sqlContext <- sparkR.session(sc) #?

setwd("~/Dropbox/KSULasso/R/SparkSandbox/")
x = read.csv("x")[,-1]
dependent = read.csv("y")[,-1]
all <- cbind(dependent, x)
features <- ncol(x)
samples <- nrow(x)
bootstraps <- round(features / samples) * 80

df <- createDataFrame(all) # SPARK DATAFRAME

# >>>>>>>>>> Sandbox <<<<<<<<<<
head(df)
head(select(df, df$dependent, df$P1))
head(select(df, columns(df)))
# >>>>>>>>>> Sandbox <<<<<<<<<<

random.features <- c(1, base::sample(features, samples, replace = FALSE))
df_sample <- sample(df, withReplacement = TRUE, fraction = 1, seed = 1)
df_sample <- subset(df, select = random.features)
df_coef <- spark.glm(dependent ~ ., data = df_sample)
summary(df_coef)
