source("CUMSUM.R")
source("LocalFunctions.R")

#Import the Excel file
file_path <- "C:/Users/Antonio/Desktop/log_yields.csv"
df <- read.csv(file_path,sep=",")
str(df)


#take the matrix of df 
X_matrix <- as.matrix(df[, -1])
rownames(X_matrix) <- df[, 1]

#X_matrix <- X_matrix[,5:10]
# Check the result
str(X_matrix)
head(X_matrix)

CUSUM_TEST(X_matrix,MC=1e3,ratio = 2/3 , plotting = TRUE)
# plot_one(X_matrix)
# plot_cusum(X_matrix)

