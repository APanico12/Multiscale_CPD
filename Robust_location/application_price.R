#import dataset price
df <-read.csv("eurchf_volatility_2014_2015.csv")
head(df)
plot(df$abs_log_return ,type="l")
price <- df$abs_log_return
#run our test
source("CUSUM.R")
res <- CUSUM.mean(x = price, k = 0.45, loss = "L2", c = 2.985, MC = 1000, linearized = TRUE, plotting = TRUE)
