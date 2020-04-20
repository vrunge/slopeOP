
######## get the code in package slopeOP

#devtools::install_github("vrunge/slopeOP", force = TRUE)
library(slopeOP)
library(ggplot2)
library(reshape2)

nb1 <- NULL
nb2 <- NULL
p <- 10 ### nb loops

mySigmas <- seq(from = 1, to = 30, by = 1)

###########################################
####### evaluate pruning efficiency #######
###########################################

for(sigma in mySigmas)
{
  n <- 500 #data length
  penalty <- 2 * sigma^2 * log(n)
  res1 <- NULL
  res2 <- NULL
  for(i in 1:p)
  {
    ## generate new hat-shaped data with given sd sigma
    y <- c(seq(from = 10, to = 50, length.out = n/2), seq(from = 50, to = 10, length.out = n/2)) + rnorm(n, 0, sigma)

    ## run with channel and pruning options
    s1 <- slopeOP(data = y, states = 0:60, penalty = penalty, type = "channel",  testMode = TRUE)
    res1 <- c(res1, s1$pruning)
    s2 <- slopeOP(data = y, states = 0:60, penalty = penalty, type = "pruning", testMode = TRUE)
    res2 <- c(res2, s2$pruning)
  }
  nb1 <- c(nb1, mean(res1))
  nb2 <- c(nb2, mean(res2))
}
df <- data.frame(mySigmas,nb1,nb2)
write.csv(df, "pruning.csv")

###############################
####### plot the result #######
###############################

df <- read.csv("pruning.csv")
df <- df[,2:4]
colnames(df) <- c("sigma", "channel pruning", "Inequality-based pruning")
df <- melt(df, id.vars = "sigma")


# Everything on the same plot
ggplot(df, aes(sigma, value, col=variable)) +
  geom_point(size = 2) +
  labs(y = "proportion of the m*m*n*(n-1) elements to scan for cost matrix Q") +
  labs(x = "sigma values") +
        theme(axis.text.x = element_text(size=18),
              axis.text.y = element_text(size=18),
              legend.text=element_text(size=18),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=16),
              legend.position = c(0.7, 0.1),
              legend.title = element_blank())


