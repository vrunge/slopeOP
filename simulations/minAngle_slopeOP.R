######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ########

######## get the code in package slopeOP
#devtools::install_github("vrunge/slopeOP", force = TRUE)
library(slopeOP)

############## parallel computing ##############
library(parallel)
library(fields)
cores <- detectCores()
cores <- 40

################################################

######## building the model scenario 4
n <- 500 #data length
y1 <- seq(from = 10, to = 70, length.out = 250)
y2 <- seq(from = 70, to = 10, length.out = 251)
y <- c(y1, y2[-1])

################################################
### one.simu function

one.simu <- function(i, data, penalty, sigma, states, angle, type = "gaussian", deg = 3)
{
  #generate data with outliers
  n <- length(data)
  if(type == "gaussian")
  {
    z <- data + rnorm(n = n, sd = sigma)
  }
  if(type == "student")
  {
    coeff <- sigma*sqrt((deg - 2)/deg)
    z <- data + coeff*rt(n = 1000, df = 3)
  }

  ## run with smoothing and unconstrained optionsga
  s1 <- slopeOP(data = z, states = states, penalty = penalty, constraint = "smoothing", minAngle = angle)
  s2 <- slopeOP(data = z, states = states, penalty = penalty, constraint = "null")

  steps1 <- diff(s1$parameters)/diff(s1$changepoints)
  response1 <- rep(steps1, diff(s1$changepoints))
  response1 <- c(s1$parameters[1], cumsum(response1) + s1$parameters[1])
  steps2 <- diff(s2$parameters)/diff(s2$changepoints)
  response2 <- rep(steps2, diff(s2$changepoints))
  response2 <- c(s2$parameters[1], cumsum(response2) + s2$parameters[1])
  MSE1 <- sum((data-response1)^2)/n
  MSE2 <- sum((data-response2)^2)/n

  df <- data.frame(numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sigma", "angle",
                    "penalty", "nbSeg_minAngle", "nbSeg_std", "MSE_minAngle", "MSE_std")
  df[1,] <- c(i, n, sigma, angle, penalty/(sigma^2*log(n)),
              length(s1$changepoints) - 1, length(s2$changepoints) - 1, MSE1, MSE2)
  return(df)
}


#############################################
#######  GRAPH 1 : MSE and NbSegment  #######
#############################################

nbSimus <- 10 ### nb simu for each experiment

sigma <- 24 #initial standard deviation
myAngle <- 130
myBetas <- seq(from = (1/2)*sigma^2*log(500), to = 2.5*sigma^2*log(500), length.out = 40)

res <- NULL
for(beta in myBetas)
{
  res <- c(res, mclapply(1:nbSimus, FUN = one.simu,
                         data = y,
                         penalty = beta,
                         sigma = sigma,
                         states = 0:80,
                         angle = myAngle,
                         type = "student",
                         mc.cores = cores))
  print(beta)
}

df <- do.call(rbind, res)
df
write.csv(df, file="minAngles.csv")


##################################################
#################################################

library(ggplot2)
library(reshape2)


###############################
####### plot the result #######
###############################

df <- read.csv("minAngles.csv")
df <- df[,-1]


dfmean <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimus+1),each=nbSimus,len=nrow(df))),mean)[-1]

#####



df_plot <- dfmean[,c(5,8,9)]
colnames(df_plot) <- c("penalty",  "nb segments minAngle", "nb segments std", "MSE minAngle", "MSE std")
df_plot <- melt(df_plot, id.vars = "penalty")

# Everything on the same plot
plot1 <- ggplot(df_plot, aes(penalty, value, col=variable)) +
  geom_point(size = 2) + geom_line() +
  labs(y = "MSE") +
  labs(x = "penalty /"~sigma^2~"log(n)") +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=16),
        legend.position = c(0.7, 0.7),
        legend.title = element_blank())



df_plot <- dfmean[,c(5,6,7)]
colnames(df_plot) <- c("penalty",  "nb segments minAngle", "nb segments std", "MSE minAngle", "MSE std")
df_plot <- melt(df_plot, id.vars = "penalty")

# Everything on the same plot
plot2 <- ggplot(df_plot, aes(penalty, value, col=variable)) +
  geom_point(size = 2) + geom_line() +
  scale_color_manual(values = c("nbSeg_minAngle" = "#ff9c33", "nbSeg_std" = "#33beff")) +
  labs(y = "number of segments") +
  labs(x = "penalty /"~sigma^2~"log(n)") +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=16),
        legend.position = c(0.7, 0.7),
        legend.title = element_blank())

library(cowplot)

plot_grid(plot1, plot2, labels = "AUTO")
