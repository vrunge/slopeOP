

########           ########
######## packages  ########
########           ########

#devtools::install_github("vrunge/slopeOP", force = TRUE)
library(slopeOP)
library(parallel)
library(ggplot2)
library(reshape2)
library(cowplot)

#cores <- detectCores()
cores <- 95


########                       ########
######## One hat-shaped signal ########
########                       ########

######## building the model scenario 4
n <- 500 #data length
y1 <- seq(from = 10, to = 70, length.out = 250)
y2 <- seq(from = 70, to = 10, length.out = 251)
y <- c(y1, y2[-1])



########                   ########
######## one.simu function ########
########                   ########

one.simu <- function(i,
                     data,
                     penalty,
                     sigma,
                     states,
                     angle,
                     type = "gaussian",
                     deg = 3) #degree pf freedom for student type
{
  #generate data with outliers
  n <- length(data)
  if(type == "gaussian")
  {
    z <- data + rnorm(n = n, sd = sigma)
  }
  if(type == "student")
  {
    coeff <- sigma * sqrt((deg - 2)/deg)
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


####################################
####### : MSE and NbSegment  #######
####################################

nbSimus <- 100 ### nb simu for each experiment

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

write.csv(df, file="minAngles.csv", row.names = FALSE)



###############################
####### read the result #######
###############################

df <- read.csv("minAngles.csv")


data_summary <- function(data, varname, groupnames)
{
  require(plyr)
  summary_func <- function(x, col)
  {
    c(mean = mean(x[[col]], na.rm=TRUE),
      q1 = quantile(x[[col]], 0.025), q3 = quantile(x[[col]], 0.975))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


########                      ########
######## data transformation  ########
########                      ########


df1 <- data_summary(df, varname="MSE_minAngle",
             groupnames=c("penalty"))
df1 <- cbind(type = "minAngle", df1)
colnames(df1)[3] <- "MSE"

df2 <- data_summary(df, varname="MSE_std",
             groupnames=c("penalty"))
df2 <- cbind(type = "std", df2)
colnames(df2)[3] <- "MSE"

df12 <- rbind(df1,df2)

df3 <- data_summary(df, varname="nbSeg_minAngle",
             groupnames=c("penalty"))
df3 <- cbind(type = "minAngle", df3)
colnames(df3)[3] <- "nbSeg"

df4 <- data_summary(df, varname="nbSeg_std",
                    groupnames=c("penalty"))
df4 <- cbind(type = "std", df4)
colnames(df4)[3] <- "nbSeg"

df34 <- rbind(df3,df4)




########                    ########
######## plot with ggplot2  ########
########                    ########


plot1 <- ggplot(df12, aes(x = penalty, y = MSE, col=type))   +
  labs(y = "MSE") +  labs(x = "penalty /"~sigma^2~"log(n)") +
  geom_point(size = 2, aes(shape = type)) +
  geom_errorbar(aes(ymin=`q1.2.5%`, ymax=`q3.97.5%`), width=.05) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=16),
        legend.position = c(0.7, 0.7),
        legend.title = element_blank())


plot2 <- ggplot(df34, aes(x = penalty, y = nbSeg, col=type))   +
  labs(y = "number of segments") +  labs(x = "penalty /"~sigma^2~"log(n)") +
  geom_point(size = 2, aes(shape = type)) +
  geom_errorbar(aes(ymin=`q1.2.5%`, ymax=`q3.97.5%`), width=.05) +
  scale_colour_manual(values = c("minAngle" = "#ff9c33", "std" = "#33beff")) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=16),
        legend.position = c(0.7, 0.7),
        legend.title = element_blank())


plot_grid(plot1, plot2, labels = "AUTO")
