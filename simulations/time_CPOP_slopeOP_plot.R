
######################
###### packages ######
######################

library(ggplot2)
library(cowplot)
library(reshape2)

############################################
###### summary function for quantiles ######
############################################

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



############################################
###### read and transform dataframe ########
############################################

MyData <- read.csv(file="timeSIGNAL.csv", header = TRUE)
df <- melt(MyData, id.vars = c("type","n","sigma"))

dfnew <- df[,-4]
df3new <- dfnew[dfnew$sigma == 3,]
df24new <- dfnew[dfnew$sigma == 24,]

df3newBIS <- data_summary(df3new, varname="value",
                          groupnames=c("n", "type"))
df3newBIS[df3newBIS$type == "CPOP",2] <- "CPOP (sigma = 3)"
df3newBIS[df3newBIS$type == "slopeOP",2] <- "slopeOP (sigma = 3)"


df24newBIS <- data_summary(df24new, varname="value",
                           groupnames=c("n", "type"))
df24newBIS[df24newBIS$type == "CPOP",2] <- "CPOP (sigma = 24)"
df24newBIS[df24newBIS$type == "slopeOP",2] <- "slopeOP (sigma = 24)"


theMin <- min(df3newBIS[,3:5],df24newBIS[,3:5])
theMax <- max(df3newBIS[,3:5],df24newBIS[,3:5])

################################
###### PLOT with ggplot2 #######
################################

# Everything on the same plot
p3 <- ggplot(df3newBIS, aes(x = n, y = value, col=type)) +  scale_x_log10()+ scale_y_log10(limits = c(theMin, theMax))  +
  labs(y = "time in seconds") +  labs(x = "length of the time series") +
  geom_point(size = 2, aes(shape = type)) +
  geom_errorbar(aes(ymin=`q1.2.5%`, ymax=`q3.97.5%`), width=.01) +
  scale_colour_manual(values = c("CPOP (sigma = 3)" = "#0080FF",
                                 "slopeOP (sigma = 3)" = " dark blue")) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = c(0.7, 0.1),
        legend.title = element_blank())

# Everything on the same plot
p24 <- ggplot(df24newBIS, aes(x = n, y = value, col=type)) +  scale_x_log10()+ scale_y_log10(limits = c(theMin, theMax))  +
  labs(y = "time in seconds") +  labs(x = "length of the time series") +
  geom_point(size = 2, aes(shape = type)) +
  geom_errorbar(aes(ymin=`q1.2.5%`, ymax=`q3.97.5%`), width=.01) +
  scale_colour_manual(values = c("CPOP (sigma = 24)" = "#0080FF",
                                 "slopeOP (sigma = 24)" = " dark blue")) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = c(0.7, 0.1),
        legend.title = element_blank())


plot_grid(p3, p24, labels = c("A", "B"))


################# coefficients and analysis #################
################# coefficients and analysis #################
################# coefficients and analysis #################

R1 <- df3newBIS[df3newBIS$type == "CPOP (sigma = 3)",c(1,3)]
R2 <- df3newBIS[df3newBIS$type == "slopeOP (sigma = 3)",c(1,3)]
R3 <- df24newBIS[df24newBIS$type == "CPOP (sigma = 24)",c(1,3)]
R4 <- df24newBIS[df24newBIS$type == "slopeOP (sigma = 24)",c(1,3)]

l1_CPOP <- lm(log(value) ~ log(n), data = R1, )
l2_slopeOP <- lm(log(value) ~ log(n), data = R2, )
l3_CPOP <- lm(log(value) ~ log(n), data = R3, )
l4_slopeOP <- lm(log(value) ~ log(n), data = R4, )

summary(l1_CPOP)
summary(l2_slopeOP)
summary(l3_CPOP)
summary(l4_slopeOP)

l1_CPOP$coefficients
l2_slopeOP$coefficients
l3_CPOP$coefficients
l4_slopeOP$coefficients


Nscale <-  df3newBIS[df3newBIS$type == "slopeOP (sigma = 3)","n"]

V3_slopeOP <- df3newBIS[df3newBIS$type == "slopeOP (sigma = 3)","value"]
V3_CPOP <- df3newBIS[df3newBIS$type == "CPOP (sigma = 3)","value"]
Nscale[ V3_slopeOP < V3_CPOP]

V24_slopeOP <- df24newBIS[df24newBIS$type == "slopeOP (sigma = 24)","value"]
V24_CPOP <- df24newBIS[df24newBIS$type == "CPOP (sigma = 24)","value"]
Nscale[ V24_slopeOP < V24_CPOP]


