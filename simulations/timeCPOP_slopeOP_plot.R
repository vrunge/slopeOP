library(ggplot2)
library(reshape2)

#################  read.csv  #################
#################  read.csv  #################

MyData1 <- read.csv(file="time3.csv", header=TRUE)
MyData2 <- read.csv(file="time24.csv", header=TRUE)

#################  ggplot2  #################
#################  ggplot2  #################

x <- log(MyData1$mytest)
y1T1 <- log(MyData1$T1)
y1T2 <- log(MyData1$T2)
y2T1 <- log(MyData2$T1)
y2T2 <- log(MyData2$T2)

df <- data.frame(x, y1T1, y1T2, y2T1, y2T2)
colnames(df) <- c("n", "CPOP sigma = 3", "slopeOP sigma = 3",
                  "CPOP sigma = 24", "slopeOP sigma = 24")
df <- exp(df)

df2 <- melt(df, id.vars = "n")
colnames(df2) <- c("n", "variable" , "value")


# Everything on the same plot
ggplot(df2, aes(n, value, col=variable)) +  scale_x_log10()+ scale_y_log10()  +
  labs(y = "time in s") +  labs(x = "length of the time series") +
  geom_point(size = 2, aes(shape = variable)) +
  scale_colour_manual(values = c("CPOP sigma = 3" = "#0080FF",
                                 "slopeOP sigma = 3" = "#0080FF",
                                "CPOP sigma = 24" = "dark blue",
                                "slopeOP sigma = 24" = "dark blue")) +
  theme(axis.text.x = element_text(size=18),
       axis.text.y = element_text(size=18),
       legend.text=element_text(size=18),
       axis.title.x=element_text(size=18),
       axis.title.y=element_text(size=18),
       legend.position = c(0.7, 0.2),
       legend.title = element_blank())



  #################  coefficient  #################
#################  coefficient  #################

l1 <- lm(y1T1 ~ x)
l2 <- lm(y1T2 ~ x)
l3 <- lm(y2T1 ~ x)
l4 <- lm(y2T2 ~ x)

summary(l1)
summary(l2)
summary(l3)
summary(l4)
