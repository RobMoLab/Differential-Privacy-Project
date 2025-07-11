data <- read.csv("~/Desktop/Plots_PCOQS/size_proportions_NonPrivateModel_NonPrivateConformal.csv")
boxplot(data[, c(1,3,4)])
boxplot(data[, c(2,5,6)])
summary(data)
apply(data, 2, FUN = function(x) 2*sd(x))

data <- read.csv("~/Desktop/Plots_PCOQS/size_proportions_NonPrivateModel_NonPrivateConformal.csv", skip = 1)
summary(data)
apply(data, 2, FUN = function(x) 2*sd(x))