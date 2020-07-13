##### plot() #####
# MASS package contains many different datasets
data <- data.frame(x=c(1,2,3,4,5,6,7,8,9,10), y=c(2,4,6,8,10,12,14,16,18,20))

# different types of plots
par(mfrow=c(3,3))
plot(data$x, data$y, type="p", col="red",
     main="Scatter Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="l", col="orange",
     main="Line Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="o", col="yellow",
     main="Line Plot with Points", xlab="x", ylab="y")
plot(data$x, data$y, type="b", col="green",
     main="Scatter Plot with Lines", xlab="x", ylab="y")
plot(data$x, data$y, type="h", col="blue",
     main="Histogram Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="s", col="purple",
     main="Staircase (H then V) Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="S", col="black",
     main="Staircase (V then H) Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="c", col="grey",
     main="Dashed Plot", xlab="x", ylab="y")
plot(data$x, data$y, type="n",
     main="None", xlab="x", ylab="y")

# legend
data2 <- data.frame(x=c(1,2,3,4,5,6,7,8,9,10), y=c(3,6,9,12,2,4,6,8,10,12))
plot(data$x, data$y, main="Legends", xlab="x", ylab="y") #defaults to type="p"
lines(data2$x, data2$y, col="red") #defaults to type="l"
legend("topleft", c("A", "B"), fill=c("black", "red"))
# other plotting functions
b <- data.frame(x=c("A", "B", "C"), y=c(54, 79, 23))
h <- c(rep(seq(1,100,10),10), rep(seq(1,100,5),7), rep(seq(1,100,3),2), rep(seq(1,100,30),8), rep(seq(1,100,18),9))
par(mfrow=c(3,1))
barplot(height=b$y, names=b$x, main="Barplot", xlab="x", ylab="y")
hist(h, main="Histogram", xlab="Data", ylab="Frequency")
boxplot(h, horizontal=TRUE, main="Boxplot")


##### PCA #####
data("iris")
data <- iris[,1:4]
species <- iris[,5]
pca <- prcomp(data, center=TRUE, scale=TRUE)
var <- cumsum(pca$sdev)/(sum(pca$sdev))
plot(var, type="l", main="Percent of Variation Explained by Each PC", xlab="PC", ylab="Percent (%)")
pc <- as.data.frame(pca$x)
plot(pc$PC1, pc$PC2, col=species, main="PCA Plot", xlab="PC1", ylab="PC2")
legend("topright", c("setosa", "versicolor", "virginica"), fill=c("black", "red", "green"))
