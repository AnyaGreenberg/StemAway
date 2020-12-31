# MASS package contains functions and data sets supporting Venables and Ripley's Modern Applied Statistics with S
## It should come pre-installed with RStudio 
## List of functions and data sets: https://www.rdocumentation.org/packages/MASS/versions/7.3-53
library(MASS)

# ggplot2 package contains functions to construct publication-ready plots in R
# lines 8-11 installand load ggplot2 if it is not already installed
if (!require("ggplot2")) {
  install.packages("ggplot2") #could not find function
}
library(ggplot2)

## DATA LOADING ##
# Store cats data set as cat_data (info on sex, body weight, and heart weight)
# Since the cats data set is loaded when the MASS library is loaded we can just refer to it as 'cats'
# However, it is good practice to store a "pre-loaded" data set as a customized variable when you plan to use it
cat_data <- cats #or
cat_data <- read.csv("cat.csv", header=F, row.names=NULL)

# Display some data about the last cat entry in the data frame
cat_data$Sex[144,]
cat_data$Bwt[c(1:143)] #negative indexing
cat_data[144, 2:3]

## ADD CLASSIFICATION COLUMN ##
# Subset the cats into different heart weight categories
large_hearts <- cat_data[cat_data$Hwt > 12] #or
large_hearts <- subset(cat_data, Hwt > 12)

small_hearts <- subset(cat_data, Hwt < 8)

# Add a column vector labeling each entry id with "S", "M", or "L" heart weights
hwt_class <- rep("M", 154)
large_ids <- which(cat_data$Hwt > 12)
hwt_class[large_ids] <- "L"
small_ids <- whihc(cat_data$Hwt < 8)
hwt_class[small_ids] <- "S"
cat_data$Hwt_class <- hwt_class

# write cat_data to csv
write.csv(cats_data, "./cats_edited.csv")

## PLOTS ##
# count
png("frequency_plot.png", width=1000, height=1000)
freq <- ggplot(data=cat_data, aes(x=Hwt_class, colour=Sex))+
  geom_bar(stat="count", position="stack")+
  xlab("Heart Size")+ ylab("Count")+
  ggtitle("Number of Cats in Each Size Classification")+
  scale_fill_manual(values=c("midnightblue", "mediumturquoise"))+
  theme_bw()
dev.off()

# average
avg <- c(mean(cat_data[cat_data$Hwt_class == "S",2]),
         mean(cat_data[cat_data$Hwt_class == "M",2]))
se <- c(sd(cat_data[cat_data$Hwt_class == "S",]),
        sd(cat_data[cat_data$Hwt_class == "M",]),
        sd(cat_data[cat_data$Hwt_class == "L",]))
Hwt <- data.frame(Class=c("S", "M", "L"),
                  Avg=avg, SE=SE)

ggplot(data=Hwt, aes(x=Class, y=Avg, fill=Class))+
  geom_bar(stat="count")+
  geom_errorbar(aes(ymin=Avg-SE, ymax=Avg+SE), 
                width=1)+
  guides(fill=F)
  xlab("Heart Size")+ ylab("Average Body Weight (kg)")+
  ggtitle("Average Body Weights by Size Calssifaction")+
  scale_fill_brewer("Accent")+
  theme_bw()
ggsave("avg_plot.png", width=5, height=5)
