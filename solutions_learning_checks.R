#########################################
# Chapter 3: Objects an data types in R #
#########################################

fav_music <- c("Prince", "REM", "Ryan Adams", "BLOF")
num_concerts <- c(0, 3, 1, 0)
num_records <- c(2, 7, 5, 1)
my_music <- data.frame(fav_music, num_concerts, num_records)
names(my_music) <- c("artist", "concerts", "records")
summary(my_music)
my_music$records
sum(my_music$records)


##############################################
## Chapter 4: Getting started with data in R #
##############################################

# Exercise 1

path <- file.path('C:/Users/u0043788/Dropbox/IIR Machine learning en Data Science opleiding/Bookdown R Intro/data')
path.na <- file.path(path, "na.txt")
data_na <- read.table(path.na, header = TRUE)
data_na
is.na(data_na)
is.na(data_na$wage)
is.na(data_na$school)
is.na(data_na$expr)
is.na(data_na$female)
is.na(data_na$industry)
# '999', 'do not know' and 'na' should also be considered as missing values
data_na <- read.table(path.na, header = TRUE, na.strings = c("Do not know", "na", "NA", 999))
is.na(data_na)
data_na$female <- as.factor(data_na$female)
data_na
is.na(data_na$school)
sum(is.na(data_na$wage))

# Exercise 2
# install.packages("AER")
library(AER)
data("Parade2005")
? Parade2005
head(Parade2005)

attach(Parade2005)

CAData <- subset(Parade2005, state == "CA")

# more useful instructions
subset(Parade2005, state %in% c("CA", "ID"))
subset(Parade2005, state == "ID" & gender == "male")

mean(CAData$earnings)
hist(CAData$earnings)

tapply(earnings, state, mean)
aggregate(earnings ~ state, Parade2005, mean)

s <- subset(Parade2005, state == "ID")
s
summary(Parade2005$state)
d <- aggregate(earnings ~ state , Parade2005, length)
d[d$state == "ID", ]

tapply(earnings, celebrity, mean)
tapply(earnings, celebrity, median)

boxplot(log(earnings) ~ celebrity, col = "light blue")

plot(density(log(earnings), bw = "SJ"), type = "l", main = "log(earnings)")
rug(log(earnings))
detach(Parade2005)


#####################################
## Chapter 5: Visualizing data in R #
#####################################

# Exercise 1

#path <- file.path('C:/Users/u0043788/Dropbox/PE Introduction to R/data')
path <- file.path('C:/Users/u0043788/Dropbox/IIR Machine learning en Data Science opleiding/Bookdown R Intro/data')
path.danish <- file.path(path, "danish.txt")
danish <- read.table(path.danish, header = TRUE)
danish$Date <- as.Date(danish$Date, "%m/%d/%Y")
str(danish)
plot(danish$Date, danish$Loss.in.DKM, type = "l", xlab = "Date", ylab = "Loss",
     main = "Fire insurance data")

# Exercise 2

library(ggplot2)
ggplot(danish, aes(x = Date, y = Loss.in.DKM)) + theme_bw() + 
	geom_line()

# Exercise 3

#path <- file.path('C:/Users/u0043788/Dropbox/PE Introduction to R/data')
path <- file.path('C:/Users/u0043788/Dropbox/IIR Machine learning en Data Science opleiding/Bookdown R Intro/data')
path.car_price <- file.path(path, "car_price.csv")
car_price <- read.csv(path.car_price)

summary(car_price$price)
summary(car_price$income)

plot(car_price$income, car_price$price)
lines(lowess(car_price$income, car_price$price), col = "blue")

plot(price ~  income, data = car_price)

library(ggplot2)

ggplot(car_price, aes(x = income, y = price))+
  geom_point(shape = 1, alpha = 1/2)+
  geom_smooth() 

ggplot(car_price, aes(x = income, y = price))+ theme_bw()+
  geom_point(shape = 1, alpha = 1/2)+ 
  geom_smooth()

# Exercise 4

ggplot(mpg, aes (x = displ, y = hwy))+ 
  geom_point()

ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
  geom_point(mapping = aes(color = class)) + 
  geom_smooth()

##############################
## Chapter 6: Data wrangling #
##############################

library(AER)
install.packages("dplyr")
library(dplyr)
data("Parade2005")
? Parade2005
head(Parade2005)

attach(Parade2005)

Parade2005 %>% filter(state == "CA")
subset(Parade2005, state == "CA")

Parade2005 %>% filter(state == "CA") %>% summarize(mean = mean(earnings))
mean(subset(Parade2005, state == "CA")$earnings)

test <- Parade2005 %>% group_by(state) %>% summarize(mean = mean(earnings))

Parade2005 %>% filter(state == "ID") %>% summarize(number = n())
length(subset(Parade2005, state == "ID"))
d <- aggregate(earnings ~ state , Parade2005, length)
d[d$state == "ID", ]

Parade2005 %>% group_by(celebrity) %>% 
  summarize(mean = mean(earnings), median = median(earnings))

Parade2005 %>% group_by(celebrity) %>% 
  summarize(mean = mean(earnings), median = median(earnings)) %>%
  ggplot(aes(x = celebrity, y = mean)) + theme_bw() +
  geom_point(color = "blue")

Parade2005 %>% group_by(celebrity) %>% 
  ggplot(aes(x = celebrity, y = earnings)) + theme_bw() +
  geom_boxplot(color = "blue")


detach(Parade2005)

###########################################################
## Chapter 7: Working with probability distributions in R #
###########################################################

# Exercise 1

set.seed(1)
random_numbers <- runif(10)
coin_tosses_1 <- ifelse(random_numbers>.5, 'head', 'tail')
set.seed(1)
coin_tosses <- rbinom(n = 10, size = 1, prob = .5)

# Exercise 2

set.seed(1)
heights <- rnorm(n = 100, mean = 1.70, sd = .1)
summary(heights)
pnorm(1.90, mean = 1.70, sd = .1)
1 - pnorm(1.60, mean = 1.70, sd = .1)
pnorm(1.60, mean = 1.70, sd = .1, lower.tail = FALSE)

# Exercise 3

set.seed(1)
patients <- rexp(rate = 1/50, n =30)
pexp(q = 10, rate = 1/50)


######################################
## Chapter 8: Writing functions in R #
######################################

# Exercise 1

f.sum <- function (x, y) {
  r <- x + y
  r
}

f.sum(5, 10)

# Exercise 2

variance <- function(x, biased = FALSE)
{
  if (biased)
  {
    n <- length(x)
    return((n - 1)/n * var(x))
  }
  else
    return(var(x))
}

x <- c(1:10)
variance(x)
var(x)
variance(x, TRUE)

# Exercise 3

f.count <- function (v, x) {
  count <- 0
  for (i in 1:length(v)) {
    if (v[i] == x) {
      count <- count + 1
    }
  }
  count
}

f.count(c(1:9, rep(10, 100)), 10)

# Exercise 4

desi <- function(x, med = FALSE) {
  
  mean <- round(mean(x), 1)
  stdv <- round(sd(x), 1)
  cat("Mean is:", mean, ", SD is:", stdv, "\n")
  
  if(med){
    median <- median(x)
    cat("Median is:", median , "\n")
  }
}

desi(1:10, med=TRUE)

#######################################
## Chapter 10: Linear Regression in R #
#######################################

library(mlbench)
library(dplyr)
library(ggplot2)
library(reshape2)
data("BostonHousing")
housing <- BostonHousing
str(housing)

housing %>%
  ggplot(aes(x = medv)) +
  stat_density() +
  labs(x = "Median Value ($1000s)", y = "Density", title = "Density Plot of Median Value House Price in Boston") +
  theme_bw()

hist(housing$medv, probability=TRUE, nclass=25, col="light blue")
lines(density(housing$medv))

summary(housing$medv)

housing %>%
  select(c(crim, medv)) %>%
  ggplot(aes(x = crim, y = medv)) +
  geom_point(alpha = 0.7) +
  geom_smooth(color = "blue") + theme_bw() 

res <- housing %>%
  select(c(crim, rm, age, rad, tax, lstat, medv)) %>%
  melt(, id.vars = "medv")

housing %>%
  select(c(crim, rm, age, rad, tax, lstat, medv)) %>%
  melt(, id.vars = "medv") %>%
  ggplot(aes(x = value, y = medv, colour = variable)) +
  geom_point(alpha = 0.7) +
  stat_smooth(color = "black") +
  facet_wrap( ~ variable, scales = "free", ncol = 2) +
  labs(x = "Variable 
       Value", y = "Median House Price ($1000s)") +
  theme_minimal()

library("caret")
set.seed(123)
to_train <- createDataPartition(y = housing$medv, p = 0.75, list = FALSE)
train <- housing[to_train, ]
test <- housing[-to_train, ]

first_lm <- lm(medv ~ crim + rm + tax + lstat, data = train)
lm1_rsqu <- summary(first_lm)$r.squared
print(paste("First linear model has an r-squared value of ", round(lm1_rsqu, 3), sep = ""))

second_lm <- lm(log(medv) ~ crim + rm + tax + lstat, data = train)
lm2_rsqu <- summary(second_lm)$r.squared
print(paste("Our second linear model has an r-squared value of ", round(lm2_rsqu, 3), sep = ""))
mean(second_lm$residuals)

predicted <- predict(second_lm, newdata = test)
results <- data.frame(predicted = exp(predicted), original = test$medv)

results %>%
  ggplot(aes(x = predicted, y = original)) +
  geom_point() +
  stat_smooth() +
  labs(x = "Predicted Values", y = "Original Values", title = "Predicted vs. Original Values") +
  theme_minimal()

