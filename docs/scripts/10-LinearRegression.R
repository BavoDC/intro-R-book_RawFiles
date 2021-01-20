## ------------------------------------------------------------------------
path <- file.path('data')
path.car <- file.path(path, "car_price.csv")
car_price <- read.csv(path.car)

## ------------------------------------------------------------------------
attach(car_price)
summary(price)
summary(income)

# average
mean(price)
mean(income)

# standard deviation
sd(price)
sd(income)

# the 5-th and 95-th percentiles
quantile(price, c(0.05, 0.95))
quantile(income, c(0.05, 0.95))

# histograms of price and income
# density histogram for 'price'
hist(price, br = 20, xlim = c(5000, 30000), col="grey", freq=FALSE)
lines(density(price), col=4)
# frequency histogram for 'income/1000'
hist(income/1000, br=10, xlab="income (in $000's)", xlim=c(0, 120), col="grey")

# scatter plot 'income/1000' versus 'price'
plot(income/1000, price, pch=21, cex=1.2, xlab="income (in $000's)")
detach(car_price)

## ------------------------------------------------------------------------
library("ggplot2")
ggplot(car_price, aes(x = income/1000, y = price)) +
  theme_bw() +
  geom_point(shape=1, alpha = 1/2) + 
  geom_smooth() 

## ------------------------------------------------------------------------
lm1 <- lm(price ~ income, data = car_price)
summary(lm1)
# check attributes of object 'lm1'
names(lm1)
# some useful stuff: 'coefficients', 'residuals', 'fitted.values', 'model'
lm1$coef
lm1$residuals
lm1$fitted.values

## ------------------------------------------------------------------------
# use built-in plot function
# you may have noticed that we have used the function plot with all kinds of arguments: 
# one or two variables, a data frame, and now a linear model fit;
# in R jargon plot is a generic function; it checks for the kind of object that you # are plotting and then calls the appropriate (more specialized) function to do the work.
plot(lm1)

## ------------------------------------------------------------------------
# add the regression line to the scatter plot
plot(car_price$income, car_price$price, pch=21, cex=1.2, xlab = "income", main = "Simple linear regression")
# add LS line like this
abline(lm1, col="blue", lwd=2)
# or like this
abline(lm1$coefficients[1], lm1$coefficients[2])

## ------------------------------------------------------------------------
ggplot(car_price, aes(x = income, y = price)) + 
 theme_bw() +
 geom_point(shape=1, alpha = 1/2) + 
 geom_smooth()+geom_abline(intercept = lm1$coef[1], slope = lm1$coef[2],   colour="red", size=1.25) 

## ------------------------------------------------------------------------
plot(car_price$income, car_price$price, pch=21, cex=1.2, xlab = "income", main = "Simple linear regression")
abline(lm1, col = "blue", lwd=2)
segments(car_price$income, car_price$price, car_price$income, lm1$fitted.values, lty=1)

## ------------------------------------------------------------------------
summary(lm1)

## ------------------------------------------------------------------------
error.SS <- sum(lm1$resid^2)
error.SS
sqrt(error.SS/(nrow(car_price)-2))

## ------------------------------------------------------------------------
attach(car_price)
total.SS <- sum((price-mean(price))^2)
total.SS
error.SS <- sum(lm1$resid^2)
error.SS

# R^2?
(total.SS-error.SS)/total.SS
detach(car_price)

## ------------------------------------------------------------------------
attach(car_price)
# anova table in R?
anova(lm1)

# F-statistic in anova and in output lm1?
lm0 <- lm(price ~ 1)
error0.SS <- sum(lm0$resid^2)

# calculate F-statistic
F <- ((anova(lm0)$"Sum Sq")-(anova(lm1)$"Sum Sq"[2]))/(anova(lm1)$"Mean Sq"[2]) 
F
# critical values
qf(0.95, 1, 60)
1-pf(F, 1, 60)

detach(car_price)

## ------------------------------------------------------------------------
path <- file.path('data')
path.mort <- file.path(path, "pollution.csv")
mort_poll <- read.csv(path.mort)

## ------------------------------------------------------------------------
attach(mort_poll)
summary(mort_poll)
# get correlation matrix
round(cor(mort_poll), 4)

# create dataframes
# weather related vars
mort_poll_1 <- data.frame(mort, prec, jant, jult, humid)
# socio-economic vars
mort_poll_2 <- data.frame(mort, ovr65, popn, educ, hous, dens, nonw, wwdrk, poor)
# pollution effects
mort_poll_3 <- data.frame(mort, hc, nox, so2)

# matrix scatterplots
pairs(mort_poll_1, cex=1, pch=19)
pairs(mort_poll_2, cex=0.5, pch=19)
pairs(mort_poll_3, cex=1, pch=19)
detach(mort_poll)

## ------------------------------------------------------------------------
attach(mort_poll)
lm1 <- lm(mort ~ educ + so2)
summary(lm1)
detach(mort_poll)

## ------------------------------------------------------------------------
anova(lm1)
attach(mort_poll)
lm0 <- lm(mort ~ 1)
lm_educ <- lm(mort ~ educ)
anova(lm_educ)
F_educ <- ((anova(lm0)$"Sum Sq")-(anova(lm_educ)$"Sum Sq"[2]))/(anova(lm1)$"Mean Sq"[3]) 
F_educ
F_so2 <- ((anova(lm_educ)$"Sum Sq"[2])-(anova(lm1)$"Sum Sq"[3]))/(anova(lm1)$"Mean Sq"[3]) 
F_so2
detach(mort_poll)

## ------------------------------------------------------------------------
attach(mort_poll)
x0 <- data.frame(educ = 10, so2 = exp(2))
predict(lm1, x0, interval = "confidence")
predict(lm1, x0, interval = "prediction")
detach(mort_poll)

## ------------------------------------------------------------------------
attach(mort_poll)
grid <- seq(8, 15, 0.1)
x.new <- data.frame(educ = grid, so2 = exp(2))
p <- predict(lm1, x.new, se=TRUE, interval="prediction")
p1 <- predict(lm1, x.new, se=TRUE, interval="confidence")
# use `matplot` to plot the columns of one matrix against the columnsof another
matplot(grid, p$fit, lty=c(1,2,2), col=c("black", "red", "red"), type = "l", xlab = "educ", ylab = "mort", main = "Predicted mort over a range of educ, log(so2)=2")
matlines(grid, p1$fit, lty = c(1, 2, 2), col = c("black", "blue", "blue"))
rug(educ)
# for an explanation wrt different shapes, see 
# http://stats.stackexchange.com/questions/85560/shape-of-confidence-interval-for-p# redicted-values-in-linear-regression
detach(mort_poll)

## ------------------------------------------------------------------------
attach(mort_poll)
lm2 <- lm(mort ~ prec + jant + jult + humid + hc + nox + so2 + ovr65 + popn + educ + hous + dens + nonw + wwdrk + poor)
lm2$coef
detach(mort_poll)

## ------------------------------------------------------------------------
# model selection based on AIC
library(MASS)
attach(mort_poll)
lm1 <- lm(mort ~ 1)
# get AIC, mind the difference
AIC(lm1)
extractAIC(lm1)
# for linear models with unknown scale (i.e., for lm and aov), 
# -2 log L is computed from the deviance and uses a different additive constant to 
# logLik and hence AIC 

# forward search
stepAIC(lm1, list(upper = ~ prec + jant + jult + ovr65 + popn + educ + hous + dens + nonw + wwdrk + poor + hc + log(nox) + log(so2) + humid, lower = ~ 1), direction = "forward")

# backward search
lm1 <- lm(mort ~ prec + jant + jult + ovr65 + popn + educ + hous + dens + nonw +
		 wwdrk + poor + hc + log(nox) + log(so2) + humid)
stepAIC(lm1, list(upper = ~ prec + jant + jult + ovr65 + popn + educ + hous + dens + nonw + wwdrk + poor + hc + log(nox) + log(so2) + humid, lower = ~ 1), direction = "backward")

# both directions search
lm1 <- lm(mort ~ 1)
lm1 <- lm(mort ~ prec + jant + jult + ovr65 + popn + educ + hous + dens + nonw +
		 wwdrk + poor + hc + log(nox) + log(so2) + humid)
stepAIC(lm1, list(upper = ~ prec + jant + jult + ovr65 + popn + educ + hous + dens + nonw + wwdrk + poor + hc + log(nox) + log(so2) + humid, lower = ~ 1), direction = "both")
detach(mort_poll)

## ------------------------------------------------------------------------
library(mlbench)
data("BostonHousing")

