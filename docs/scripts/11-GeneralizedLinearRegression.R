## ---- eval=FALSE---------------------------------------------------------
## n <- scan(n = 54)
## 1  8 10  8  5 11 14 12 11 10  5 12 13 12 15 13 12 24
## 12 11  6  8 16 19 28 11 14  4 12  8 18  3 17  6 11 18
## 12  3 10 18 10 13 12 31 16 16 13 14  8 19 20  9 23 27
## 
## expo <- scan(n = 54) * 7
## 10 22 30 11 15 20 25 25 23 28 19 22 19 21 19 16 18 29
## 25 18 20 13 26 21 27 14 16 11 23 26 29 13 26 13 17 27
## 20 18 20 29 27 24 23 26 18 25 17 29 11 24 16 11 22 29
## 
## n
## expo

## ---- echo=FALSE---------------------------------------------------------
n <- scan(textConnection("
1  8 10  8  5 11 14 12 11 10  5 12 13 12 15 13 12 24
12 11  6  8 16 19 28 11 14  4 12  8 18  3 17  6 11 18
12  3 10 18 10 13 12 31 16 16 13 14  8 19 20  9 23 27"));

expo <- scan(textConnection(" 
10 22 30 11 15 20 25 25 23 28 19 22 19 21 19 16 18 29
25 18 20 13 26 21 27 14 16 11 23 26 29 13 26 13 17 27
20 18 20 29 27 24 23 26 18 25 17 29 11 24 16 11 22 29")) * 7;

n
expo

## ------------------------------------------------------------------------
sex <- as.factor(rep(1:2, each=27, len=54))
region <- as.factor(rep(1:3, each=9, len=54))
type <- as.factor(rep(1:3, each=3, len=54))
job <- as.factor(rep(1:3, each=1, len=54))
sex
region
type
job

## ------------------------------------------------------------------------
g1 <- glm(n ~ sex + region + type + job + offset(log(expo)), fam = poisson(link = log))

## ------------------------------------------------------------------------
summary(g1)

## ------------------------------------------------------------------------
names(g1)
g1$coef

## ------------------------------------------------------------------------
g1$fitted.values
g1$linear.predictors

plot(g1$fitted.values, n, xlab = "Fitted values", ylab = "Observed claims")
abline(lm(g1$fitted ~ n), col="light blue", lwd=2)
abline(0, 1, col = "dark blue", lwd=2)

## ------------------------------------------------------------------------
AIC(g1)

## ---- warning=FALSE, message=FALSE---------------------------------------
g2 <- glm(n/expo ~ sex+region+type+job,fam=poisson(link=log))
summary(g2)

## ---- warning=FALSE, message=FALSE---------------------------------------
g3 <- glm(n/expo ~ sex+region+type+job,weights=expo,fam=poisson(link=log))
summary(g3)

## ------------------------------------------------------------------------
g1 <- glm(n ~ 1 + region + type + job, poisson, offset = log(expo))
anova(g1, test="Chisq")

## ------------------------------------------------------------------------
# p-value for region
1 - pchisq(21.597, 2)
# or
pchisq(21.597, 2, lower.tail = FALSE)

## ----eval=FALSE----------------------------------------------------------
## anova(g1,test="Chisq")

## ------------------------------------------------------------------------
# what if we use 'F' instead of 'Chisq'?
anova(g1,test="F") 
# not appropriate for regular Poisson regression, see Warning message in the console!

## ------------------------------------------------------------------------
# Warning message:
# In anova.glm(g1, test = "F") :
#   using F test with a 'poisson' family is inappropriate

## ------------------------------------------------------------------------
(21.597/2)/1

## ------------------------------------------------------------------------
# construct an analysis-of-deviance table
g1 <- glm(n ~ 1, poisson , offset=log(expo))
g2 <- glm(n ~ sex, poisson , offset=log(expo))
g3 <- glm(n ~ sex+region, poisson, offset=log(expo))
g4 <- glm(n ~ sex+region+sex:region, poisson, offset=log(expo))
g5 <- glm(n ~ type, poisson, offset=log(expo))
g6 <- glm(n ~ region, poisson, offset=log(expo))
g7 <- glm(n ~ region+type, poisson, offset=log(expo))
g8 <- glm(n ~ region+type+region:type, poisson, offset=log(expo))
g9 <- glm(n ~ region+type+job, poisson, offset=log(expo))
g10 <- glm(n ~ region+type+sex, poisson, offset=log(expo))

## ------------------------------------------------------------------------
summary(g8)
g8$deviance

## ------------------------------------------------------------------------
anova(g1, g2, test = "Chisq")

## ------------------------------------------------------------------------
anova(g7, g8, test = "Chisq")

## ------------------------------------------------------------------------
g.poi <- glm(n ~ 1 + region + type, poisson, offset = log(expo))
summary(g.poi)

g.quasi <- glm(n ~ 1 + region + type, quasipoisson, offset = log(expo))
summary(g.quasi)

## ------------------------------------------------------------------------
# dispersion parameter in g is estimated as follows
phi <- sum(residuals(g.poi, "pearson")^2)/g.poi$df.residual
phi

## ------------------------------------------------------------------------
anova(g.quasi, test = "F")

## ------------------------------------------------------------------------
F <- (21.597/2)/phi
F

## ------------------------------------------------------------------------
pf(F, 2, 49, lower.tail = FALSE)

## ------------------------------------------------------------------------
# install.packages("MASS")
library(MASS)
g.nb <- glm.nb(n ~ 1+region+sex+offset(log(expo)))
summary(g.nb)

