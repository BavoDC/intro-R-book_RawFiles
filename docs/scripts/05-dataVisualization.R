## ------------------------------------------------------------------------
# load the 'Journals' data set in the AER package
data("Journals")
# scan the data
head(Journals)
names(Journals)
# e.g. get variable 'price' 
Journals$price
summary(Journals$price)
# focus on price of journal per citation
Journals$citeprice <- Journals$price/Journals$citations

## ------------------------------------------------------------------------
attach(Journals)
plot(log(subs), log(citeprice))
rug(log(subs))	# adds ticks, thus visualizing the marginal distributions of
			# the variables, along one or both axes of an existing plot.		
rug(log(citeprice), side = 2)
detach(Journals)
# avoid "attach()" and "detach()"
plot(log(subs) ~ log(citeprice), data = Journals)

## ------------------------------------------------------------------------
plot(log(citeprice)~log(subs), data = Journals, pch = 19, col = "blue", xlim = c(0, 8), ylim = c(-7, 4), main = "Library subscriptions")
rug(log(Journals$subs))
rug(log(Journals$citeprice), side=2)
# subset data, look at journal entitled "Econometrica"
journal <- "Econometrica"
journal_info <- subset(Journals, title==journal)
x.val <- log(journal_info$subs)
y.val <- log(journal_info$citeprice)
text(x.val, y.val, journal, pos=2)

## ---- eval=FALSE---------------------------------------------------------
## path <- file.path('C:/Users/u0043788/Dropbox/PE Introduction to R/graphs')
## graph.path <- file.path(path, "myfile.pdf")
## pdf(graph.path, height = 5, width = 6)
## 	plot(log(citeprice)~log(subs), data = Journals, pch = 19, col = "blue", xlim =
## 	       c(0, 8), ylim = c(-7, 4),
## 	main = "Library subscriptions")
## 	rug(log(Journals$subs))
## 	rug(log(Journals$citeprice),side=2)
## 	journal <- "Econometrica"
##   journal_info <- subset(Journals, title==journal)
##   x.val <- log(journal_info$subs)
##   y.val <- log(journal_info$citeprice)
##   text(x.val, y.val, journal, pos=2)
## dev.off()

## ------------------------------------------------------------------------
curve(dnorm, from = -5, to = 5, col = "slategray", lwd = 3, main = "Density of the standard normal distribution")
text(-5, 0.3, expression(f(x) == frac(1, sigma ~~ sqrt(2*pi)) ~~ e^{-frac((x - mu)^2, 2*sigma^2)}),adj=0)

## ------------------------------------------------------------------------
library(ggplot2)

# use default theme
ggplot(data = mtcars, mapping = aes(x = hp, y = mpg)) +
  geom_point(shape = 1, alpha = 1/2) +
  geom_smooth() 
# shorter
ggplot(mtcars, aes(x = hp, y = mpg)) +
  geom_point(shape = 1, alpha = 1/2) +
  geom_smooth() 
# use black and white lay-out
ggplot(mtcars, aes(x = hp, y = mpg)) + theme_bw() +
  geom_point(shape = 1, alpha = 1/2)+ 
  geom_smooth() 

## ------------------------------------------------------------------------
ggplot(mtcars, aes(x = hp, y = mpg))+
  geom_point(mapping = aes(color = gear))

## ------------------------------------------------------------------------
ggplot(mtcars, aes(x = hp, y = mpg))+
  geom_point(mapping = aes(alpha = gear))

ggplot(mtcars, aes(x = hp, y = mpg))+
  geom_point(mapping = aes(size = gear))

## ------------------------------------------------------------------------
ggplot(mtcars, aes(factor(cyl), mpg)) +
  geom_boxplot() + geom_jitter() + theme_bw()

## ------------------------------------------------------------------------
p <- ggplot(mtcars, aes(factor(cyl), mpg))
p + geom_boxplot() + geom_jitter() + theme_bw()

## ------------------------------------------------------------------------
library(corrplot)
# get correlation matrix
M <- cor(mtcars)
str(M)
M
# visualize the correlation structure
corrplot(M, method="circle")
corrplot(M, method="square")
corrplot(M, method="color")
corrplot(M, type="upper")
corrplot(M, type="upper", method="square")

