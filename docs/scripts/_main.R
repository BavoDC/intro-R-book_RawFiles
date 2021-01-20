## ---- include=FALSE------------------------------------------------------

# require(bookdown)
# render_book('index.Rmd', 'bookdown::gitbook')

# Set output options
options(width = 80, digits = 4, bookdown.clean_book = TRUE)
knitr::opts_chunk$set(
  tidy = FALSE,
  #out.width='\\textwidth',
  fig.width = 6, fig.height = 4,
  fig.align = "center",
  comment = NA
)

# Packages needed for following code in book
needed_pkgs <- c(
  # Data packages:
  "ISLR",
  # Used packages:
  "ggplot2", "tibble", "tidyr", "dplyr", "readr", "dygraphs", "rmarkdown",
  "knitr", "mosaic", "broom", "remotes", "forcats", "plotly", "moderndive",
  "janitor", "infer", "sas7bdat", "readxl", "AER", "corrplot", "car", "data.table",
  # Internally used packages:
  "devtools", "webshot", "tufte", "mvtnorm", "stringr", "gridExtra"
)
new_pkgs <- needed_pkgs[!(needed_pkgs %in% installed.packages())]
if(length(new_pkgs)) {
  install.packages(new_pkgs, repos = "http://cran.rstudio.com")
}


# Automatically create a bib database for R packages
knitr::write_bib(
  c(.packages(), "bookdown", "knitr", "rmarkdown", "nycflights13", "devtools",
    "ggplot2", "webshot", "dygraphs", "tufte", "okcupiddata", "mosaic", "dplyr",
    "ggplot2movies", "fivethirtyeight", "tibble", "readr", "tidyr"),
  "bib/packages.bib"
)

# Add all simulation results here
dir.create("rds")

dir.create("docs/scripts")



# Add all knitr::purl()'ed chapter R scripts here

#bookdown::render_book('index.Rmd', 'bookdown::gitbook')

# purl R scripts. For some reason this needs to be run manually:
#if(FALSE){
#  # Note order matters here:
#  chapter_titles <- c("objects-data-types", "started-with-data")
#  chapter_numbers <- stringr::str_pad(2:(length(chapter_titles) + 1), 2, "left", pad = "0")
#  for(i in 1:length(chapter_numbers)){
#    Rmd_file <- stringr::str_c(chapter_numbers[i], "-", chapter_titles[i], ".Rmd")
#    R_file <- stringr::str_c("docs/scripts/", chapter_numbers[i], "-", chapter_titles[i], ".R")
#    knitr::purl(Rmd_file, R_file)
#  }
#  file.exists("my_index_Arcturus.Rmd")
#  knitr::purl("my_index_Arcturus.Rmd")
#}


## ----pipeline-figure, echo=FALSE, fig.align='center', fig.cap="Data/Science Pipeline"----
knitr::include_graphics("images/tidy1.png")

## ----developers-figure, echo=FALSE, fig.align='center'-------------------
knitr::include_graphics("images/NYTimesR2009.jpg")

## ---- eval=FALSE---------------------------------------------------------
## library(ggplot2)
## library(dplyr)

## ------------------------------------------------------------------------
# use 'right click, run line or selection', of Ctrl+R
10^2+36

## ------------------------------------------------------------------------
# assign value '4' to 'a'
a <- 4
a
# now R remembers what 'a' is
# calculations with 'a'
a*5
(a+10)/2
# or give a new value to 'a'
a <- a+1
a

## ------------------------------------------------------------------------
my_numeric <- 42.5

my_character <- "some text"

my_logical <- TRUE

my_date <- as.Date("05/29/2018", "%m/%d/%Y")

## ------------------------------------------------------------------------
class(my_numeric)

# your turn to check the type of 'my_character' and 'my_logical' and 'my_date'

## ------------------------------------------------------------------------
# see all objects stored in R's memory, where 'ls()' is for 'List Objects' 
# and returns a vector of character strings
# giving the names of the objects in the specified environment
ls()
# to remove objects from R's memory, use
rm(a)
rm(my_character, my_logical)
rm(list=c('my_date', 'my_numeric'))
rm(list=ls())

## ------------------------------------------------------------------------
# To combine elements into a vector, use c():
c(1, 2, 3, 4)
# or
1:10
# or
seq(from=0, to=10, by=0.5)
# create a variable x
x <- 1:20
x
x[6:9]
x[c(2, 5, 13)]
# or
xx <- c(0, 3:5, 20, 0)
xx
xx[2:3]
length(xx)
# but c(.) can also concatenate other things than numbers
family <- c("Katrien", "Jan", "Leen")
family
family[2]
str(family) # str() displays the structure of an R object in compact way
class(family)

## ------------------------------------------------------------------------
my_vector <- c("Katrien Antonio", "teacher")
names(my_vector) <- c("Name", "Profession")
my_vector

## ------------------------------------------------------------------------
# a 3x4 matrix, filled with 1,2,..., 12
matrix(1:12, 3, 4, byrow = TRUE)
matrix(1:12, byrow = TRUE, nrow=3)
# hmmm, check help on 'matrix'
? matrix
# one way of creating matrices is to bind vectors together
cbind(1:2, 6:9)     # by columns
rbind(1:3, -(1:3))  # by rows
# create matrix object 'm'
m <- matrix(1:12, 3, 4)
m
m[1,4] # extract an element
m[,2]  # extract a column
nrow(m);ncol(m);dim(m) # useful stuff
# another example
m <- cbind(a = 1:3, b = letters[1:3])
m
# ask help, what is the built-in 'letters'?
? letters

## ---- eval=FALSE---------------------------------------------------------
## mtcars

## ------------------------------------------------------------------------
str(mtcars)
head(mtcars)
tail(mtcars)

## ------------------------------------------------------------------------
t <- data.frame(x = c(11, 12, 7), y = c(19, 20, 21), z = c(10, 9, 7))
t$x
t[["x"]]
# quick scan of the object 't'
summary(t)
str(t)
# another way to create the same data frame
x <- c(11, 12, 7)
y <- c(19, 20, 21)
z <- c(10, 9, 7)
t <- data.frame(x, y, z)

## ------------------------------------------------------------------------
mean(t$z)   
mean(z)   # does not work, why not?
attach(t) # but...
mean(z)
detach(t) # does the job
# or, avoid "attach(.)" and "detach(.)"
with(t, mean(z))

## ------------------------------------------------------------------------
# this does not work
# t <- data.frame(x = c(11,12), y = c(19,20,21), z = c(10,9,7)) 
# but you _can_ do
t <- data.frame(x = c(11, 12, NA), y = c(19, 20, 21), z = c(10, 9, 7))
# data frame with different types of information
b <- data.frame(x = c(11, 12, NA), y = c("me", "you", "everyone"))
str(b)
# hey there! 'y' should not be factor, but character variable
b$y <- as.character(b$y)
str(b)

## ------------------------------------------------------------------------
# a first example of a list
L <- list(one = 1, two = c(1, 2), five = seq(1, 4, length=5),
          six = c("Katrien", "Jan"))
names(L)
summary(L)
class(L)
str(L)

# list within a list
# a list containing: a sample from a N(0,1), plus some markup
# list within list
mylist <- list(sample = rnorm(5), family = "normal distribution", parameters = list(mean = 0, sd = 1))
mylist
str(mylist)

# now check
mylist[[1]]
mylist$sample
mylist$parameters
mylist$parameters$mean

## ---- eval=FALSE---------------------------------------------------------
## # what is the current working directory?
## getwd()
## # which files are currently stored in my working directory?
## dir()

## ------------------------------------------------------------------------
# where are my data files?
path <- file.path('data')
# how to find a path name on your computer?
# file.choose()

## ------------------------------------------------------------------------
path.pools <- file.path(path, "swimming_pools.csv")
pools <- read.csv(path.pools)
str(pools)

## ------------------------------------------------------------------------
pools <- read.csv(path.pools, stringsAsFactors = FALSE)
str(pools)

## ------------------------------------------------------------------------
path.fire <- file.path(path, "danish.txt")
danish <- read.table(path.fire, header = TRUE)
head(danish) # use the argument 'n' to display less/more records
tail(danish)
str(danish)
names(danish)
dim(danish)

## ------------------------------------------------------------------------
path.hotdogs <- file.path(path, "hotdogs.txt")
hotdogs <- read.table(path.hotdogs, header = FALSE, col.names = c("type", "calories", "sodium"))
# display structure of hotdogs
str(hotdogs)
# edit the colClasses argument to import the data correctly: hotdogs2
hotdogs2 <- read.table(path.hotdogs, header = FALSE, 
                       col.names = c("type", "calories", "sodium"),
                       colClasses = c("factor", "NULL", "numeric"))
# display structure of hotdogs2
str(hotdogs2)

## ------------------------------------------------------------------------
danish$Date <- as.Date(danish$Date, "%m/%d/%Y")
str(danish)

## ------------------------------------------------------------------------
path.fire <- file.path(path, "danish.txt")
danish <- read.table(path.fire, header = TRUE, colClasses = c("Date", "numeric"))

## ------------------------------------------------------------------------
setAs("character","myDate", function(from) as.Date(from, format="%m/%d/%Y") )
danish2 <- read.table(path.fire, header = TRUE, colClasses = c("myDate", "numeric"))
str(danish2)

## ------------------------------------------------------------------------
policy.path <- file.path(path, "policy.csv")
policy <- read.table(policy.path, header=TRUE, sep=";")
head(policy)
tail(policy)
str(policy)
names(policy)
dim(policy)

## ------------------------------------------------------------------------
#install.packages('sas7bdat')
library(sas7bdat)
path.severity <- file.path(path, "severity.sas7bdat")
severity <- read.sas7bdat(path.severity)
head(severity)
tail(severity)
str(severity)
names(severity)
dim(severity)

## ------------------------------------------------------------------------
# load the readxl package
library(readxl)
path.urbanpop <- file.path(path, "urbanpop.xlsx")
excel_sheets(path.urbanpop)

## ------------------------------------------------------------------------
pop_1 <- read_excel(path.urbanpop, sheet = 1)
pop_2 <- read_excel(path.urbanpop, sheet = 2)
pop_3 <- read_excel(path.urbanpop, sheet = 3)
str(pop_1)
# put pop_1, pop_2 and pop_3 in a list: pop_list
pop_list <- list(pop_1, pop_2, pop_3)

## ------------------------------------------------------------------------
pop_1_df <- as.data.frame(pop_1)
str(pop_1_df)

## ------------------------------------------------------------------------
pop_list <- lapply(excel_sheets(path.urbanpop), read_excel, path = path.urbanpop)
str(pop_list)

## ------------------------------------------------------------------------
path.urbanpop_nonames <- file.path(path, "urbanpop_nonames.xlsx")
# Import the the first Excel sheet of urbanpop_nonames.xlsx (R gives names): pop_a
pop_a <- read_excel(path.urbanpop_nonames, col_names = FALSE)
# Import the the first Excel sheet of urbanpop_nonames.xlsx (specify col_names): pop_b
cols <- c("country", paste0("year_", 1960:1966))
pop_b <- read_excel(path.urbanpop_nonames, col_names = cols)
# Print the summary of pop_a
summary(pop_a)
# Print the summary of pop_b
summary(pop_b)

## ------------------------------------------------------------------------
states <- data.frame(state.x77)
str(states) 
names(states)
dim(states)
head(states)
states[14, ]
states[3, 6]
states[, 'Frost'] 
states$Frost

## ------------------------------------------------------------------------
state.region
length(state.region)
# select those states that are in the south of the US 
mysubset <- subset(states, state.region == "South")
# subset a selection of variables
str(states)
mysubset <- states[, c(1:2, 7:8)]
mysubset <- states[, c("Population", "Income", "Frost", "Area")]

## ------------------------------------------------------------------------
least_pop <- which.min(states$Population)
states[least_pop, ]
most_pop <- which.max(states$Population)
states[most_pop, ]

## ------------------------------------------------------------------------
sort(states$Population)
# not what we want, thus
sort1.states <- states[order(states$Population), ]
head(sort1.states)
# sort by two variables
sort2.states <- states[order(states$Illiteracy, states$Income), ]
head(sort2.states)
# sort in reverse order
sort3.states <- states[order(-states$Life.Exp), ]
head(sort3.states)

## ------------------------------------------------------------------------
mydat <- data.frame(id = factor(1:12), 
                    group = factor(rep(1:2, each = 3)))
str(mydat)
head(mydat)
x <- rnorm(12)
y <- sample(70:100, 12)
x2 <- rnorm(12)
# add a column
mydat$grade <- y  
head(mydat)

## ------------------------------------------------------------------------
df <- data.frame(id = mydat$id, y)
head(df)
mydat2 <- merge(mydat, df, by = "id", sort = F) # using merge
head(mydat2)
mydat3 <- cbind(mydat, x) # using cbind()
head(mydat3)

## ------------------------------------------------------------------------
# add rows
df <- data.frame(id = factor(13:24), 
                 group = factor(rep(1:2, e = 3)), grade = sample(y))
df
mydat2 <- rbind(mydat, df)
mydat2

## ------------------------------------------------------------------------
library(ggplot2)
head(diamonds)
# average price for each type of cut
aggregate(price ~ cut, diamonds, mean)
# do a manual check, using `subset()`
s <- subset(diamonds, cut == 'Fair')
mean(s$price)
# add arguments to the function called
aggregate(price ~ cut, diamonds, mean, na.rm=TRUE)

## ------------------------------------------------------------------------
s <- aggregate(price ~ cut, diamonds, mean)
s
dd <- merge(diamonds, s, by="cut", sort = "FALSE")
head(dd)
head(diamonds)
head(subset(dd, cut == "Very Good"))
# change name of the column
names(dd)[names(dd) == 'price.y'] <- 'average price'
# add additional grouping variable
aggregate(price ~ cut + color, diamonds, mean, na.rm=TRUE)
# store results in an object
res <- aggregate(price ~ cut + color, diamonds, mean, na.rm=TRUE)
str(res)
head(res)
# aggregate two variables, combine with 'cbind'
aggregate(cbind(price, carat) ~ cut, diamonds, mean)
aggregate(cbind(price, carat) ~ cut + color, diamonds, mean)

## ---- warning = FALSE, message = FALSE-----------------------------------
library(AER)
data("CPS1985")
str(CPS1985)
head(CPS1985) 
summary(CPS1985$wage)

## ------------------------------------------------------------------------
# attach the data set; R knows where to find the variables
attach(CPS1985)
summary(wage)
is.numeric(wage)
mean(wage)
median(wage)
fivenum(wage)	# Tukey's five number summary
min(wage)
max(wage)
var(wage)
sd(wage)
hist(wage, freq = FALSE)
hist(log(wage), freq=FALSE, nclass=20, col="pink")
lines(density(log(wage)), col=4)
detach(CPS1985)

## ------------------------------------------------------------------------
attach(CPS1985)
str(occupation) # factor variable with 6 levels
summary(occupation)
detach(CPS1985)

## ------------------------------------------------------------------------
levels(CPS1985$occupation)[c(2, 6)] <- c("techn", "mgmt")
summary(CPS1985$occupation)

## ------------------------------------------------------------------------
attach(CPS1985)
tab <- table(occupation)
tab
prop.table(tab)
barplot(tab)
pie(tab)
pie(tab,col = gray(seq(0.4, 1.0, length = 6)))
detach(CPS1985)

## ------------------------------------------------------------------------
attach(CPS1985)
table(gender, occupation)
prop.table(table(gender, occupation))
prop.table(table(gender, occupation), 2)
# use mosaic plot 
plot(gender ~ occupation, data = CPS1985)
detach(CPS1985)

## ------------------------------------------------------------------------
attach(CPS1985)
# here: apply 'mean(.)' to 'log(wage)' by 'gender'
tapply(wage, gender, mean)
options(digits=5)
tapply(log(wage), list(gender, occupation), mean)
detach(CPS1985)
# let's check these results
# use subset(.) to extract part of the data
s <- subset(CPS1985, select=c(gender, occupation, wage))
s1 <- subset(s, gender == "female" & occupation == "techn")
mean(log(s1$wage))

## ------------------------------------------------------------------------
attach(CPS1985)
# see e.g. http://www.r-bloggers.com/box-plot-with-r-tutorial/
boxplot(log(wage) ~ gender)
boxplot(log(wage) ~ gender + occupation, col="light blue")
boxplot(log(wage) ~ gender + occupation, col="light blue", las=2)
# make it a nice graph
.pardefault <- par(no.readonly = T) # to store the default settings of par(.)
boxplot(log(wage) ~ gender + occupation, col="light blue", las=2, par(mar = c(12, 5, 4, 2) + 0.1))
par(.pardefault)
detach(CPS1985)

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
  geom_point(shape=1, alpha = 1/2)+
  geom_smooth() 
# shorter
ggplot(mtcars, aes(x = hp, y = mpg)) +
  geom_point(shape=1, alpha = 1/2)+
  geom_smooth() 
# use black and white lay-out
ggplot(mtcars, aes(x = hp, y = mpg)) + theme_bw() +
  geom_point(shape=1, alpha = 1/2)+ 
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
ggplot(mtcars, aes(factor(cyl), mpg))+
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

## ------------------------------------------------------------------------
library(ggplot2)
diamonds

## ------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
diamonds %>% filter(cut == "Ideal")

## ------------------------------------------------------------------------
diamonds %>% filter(cut == "Ideal" & color == "E")
diamonds %>% filter(cut == "Ideal" & color == c("E", "D"))

## ------------------------------------------------------------------------
diamonds %>% filter(cut == "Ideal") %>% summarize(mean = mean(price), std_dev = sd(price))

## ------------------------------------------------------------------------
# base R way with aggregate
aggregate(price ~ cut, diamonds, mean)

## ------------------------------------------------------------------------
diamonds %>% group_by(cut) %>% summarize(mean = mean(price))

## ------------------------------------------------------------------------
diamonds %>% group_by(cut, color) %>% summarize(price = mean(price))

## ------------------------------------------------------------------------
diamonds %>% group_by(cut) %>% summarize(price = mean(price), carat = mean(carat))

## ------------------------------------------------------------------------
diamonds %>% group_by(cut, color) %>% summarize(price = mean(price), carat = mean(carat))

## ------------------------------------------------------------------------
d <- diamonds %>% group_by(cut) %>% summarize(price = mean(price), carat = mean(carat))
new_diamonds <- diamonds %>% inner_join(d, by = "cut")
View(diamonds)
View(new_diamonds)

## ------------------------------------------------------------------------
library(data.table)
library(ggplot2)
str(diamonds)
diamonds_DT <- data.table(diamonds)
diamonds_DT # notice intelligent printing of this DT
summary(diamonds_DT$cut)

## ------------------------------------------------------------------------
# key is used to index the data.table and will provide the extra speed
setkey(diamonds_DT, cut)
tables()
diamonds_DT[J("Ideal"), ]
# more than one column can be set as key
setkey(diamonds_DT, cut, color)
tables()
# access rows according to both keys, use function 'J'
diamonds_DT[J("Ideal", "E"), ]
diamonds_DT[J("Ideal", c("E", "D")), ]
# what would be the alternative with base R?
subset(diamonds, diamonds$cut=="Ideal" && diamonds$color==c("E", "D"))

## ------------------------------------------------------------------------
# base R way with aggregate
aggregate(price ~ cut, diamonds, mean)
system.time(aggregate(price ~ cut, diamonds, mean))
# aggregation with data.table
# will go faster thanks to indexing
diamonds_DT[ , mean(price), by=cut]
system.time(diamonds_DT[ , mean(price), by=cut])
# give variable names in the create date.table
diamonds_DT[ , list(price = mean(price)), by=cut]
# aggregate on multiple columns
diamonds_DT[ , mean(price), by=list(cut,color)]
# aggregate multiple arguments
diamonds_DT[ , list(price = mean(price), carat = mean(carat)), by = cut]
diamonds_DT[ , list(price = mean(price), carat = mean(carat), caratSum = sum(carat)), by=cut]
# multiple metrics and multiple grouping variables
diamonds_DT[ , list(price = mean(price), carat = mean(carat)), by = list(cut, color)]

## ------------------------------------------------------------------------
# join two data.tables
d <- diamonds_DT[ , list(price = mean(price), carat = mean(carat)), by = cut]
d
setkey(diamonds_DT, cut)
dmerge <- diamonds_DT[d]
dmerge

## ------------------------------------------------------------------------
nSim 	   <- 100
p        <- 0.3
n	       <- 6

# generate 'nSim' obs. from Bin(n,p) distribution 
data_binom <- rbinom(nSim, n, p)

# calculate mean and variance
mean(data_binom) # empirical mean
var(data_binom)  # empirical variance

n*p 		      # theoretical mean
n*p*(1-p)	    # theoretical variance

# visualize
range <- seq(-1,n,1/1000)
plot(ecdf(data_binom))   # ecdf
lines(range,pbinom(range, n, p), col = 'red') # cdf

par(mfrow=c(1,2))
plot(0:n, dbinom(0:n, n, p), type = 'h') # pdf
plot(prop.table(table(data_binom)))
par(mfrow=c(1,1))

## ------------------------------------------------------------------------
nSim 	 <- 100
lambda <- 1

# generate 'nSim' observations from Poisson(\lambda) distribution
data_pois <- rpois(nSim, lambda)

# calculate mean and variance
mean(data_pois) # empirical mean
var(data_pois)  # empirical variance

lambda	    # theoretical mean
lambda	    # theoretical variance

# visualize
range  <- seq(0,8, 1/1000)
plot(ecdf(data_pois))   # ecdf
lines(range,ppois(range, lambda), col = 'red') # cdf

par(mfrow=c(1,2))
plot(0:8, dpois(0:8, lambda), type = 'h') # pdf
plot(prop.table(table(data_pois)))
par(mfrow=c(1,1))

## ------------------------------------------------------------------------
# evaluate cdf of N(0,1) in 0
pnorm(0, mean=0, sd=1)
# or shorter
pnorm(0, 0, 1)
# 95% quantile of N(0,1) 
qnorm(0.95, mean=0, sd=1)
# a set of quantiles
qnorm(c(0.025, 0.05, 0.5, 0.95, 0.975), 0, 1)
# generate observations from N(0,1)
x <- rnorm(10000, mean=10, sd=1)
# visualize
hist(x, probability=TRUE, nclass=55, col="pink")
curve(dnorm(x, mean=10, sd=1), xlim=range(x), col="black",add=TRUE)

## ------------------------------------------------------------------------
# check parametrization of gamma density in R
? dgamma
# grid of points to evaluate the gamma density
x <- seq(from = 0, to = 20, by = 0.001)
# choose a color palette
colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# shape and rate parameter combinations shown in the plot
shape <- c(1, 2, 3)
rate <- c(0.5, 0.5, 0.5)
plot(x, dgamma(x, shape = shape[1], rate = rate[1]), type='l', xlab ='x', ylab='Gamma density', main='Effect of the shape parameter on the Gamma density')
for(i in 2:length(shape)){
    lines(x, dgamma(x, shape = shape[i], rate = rate[i]), col=colors[i])
}
# add a legend  
legend("topright", paste("shape = ", shape, ", rate = ", rate, sep=""), col = colors, lty=1)

## ------------------------------------------------------------------------
3 == (2 + 1)
"intermediate" != "r"
TRUE != FALSE
"Rchitect" != "rchitect"

## ------------------------------------------------------------------------
(1 + 2) > 4
"dog" < "Cats"
TRUE <= FALSE

## ------------------------------------------------------------------------
katrien <- c(19, 22, 4, 5, 7)
katrien > 5
jan <- c(34, 55, 76, 25, 4)
jan <= 30

## ------------------------------------------------------------------------
TRUE & TRUE
FALSE | TRUE
5 <= 5 & 2 < 3
3 < 4 | 7 < 6

## ------------------------------------------------------------------------
katrien > 5 & jan <= 30

## ------------------------------------------------------------------------
!TRUE
!(5 > 3)
!!FALSE

## ------------------------------------------------------------------------
num_attendees <- 30
if (num_attendees > 5) {
  print("You're popular!")
}

## ------------------------------------------------------------------------
num_attendees <- 5
if (num_attendees > 5) {
  print("You're popular!")
}else{
  print("You are not so popular!")
}

## ------------------------------------------------------------------------
num_attendees <- 5
num_questions <- 3
if (num_attendees <= 5 & num_questions <=2) {
  print("Easy workshop")
}else{
  print("Work to do!")
}

## ------------------------------------------------------------------------
todo <- 64

while (todo > 30) {
  print("Work harder")
  todo <- todo - 7
}

todo

## ------------------------------------------------------------------------
i <- 1

while (i <= 10) {
  print(3 * i)
  if ( (3 * i) %% 8 == 0) {
    break
  }
  i <- i + 1
}

## ------------------------------------------------------------------------
primes <- c(2, 3, 5, 7, 11, 13)

# loop version 1
for (p in primes) {
  print(p)
}

# loop version 2
for (i in 1:length(primes)) {
  print(primes[i])
}

## ------------------------------------------------------------------------
rquote <- "r's internals are irrefutably intriguing"
chars <- strsplit(rquote, split = "")[[1]]
chars
# Initialize rcount
rcount <- 0

# Finish the for loop
for (char in chars) {
  if (char == "r") {
    rcount <- rcount + 1
  }
  if (char == "u") {
    break
  }
}

# Print out rcount
rcount

## ------------------------------------------------------------------------
? mean
help(mean)
args(mean)

## ------------------------------------------------------------------------
katrien <- c(2, 9, 6, 8, NA)

mean(katrien)

mean(katrien, na.rm = TRUE)

## ------------------------------------------------------------------------
katrien <- c(2, 9, 6, 8, NA)
jan <- c(0, 3, 2, NA, 5)
katrien - jan
mean(abs(katrien - jan), na.rm = TRUE)

## ------------------------------------------------------------------------
my_sqrt <- function(x) {
  sqrt(x)
}

# Use the function
my_sqrt(12)
my_sqrt(16)

sum_abs <- function(x, y) {
  abs(x) + abs(y)
}

# Use the function
sum_abs(-2, 3)

## ------------------------------------------------------------------------
my_sqrt <- function(x, print_info = TRUE) {
  y <- sqrt(x)
  if (print_info) {
    print(paste("sqrt", x, "equals", y))
  }
  return(y)
}

# some calls of the function
my_sqrt(16)
my_sqrt(16, FALSE)
my_sqrt(16, TRUE)

## ------------------------------------------------------------------------
v <- c(16, 25, 36)
my_sqrt(v)

## ------------------------------------------------------------------------
my_matrix <- matrix(1:9, nrow = 3)
# sum the rows
apply(my_matrix, 1, sum)
# sum the columns
apply(my_matrix, 2, sum)
# impute a missing observation in my_matrix
my_matrix[2,1] <- NA
apply(my_matrix, 1, sum)
apply(my_matrix, 1, sum, na.rm = TRUE)

## ------------------------------------------------------------------------
wages <- c(5500, 3500, 6500, 7500)
gender <- c("F", "F", "M", "M")
region <- c("North", "South", "North", "South")
salary <- data.frame(wages, gender, region)
tapply(salary$wages, salary$gender, mean)
tapply(salary$wages, list(salary$gender, salary$region), mean)

## ------------------------------------------------------------------------
my_list <- list(A = matrix(1:9, 3), B = 1:5, C = matrix(1:4, 2), D = 2)
my_list
lapply(my_list, sum)
sapply(my_list, sum)
my_names <- c("Katrien", "Jan", "Leen")
lapply(my_names, nchar)
sapply(my_names, nchar)

## ------------------------------------------------------------------------
# in one line of code
uniroot(function(x) x^2-3^(-x), lower=0, upper=1)
? uniroot
# in more lines of code
f <- function(x){
	x^2-3^(-x)
}
# calculate root
opt <- uniroot(f, lower=0, upper=1)
# check arguments
names(opt)
# evaluate 'f(.)' in the root
f(opt$root)
# visualize the function
range <- seq(-2, 2, by=0.2)
plot(range, f(range), type="l")
points(opt$root, f(opt$root), pch=20)
segments(opt$root, -7, opt$root, 0, lty=2)
segments(-3, 0, opt$root, 0, lty=2)

## ------------------------------------------------------------------------
# visualize the density
shape1 <- 3
shape2 <- 2
x <- seq(from=0, to=1, by=0.01)
curve(dbeta(x,shape1,shape2), xlim=range(x))

opt_beta <- optimize(dbeta, interval=c(0,1), maximum=TRUE, shape1, shape2)
points(opt_beta$maximum, opt_beta$objective, pch=20, cex=1.5)
segments(opt_beta$maximum, 0, opt_beta$maximum, opt_beta$objective, lty=2)

## ------------------------------------------------------------------------
nsim <- 10000
x <- rgamma(nsim, shape=3, rate=1.5)
# calculate log-likelihood
f <- function(p,x){
	-sum(dgamma(x, shape=p[1], rate=p[2], log=TRUE))
}
nlm(f, c(1, 1), x=x)

# same example, now use 'optim'
optim(c(1, 1), f, x=x)

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
# use built-in plot function
# you may have noticed that we have used the function plot with all kinds of arguments: 
# one or two variables, a data frame, and now a linear model fit;
# in R jargon plot is a generic function; it checks for the kind of object that you # are plotting and then calls the appropriate (more specialized) function to do the work.
plot(lm1)

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
 geom_smooth()+geom_abline(intercept = lm1$coef[1], slope = lm1$coef[2],   colour="red", size=2) 

## ------------------------------------------------------------------------
plot(car_price$income, car_price$price, pch=21, cex=1.2, xlab = "income", main = "Simple linear regression")
abline(lm1, col = "blue", lwd=2)
segments(car_price$income, car_price$price, car_price$income, lm1$fitted.values, lty=1)

## ------------------------------------------------------------------------
# anova table
anova(lm1)
attach(car_price)
total.SS <- sum((price-mean(price))^2)
total.SS
error.SS <- sum(lm1$resid^2)
error.SS

# R^2?
(total.SS-error.SS)/total.SS

# F-statistic in anova?
lm0 <- lm(price ~ 1)
error0.SS <- sum(lm0$resid^2)

# calculate F-statistic
F <- ((anova(lm0)$"Sum Sq")-(anova(lm1)$"Sum Sq"[2]))/(anova(lm1)$"Mean Sq"[2]) 
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
anova(lm1)

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

