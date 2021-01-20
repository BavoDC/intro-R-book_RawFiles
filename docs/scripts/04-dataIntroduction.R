## ---- eval=FALSE---------------------------------------------------------
## # what is the current working directory?
## getwd()
## # which files are currently stored in my working directory?
## dir()

## ------------------------------------------------------------------------
# where are my data files?
path <- file.path('C:\\Users\\u0043788\\Dropbox\\IIR Machine learning en Data Science opleiding\\Bookdown R Intro\\data')
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
mysubset_2 <- subset(states, states$Life.Exp > 70)
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
install.packages("ggplot2")
library(ggplot2)
head(diamonds)
# average price for each type of cut
aggregate(price ~ cut, diamonds, mean)
# do a manual check, using `subset()`
s <- subset(diamonds, cut == 'Fair')
mean(s$price)
# add arguments to the function called
aggregate(price ~ cut, diamonds, mean, na.rm = TRUE)

## ------------------------------------------------------------------------
s <- aggregate(price ~ cut, diamonds, mean)
s
dd <- merge(diamonds, s, by = "cut", sort = "FALSE")
head(dd)
head(diamonds)
head(subset(dd, cut == "Very Good"))
# change name of the column
names(dd)[names(dd) == 'price.y'] <- 'average price'
# add additional grouping variable
aggregate(price ~ cut + color, diamonds, mean, na.rm = TRUE)
# store results in an object
res <- aggregate(price ~ cut + color, diamonds, mean, na.rm = TRUE)
str(res)
head(res)
# aggregate two variables, combine with 'cbind'
aggregate(cbind(price, carat) ~ cut, diamonds, mean)
aggregate(cbind(price, carat) ~ cut + color, diamonds, mean)

## ---- warning = FALSE, message = FALSE-----------------------------------
install.packages("AER")
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
aggregate(wage ~ gender, CPS1985, mean)
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

