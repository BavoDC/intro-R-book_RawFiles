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
mtcars

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


