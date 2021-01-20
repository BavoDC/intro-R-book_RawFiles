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
if (num_attendees <= 5 & num_questions <= 2) {
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

