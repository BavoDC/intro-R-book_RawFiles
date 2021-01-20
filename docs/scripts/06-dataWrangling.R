## ------------------------------------------------------------------------
library(ggplot2)
diamonds

## ------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
diamonds %>% filter(cut == "Ideal")

## ------------------------------------------------------------------------
diamonds %>% filter(cut == "Ideal" & color == "E")
diamonds %>% filter(cut == "Ideal" & color %in% c("E", "D"))

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

