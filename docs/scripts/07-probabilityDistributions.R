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

