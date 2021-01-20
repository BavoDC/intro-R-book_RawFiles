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

