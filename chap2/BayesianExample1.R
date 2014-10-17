#-----------------------#
#
# Author: Napo
#
#-----------------------#

library(TeachingDemos)
set.seed(123)
mCows = rep(1, 100)
hCows = rep(0, 400)
Data = c(mCows, hCows)
Y = sample(Data)
hist(Y, main = "Mastitis distribution", xlab = "Infected/Healthy", col = 'blue')

# Classical inference
# ML approach

BMLE = function(rho, Y, n){
  Bern.likelihood = rho^(sum(Y))*(1 - rho)^(n - sum(Y))
}

rho.prop = seq(0, 1, by = 0.001)
lh = apply(matrix(rho.prop), 1, BMLE, Y, length(Y))
D = data.frame(cbind(rho.prop, lh))
D = D[order(D$lh),]

plot(rho.prop, lh, ylab = "Likelihood", type = 'b')
abline(v = tail(D$rho.prop, n = 1), col = 'red')
mean(Y)


# Bayesian inference
a = 10
b = 10
n = length(Y)
Prior.dist = rbeta(10000, a, b)
hist(Prior.dist, col = "gray")
aStar = sum(Y) + a
bStar = n - sum(Y) + b


Post.samples = rbeta(10000, aStar, bStar)
mean(Post.samples)
hist(Post.samples, col = 'grey', xlab = 'Posterior distribution for rho'); abline(v = mean(Post.samples), col = 'red')
plot.ts(Post.samples, xlab = 'Iteration'); abline(h = mean(Post.samples), col = 'red')

#Credible interval 95%
lwr = qbeta(0.025, aStar, bStar); lwr
upr = qbeta(0.975, aStar, bStar); upr

#HPD
hpd(qbeta, shape1 = aStar , shape2 = bStar, .95, tol = 0.0001)


#----------------------------------------#
# Normal Example on Rat weight at day 36 #
#----------------------------------------#
Y = c(320, 354, 328, 297, 323, 331, 305, 338, 376, 296, 
      352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
      334, 302, 302, 323, 331, 345, 333, 316, 291, 324)
n = length(Y)
hist(Y, col = 'gray'); abline(v = mean(Y), col  = 'red')
sum(Y)
mean(Y)
v = round(sd(Y), 1); v

w = 0.015
m = 275
prior.dist = rnorm(10000, m, 1/w)
hist(prior.dist, col = 'gray')

mStar = (v*sum(Y) + w*m)/(n*v + w)
vStar = n*v + w
post.dist = rnorm(10000, mStar, sqrt(vStar))
mean(post.dist)
hist(post.dist, col = 'gray', xlab = 'Posterior distribution for mu'); abline(v = mean(post.dist), col = 'red')
plot.ts(post.dist, xlab = 'Iteration'); abline(h = mean(post.dist), col = 'red')

# Credible intervals
lwr = qnorm(0.025, mStar, sqrt(vStar)); lwr
upr = qnorm(0.975, mStar, sqrt(vStar)); upr

# HPD regions
emp.hpd(post.dist, conf = 0.95)
