###############2222
load("C:/Users/dell/Downloads/dataex2(5).Rdata")
install.packages("maxLik")
require(maxLik)
#loglikelihood function
log_like <- function(param, data){
  #Define variable
  x <- data[,1]; r <- data[,2]
  μ <- param[1]
  #formula
  sum(r*log(dnorm(x, mean=μ, sd=1.5)) + (1-r)*log(pnorm(x, mean=μ, sd=1.5)))
}
# maximum likelihood estimate
mle <- maxLik(logLik = log_like, data = dataex2, start = c(μ=0))
summary(mle)

###################4
load("C:/Users/dell/Downloads/dataex4(5) (1).Rdata")
#EM algorithm
em.two <- function(data, theta0, eps){
  theta <- theta0
  beta0t <- theta[1];  beta1t <- theta[2]
  diff <- 1
  
  while(diff > eps){
    theta.old <- theta
    #E-step
    log_likelihood_exp <- function(param, x, y){
      #Distinguish between two types of variables
      obs <- which(is.na(y) == FALSE)
      xobs <- x[obs]; xmis <- x[-obs]
      yobs <- y[obs]; 
      beta0 <- param[1]
      beta1 <- param[2]
      #loglikelihood-E step
      sum(yobs*(beta0+beta1*xobs))-sum(log(1+exp(beta0+beta1*x)))+sum((beta0+beta1*xmis)*exp(beta0t+xmis*beta1t)/(1+exp(beta0t+xmis*beta1t)))
    }
    #maximum likelihood estimates
    mle_optim <- optim(c(beta0t,beta1t), log_likelihood_exp, x=dataex4[ ,1], y=dataex4[ ,2],
                       control=list("fnscale"=-1),hessian=TRUE)
    theta <- mle_optim$par
    #several updates
    beta0t <- theta[1]
    beta1t <- theta[2]
    diff <- sum(abs(theta - theta.old))
  }
  return(theta)
}
#compute the maximum likelihood estimate of β
res <- em.two(dataex4, c(1,1), eps = 0.0001)
beta0 <- res[1]; beta1 <- res[2]
beta0;beta1

############################5
load("C:/Users/dell/Downloads/dataex5(3).Rdata")
# EM algorithm
em.mixture.two <- function(y, theta0, eps){
  theta <- theta0
  p <- theta[1];  lambda1 <- theta[2]; μ <- theta[3]
  diff <- 1
  while(diff > eps){
    theta.old <- theta
    #E-step
    ptilde1 <- p*lambda1/y^(lambda1+1)
    ptilde2 <- (1-p)*μ/y^(μ+1)
    ptilde <- ptilde1/(ptilde1 + ptilde2)
    #M-step
    p <- mean(ptilde)
    lambda1 <- sum(ptilde)/sum(log(y)*ptilde)
    μ <- sum(1-ptilde)/sum((1-ptilde)*log(y))
    theta <- c(p, lambda1, μ)
    diff <- sum(abs(theta-theta.old))
  }
  return(theta)
}
#implement the algorithm and find the maximum likelihood estimates
res <- em.mixture.two(y = dataex5, theta0 = c(0.3, 0.3, 0.4), eps = 0.0001)
pest <- res[1]; lambdaest <- res[2]; μest <- res[3]
pest; lambdaest; μest

#Draw the histogram of the data with the estimated density superimposed
hist(dataex5, breaks= 'freedman-diaconis',  xlab = "Y", xlim = c(0, 10), 
     main = "Mixture distribution Y" , prob = TRUE)
#estimated density
d <- function(y){
  return(pest*lambdaest*y^(-lambdaest-1)+(1-pest)*μest*y^(-μest-1))
}
curve(d , col = 4, lwd = 2, add = TRUE)
