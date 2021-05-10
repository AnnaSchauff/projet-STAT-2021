# Generate Bivariate Gaussian Observations
gen_normal <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0){
  x <- rnorm(n, muX, sX)
  z <- rnorm(n, muY, sY * sqrt(1-rho^2))
  y <- sY/sX * rho * (x-muX) + z
  df <- data.frame(x=x, y=y)
  return(df)
}
hi gorgeous!!!

# Generate Discretized Bivariate Gaussian Observations
gen_discrete <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0){
  x <- rnorm(n, muX, sX)
  z <- rnorm(n, muY, sY * sqrt(1-rho^2))
  y <- sY/sX * rho * (x-muX) + z
  df <- data.frame(x=trunc(x*5)/5, y=trunc(y*5)/5)
  return(df)
}


# Generate data with outliers
gen_outliers <- function(n, muX=0, muY=0, sX=1, sY=1, rho=0, out=0.05, dev=1){
  df <- gen_normal(n, muX, muY, sX, sY, rho)
  
  if(exists(".Random.seed")){
    INITIALSEED <- .Random.seed
    lockBinding("INITIALSEED", environment())
    on.exit(.Random.seed <<- INITIALSEED)
  }
  outliers=sample(n,ceiling(out*n))
  values_from <- sample(n, ceiling(out*n))
  
  df$y[outliers] <- df$y[values_from] + (df$y[values_from]-mean(df$y)) * dev
  return(df)
}

# Generate Bivariate Data with Nonlinear Relationship
gen_nonlinear <- function(n, muX, muY, sX, sY, angle){
  rotmat <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),2,2)
  x <- rnorm(n, muX, sX)
  y <- x**2 + rnorm(n, muY, sY)
  df <- cbind(x,y)
  df <- df %*% rotmat
  df <- as.data.frame(df)
  names(df) <- c("x", "y")
  return(df)
}

# True correlation value for the nonlinear setting
cor_nonlinear <- function(muX, muY, sX, sY, angle){
  # muX and muY have no effect
  sQ <- sqrt(cos(angle)**2 * (2*sX**4 + 4*sX**2*muX**2 +  sY**2) + sin(angle)**2 * sX**2 + 2 * cos(angle)*sin(angle)*(2*sX**2*muX))
  sR <- sqrt(sin(angle)**2 * (2*sX**4 + 4*sX**2*muX**2 + sY**2) + cos(angle)**2 * sX**2 - 2 * cos(angle)*sin(angle)*(2*sX**2*muX))
  cov <- sin(angle)*cos(angle) * (sX**2 - (2*sX**4 + 4*sX**2*muX**2 + sY**2)) +
    (cos(angle)**2 - sin(angle)**2) * (2*sX**2*muX)
  cor <- round(cov/(sQ*sR), 4)
  return(cor)
}

# Bootstrap Confindence Interval
boot.CI <- function(x,y, B, method){
  n <- length(x)
  muX.est <- mean(x)
  muY.est <- mean(y)
  sX.est  <- sd(x)
  sY.est  <- sd(y)
  rho.est <- cor(x, y)
  
  rXY.boot <- rep(0, B)
  # set seed here!
  for(i in 1:B){
    data.boot <-   gen_normal(n=n, 
                          muX=muX.est, 
                          muY=muY.est, 
                          sX=sX.est, 
                          sY=sY.est, 
                          rho=rho.est)
    rXY.boot[i] <- cor(data.boot$x, data.boot$y, method=method)
  }
  
  # Get the .95 bootstrap percentile confidence interval
  boot.CI <- quantile(rXY.boot, probs=c(0.025, 0.975))
  boot.CI
}

# Nonparametric Bootstrap Confindence Interval
np.boot.CI <- function(x,y, B, method){
  n <- length(x)
  
  rXY.boot <- rep(0, B)
  # set seed here!
  for(i in 1:B){
    s <- sample(n, replace = T)
    rXY.boot[i] <- cor(x[s], y[s], method=method)
  }
  
  # Get the .95 bootstrap percentile confidence interval
  np.boot.CI <- quantile(rXY.boot, probs=c(0.025, 0.975))
  np.boot.CI
}


# The following function computes the confidence intervals of different types (type = "boot" or "npboot"),
# using different estimators (method = "pearson" or "spearman"), in different data
# settings (setting= "normal", "outliers", "discrete", or "nonlinear")
CI.sims <- function(n, setting, method, type, B, M, 
                   muX, muY, sX, sY, rho=NULL, out=NULL, dev=NULL, angle=NULL, level=0.8){
  if(!setting %in% c("normal", "discrete", "outliers", "nonlinear")) stop("Wrong argument `setting`")
  if(!method %in% c("pearson", "spearman")) stop("Wrong argument `method`")
  if(!type %in% c("boot", "npboot")) stop("Wrong argument `type`")
  
  if(setting=="normal"){
    if(is.null(rho)) stop("Set a value for rho.")
    gen_data <- function() gen_normal(n, muX, muY, sX, sY, rho)
  } else if(setting=="outliers"){
    if(is.null(rho)) stop("Set a value for rho.")
    if(is.null(out)) stop("Set a value for out.")
    if(is.null(dev)) stop("Set a value for dev.")
    if(out>1) stop("The value for out must be between 0 and 1.")
    gen_data <- function() gen_outliers(n, muX, muY, sX, sY, rho, out=out, dev=dev)
  } else if(setting=="discrete"){
    gen_data <- function() gen_discrete(n, muX, muY, sX, sY, rho)
  } else if(setting=="nonlinear"){
    if(is.null(angle)) stop("Set a value for angle.")
    gen_data <- function() gen_nonlinear(n, muX, muY, sX, sY, angle)
  }
  probs <- c((1-level)/2, level + (1-level )/2)
  
  if(type=="boot"){
    CI <- t(sapply(1:M, function(ni){
      df <- gen_data()
      muX <- mean(df$x)
      muY <- mean(df$y)
      sX <- sd(df$x)
      sY <- sd(df$y)
      rho <- cor(df$x,df$y)
      quantile(sapply(1:B, function(na){
      cor(gen_normal(n, muX, muY, sX, sY, rho), method=method)[2]
    }), probs=probs)}
    ))
  }else if (type=="npboot"){
    CI <- t(sapply(1:M, function(ni){
      df <- gen_data()
      quantile(sapply(1:B, function(na){
      s <- sample(n, replace=T)
      cor(df[s,], method=method)[2]
    }), probs=probs)}
    ))
  } else {
    stop("type must be one of ('boot', 'npboot'")
  }
  CI
}

# Examples
if(0){

# We will look at the following functions:
gen_normal
gen_outliers
gen_discrete
gen_nonlinear
cor_nonlinear
CI.sims

# Parameters:
SAMPLESEED <- 91637
SIMULATIONSEED <- 28528
n <- 130
muX <- -0.1
muY <- -1.7
sX  <- 1.1
sY  <- 0.6
rho <- 0.84
out <- 0.03
dev <- 3
angle <- 0.3


# Generate data, normal setting
set.seed(SAMPLESEED)
df <- gen_normal(n=n, muX=muX, muY=muY, sX=sX, sY=sY, rho=rho)
# true correlation coefficient in this setting
rho <- rho
plot(df, main=paste("Population correlation = ", round(rho, 3)))



# Generate data, outliers setting
set.seed(SAMPLESEED)
df <- gen_outliers(n=n, muX=muX, muY=muY, sX=sX, sY=sY, rho=rho, out=out, dev=dev)
# population correlation in this setting (ignoring outliers)
rho <- rho
plot(df, main=paste("Population correlation = ", round(rho, 3)))

# Generate data, discrete setting
set.seed(SAMPLESEED)
df <- gen_discrete(n=n, muX=muX, muY=muY, sX=sX, sY=sY, rho=rho)
# population correlation in this setting (ignoring discretization)
rho <- rho
plot(df, main=paste("Population correlation = ", round(rho, 3)))

# Generate data, nonlinear setting
n <- 1000
set.seed(SAMPLESEED)
df <- gen_nonlinear(n=n, muX=muX, muY=muY, sX=sX, sY=sY, angle=pi)
# population correlation in this setting 
rho <- cor_nonlinear(muX=muX, muY=muY, sX=sX, sY=sY, angle=angle)
plot(df, main=paste("Population correlation = ", round(rho, 3)))

# Compute simulated confidence interval, for normal setting, nonparametric bootstrap
# , spearman correlation coefficient
B <- 100  # number of bootstrap samples
M <- 1000  # number of monte carlo simulations

set.seed(SIMULATIONSEED)
CI.sim <- CI.sims(n = n,
              setting = "nonlinear", # one of "normal", "outliers", "discrete", or "nonlinear"
              method = "spearman",  # one of "pearson", "spearman"
              type = "boot",  # one of "boot", "npboot"
              B = B,
              M = M,
              muX = muX,
              muY = muY,
              sX = sX,
              sY = sY,
              rho = rho,
              out = out, # this parameter is ignored in the normal setting
              dev = dev, # this parameter is ignored in the normal setting
              angle = angle) # this parameter is ignored in the normal setting

# Compute coverage:
population_correlation <- rho  # this is true in the normal setting
coverage <- mean((CI.sim[,1]< population_correlation) * (CI.sim[,2]> population_correlation))
}
