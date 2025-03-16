library(dirichletprocess)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(2)

library(dirichletprocess)


# Likelihood function
Likelihood.mvlogistic <- function(mdObj, x, theta) {
  nComp <- dim(theta[[1]])[3]  
  nObs <- nrow(x)
  d <- ncol(x)
  likelihoods <- matrix(0, nrow = nObs, ncol = nComp)
  
  for (i in 1:nComp) {
    mu_i <- theta[[1]][, , i] 
    s_i <- theta[[2]][, , i]  
    likelihoods[, i] <- apply(x, 1, function(row) {
      prod(exp(-(row - as.vector(mu_i)) / as.vector(s_i)) /
             (as.vector(s_i) * (1 + exp(-(row - as.vector(mu_i)) / as.vector(s_i)))^2))
    })
  }
  return(likelihoods)
}

# Prior draw function
PriorDraw.mvlogistic <- function(mdObj, n = 1) {
  priorParameters <- mdObj$priorParameters
  d <- length(priorParameters$mu0)
  
  # Draw mu: multinorm
  mu <- array(0, dim = c(d, 1, n))
  for (i in 1:n) {
    mu[, , i] <- matrix(mvtnorm::rmvnorm(1, mean = priorParameters$mu0, sigma = priorParameters$S0), ncol = 1)
  }
  
  # Draw s: inverse gamma
  s <- array(0, dim = c(d, 1, n))
  for (i in 1:n) {
    s_vec <- numeric(d)
    for (j in 1:d) {
      s_vec[j] <- 1 / rgamma(1, shape = priorParameters$a0, rate = priorParameters$b0)
    }
    s[, , i] <- matrix(s_vec, ncol = 1)
  }
  
  theta <- list(mu = mu, s = s)
  return(theta)
}

# Posterior draw 
PosteriorDraw.mvlogistic <- function(mdObj, x, n = 1, n_iter = 1000, burn_in = 500, ...) {
  priorParameters <- mdObj$priorParameters
  d <- ncol(x)
  
  # Initialize parameters
  # set mu to sample mean
  # Set s to sample standard deviation
  current_mu <- colMeans(x)
  current_s <- apply(x, 2, sd)
  current_s[current_s == 0] <- 1
  current_log_s <- log(current_s)
  
  log_posterior <- function(mu, log_s) {
    s <- exp(log_s)
    # Compute log likelihood
    log_lik <- sum(apply(x, 1, function(row) {
      val <- sum( -(row - mu) / s - log(s) - 2 * log(1 + exp(-(row - mu) / s)) )
      if (!is.finite(val)) return(-Inf)
      return(val)
    }))
      
    # Priors, assume iid
    log_prior_mu <- sum(dnorm(mu, mean = priorParameters$mu0, 
                              sd = sqrt(diag(priorParameters$S0)), log = TRUE))
    a0 <- priorParameters$a0; b0 <- priorParameters$b0
    log_prior_s <- sum(a0 * log(b0) - lgamma(a0) - a0 * log_s - b0 / exp(log_s))
    
    total <- log_lik + log_prior_mu + log_prior_s
    if (!is.finite(total)) total <- -Inf
    return(total)}
  
  chain_mu <- matrix(0, nrow = n_iter, ncol = d)
  chain_log_s <- matrix(0, nrow = n_iter, ncol = d)
  chain_mu[1, ] <- current_mu
  chain_log_s[1, ] <- current_log_s
  current_lp <- log_posterior(current_mu, current_log_s)
  
  # Proposal standard deviations 
  prop_sd_mu <- rep(0.1, d)
  prop_sd_log_s <- rep(0.1, d)
  
  for (iter in 2:n_iter) {
    # tuning parameters – adjust
    proposed_mu <- current_mu + rnorm(d, 0, prop_sd_mu)
    proposed_log_s <- current_log_s + rnorm(d, 0, prop_sd_log_s)
    proposed_lp <- log_posterior(proposed_mu, proposed_log_s)
    
    # reject proposals with infinite proposed_lp  
    if (!is.finite(proposed_lp)) {
      accept <- FALSE
    } else {
      accept_prob <- exp(proposed_lp - current_lp)
      if (is.na(accept_prob)) accept_prob <- 0  # safeguard if NA
      accept <- runif(1) < accept_prob
    }
    
    if (accept) {
      current_mu <- proposed_mu
      current_log_s <- proposed_log_s
      current_lp <- proposed_lp
    }
    
    chain_mu[iter, ] <- current_mu
    chain_log_s[iter, ] <- current_log_s
  }
  
  
  # Compute the posterior parameters: average of the post-burn-in
  final_mu <- colMeans(chain_mu[(burn_in + 1):n_iter, , drop = FALSE])
  final_s  <- exp(colMeans(chain_log_s[(burn_in + 1):n_iter, , drop = FALSE]))
  
  # Return the parameters 3D array
  mu_out <- array(matrix(final_mu, ncol = 1), dim = c(d, 1, n))
  s_out  <- array(matrix(final_s,  ncol = 1), dim = c(d, 1, n))
  return(list(mu = mu_out, s = s_out))
}

Predictive.mvlogistic <- function(mdObj, x, n_draws = 1000, ...) {
  theta_draws <- PosteriorDraw.mvlogistic(mdObj, x, n = n_draws, ...)
  lik <- Likelihood.mvlogistic(mdObj, x, theta_draws)
  pred <- rowMeans(lik)
  return(pred)
}

MvlogisticCreate <- function(priorParameters) {
  mdObj <- MixingDistribution("mvlogistic", priorParameters, "nonconjugate")
  return(mdObj)
}

DirichletProcessMvlogistic <- function(y,
                                       g0Priors,
                                       alphaPriors = c(2, 4),
                                       numInitialClusters = 1) {
  if (!is.matrix(y)) {
    y <- matrix(y, ncol = length(y))
  }
  
  # set defaults base priors
  if (missing(g0Priors)) {
    d <- ncol(y)
    g0Priors <- list(mu0 = rep(0, d),
                     S0 = diag(d),
                     a0 = 2,   
                     b0 = 2)  
  } 
  
  mdobj <- MvlogisticCreate(g0Priors)
  dpobj <- DirichletProcessCreate(y, mdobj, alphaPriors)
  dpobj <- Initialise(dpobj, numInitialClusters = numInitialClusters)
  
  return(dpobj)
}



data("iris")

# Standardize 
df_iris = scale(iris[,1:4])
dp_iris <- DirichletProcessMvlogistic(df_iris, numInitialClusters = nrow(iris))
dp_iris = Fit(dp_iris, 1000)
# Plot
plot(dp_iris)

