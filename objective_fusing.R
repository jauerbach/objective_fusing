library("MASS")
library("CVXR")

#########################
# Simulation Parameters #
#########################

n <- 50
p <- 20
lambdas <- c(0, .01, .1, 1, 10, 100)
B <- 1
n_sim <- 500

##############################
# Simulation when IV correct #
##############################

results <- matrix(ncol = 7, nrow = n_sim * length(lambdas))
colnames(results) <- c("lambda", "b_ols", "b_iv", "b_fac", 
                       "b_iv_1", "b_iv_2", "b_iv_3")

i <- 0
pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)
for(lambda in lambdas) {
  for(sim in 1:n_sim) {
    i <- i + 1
    set.seed(i)
    setTxtProgressBar(pb, i)
    
    a_1 <- rnorm(n) 
    b_1 <- rnorm(n)
    c_1 <- rnorm(n)
    d_1 <- rnorm(n)
    e_1 <- rnorm(n)
    a_2 <- rnorm(p) 
    b_2 <- rnorm(p)
    c_2 <- rnorm(p)
    d_2 <- rnorm(p)
    e_2 <- rnorm(p)
    
    u <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    v <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    w <- matrix(rnorm(n * p, sd = 1), ncol = p)
    e <- matrix(rnorm(n * p, sd = 1), ncol = p)
    
    X <- a_1 %*% t(a_2) + b_1 %*% t(b_2) + u
    A <- b_1 %*% t(b_2) + c_1 %*% t(c_2) + v
    Z <- a_1 %*% t(a_2) + d_1 %*% t(d_2) + w
    Y <- B * X + A + e
    
    b_ols <- cov(c(X), c(Y)) / cov(c(X), c(X))
    b_iv  <- cov(c(Z), c(Y)) / cov(c(Z), c(X))
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 +
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_fac <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 + 
                            lambda * norm_nuc(t(cbind(1, Z)) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_1 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(cbind(1, Z)) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_2 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(Z) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(t(Z) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_3 <- result[[1]][2]
    
    results[i, ] <- cbind(lambda, b_ols, b_iv, b_fac, b_iv_1, b_iv_2, b_iv_3)
  }
}

mse_summary_iv_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) sqrt(sum((x[results[,"lambda"] == lambda] - 1)^2)))
  )

bias_summary_iv_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) mean((x[results[,"lambda"] == lambda] - 1)))
  )

colnames(bias_summary_iv_correct)  <- colnames(mse_summary_iv_correct) <- lambdas
round(bias_summary_iv_correct, 3)
round(mse_summary_iv_correct, 3)

###############################
# Simulation when OLS correct #
###############################

results <- matrix(ncol = 7, nrow = n_sim * length(lambdas))
colnames(results) <- c("lambda", "b_ols", "b_iv", "b_fac", "b_iv_1", "b_iv_2", "b_iv_3")

i <- 0
pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)

for(lambda in lambdas) {
  for(sim in 1:n_sim) {
    i <- i + 1
    set.seed(i)
    setTxtProgressBar(pb, i)
    
    a_1 <- rnorm(n) 
    b_1 <- rnorm(n)
    c_1 <- rnorm(n)
    d_1 <- rnorm(n)
    e_1 <- rnorm(n)
    a_2 <- rnorm(p) 
    b_2 <- rnorm(p)
    c_2 <- rnorm(p)
    d_2 <- rnorm(p)
    e_2 <- rnorm(p)
    
    u <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    v <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    w <- matrix(rnorm(n * p, sd = 1), ncol = p)
    e <- matrix(rnorm(n * p, sd = 1), ncol = p)
    
    X <- a_1 %*% t(a_2) + b_1 %*% t(b_2) + u
    A <- c_1 %*% t(c_2) + d_1 %*% t(d_2) + v
    Z <- d_1 %*% t(d_2) + e_1 %*% t(e_2) + w
    Y <- B * X + A + e
    
    b_ols <- cov(c(X), c(Y)) / cov(c(X), c(X))
    b_iv  <- cov(c(Z), c(Y)) / cov(c(Z), c(X))
 
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 +
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_fac <- result[[1]][2]

    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 + 
                            lambda * norm_nuc(t(cbind(1, Z)) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_1 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(cbind(1, Z)) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_2 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(Z) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(t(Z) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_3 <- result[[1]][2]

    results[i, ] <- cbind(lambda, b_ols, b_iv, b_fac, b_iv_1, b_iv_2, b_iv_3)
  }
}

mse_summary_ols_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) sqrt(sum((x[results[,"lambda"] == lambda] - 1)^2)))
  )

bias_summary_ols_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) mean((x[results[,"lambda"] == lambda] - 1)))
  )

colnames(bias_summary_ols_correct)  <- colnames(mse_summary_ols_correct) <- lambdas
round(bias_summary_ols_correct, 3)
round(mse_summary_ols_correct, 3)

###################################
# Simulation when neither correct #
###################################

results <- matrix(ncol = 7, nrow = n_sim * length(lambdas))
colnames(results) <- c("lambda", "b_ols", "b_iv", "b_fac", "b_iv_1", "b_iv_2", "b_iv_3")

i <- 0
pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)
for(lambda in lambdas) {
  for(sim in 1:n_sim) {
    i <- i + 1
    set.seed(i)
    setTxtProgressBar(pb, i)

    a_1 <- rnorm(n) 
    b_1 <- rnorm(n)
    c_1 <- rnorm(n)
    d_1 <- rnorm(n)
    e_1 <- rnorm(n)
    a_2 <- rnorm(p) 
    b_2 <- rnorm(p)
    c_2 <- rnorm(p)
    d_2 <- rnorm(p)
    e_2 <- rnorm(p)
    
    u <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    v <- matrix(rnorm(n * p, sd = 1), ncol = p) 
    w <- matrix(rnorm(n * p, sd = 1), ncol = p)
    e <- matrix(rnorm(n * p, sd = 1), ncol = p)
    
    X <- a_1 %*% t(a_2) + b_1 %*% t(b_2) + u
    A <- b_1 %*% t(b_2) + c_1 %*% t(c_2) + v
    Z <- a_1 %*% t(a_2) + c_1 %*% t(c_2) + w
    Y <- B * X + A + e
    
    b_ols <- cov(c(X), c(Y)) / cov(c(X), c(X))
    b_iv  <- cov(c(Z), c(Y)) / cov(c(Z), c(X))
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 +
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_fac <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(Y - beta_hat[1] - X * beta_hat[2] - A_hat, "F")^2 + 
                            lambda * norm_nuc(t(Z) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_1 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(Z) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_2 <- result[[1]][2]
    
    beta_hat <- Variable(2, 1)
    A_hat <- Variable(n, p)
    objective <- Minimize(norm(t(Z) %*% (Y - beta_hat[1] - X * beta_hat[2] - A_hat), "F")^2 + 
                            lambda * norm_nuc(t(Z) %*% A_hat))
    problem <- Problem(objective)
    result <- solve(problem)
    b_iv_3 <- result[[1]][2]
    
    results[i, ] <- cbind(lambda, b_ols, b_iv, b_fac, b_iv_1, b_iv_2, b_iv_3)
  }
}

mse_summary_neither_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) sqrt(sum((x[results[,"lambda"] == lambda] - 1)^2)))
  )

bias_summary_neither_correct <-
  sapply(lambdas, 
         function(lambda) 
           apply(results[,-1], 2, 
                 function(x) mean((x[results[,"lambda"] == lambda] - 1)))
  )

colnames(bias_summary_neither_correct)  <- colnames(mse_summary_neither_correct) <- lambdas
round(bias_summary_neither_correct, 3)
round(mse_summary_neither_correct, 3)

#####################
# Make Latex Tables #
#####################

rownames(mse_summary_iv_correct)

library("kable")
library("kableExtra")
kable(mse_summary_iv_correct, format="latex", digits = 1)
kable(mse_summary_ols_correct, format="latex", digits = 1)
kable(mse_summary_neither_correct, format="latex", digits = 1)

