library(plant)
library(greta)
library(greta.dynamics)



# One species -------------------------------------------------------------
# Following https://traitecoevo.github.io/plant/articles/individuals.html
ind    <- FF16_Individual()
env    <- FF16_fixed_environment(1.0)
times  <- seq(0, 50, length.out = 101)
result <- grow_individual_to_time(ind, times, env)

# plot(times, result$state[, "height"], 
#      xlab="Time (years)", ylab="Height (m)")
# rates <- sapply(result$individual, function(x) x$rate("height"))
# plot(times, rates)
# plot(result$state[, "height"], rates)

# generate noisy data
y <- rlnorm(length(times),
            log(result$state[, "height"]),
            0.1)

# fit ODE model using greta
# the Zeide function for ODE
zeide <- function(y, t, a, b, c) {
    dy <- a * y^b * exp(-c*y)
    return(dy)
}

# priors
a <- normal(0, 2.5, truncation = c(0, Inf))
b <- normal(0, 2.5, truncation = c(0, Inf))
c <- normal(0, 2.5, truncation = c(0, Inf))

# initial values and observation error
y0 <- lognormal(log(0.5), 0.1)
sigma <- exponential(1)

# solution to the ODE
log_mu <- log(ode_solve(zeide, y0, times, a, b, c))
mu <- exp(log_mu)

# likelihood
distribution(y) <- lognormal(log_mu, sigma)
mod <- model(a, b, c, sigma, y0)

# inference using maximum a posteriori (we can switch to MCMC later)
est <- opt(mod)

# calculate fitted values
y_hat <- calculate(mu, values = est$par)

# compare (looks okay)
# we used the "wrong" phenomenological model and potentially inappropriate priors
plot(times, y)
lines(times, y_hat$mu)




# Multiple species varying in LMA -----------------------------------------

n_sp   <- 10   # number of species
params <- scm_base_parameters("FF16")
s1     <- strategy_list(trait_matrix(seq(0.05, 1, length.out = n_sp), "lma"), 
                        params, 
                        birth_rate_list = 0.1)  # birth_rate_list?
ind_list    <- lapply(s1, FF16_Individual)
result_list <- lapply(ind_list, function(x) grow_individual_to_time(x, times, env))
result_mat  <- sapply(result_list, function(x) x$state[, "height"])
# 

# for internal reasons I needed to loop across columns (species) rather than 
# fitting a single multispecies model, for now
est_list   <- list()
y_hat_list <- list()
y_list     <- list()
for (i in seq_len(n_sp)) {
    # generate noisy data
    y_list[[i]] <- rlnorm(length(times),
                          log(result_mat[, i]),
                          0.1)

    # priors
    a <- normal(0, 2.5, truncation = c(0, Inf))
    b <- normal(0, 2.5, truncation = c(0, Inf))
    c <- normal(0, 2.5, truncation = c(0, Inf))
    
    # initial values and observation error
    y0 <- lognormal(log(0.5), 0.1)
    sigma <- exponential(1)
    
    # solution to the ODE
    log_mu <- log(ode_solve(zeide, y0, times, a, b, c))
    mu <- exp(log_mu)
    
    # likelihood
    distribution(y_list[[i]]) <- lognormal(log_mu, sigma)
    mod <- model(a, b, c, sigma, y0)
    
    # inference using maximum a posteriori (we can switch to MCMC later)
    est_list[[i]] <- opt(mod)
    
    # calculate fitted values
    y_hat_list[[i]] <- calculate(mu, values = est_list[[i]]$par)
}

# examine the relationship between LMA and the Zeide parameters
par_list <- t(sapply(est_list, function(x) unlist(x$par)))
par_list <- cbind(LMA = sapply(s1, function(x) x$lma), par_list)
pairs(par_list[, 1:4])

# 
par(mfrow = c(1, 2))
matplot(result_mat, type = "l", ylim = c(0, 25))
matplot(simplify2array(y_list), ylim = c(0, 25), pch = 16)
matlines(simplify2array(lapply(y_hat_list, unlist)))

         