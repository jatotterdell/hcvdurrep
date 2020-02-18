#' Simulate data
#' 
#' @param m The number of participants to sample
#' @param p The response probability vector
#' @param d The doses
#' @param w The allocation weights
#' @return A `list` of data
sim_dat <- function(m, p, w = rep(1 / length(p), length(p))) {
  if(!all.equal(length(p), length(w))) stop("w and d must have same length.")
  N <- length(p)
  n <- rmultinom(1, m, w)[, 1]
  y <- rbinom(N, n, p)
  dat <- list(N = N, d = 1:N, n = n, y = y) 
  return(dat)
}

#' Normalise a vector to sum to one
#' 
#' @param x A positive numeric vector
#' @return `x` normalised to sum to one
normalise <- function(x) x / sum(x)

#' Compile Stan model
#' 
#' @param mod The model to compile
#' @return A `stanmodel` object
create_mod <- function(mod) {
  rstan::stan_model(file = paste0("stan/", mod, ".stan"))
}

#' Simulate a duration-response trial
#' 
#' @param interim_seq The sample size at which each interim occurs.
#' @param model A list containing a compiled Stan model in `[[1]]` which takes data arguments: `N`, `n`, `y`, `d`
#' and corresponding prior data list in `[[2]]` which will be composed with the simulated trial data for model fitting.
#' The model should generate quantities `theta[d]` corresponding to the probability of response at duration `d`.
#' @return The simulated trial data and results
#' 
#' Alternatively, the number of failures at which each interim occurs if `use_failures = TRUE`. 
sim_trial <- function(
  interim_seq, 
  p, 
  model, 
  epsilon = 0.9, 
  delta = 0.05,
  rar = FALSE,
  use_failures = FALSE, 
  ...) {
  
  require(rstan)
  require(dplyr)
  require(tidybayes)
  
  n_max <- max(interim_seq)
  n_int <- diff(c(0, interim_seq))
  K <- length(interim_seq)
  N <- length(p)
  
  alloc_ratio <- vector("list", K)
  new_dat <- vector("list", K)
  dat <- vector("list", K)
  theta <- vector("list", K)
  prob_noninf <- vector("list", K)
  prob_eff <- vector("list", K)
  prob_med <- vector("list", K)
  prob_bnd <- vector("list", K)
  
  alloc_ratio[[1]] <- rep(1 / N, N)
  
  for(k in 1:K) {
    
    # Generate data
    new_dat[[k]] <- sim_dat(n_int, p, alloc_ratio[[k]])
    if(k == 1) {run_dat <- new_dat[[k]]} else {
      run_dat$n <- run_dat$n + new_dat[[k]]$n
      run_dat$y <- run_dat$y + new_dat[[k]]$y
    }
    
    # Fit model
    moddat <- tidybayes::compose_data(run_dat, model[[2]])
    fit <- rstan::sampling(model[[1]], data = moddat, pars = c("theta", "zeta"), ...)
    draws <- tidybayes::gather_draws(fit, theta[d], zeta[d]) 
    dat[[k]] <- as_tibble(run_dat)
    
    
    # Calculate quantities of interest #
    #----------------------------------#
    
    # Estimate summary
    theta[[k]] <- draws %>% 
      median_hdci()
    
    # Relative to maximum
    prob_noninf[[k]] <- draws %>%
      filter(.variable == "zeta") %>%
      summarise(prob_noninf = mean(.value >= -delta))
    
    # Efficacy of duration wrt epsilon
    prob_eff[[k]] <- draws %>%
      summarise(prob_eff = mean(.value >= epsilon))
    
    # MED wrt epsilon
    prob_med[[k]] <- draws %>%
      group_by(.variable, .draw) %>%
      mutate(med1 = .value >= epsilon,
             med2 = lag(.value, 1, 0) < epsilon,
             med3 = med1 & med2) %>%
      group_by(d) %>%
      summarise(prob_eff = mean(med1),
                prob_lageff = mean(med2),
                prob_med = mean(med3))
    
    # Probability bounded
    prob_bnd[[k]] <- draws %>%
      summarise(prob_bnd = mean(.value >= epsilon & .value < epsilon + 0.05))
    
    # Update allocation ratios according to MED prob
    # (if ALL prob_med = 0 or ALL prob_eff < alpha, stop?)
    if(k < K) {
      if(!rar) {
        alloc_ratio[[k + 1]] <- rep(1 / N, N)
      } else {
        alloc_ratio[[k + 1]] <- sqrt(prob_med[[k]]$prob_med)  
      }
    }
  }
  
  return(list(
    dat = bind_rows(dat, .id = "Interim"),
    theta = bind_rows(theta, .id = "Interim"),
    prob_noninf = bind_rows(prob_noninf, .id = "Interim"),
    prob_eff = bind_rows(prob_eff, .id = "Interim"),
    prob_med = bind_rows(prob_med, .id = "Interim"),
    prob_bnd = bind_rows(prob_bnd, .id = "Interim"),
    alloc_prob = normalise(alloc_ratio[[K]])
  ))
}

#' Simulate a duration-response trial scenario
#' 
#' @param sce_dat The scenario tibble
#' @param reps The number of trial replications
#' @param ... Additional arguments to functions
#' @return A tibble with the scenario data and simulation results
sim_scenario <- function(sce_dat, reps, ...) {
  require(parallel)
  
  res <- mclapply(1:nrow(sce_dat), function(i, ...) {
    tibble::enframe(
      sapply(1:reps, function(j)
                sim_trial(
                  interim_seq = sce_dat$seq_ss[[i]], 
                  p = sce_dat$p[[i]], 
                  model = sce_dat$mod[[i]], ...), simplify = F), 
      "Trial", "Result")
  }, ..., mc.cores = parallel::detectCores() - 1) 
  
  return(bind_cols(sce_dat, enframe(res, "ID", "Sim") %>% select(Sim)))
}


med_stopping_rule <- function(sce, eff_thres, noneff_thres) {
  
}

noninf_stopping_rule <- function(sce, noninf_thres) {
  
}