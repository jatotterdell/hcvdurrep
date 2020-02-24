source("r/rar_duration_functions.R")

library(parallel)
library(plyr)
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(rstan)
library(tictoc)

# Compile models ----

mod0 <- create_mod("ind_dose")
prdat0 <- list(mu = 0, sigma = 5)
model0 <- list(mod0, prdat0)

mod1 <- create_mod("sigmoidal_emax")
prdat1 <- list(mu = rep(0, 4), sigma = rep(2.5, 4))
model1 <- list(mod1, prdat1)

mod2 <- create_mod("rw1mono")
prdat2 <- list(mu = 0, sigma = 2.5, tau_n = 1e-2, tau_m = 1e-2)
model2 <- list(mod2, prdat2)

mod3 <- create_mod("rw1monohs")
prdat3 <- list(mu = 0, sigma = 2.5, lambda = 1)
model3 <- list(mod3, prdat3)

# Setup scenarios ----

p1 <- setNames(rep(0.9, 5), paste0("p", 1:5))
p2 <- setNames(c(rep(0.80, 3), rep(0.95, 2)), paste0("p", 1:5))
p3 <- setNames(plogis(-1.5 + 2.1*log(4:8)), paste0("p", 1:5))
p4 <- setNames(c(0.825, 0.875, 0.925, 0.95, 0.975), paste0("p", 1:5))
p5 <- setNames(c(0.85, 0.875, 0.9, 0.925, 0.95), paste0("p", 1:5))

ptib <- tibble(
  Scenario = 1:5,
  p = list(p1, p2, p3, p4, p5))
modtib <- tibble(
  Model = c("Independent duration", "Sigmoidal Emax", "RW(1) Mono", "Locally Adaptive RW(1) Mono"),
  mod = list(model0, model1, model2, model3))
sstib <- tibble(
  max_ss = c(500, 1000, 2500), 
  seq_ss = list(seq(100, 500, 100), seq(200, 1000, 200), seq(500, 2500, 500)))
scenario_dat <- expand_grid(ptib, modtib, sstib, tibble(rar = c("none", "med", "mnd")))

reps <- 1000

# Independent duration model ----

# scenario_dat_sub <- scenario_dat %>% filter(Model == "Independent duration")
# res <- sim_scenario_alt(scenario_dat_sub, reps,
#                         refresh = 0, chains = 1, warmup = 500, iter = 3000)
# saveRDS(res, "out/rar_inddur2.rds")

# Locally adaptive mono RW(1) ----

scenario_dat_sub <- scenario_dat %>% filter(Model == "Locally Adaptive RW(1) Mono")
res <- sim_scenario_alt(scenario_dat_sub, reps,
                        refresh = 0, chains = 1, warmup = 500, iter = 3000)
saveRDS(res, "out/rar_rw1monohs2.rds")

