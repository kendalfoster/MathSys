##############################
###### PS918 Workshop 1 ######
####### Kendal Foster ########
##############################



# Preamble ----------------------------------------------------------------
## Need to install packages "rtdists" and "tidyverse" ##
library("rtdists")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank(), 
                  panel.grid.minor.x = element_blank()))




# Preliminaries -----------------------------------------------------------
## Simulate data ##
set.seed(23)
N <- 500 # number of simulated responses
v1 <- 1.5 # drift rate
a1 <- 1.2 # boundary separation
t01 <- 0.3 # no-decision time
d1 <- rdiffusion(N, a=a1, v=v1, t0=t01)
str(d1)

## Plot the simulated data as a histogram
ggplot(d1, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")

## Function to generate the theoretical density
get_diffusion_df_old <- function(a, v, t0, z = a/2, sz = 0, sv = 0, st0 = 0, n = 401, start = 0, end = 3) {
  rt <- seq(start, end, length.out = n)
  tibble(
    rt = rt,
    upper = ddiffusion(rt, response = "upper", a = a, v = v, t0 = t0, z = z, sz = sz, sv = sv, st0 = st0),
    lower = ddiffusion(rt, response = "lower", a = a, v = v, t0 = t0, z = z, sz = sz, sv = sv, st0 = st0)
  ) %>%
    gather(key = "response", value = "Density", upper, lower) %>%
    mutate(response = factor(response, levels = c("lower", "upper")))
}

## Plot the theoretical PDF over the simulated data
pred_d1 <- get_diffusion_df_old(a = a1, v = v1, t0 = t01)
ggplot(d1, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  geom_line(data = pred_d1, aes(y = Density, group = 1)) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")

## Empirical statistics
d1 %>%
  group_by(response) %>%
  summarise(acc = n()/nrow(d1),
            q10 = quantile(rt, probs = 0.1),
            q50 = quantile(rt, probs = 0.5),
            q90 = quantile(rt, probs = 0.9))

## Theoretical statistics
pdiffusion(Inf, response = "upper", a = a1, v = v1, t0 = t01)

## Theoretical response time quantiles
qdiffusion(c(0.1, 0.5, 0.9), response = "upper", a = a1, v = v1, t0 = t01, scale_p = TRUE)
qdiffusion(c(0.1, 0.5, 0.9), response = "lower", a = a1, v = v1, t0 = t01, scale_p = TRUE)





# Task 1: Vary Parameters -------------------------------------------------
## Redefine function to generate the theoretical density
get_diffusion_df <- function(a, v, t0, z, sz = 0, sv = 0, st0 = 0, n = 401, start = 0, end = 3) {
  rt <- seq(start, end, length.out = n)
  tibble(
    rt = rt,
    upper = ddiffusion(rt, response = "upper", a = a, v = v, t0 = t0, z = z, sz = sz, sv = sv, st0 = st0),
    lower = ddiffusion(rt, response = "lower", a = a, v = v, t0 = t0, z = z, sz = sz, sv = sv, st0 = st0)
  ) %>%
    gather(key = "response", value = "Density", upper, lower) %>%
    mutate(response = factor(response, levels = c("lower", "upper")))
}

## Experiment 1: Increase Drift Rate, v
v2 <- 2*v1
d2 <- rdiffusion(N, a=a1, v=v2, t0=t01)
pred_d2 <- get_diffusion_df(a = a1, v = v2, t0 = t01, z = a1/2)
ggplot(d2, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  geom_line(data = pred_d2, aes(y = Density, group = 1)) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")

## Experiment 2: Increase Boundary Separation, a
a2 <- 2*a1
d3 <- rdiffusion(N, a=a2, v=v1, t0=t01)
pred_d3 <- get_diffusion_df(a = a2, v = v1, t0 = t01, z = a2/2)
ggplot(d3, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  geom_line(data = pred_d3, aes(y = Density, group = 1)) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")

## Experiment 3: Increase Non-Decision Time, t0
t02 <- 2*t01
d4 <- rdiffusion(N, a=a1, v=v2, t0=t02)
pred_d4 <- get_diffusion_df(a = a1, v = v1, t0 = t02, z = a1/2)
ggplot(d4, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  geom_line(data = pred_d4, aes(y = Density, group = 1)) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")

## Experiment 4: Increase Bias, z
z2 <- a1/1.25
d5 <- rdiffusion(N, a=a1, v=v2, t0=t02) # may have to change z? Look in documentation
pred_d5 <- get_diffusion_df(a = a1, v = v1, t0 = t02, z = z2)
ggplot(d5, aes(x = rt)) +
  geom_histogram(aes(y=stat(count/sum(count)*1/width)), binwidth = 0.1) +
  geom_line(data = pred_d5, aes(y = Density, group = 1)) +
  facet_wrap(~response) +
  ylab("Density") + xlab("Response time (s)")



# Task 2: Maximum Likelihood Estimation -----------------------------------
## Set up parameters
par1 <- c(a = 0.8, v = 1, t0 = 0.1)

## Define likelihood function
ll_diffusion <- function(pars, rt, response) {
  densities <- tryCatch(
    ddiffusion(rt, response=response,
               a=pars[1], v=pars[2], t0=pars[3]),
    error = function(e) 0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

## Test likelihood function, using par1 and d1 from above; should get 1063
ll_diffusion(pars = par1, rt = d1$rt, response = d1$response)

## Finding the ML Estimate
res1 <- nlminb(par1, ll_diffusion,
               rt = d1$rt, response = d1$response,
               lower = c(0, 0, -Inf))
res1

## Avoiding Local Optima
get_start_values <- function() {
  c(
    a = runif(1, 0.5, 3),
    v = rnorm(1, 0, 2),
    t0 = runif(1, 0, 0.2)
  )
}

res2 <- rerun(5, nlminb(get_start_values(), ll_diffusion,
                        rt = d1$rt, response = d1$response,
                        lower = c(0, 0, -Inf))) %>%
  map_dfr(~as_tibble(cbind(t(.$par),
                           logLik = -.$objective,
                           convergence = .$convergence)))
res2

mle <- res2[which.max(res2$logLik), 1:3 ]
mle

## Detecting Identifiability Issues
round(max(res2$logLik), 3) == round(res2$logLik, 3)

which_max <- which(round(max(res2$logLik), 3) == round(res2$logLik, 3))
## exclude actual ML estimate:
which_max <- which_max[which_max != which.max(res2$logLik)]

## copy ML estimates
mle2 <- mle
## remove all estimates that are not identifiable:
mle2[abs(mle - res2[which_max[1], 1:3 ]) > 0.01] <- NA
mle2
