##############################
##### PS918 Assingment 1 #####
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

## Set Working Directory ##
#setwd("G:/My Drive/Scholastic/Grad School/Warwick/Term 2/PS918/Workshops/Week 3-4")
setwd("/home/kfoster/Documents/PS918/Workshops/Week 3-4")


# Data Processing ---------------------------------------------------------

## Import data ##
med <- read_csv("medical_dm.csv")

## Filter out bad response times ##
med <- filter(med, rt > 0.25 & rt < 2.5)

## Covert blast/non-blast response to upper/lower (for use with ddiffusion) ##
med$ul <- with(med, "filler")

for (i in 1:length(med$response)) {
  if (med$response[i] == "blast") {
    med$ul[i] <- "upper"
  } else {
    med$ul[i] <- "lower"
  }
}

## Make sure the conversion worked; should be 0 obs
med_err <- filter(med, ul != "upper" & ul != "lower")

# Split individuals into Novice and Expert categories
med_novice <- filter(med, group == "novice")
med_expert <- filter(med, group == "experienced" | group == "inexperienced")


# Maximum Likelihood Estimation -------------------------------------------
# Define likelihood function
ll_diffusion <- function(pars, rt, response) {
  densities <- tryCatch(
    ddiffusion(rt, response=response,
               a=pars[1], v=pars[2], t0=pars[3],
               z=pars[4], sv=pars[5], st0=pars[6], sz=0),
    error = function(e) 0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# Testing
par1 <- c(
  a = runif(1, 0.5, 3),
  v = rnorm(1, 0, 2),
  t0 = runif(1, 0, 0.5),
  z = runif(1, 0.4, 0.6),
  sv = runif(1, 0, 0.5),
  st0 = runif(1, 0, 0.5)
)
ll_diffusion(pars=par1, rt=med_expert$rt, response=med_expert$ul)

nlminb(par1, ll_diffusion,
       rt = med_expert$rt, response = med_expert$ul,
       lower = c(0, 0, -Inf))


# Model Fitting -----------------------------------------------------------

# Avoiding Local Optima
get_start_values <- function() {
  c(
    a = runif(1, 0.5, 3),
    v = rnorm(1, 0, 2),
    t0 = runif(1, 0, 0.5),
    z = runif(1, 0.4, 0.6),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5)
  )
}

# Experts
nruns <- 5
unique_experts <- sort(unique(med_expert$id))
res_expert <- matrix()
for (i in 1:length(unique_experts)) {
  med_temp <- filter(med_expert, id==unique_experts[i])
  res_temp <- rerun(nruns, nlminb(get_start_values(), ll_diffusion,
                                  rt = med_temp$rt, response = med_temp$ul,
                                  lower = c(0, 0, -Inf))) %>%
    map_dfr(~as_tibble(cbind(t(.$par),
                             logLik = -.$objective,
                             convergence = .$convergence)))
  res_expert <- cbind.data.frame(res_expert, res_temp)
}
res_expert <- subset(res_expert, select=-c(res_expert))

# Novices
nruns <- 5 # redefined here
unique_novices <- sort(unique(med_novice$id))
res_novice <- matrix()
for (i in 1:length(unique_novices)) {
  med_temp <- filter(med_novice, id==unique_novices[i])
  res_temp <- rerun(nruns, nlminb(get_start_values(), ll_diffusion,
                                  rt = med_temp$rt, response = med_temp$ul,
                                  lower = c(0, 0, -Inf))) %>%
    map_dfr(~as_tibble(cbind(t(.$par),
                             logLik = -.$objective,
                             convergence = .$convergence)))
  res_novice <- cbind.data.frame(res_novice, res_temp)
}
res_novice <- subset(res_novice, select=-c(res_novice))







# Testing Space -----------------------------------------------------------




filtered_diffusion <- function(pars, rt, response, classification) {
  idx <- which(classification == "blast")
  rt <- rt[c(idx)]
  response <- response[c(idx)]
  densities <- tryCatch(
    ddiffusion(rt, response=response,
               a=pars[1], v=pars[2], t0=pars[3],
               z=pars[4], sv=pars[5], st0=pars[6], sz=0),
    error = function(e) 0)
  if (any(densities == 0)) return(1e6)
  return(densities)
}


bnb_ll <- function(pars, rt, response, classification) {
  idx <- which(classification == "blast")
  rt_b <- rt[c(idx)]
  response_b <- response[c(idx)]
  densities_b <- tryCatch(
    ddiffusion(rt_b, response=response_b,
               a=pars[1], v=pars[2], t0=pars[4],
               z=pars[5], sv=pars[6], st0=pars[7], sz=0),
    error = function(e) 0)
  if (any(densities_b == 0)) return(1e6)
  rt_n <- rt[-c(idx)]
  response_n <- response[-c(idx)]
  densities_n <- tryCatch(
    ddiffusion(rt_n, response=response_n,
               a=pars[1], v=pars[3], t0=pars[4],
               z=pars[5], sv=pars[6], st0=pars[7], sz=0),
    error = function(e) 0)
  if (any(densities_n == 0)) return(1e6)
  dens <- append(densities_b, densities_n)
  return(-sum(log(dens)))
}

bnb_ll(pars=get_start_values_bnb(), rt=med_expert$rt, response=med_expert$ul, classification=med_expert$classification)


get_start_values_bnb <- function() {
  c(
    a = runif(1, 0.5, 3),
    v_b = rnorm(1, 0, 2),
    v_n = rnorm(1, 0, 2),
    t0 = runif(1, 0, 0.5),
    z = runif(1, 0.4, 0.6),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5)
  )
}


nlminb(get_start_values_bnb(), bnb_ll,
       rt = med_temp$rt, response = med_temp$ul, classification = med_temp$classification,
       lower = c(0, 0, -Inf))







