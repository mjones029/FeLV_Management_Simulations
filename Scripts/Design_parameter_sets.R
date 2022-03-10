## Design_parameter_sets.R
# 
#========================================================	
# ---
### title: Design parameter sets
# author: Marie Gilbertson
# date: "10/09/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Parameter set design for FeLV management simulations

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(AlgDesign)


#### set seed ####
set.seed(65213)



#### BASELINE POPULATION SIMULATION ####

### pop size ###
# based on McClintock pop estimates for ~2012
pop.size <- 150

### network density ###
# based on ideal FIV/FeLV manuscript random forest (doi: 10.1101/2021.01.09.426055)
net.den <- 0.08

### adult:subadult proportion ###
# min from FIV/FeLV manuscript random forest = 0.86 (as of 10/13/20)
# mean from FIV/FeLV manuscript = 0.91
a_prop <- 0.89


net.baselines <- c(pop.size, net.den, a_prop)


#### BASELINE FELV TRANSMISSION ####
# based on FIV/FeLV manuscript random forest (as of 10/13/20)
beta <- 0.27                          # base transmission rate (applies to all progressive individuals)
c.beta <- 0.1                         # constant applied to beta to determine transmission rate of latent/regressive individuals 
c.rate <- 0.3                         # "contact rate" which is just a lower probability of transmission occuring
prop.outcomes <- c(0.25, 0.25, 0.5)   # proportion per each infected state
progr.dur <- 1/18                     # weekly probability of mortality for progressives (=1/progressive infection duration in weeks)
regr.dur.c <- 1                       # regressive infection duration multiplier
terr.repop <- 0.11                    # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)



trans.baselines <- c(beta, c.beta, c.rate, prop.outcomes, progr.dur, regr.dur.c, terr.repop)



#### PROACTIVE VACCINATION ####

### proportion vaccinated at onset ###
prop.pro.vax <- seq(0.1, 0.8, by = 0.1)


### ratio of unboosted to boosted ###
vax.ratio <- c(1, 0, 0.5) # all unboosted, all boosted, 50/50 boosted


### vaccine efficacy (unboosted, boosted) ###
unboost.ve <- 0.4
boost.ve <- 0.8


### generate factorial design for proactive vaccination simulations ###
fact.provax.1 <- gen.factorial(levels = c(
                                     length(prop.pro.vax),
                                     length(vax.ratio)
                                      ), 
                        nVars = 2, center = F, 
                        varNames = c("prop.pro.vax", "vax.ratio"))

fact.provax.2 <- data.frame(prop.pro.vax = prop.pro.vax[fact.provax.1$prop.pro.vax],
                            vax.ratio = vax.ratio[fact.provax.1$vax.ratio]
)



#### REACTIVE VACCINATION ####

### proportion vaccinated at onsest ###
# 0, 0.2, 0.4, 0.6
prop.pro.vax_r <- seq(0.0, 0.6, by = 0.2)


### for proactive: ratio of unboosted to boosted ###
vax.ratio_r <- 0.5 # 50/50 boosted


### vaccine efficacy (unboosted, boosted) ###
unboost.ve_r <- 0.4
boost.ve_r <- 0.8


### onset of reactive vaccination ###
# 6mo or 1yr (in weeks)
onset_r <- c(26, 52) 


### distribution type ###
# random or spatial
dist.types <- c("random", "spatial")


### use vaccination window? (i.e. is time the rate-limiting factor?) ###
# 6 month window or no window
window.dur <- c(T, F)


### generate factorial design for proactive vaccination simulations ###
fact.reactvax.1 <- gen.factorial(levels = c(
  length(prop.pro.vax_r),
  length(vax.ratio_r),
  length(unboost.ve_r),
  length(boost.ve_r),
  length(onset_r),
  length(dist.types),
  length(window.dur)
), 
nVars = 7, center = F, 
varNames = c("prop.pro.vax", "vax.ratio", "unboost.ve", "boost.ve", "onset", "dist.type",
             "window.dur"))

fact.reactvax.2 <- data.frame(prop.pro.vax = prop.pro.vax_r[fact.reactvax.1$prop.pro.vax],
                            vax.ratio = vax.ratio_r[fact.reactvax.1$vax.ratio],
                            unboost.ve = unboost.ve_r[fact.reactvax.1$unboost.ve],
                            boost.ve = boost.ve_r[fact.reactvax.1$boost.ve],
                            onset = onset_r[fact.reactvax.1$onset],
                            dist.type = dist.types[fact.reactvax.1$dist.type],
                            window.dur = window.dur[fact.reactvax.1$window.dur]
)



#### TEST AND REMOVAL ####

### proportion vaccinated at onsest ###
# 0, 0.2, 0.4, 0.6
prop.pro.vax_t <- seq(0.0, 0.6, by = 0.2)


### for proactive: ratio of unboosted to boosted ###
vax.ratio_t <- 0.5 # 50/50 boosted


### vaccine efficacy (unboosted, boosted) ###
unboost.ve_t <- 0.4
boost.ve_t <- 0.8


### onset of reactive vaccination ###
# 6mo or 1yr (in weeks)
onset_t <- c(26, 52) 


### distribution type ###
# random or spatial
dist.types_t <- c("random", "spatial")


### maximum infected individuals that can be "pulled" at one time ###
max.occupancy <- c(5)


### generate factorial design for proactive vaccination simulations ###
fact.tandr.1 <- gen.factorial(levels = c(
  length(prop.pro.vax_t),
  length(vax.ratio_t),
  length(unboost.ve_t),
  length(boost.ve_t),
  length(onset_t),
  length(dist.types_t),
  length(max.occupancy)
), 
nVars = 7, center = F, 
varNames = c("prop.pro.vax", "vax.ratio", "unboost.ve", "boost.ve", "onset", "dist.type",
             "max.occupancy"))

fact.tandr.2 <- data.frame(prop.pro.vax = prop.pro.vax_t[fact.tandr.1$prop.pro.vax],
                              vax.ratio = vax.ratio_t[fact.tandr.1$vax.ratio],
                              unboost.ve = unboost.ve_t[fact.tandr.1$unboost.ve],
                              boost.ve = boost.ve_t[fact.tandr.1$boost.ve],
                              onset = onset_t[fact.tandr.1$onset],
                              dist.type = dist.types_t[fact.tandr.1$dist.type],
                              max.occupancy = max.occupancy[fact.tandr.1$max.occupancy]
)



#### UNDERPASS CLOSURES ####

### proportion vaccinated at onsest ###
# 0, 0.2, 0.4, 0.6
prop.pro.vax_u <- seq(0.0, 0.6, by = 0.2)


### for proactive: ratio of unboosted to boosted ###
vax.ratio_u <- 0.5 # 50/50 boosted


### vaccine efficacy (unboosted, boosted) ###
unboost.ve_u <- 0.4
boost.ve_u <- 0.8


### onset of reactive vaccination ###
# 6mo or 1yr (in weeks)
onset_u <- c(26, 52) 


### duration of closure ###
# 4 weeks, 3 mo, 6 mo, 1 yr (1 year is too long, but including for purposes of comparison)
upc.dur <- c(4, 13, 26, 52)



### generate factorial design for proactive vaccination simulations ###
fact.upc.1 <- gen.factorial(levels = c(
  length(prop.pro.vax_u),
  length(vax.ratio_u),
  length(unboost.ve_u),
  length(boost.ve_u),
  length(onset_u),
  length(upc.dur)
), 
nVars = 6, center = F, 
varNames = c("prop.pro.vax", "vax.ratio", "unboost.ve", "boost.ve", "onset", "upc.dur"))

fact.upc.2 <- data.frame(prop.pro.vax = prop.pro.vax_u[fact.upc.1$prop.pro.vax],
                           vax.ratio = vax.ratio_u[fact.upc.1$vax.ratio],
                           unboost.ve = unboost.ve_u[fact.upc.1$unboost.ve],
                           boost.ve = boost.ve_u[fact.upc.1$boost.ve],
                           onset = onset_u[fact.upc.1$onset],
                           upc.dur = upc.dur[fact.upc.1$upc.dur]
)




#### SAVE PARAMETER SETS ####

### baseline tranmsission simulation parameters ###
# (baseline parameters for no intervention scenarios)
# (will also be used in main intervention scenarios)
nointer.baseline <- data.frame(matrix(c(net.baselines, trans.baselines), nrow = 1))
colnames(nointer.baseline) <- c("pop.size", "net.den", "a_prop", "beta", "c.beta", "c.rate",
                                "prop.outcomes_p", "prop.outcomes_r", "prop.outcomes_a", 
                                "progr.dur", "regr.dur.c", "terr.repop")                    

# save(nointer.baseline, file = "Parameter_sets/nointer_baseline_params.Rdata")


### proactive vaccination simulation parameters ###

fact.provax.2$unboost.ve <- unboost.ve
fact.provax.2$boost.ve <- boost.ve
  
net_trans_baseline <- nointer.baseline[rep(seq_len(nrow(nointer.baseline)), nrow(fact.provax.2)),]

pro.vax.params <- cbind(fact.provax.2, net_trans_baseline)

# save(pro.vax.params, file = "Parameter_sets/proactivevax_baseline_params.Rdata")


### reactive vaccination simulation parameters ###

net_trans_baseline <- nointer.baseline[rep(seq_len(nrow(nointer.baseline)), nrow(fact.reactvax.2)),]

react.vax.params <- cbind(fact.reactvax.2, net_trans_baseline)

# save(react.vax.params, file = "Parameter_sets/reactivevax_baseline_params.Rdata")


### test and removal simulation parameters ###

net_trans_baseline <- nointer.baseline[rep(seq_len(nrow(nointer.baseline)), nrow(fact.tandr.2)),]

react.tandr.params <- cbind(fact.tandr.2, net_trans_baseline)

# save(react.tandr.params, file = "Parameter_sets/reactivetandr_baseline_params.Rdata")


### underpass closure simulation parameters ###

net_trans_baseline <- nointer.baseline[rep(seq_len(nrow(nointer.baseline)), nrow(fact.upc.2)),]

react.upc.params <- cbind(fact.upc.2, net_trans_baseline)

# save(react.upc.params, file = "Parameter_sets/reactiveupc_baseline_params.Rdata")




