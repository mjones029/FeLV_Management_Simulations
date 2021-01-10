## Nointer_sims.R
# 
#========================================================	
# ---
### title: Simuations of FeLV transmission in the absence of interventions
# author: Marie Gilbertson
# date: "10/10/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates FeLV transmission on FIV-based networks without management interventions

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(adehabitatHR) # kerneloverlap, mcp
library(sp) # coordinates() function
library(rgeos) # centroids
library(geosphere) # geographic distances between points (distm() function)
library(dplyr) # sample_n function
library(spatstat)
library(DCG) # as.symmetricAdjacencyMatrix function
library(ergm)
library(intergraph)





#### load external functions ####
source("Scripts/simulate_pop_ergm.R")
source("Scripts/trans_sim_basic.R")
source("Scripts/post_process_outbreak_data_basic.R")
source("Scripts/props affected_births included_basic.R")


#### load data ####
# load parameter sets
param.data <- get(load("Parameter_sets/nointer_baseline_params.Rdata"))

# load telemetry
study.telemetry <- get(load("Attribute_data/2010_2012_telemetry_deid.Rdata"))

# load ergm (main FIV ERGM from FIV/FeLV project)
best.ergm <- get(load("Attribute_data/bestERGM.Rdata"))



#### SET SIMULATION BOUNDS ####
# number of sims
nsims <- seq(1, 100)

# which parameter sets
# (not particularly relevant for baseline no-intervention scenarios)
param.sets <- seq(1, nrow(param.data))

# define scenario type for file naming
scenario.type <- "nointer_baseline" 


for(l in 1:length(param.sets)){
  
  # set parameters for this round of simulations
  temp.param.num <- param.sets[l]
  temp.params <- param.data[temp.param.num,]
  
  
  # prep for storing z-loop results
  full.sims.results <- data.frame(sim.num = numeric(),
                                  dur.time = numeric(),
                                  total.prog = numeric(),
                                  total.lr = numeric(),
                                  num.failed = numeric()
  )

  
  for(z in 1:length(nsims)){
    
    sim.num <- nsims[z]
    print(paste("Simulation number ", sim.num, sep = ""))
   
    #### set seed ####
    set.seed(1453+temp.param.num+sim.num) # unique but reproducible seeding
    
    
    #### SIMULATE PANTHER POPULATION ####
    net3 <- simulate_pop_ergm(telemetry = study.telemetry,           # telemetry data to use for simulating HR centroids
                              net.den = temp.params$net.den,         # network density to constrain simulated network to
                              a_prop = temp.params$a_prop,           # approximate proportion of population that should be categorized as adults (vs juveniles; ignore kittens) 
                              ergm.mod = best.ergm,                  # ergm model results to use
                              pop.size = temp.params$pop.size,       # size of simulated population
                              net.style = "bestFIV"                  # (don't change from "bestFIV")
                              )
    
    
    ### load libraries for transmission simulations ###
    suppressMessages(library(igraph)) # have to detach for other function, so have to continually detach and reattach
    
    
    
    # convert network object to igraph object
    g <- asIgraph(net3)
    
    # set vaccine status to zero for all (just to keep consistent when using proactive vax)
    V(g)$status <- 0
    
    
    #### SIMULATE FELV TRANSMISSION ####
    m.list <- trans_sim_basic(g = g,                             # network on which to simulate (should be an igraph object)
                    beta = temp.params$beta,                     # base transmission rate (applies to all progressive individuals)
                    c.beta = temp.params$c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                    c.rate = temp.params$c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                    duration.yrs = 5,                            # pick large number so have enough time for epidemic to run out
                    prop.outcomes = c(temp.params$prop.outcomes_p, temp.params$prop.outcomes_r, temp.params$prop.outcomes_a),   # proportion per each infected state
                    progr.dur = temp.params$progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                    regr.dur.c = temp.params$regr.dur.c,         # regressive infection duration multiplier
                    terr.repop = temp.params$terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                    use.react.vax = F,                           # True/False to use REACTIVE vaccination
                    initiate.vax = 10*52,                        # time point to initiate reactive vaccination (in weeks); arbitrarily large here
                    vax.rate = 0                                 # vaccine rate (1/avg weeks to vaccinate) - this is actually the probability for binomial trials, but functioning as a rate here; arbitrarily low
                    )
     
    
    
    #### EXTRACT OUTBREAK RESULTS ####
    m.new <- m.list[[1]]
    
    
    # number progressively infected
    pi <- sum(apply(m.new, 2, function(r) any(r %in% c(1))))
    
    # number regressively infected
    ri <- sum(apply(m.new, 2, function(r) any(r %in% c(2))))
    
    # duration of outbreak
    dur <- nrow(m.new)-1 # row 1 = time 0
    
    # number of "failed" outbreaks (fewer than 5 progressively or regressively infected)
    fails <- m.list[[3]]
    
    total.fails <-length(fails) 

    ## assemble results
    temp.sim.results <- data.frame(sim.num = sim.num,
                              dur.time = dur,
                              total.prog = pi,
                              total.lr = ri,
                              num.failed = total.fails
                              )
    
    
    #### PLOT FOR ERROR CHECKING ####
    p.results <- props_affected_births_included(m.new)
    
    eplot.name <- paste("Simulation_figures/", scenario.type, "/epiplot_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".jpg", sep = "")
    jpeg(eplot.name)
    
    plot(p.results$prop.i, col="red", type='l', main="Simulation of FeLV-like spread on a contact network",sub="",xlab="Time (weeks)",ylab="Population numbers",ylim=c(0,1),xlim=c(0, nrow(p.results)))   
    lines(p.results$prop.s, col="green")   #green is sum of #0  SUSCEPTIBLE
    lines(p.results$prop.lr, col="blue")     #blue is sum of #2  INFECTIOUS LATENT/REGRESSIVE
    lines(p.results$prop.r, col = "purple") #purple is sum of #3 ABORTIVE
    lines(p.results$prop.v, col = "yellow") # yellow is #6 VACCINATED
 
    dev.off()
    
    
    #### SAVE SIM RESULTS ####
    # store sim results
    full.sims.results[z,] <- temp.sim.results
    
    # save network
    g.name <- paste("Simulation_results/", scenario.type, "/gdata_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".Rdata", sep = "")
    save(g, file = g.name)
    
    # save outbreak results
    mlist.name <- paste("Simulation_results/", scenario.type, "/mlistdata_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".Rdata", sep = "")
    save(m.list, file = mlist.name)
    
    
  }
  
  #### SAVE PARAMETER SET RESULTS ####
  full.sims.results$param.set <- paste(scenario.type, "_paramset", temp.param.num, sep = "")
  full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
  save(full.sims.results, file = full.name)
  
}

