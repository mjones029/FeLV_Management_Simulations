## Proactive_vax_sims.R
# 
#========================================================	
# ---
### title: Simuations of FeLV transmission under proactive vaccination scenarios
# author: Marie Gilbertson
# date: "10/09/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates FeLV transmission on panther networks with the inclusion of proactive vaccination.

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(stringr) # for proactive_vax function
library(adehabitatHR) # kerneloverlap, mcp
library(sp) # coordinates() function
library(rgeos) # centroids
library(geosphere) # geographic distances between points (distm() function)
library(plyr)
library(dplyr) # sample_n function
library(spatstat)
library(DCG) # as.symmetricAdjacencyMatrix function
library(ergm)
library(intergraph)



#### load external functions ####
source("Scripts/simulate_pop_ergm.R")
source("Scripts/proactive_vax.R")
source("Scripts/trans_sim_basic.R")
source("Scripts/post_process_outbreak_data_basic.R")
source("Scripts/props affected_births included_basic.R")



#### load data ####
# load parameter sets
param.data <- get(load("Parameter_sets/proactivevax_baseline_params.Rdata"))

# load telemetry
study.telemetry <- get(load("Attribute_data/2010_2012_telemetry_deid.Rdata"))

# load ergm (main FIV ERGM from FIV/FeLV project)
best.ergm <- get(load("Attribute_data/bestERGM.Rdata"))



#### SET SIMULATION BOUNDS ####
# number of sims
nsims <- seq(1, 100)


# which parameter sets
param.sets <- seq(1, nrow(param.data))


# define scenario type for file naming
scenario.type <- "provax_baseline"


# set maximum duration of simulations in years (5 years for baseline; 8 years for sensitivity)
max.duration <- 5


#### START OF MAIN LOOP ####
for(l in 1:length(param.sets)){
  
  # set parameters for this round of simulations
  temp.param.num <- param.sets[l]
  temp.params <- param.data[temp.param.num,]
  
  
  # prep for storing z-loop results
  full.sims.results <- data.frame(sim.num = numeric(),
                                  dur.time = numeric(),
                                  total.prog = numeric(),
                                  total.lr = numeric(),
                                  num.failed = numeric(),
                                  total.vaxed = numeric()
  )
  
  
  for(z in 1:length(nsims)){
    
    sim.num <- nsims[z]
    print(paste("Simulation number ", sim.num, sep = ""))
    
    #### set seed ####
    set.seed(2845+temp.param.num+sim.num) # unique, reproducible seeding
    
    
    #### SIMULATE PANTHER POPULATION ####
    net3 <- simulate_pop_ergm(telemetry = study.telemetry,           # telemetry data to use for simulating HR centroids
                              net.den = temp.params$net.den,         # network density to constrain simulated network to
                              a_prop = temp.params$a_prop,           # approximate proportion of population that should be categorized as adults (vs juveniles; ignore kittens) 
                              ergm.mod = best.ergm,                  # ergm model results to use
                              pop.size = temp.params$pop.size,       # size of simulated population
                              net.style = "bestFIV"                    # (don't change from "bestFIV")
    )
    
    
    ### load libraries for transmission simulations ###
    suppressMessages(library(igraph)) # have to detach for other function, so have to continually detach and reattach 
    
    
    
    # convert network object to igraph object
    g <- asIgraph(net3)
    
    
    #### PROACTIVELY VACCINATE ####
    new.g <- proactive_vax(g = g,
                          prop.vax = temp.params$prop.pro.vax,
                          single.vax.eff = temp.params$unboost.ve,
                          boosted.vax.eff = temp.params$boost.ve,
                          single.prob = temp.params$vax.ratio
                          )
    
    
    #### SIMULATE FELV TRANSMISSION ####
    m.list <- trans_sim_basic(g = new.g,                                   # network on which to simulate (should be an igraph object)
                              beta = temp.params$beta,                     # base transmission rate (applies to all progressive individuals)
                              c.beta = temp.params$c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                              c.rate = temp.params$c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                              duration.yrs = max.duration,                 # pick large number so have enough time for epidemic to die out
                              prop.outcomes = c(temp.params$prop.outcomes_p, temp.params$prop.outcomes_r, temp.params$prop.outcomes_a),   # proportion per each infected state
                              progr.dur = temp.params$progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                              regr.dur.c = temp.params$regr.dur.c,         # regressive infection duration multiplier 
                              terr.repop = temp.params$terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
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
    
    # number of individuals ever vaccinated (number proactively vax + reactively vax)
    ever.vax <- length(V(new.g)[V(new.g)$status==6]) + length(m.list[[4]])
    
    ## assemble results
    temp.sim.results <- data.frame(sim.num = sim.num,
                                   dur.time = dur,
                                   total.prog = pi,
                                   total.lr = ri,
                                   num.failed = total.fails,
                                   total.vaxed = ever.vax
    )
    
    
    #### PLOT FOR ERROR CHECKING ####
    p.results <- props_affected_births_included(m.new)
    
    eplot.name <- paste("Simulation_figures/", scenario.type, "/epiplot_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".jpg", sep = "")
    jpeg(eplot.name)
    
    plot(p.results$prop.i, col="red", type='l', main="Simulation of FeLV-like spread on a contact network",sub="",xlab="Time (weeks)",ylab="Population numbers",ylim=c(0,1),xlim=c(0, nrow(p.results)))   
    lines(p.results$prop.s, col="green")   #green is sum of #0  SUSCEPTIBLE
    lines(p.results$prop.lr, col="blue")     #blue is sum of #2  LATENT/REGRESSIVE
    lines(p.results$prop.r, col = "purple") #purple is sum of #3 IMMUNE
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

 



#### ANALYZE SIMULATION RESULTS ####
library(ggplot2)
# load function for extracting summary data for infections
source("Scripts/analyze_infections.R")

# load proactive vax parameter sets (if not already loaded)
param.data <- get(load("Parameter_sets/proactivevax_baseline_params.Rdata"))
param.data$param.set <- paste("provax_baseline_paramset", seq(1, 24), sep = "")

# load all "full simulation results" files and combine for plotting
all.results <- NULL

for(i in 1:length(param.sets)){
  
  temp.param.num <- param.sets[i]

  full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
  full.sims.results <- get(load(full.name))
  
  all.results <- rbind(all.results, full.sims.results)
}

# join simulation results with corresponding parameter data
all.results <- left_join(all.results, param.data, by = "param.set")



#### progressives/mortalities plot ####
# First, get summary values from no-intervention scenarios to compare to proactive vax results
no.inter.results <- analyze_infections(scenario.type = "nointer_baseline", param.sets = seq(1,1))

# Then calculate summary values for proactive vax scenarios
provax.results <- analyze_infections(scenario.type = "provax_baseline", param.sets = seq(1,24))
param.data$param.set2 <- seq(1, 24)
provax.results <- left_join(provax.results, param.data, by = c("param.set" = "param.set2"))


# create factor level for proportion vax so can control ordering in facets
all.results$prop.pro.vax_f <- factor(all.results$prop.pro.vax, levels = seq(0.8, 0.1, -0.1))
all.results$vax.ratio_f <- factor(all.results$vax.ratio, levels = c(1, 0.5, 0))

# edit facet labels for clarity
vax.ratio.labs <- c("0% boosted", "100% boosted", "50% boosted")
names(vax.ratio.labs) <- c(1.0, 0.0, 0.5)

prop.labs <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%")
names(prop.labs) <- seq(0.1, 0.8, by = 0.1)


p <- ggplot(all.results, aes(x = total.prog)) + geom_histogram(binwidth = 2)+
  facet_grid(prop.pro.vax_f ~ vax.ratio_f, labeller = labeller(prop.pro.vax_f = prop.labs, vax.ratio_f = vax.ratio.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")
  
p <- p + 
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==0), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==1), aes(xintercept = provax.results$med.pi[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==1]), colour = "#377eb8") 
  
  
p




#### duration of epidemics ####

# First, get summary values from no-intervention scenarios to compare to proactive vax results
no.inter.full <- get(load("Simulation_results/nointer_baseline/fulldata_nointer_baseline_paramset1.Rdata"))


p3 <- ggplot(all.results, aes(x = dur.time)) + geom_histogram(binwidth = 4.5)+
  facet_grid(prop.pro.vax_f ~ vax.ratio_f, labeller = labeller(prop.pro.vax_f = prop.labs, vax.ratio_f = vax.ratio.labs)) +
  geom_vline(xintercept = median(no.inter.full$dur.time), colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3 <- p3 + 
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==0), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==0]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==0.5), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==0.5]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.8 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.8 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.7" & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.7" & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.6 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.6 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.5 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.5 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.4 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.4 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.3" & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax=="0.3" & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.2 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.2 & provax.results$vax.ratio==1]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f==0.1 & vax.ratio_f==1), aes(xintercept = provax.results$med.dur[provax.results$prop.pro.vax==0.1 & provax.results$vax.ratio==1]), colour = "#377eb8") 


p3


