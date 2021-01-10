## Reactive_upc_sims.R
# 
#========================================================	
# ---
### title: Simuations of FeLV transmission under reactive underpass closure scenarios
# author: Marie Gilbertson
# date: "10/09/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates FeLV transmission on panther networks with the inclusion of reactive wildlife highway underpass closures.

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
source("Scripts/up_edge_id.R")
source("Scripts/trans_sim_upc.R")
source("Scripts/post_process_outbreak_data_basic.R")
source("Scripts/props affected_births included_basic.R")



#### load data ####
# load parameter sets
param.data <- get(load("Parameter_sets/reactiveupc_baseline_params.Rdata"))

# load telemetry
study.telemetry <- get(load("Attribute_data/2010_2012_telemetry_deid.Rdata"))

# load ergm (main FIV ERGM from FIV/FeLV project)
best.ergm <- get(load("Attribute_data/bestERGM.Rdata"))



#### SET SIMULATION BOUNDS ####
# number of sims
nsims <- seq(1, 100)

# which parameter sets
param.sets <- seq(1, 32)

# define scenario type for file naming
scenario.type <- "reactupc_baseline"



#### START OF MAIN LOOP ####
for(l in 1:length(param.sets)){
  
  # set parameters for this round of simulations
  temp.param.num <- param.sets[l]
  temp.params <- param.data[temp.param.num,]
  print(paste("Parameter set", temp.param.num, sep = " "))
  
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
    set.seed(9762+temp.param.num+sim.num) # unique, reproducible seeding
    
    
    #### SIMULATE PANTHER POPULATION ####
    net3 <- simulate_pop_ergm(telemetry = study.telemetry,           # telemetry data to use for simulating HR centroids
                              net.den = temp.params$net.den,         # network density to constrain simulated network to
                              a_prop = temp.params$a_prop,           # approximate proportion of population that should be categorized as adults (vs juveniles; ignore kittens) 
                              ergm.mod = best.ergm,                  # ergm model results to use
                              pop.size = temp.params$pop.size,       # size of simulated population
                              net.style = "bestFIV"                  # (don't change from "bestFIV")
    )
    
    
    ### load libraries for transmission simulations ###
    suppressMessages(library(igraph)) # have to detach for other function, so have to continually detach and reattach :(
    
    
    
    # convert network object to igraph object
    g <- asIgraph(net3)
    
    
    #### PROACTIVELY VACCINATE ####
    if(temp.params$prop.pro.vax>0){
      new.g <- proactive_vax(g = g,
                             prop.vax = temp.params$prop.pro.vax,
                             single.vax.eff = temp.params$unboost.ve,
                             boosted.vax.eff = temp.params$boost.ve,
                             single.prob = temp.params$vax.ratio
      )
    }else{
      # if no proactive vaccination, assign status and vax efficacy of 0 to all
      new.g <- g
      # add status = 0; veff = 0 to all
      new.g <- set_vertex_attr(new.g, "status", V(new.g), value = 0)
      new.g <- set_vertex_attr(new.g, "veff", V(new.g), value = 0)
    }
    
    
    
    #### DETERMINE EDGES ACROSS I-75 FREEWAY ####
    new.g2 <- up_edge_id(net = new.g)
    
    
    
    #### SIMULATE FELV TRANSMISSION ####
    m.list <- trans_sim_upc(g = new.g2,                                   # network on which to simulate (should be an igraph object)
                             beta = temp.params$beta,                     # base transmission rate (applies to all progressive individuals)
                             c.beta = temp.params$c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals
                             c.rate = temp.params$c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                             duration.yrs = 5,                            # pick large number so have enough time for epidemics to fade out
                             prop.outcomes = c(temp.params$prop.outcomes_p, temp.params$prop.outcomes_r, temp.params$prop.outcomes_a),   # proportion per each infected state
                             progr.dur = temp.params$progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                             regr.dur.c = temp.params$regr.dur.c,         # regressive infection duration multiplier
                             terr.repop = temp.params$terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                             initiate.upc = temp.params$onset,            # time point to initiate reactive underpass closures (in weeks); for simplicity, happens instantaneously
                             upc.dur = temp.params$upc.dur                # duration of underpass closures (re-opening also happens instantaneously)
    )
    
    
    #### PROCESS OUTBREAK RESULTS ####
    m.new <- m.list[[1]]
    
    
    # number progressively infected
    pi <- sum(apply(m.new, 2, function(r) any(r %in% c(1))))
    
    # number regressively infected
    ri <- sum(apply(m.new, 2, function(r) any(r %in% c(2))))
    
    # duration of outbreak
    dur <- nrow(m.new)-1 # row 1 = time 0
    
    # number of "failed" outbreaks (fewer than 5 progressively or regressively infected)
    fails <- m.list$fails
    
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
    lines(p.results$prop.lr, col="blue")     #blue is sum of #2  LATENT/REGRESSIVE
    lines(p.results$prop.r, col = "purple") #purple is sum of #3 IMMUNE
    # lines(p.results$prop.rlr, col = "lightblue") #lightblue is sum of #5 IMMUNE LATENTS/REGRESSIVES
    lines(p.results$prop.v, col = "yellow") # yellow is #6 VACCINATED
    # legend("topright", inset=.05,c("Susceptible","Progressive","Latent/Regressive", "Immune"),col=c("green","red","blue", "purple"),lty=c(1,1,1,1))
    
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

# load parameter sets
param.data <- get(load("Parameter_sets/reactiveupc_baseline_params.Rdata"))
param.data$param.set <- paste("reactupc_baseline_paramset", seq(1, 32), sep = "")

# load all "full simulation results" files and combine for plotting
all.results <- NULL

for(i in 1:length(param.sets)){
  
  temp.param.num <- param.sets[i]
  
  full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
  full.sims.results <- get(load(full.name))
  
  all.results <- rbind(all.results, full.sims.results)
}

all.results <- left_join(all.results, param.data, by = "param.set")


#### pull results from other "control" scenarios ####
# get summary values from no-intervention scenarios
no.inter.results <- analyze_infections(scenario.type = "nointer_baseline", param.sets = seq(1,1))

# get summary values from proactive vaccination scenarios
provax.results.summ <- analyze_infections(scenario.type =  "provax_baseline", param.sets = seq(1,24))

# load and join proactive vax parameter sets with corresponding results
provax.param.data <- get(load("Parameter_sets/proactivevax_baseline_params_101320.Rdata"))
provax.param.data$param.set <- seq(1, nrow(provax.param.data))
provax.results.summ <- left_join(provax.results.summ, provax.param.data, by = "param.set")

# keep only those results that match the proactive vax levels used within the reactive scenarios
provax.summ.sub <- provax.results.summ[provax.results.summ$prop.pro.vax %in% c(0.2, 0.4, 0.6) & provax.results.summ$vax.ratio==0.5,]


#### extract results from test and removal scenarios ####

upc.results.summ <- analyze_infections(scenario.type = "reactupc_baseline", param.sets = seq(1,32))
param.data$param.set2 <- seq(1, 32)
# join with corresponding parameter data
upc.results.summ <- left_join(upc.results.summ, param.data, by = c("param.set" = "param.set2"))



#### progressives/mortalities plot ####

# create factor level for proportion vax so can control ordering in facets
all.results$prop.pro.vax_f <- factor(all.results$prop.pro.vax, levels = seq(0.6, 0, -0.2))

prop.labs <- c("0%", "20%", "40%", "60%")
names(prop.labs) <- seq(0, 0.6, by = 0.2)

dur.labs <- c("4 weeks", "13 weeks", "26 weeks", "52 weeks")
names(dur.labs) <- unique(all.results$upc.dur)

# split into two datasets for two different onsets (at 26 vs 52 weeks); for visibility in plotting
all26 <- all.results[all.results$onset==26,]
all52 <- all.results[all.results$onset==52,]


# plot for onset = 26 weeks
pa <- ggplot(all26, aes(x = total.prog)) + geom_histogram(binwidth = 3)+
  facet_grid(prop.pro.vax_f ~ upc.dur, labeller = labeller(prop.pro.vax_f = prop.labs, upc.dur = dur.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")

pa <- pa + 
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +

  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +

  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") 


pa


# plot for onset = 52 weeks
pb <- ggplot(all52, aes(x = total.prog)) + geom_histogram(binwidth = 3)+
  facet_grid(prop.pro.vax_f ~ upc.dur, labeller = labeller(prop.pro.vax_f = prop.labs, upc.dur = dur.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")

pb <- pb + 
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.pi[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") 


pb

# arrange into one plot
library(ggpubr)
pc <- ggarrange(pa, pb, labels = c("A", "B"), nrow = 2)
pc



#### duration of epidemics ####

# get summary values from no-intervention scenarios
no.inter.full <- get(load("Simulation_results/nointer_baseline/fulldata_nointer_baseline_paramset1.Rdata"))


# check how many lasted the "full" 5 years
length(which(all.results$dur.time==260))

# again, plot just for onset = 26 weeks first
p3a <- ggplot(all26, aes(x = dur.time)) + geom_histogram(binwidth = 10)+
  facet_grid(prop.pro.vax_f ~ upc.dur, labeller = labeller(prop.pro.vax_f = prop.labs, upc.dur = dur.labs)) +
  geom_vline(xintercept = no.inter.results$med.dur, colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3a <- p3a + 
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.6" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.4" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0.2" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="4" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="13" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="26" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="26"]), colour = "#ff7f00") +
  geom_vline(data=filter(all26, prop.pro.vax_f=="0" & upc.dur=="52" & onset=="26"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="26"]), colour = "#ff7f00") 


p3a


# then plot for onset = 52 weeks
p3b <- ggplot(all52, aes(x = dur.time)) + geom_histogram(binwidth = 10)+
  facet_grid(prop.pro.vax_f ~ upc.dur, labeller = labeller(prop.pro.vax_f = prop.labs, upc.dur = dur.labs)) +
  geom_vline(xintercept = no.inter.results$med.dur, colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3b <- p3b + 
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.6" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.6 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.4" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.4 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0.2" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0.2 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="4" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="4" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="13" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="13" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="26" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="26" & upc.results.summ$onset=="52"]), colour = "#ff7f00") +
  geom_vline(data=filter(all52, prop.pro.vax_f=="0" & upc.dur=="52" & onset=="52"), aes(xintercept = upc.results.summ$med.dur[upc.results.summ$prop.pro.vax==0 & upc.results.summ$upc.dur=="52" & upc.results.summ$onset=="52"]), colour = "#ff7f00") 


p3b

# now can arrange both plots together
library(ggpubr)
p3c <- ggarrange(p3a, p3b, labels = c("A", "B"), nrow = 2)
p3c



#### EPI SUMMARY PLOT ####
# decide which parameter sets to plot
param.data2 <- subset(param.data, param.data$onset==26 & param.data$upc.dur==26 & param.data$prop.pro.vax==0.6)


# extract parameter sets to plot
sets.tp <- param.data2$param.set2

### set up loops and empty results data frame ###
full.results <- NULL

nsims <- 100

### loop through parameter sets
for(i in 1:length(sets.tp)){
  print(i)
  temp.param.num <- sets.tp[i]
  ### within each parameter set, loop through simulations
  for(j in 1:nsims){
    
    ### for each simulation, calculate proportions per state at each time point
    sim.num <- j
    
    mlist.name <- paste("Simulation_results/", scenario.type, "/mlistdata_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".Rdata", sep = "")
    m.list <- get(load(mlist.name))
    
    m.new <- m.list$m.new
    
    p.results <- props_affected_births_included(m.new)
    
    # if outbreak didn't last the max time, append the last observation out to the "max time" so means 
    # across ALL simulations are not biased by just a few long simulations in later weeks of outbreaks
    if(nrow(p.results)<261){
      append <- tail(p.results, 1)
      no.app <- 261-nrow(p.results)
      append.no <- append[rep(seq_len(nrow(append)), no.app), ]
      append.no$time <- seq(append$time+1, 260)
      
      p.results <- rbind(p.results, append.no)
    }
    
    p.results$param.set <- temp.param.num
    p.results$sim.num <- sim.num
    
    full.results <- rbind(full.results, p.results)
    
  } # end of j loop
  
  
} # end of i loop


full.results$plot.id <- paste(full.results$param.set, full.results$sim.num, sep = "_")


# mean by time point and parameter set
means.per.time <- ddply(full.results[,c("time", "param.set", "prop.i")], .(time, param.set), summarize, mean = mean(prop.i))

# plot individual simulations and mean per time
g <- ggplot() + 
  geom_line(data = full.results, aes(x=time, y=prop.i, group=plot.id), alpha = 0.2, colour = "#ca0020") 

g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

g <- g + geom_line(data = means.per.time, aes(x = time, y = mean), size = 2, colour = "#ca0020")

g <- g +  xlab("Time (weeks)") + ylab("Proportion progressively infected") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")

g <- g + geom_vline(xintercept = param.data2$onset, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = (param.data2$onset + param.data2$upc.dur), colour = "black", linetype = "dashed")
g


### repeat for second set of parameters ###
# decide which parameter sets to plot
param.data2 <- subset(param.data, param.data$onset==52 & param.data$upc.dur==26 & param.data$prop.pro.vax==0.6)


# extract parameter sets to plot
sets.tp <- param.data2$param.set2

### set up loops and empty results data frame ###
full.results <- NULL

nsims <- 100

### loop through parameter sets
for(i in 1:length(sets.tp)){
  print(i)
  temp.param.num <- sets.tp[i]
  ### within each parameter set, loop through simulations
  for(j in 1:nsims){
    
    ### for each simulation, calculate proportions per state at each time point
    sim.num <- j
    
    mlist.name <- paste("Simulation_results/", scenario.type, "/mlistdata_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".Rdata", sep = "")
    m.list <- get(load(mlist.name))
    
    m.new <- m.list$m.new
    
    p.results <- props_affected_births_included(m.new)
    
    # if outbreak didn't last the max time, append the last observation out to the "max time" so means 
    # across ALL simulations are not biased by just a few long simulations in later weeks of outbreaks
    if(nrow(p.results)<261){
      append <- tail(p.results, 1)
      no.app <- 261-nrow(p.results)
      append.no <- append[rep(seq_len(nrow(append)), no.app), ]
      append.no$time <- seq(append$time+1, 260)
      
      p.results <- rbind(p.results, append.no)
    }
    
    p.results$param.set <- temp.param.num
    p.results$sim.num <- sim.num
    
    full.results <- rbind(full.results, p.results)
    
  } # end of j loop
  
  
} # end of i loop


full.results$plot.id <- paste(full.results$param.set, full.results$sim.num, sep = "_")


# mean by time point and parameter set
means.per.time <- ddply(full.results[,c("time", "param.set", "prop.i")], .(time, param.set), summarize, mean = mean(prop.i))

# plot individual simulations and mean per time
g2 <- ggplot() + 
  geom_line(data = full.results, aes(x=time, y=prop.i, group=plot.id), alpha = 0.2, colour = "#0571b0") 

g2 <- g2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

g2 <- g2 + geom_line(data = means.per.time, aes(x = time, y = mean), size = 2, colour = "#0571b0")

g2 <- g2 +  xlab("Time (weeks)") + ylab("Proportion progressively infected") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")

g2 <- g2 + geom_vline(xintercept = param.data2$onset, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = (param.data2$onset + param.data2$upc.dur), colour = "black", linetype = "dashed")

g2

# arrange both plot together
library(ggpubr)
g3 <- ggarrange(g, g2, labels = c("A", "B"), nrow = 2)
g3
