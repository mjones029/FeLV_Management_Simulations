## Reactive_TandR_sims.R
# 
#========================================================	
# ---
### title: Simuations of FeLV transmission under reactive test and removal scenarios
# author: Marie Gilbertson
# date: "10/23/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates FeLV transmission on panther networks with the inclusion of reactive test-and-removal.


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
source("Scripts/trans_sim_tandr.R")
source("Scripts/post_process_outbreak_data_tandr.R")
source("Scripts/props affected_births included_tandr.R")



#### load data ####
# load parameter sets
param.data <- get(load("Parameter_sets/reactivetandr_baseline_params.Rdata"))

# load telemetry
study.telemetry <- get(load("Attribute_data/2010_2012_telemetry_deid.Rdata"))

# load ergm (main FIV ERGM from FIV/FeLV project)
best.ergm <- get(load("Attribute_data/bestERGM.Rdata"))



#### SET SIMULATION BOUNDS ####
# number of sims
nsims <- seq(1, 100)

# which parameter sets
param.sets <- seq(1, 16)

# define scenario type for file naming
scenario.type <- "reacttandr_baseline"



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
                                  num.failed = numeric(),
                                  pro.vax = numeric(),
                                  total.caps = numeric(),
                                  total.rems = numeric(),
                                  total.euth = numeric()
  )
  

  for(z in 1:length(nsims)){
    
    sim.num <- nsims[z]
    print(paste("Simulation number ", sim.num, sep = ""))
    
    #### set seed ####
    set.seed(8924+temp.param.num+sim.num) # unique, reproducible seeding
    
    
    #### SIMULATE PANTHER POPULATION ####
    net3 <- simulate_pop_ergm(telemetry = study.telemetry,           # telemetry data to use for simulating HR centroids
                              net.den = temp.params$net.den,         # network density to constrain simulated network to
                              a_prop = temp.params$a_prop,           # approximate proportion of population that should be categorized as adults (vs juveniles; ignore kittens) 
                              ergm.mod = best.ergm,                  # ergm model results to use
                              pop.size = temp.params$pop.size,       # size of simulated population
                              net.style = "bestFIV"                    # (don't change from "bestFIV")
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
    
    
    #### SIMULATE FELV TRANSMISSION ####
    m.list <- trans_sim_tandr(g = new.g,                                   # network on which to simulate (should be an igraph object)
                              beta = temp.params$beta,                     # base transmission rate (applies to all progressive individuals)
                              c.beta = temp.params$c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                              c.rate = temp.params$c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                              duration.yrs = 5,                            # pick large number so have enough time for epidemics to die out
                              prop.outcomes = c(temp.params$prop.outcomes_p, temp.params$prop.outcomes_r, temp.params$prop.outcomes_a),   # proportion per each infected state
                              progr.dur = temp.params$progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                              regr.dur.c = temp.params$regr.dur.c,         # regressive infection duration multiplier
                              terr.repop = temp.params$terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                              initiate.int = temp.params$onset,            # time point to initiate intervention (in weeks); in this case, to initiate test and removal
                              max.occupancy = temp.params$max.occupancy,   # maximum number of temporary removals possible at one time
                              dist.type = temp.params$dist.type,           # distribution of captures; "random" or "spatial"
                              int.dur = 17                                 # duration of capture season each year (17 weeks ~ 4 months)
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
    
    # proactive vax
    total.vax <- length(V(new.g)[V(new.g)$status==6])
    
    # number of individuals ever captured
    total.tr <- m.list$total.tr
    
    total.caps <- length(total.tr)
    
    
    # number ever removed
    total.rems <- length(total.tr[total.tr=="r"])
    
    
    # number ever euthanized 
    total.euth <- length(total.tr[total.tr=="e"])
    
    
    ## assemble results
    temp.sim.results <- data.frame(sim.num = sim.num,
                                   dur.time = dur,
                                   total.prog = pi,
                                   total.lr = ri,
                                   num.failed = total.fails,
                                   pro.vax = total.vax,
                                   total.caps = total.caps,
                                   total.rems = total.rems,
                                   total.euth = total.euth
    )
    
    
    #### PLOT FOR ERROR CHECKING ####
    p.results <- props_affected_births_included_tandr(m.new)
    
    eplot.name <- paste("Simulation_figures/", scenario.type, "/epiplot_", scenario.type, "_paramset", temp.param.num, "_sim", sim.num, ".jpg", sep = "")
    jpeg(eplot.name)
    
    plot(p.results$prop.i, col="red", type='l', main="Simulation of FeLV-like spread on a contact network",sub="",xlab="Time (weeks)",ylab="Population numbers",ylim=c(0,1),xlim=c(0, nrow(p.results)))   
    lines(p.results$prop.s, col="green")   #green is sum of #0  SUSCEPTIBLE
    lines(p.results$prop.lr, col="blue")     #blue is sum of #2  LATENT/REGRESSIVE
    lines(p.results$prop.r, col = "purple") #purple is sum of #3 IMMUNE
    lines(p.results$prop.v, col = "yellow") # yellow is #6 VACCINATED
    par(new = TRUE)
    plot(p.results$num.rem, type = "l", axes = FALSE, xlab = "", ylab = "") # can add documentation of number temporarily removed over time
    axis(side = 4, at = pretty(range(p.results$num.rem)))      # Add second axis
    mtext("number removed", side = 4, line = 0)             # Add second axis label
    
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





# #### ANALYZE SIMULATION RESULTS ####
library(ggplot2)
# load function for extracting summary data for infections
source("Scripts/analyze_infections.R")
source("Scripts/analyze_reactive_tandr.R")

# load parameter sets
param.data <- get(load("Parameter_sets/reactivetandr_baseline_params.Rdata"))
param.data$param.set <- paste("reacttandr_baseline_paramset", seq(1, 16), sep = "")

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

# load and join proactive vax parameter sets with corresponding results data
provax.param.data <- get(load("Parameter_sets/proactivevax_baseline_params.Rdata"))
provax.param.data$param.set <- seq(1, nrow(provax.param.data))
provax.results.summ <- left_join(provax.results.summ, provax.param.data, by = "param.set")

# keep only those results that match the proactive vax levels used within the reactive scenarios
provax.summ.sub <- provax.results.summ[provax.results.summ$prop.pro.vax %in% c(0.2, 0.4, 0.6) & provax.results.summ$vax.ratio==0.5,]


#### extract results from test and removal scenarios ####
# get infection summary values
tandr.results.summ <- analyze_infections(scenario.type = "reacttandr_baseline", param.sets = seq(1,16))
param.data$param.set2 <- seq(1, 16)
# join with test-and-removal parameter information
tandr.results.summ <- left_join(tandr.results.summ, param.data, by = c("param.set" = "param.set2"))

# also calculate summary values specific to test-and-removal scenarios (e.g. proportion of captures that were "successful")
tandr.rem.results <- analyze_reactive_tandr(scenario.type = "reacttandr_baseline", param.sets = seq(1, 16))
# join with parameter information
tandr.rem.results <- left_join(tandr.rem.results, param.data, by = c("param.set" = "param.set2"))
# relevel for plotting
tandr.rem.results$prop.pro.vax_f <- factor(tandr.rem.results$prop.pro.vax, levels = seq(0.6, 0, -0.2))


#### progressives/mortalities plot ####

# create factor level for proportion proactively vax'ed so can control ordering in facets
all.results$prop.pro.vax_f <- factor(all.results$prop.pro.vax, levels = seq(0.6, 0, -0.2))

# update facet labels for interpretability
prop.labs <- c("0%", "20%", "40%", "60%")
names(prop.labs) <- seq(0, 0.6, by = 0.2)

dist.labs <- c("Random", "Spatial")
names(dist.labs) <- unique(all.results$dist.type)

on.labs <- c("Onset at 26 weeks", "Onset at 52 weeks")
names(on.labs) <- unique(all.results$onset)



p <- ggplot(all.results, aes(x = total.prog)) + geom_histogram(binwidth = 3)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")

p <- p + 
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.pi[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") 

# (optional): add test showing median proportion of successful captures
p <- p + 
  geom_text(data = tandr.rem.results,
            mapping = aes(x = Inf, y = Inf, label = signif(med.psc, 1)),
            vjust = 1.5,
            hjust = 1.15
  )


p




#### duration of epidemics ####

# get summary values from no-intervention scenarios
no.inter.full <- get(load("Simulation_results/nointer_baseline/fulldata_nointer_baseline_paramset1.Rdata"))


# check how many lasted the "full" 5 years
length(which(all.results$dur.time==260))



p3 <- ggplot(all.results, aes(x = dur.time)) + geom_histogram(binwidth = 10)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.dur, colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3 <- p3 + 
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.6 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.4 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0.2 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="random" & tandr.results.summ$onset=="52"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="26"]), colour = "#984ea3") +
  geom_vline(data=filter(all.results, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = tandr.results.summ$med.dur[tandr.results.summ$prop.pro.vax==0 & tandr.results.summ$dist.type=="spatial" & tandr.results.summ$onset=="52"]), colour = "#984ea3") 


p3
