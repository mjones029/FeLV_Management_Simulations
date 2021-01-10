## Reactive_vax_sims.R
# 
#========================================================	
# ---
### title: Simuations of FeLV transmission under reactive vaccination scenarios
# author: Marie Gilbertson
# date: "10/09/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Simulates FeLV transmission on panther networks with the inclusion of reactive vaccination.

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
source("Scripts/trans_sim_reactvax.R")
source("Scripts/post_process_outbreak_data_basic.R")
source("Scripts/props affected_births included_basic.R")



#### load data ####
# load parameter sets
param.data <- get(load("Parameter_sets/reactivevax_baseline_params.Rdata"))

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
scenario.type <- "reactvax_baseline"



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
                                  total.vaxed = numeric(),
                                  react.boosts = numeric(),
                                  total.vax.used = numeric()
  )

  
  for(z in 1:length(nsims)){
    
    sim.num <- nsims[z]
    print(paste("Simulation number ", sim.num, sep = ""))
    
    #### set seed ####
    set.seed(2657+temp.param.num+sim.num) # unique, reproducible seeding
    
    
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
    m.list <- trans_sim_reactvax(g = new.g,                                  # network on which to simulate (should be an igraph object)
                                beta = temp.params$beta,                     # base transmission rate (applies to all progressive individuals)
                                c.beta = temp.params$c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                                c.rate = temp.params$c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                                duration.yrs = 5,                            # pick large number so have enough time for epidemics to die out
                                prop.outcomes = c(temp.params$prop.outcomes_p, temp.params$prop.outcomes_r, temp.params$prop.outcomes_a),   # proportion per each infected state
                                progr.dur = temp.params$progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                                regr.dur.c = temp.params$regr.dur.c,         # regressive infection duration multiplier
                                terr.repop = temp.params$terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                                single.vax.eff = temp.params$unboost.ve,     # unboosted vax efficacy
                                boosted.vax.eff = temp.params$boost.ve,      # boosted vax efficacy
                                initiate.vax = temp.params$onset,            # time point to initiate reactive vaccination (in weeks)
                                window.dur = temp.params$window.dur,         # use a vaccination window? If so, currently set for 26 week window; otherwise, reactive vax year-round after onset
                                dist.type = temp.params$dist.type            # reactive vax distribution type (i.e. random vs spatially-targeted)
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
    
    # number of individuals ever **successfully** vaccinated (number proactively vax + single reactive vax)
    ever.vax <- length(V(new.g)[V(new.g)$status==6]) + length(m.list$total.rvax[m.list$total.rvax=="i"])
    
    # number **successfully** reactively boosted
    react.boosts <- length(m.list$total.rvax[m.list$total.rvax=="b"])
    
    # total vaccines used (for proactive, veff=0.8 --> 2 vaccines)
    total.reacts <- length(m.list$total.rvax[m.list$total.rvax!="0"])
    total.proacts1 <- length(V(new.g)$veff[V(new.g)$veff==temp.params$unboost.ve])
    total.proacts2 <- (length(V(new.g)$veff[V(new.g)$veff==temp.params$boost.ve]))*2
    
    total.vax.used <- total.reacts + total.proacts1 + total.proacts2
    
    # number of attempts at vaccination
    vax.attempts <- length(m.list$attempted.vax)
    
    
    
    ## assemble results
    temp.sim.results <- data.frame(sim.num = sim.num,
                                   dur.time = dur,
                                   total.prog = pi,
                                   total.lr = ri,
                                   num.failed = total.fails,
                                   total.vaxed = ever.vax,
                                   react.boosts = react.boosts,
                                   total.vax.used = total.vax.used
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
source("Scripts/analyze_reactive_vaccinations.R")

# load parameter data
param.data <- get(load("Parameter_sets/reactivevax_baseline_params.Rdata"))
param.data$param.set <- paste("reactvax_baseline_paramset", seq(1, nrow(param.data)), sep = "")



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

# load and join proactive vax parameter sets with their corresponding set of results
provax.param.data <- get(load("Parameter_sets/proactivevax_baseline_params.Rdata"))
provax.param.data$param.set <- seq(1, nrow(provax.param.data))
provax.results.summ <- left_join(provax.results.summ, provax.param.data, by = "param.set")

# keep only those results that match the proactive vax levels used within the reactive scenarios
provax.summ.sub <- provax.results.summ[provax.results.summ$prop.pro.vax %in% c(0.2, 0.4, 0.6) & provax.results.summ$vax.ratio==0.5,]


#### summary results from reactive scenarios ####
# get summary values from proactive vaccination scenarios
reactvax.results.summ <- analyze_infections(scenario.type =  "reactvax_baseline", param.sets = seq(1,32))

# join reactive vax parameter sets with corresponding results
reactvax.param.data <- param.data
reactvax.param.data$param.set <- seq(1, nrow(reactvax.param.data))
reactvax.results.summ <- left_join(reactvax.results.summ, reactvax.param.data, by = "param.set")
# relevel for plotting
reactvax.results.summ$prop.pro.vax_f <- factor(reactvax.results.summ$prop.pro.vax, levels = seq(0.6, 0, -0.2))


#### vaccination summaries from reactive scenarios ####
reactvax.vax.results <- analyze_reactive_vaccinations(scenario.type = "reactvax_baseline",
                                                      param.sets = seq(1,32))

# join with corresponding parameter data
reactvax.vax.results <- left_join(reactvax.vax.results, reactvax.param.data, by = "param.set")

# relevel for plotting
reactvax.vax.results$prop.pro.vax_f <- factor(reactvax.vax.results$prop.pro.vax, levels = seq(0.6, 0, -0.2))

# split into subsets for better viewing
vax.results.allyear <- reactvax.vax.results[reactvax.vax.results$window.dur==F,]
vax.results.window <- reactvax.vax.results[reactvax.vax.results$window.dur==T,]



#### prep data for plotting ####

# create factor level for proportion vax so can control ordering in facets
all.results$prop.pro.vax_f <- factor(all.results$prop.pro.vax, levels = seq(0.6, 0, -0.2))

# update facet labelling
prop.labs <- c("0%", "20%", "40%", "60%")
names(prop.labs) <- seq(0, 0.6, by = 0.2)

dist.labs <- c("Random", "Spatial")
names(dist.labs) <- unique(all.results$dist.type)

on.labs <- c("Onset at 26 weeks", "Onset at 52 weeks")
names(on.labs) <- unique(all.results$onset)


#### subset data into chunks of primary interest ####

# no vax duration limit
window <- subset(all.results, all.results$window.dur==T)
allyear <- subset(all.results, all.results$window.dur==F)

### start by plotting results when reactive vaccination occurs year-round
reactvax.results.summ <- subset(reactvax.results.summ, reactvax.results.summ$window.dur==F)



#### progressives/mortalities plot ####
p <- ggplot(allyear, aes(x = total.prog)) + geom_histogram(binwidth = 3)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")

p <- p + 
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") 
  
## (optional): add test regarding median total individuals vaccianted and median vaccines per individuals
# p <- p + geom_text(data = vax.results.allyear,
#                    mapping = aes(x = Inf, y = Inf, label = med.tv),
#                    vjust = 1.5,
#                    hjust = 1.15
#                   ) +
#           geom_text(data = vax.results.allyear,
#                   mapping = aes(x = Inf, y = Inf, label = signif(med.vpi, 3)),
#                   vjust = 3.5,
#                   hjust = 1.15,
#                   colour = "purple"
#                   )

p



### duration of epidemics ####


p3 <- ggplot(allyear, aes(x = dur.time)) + geom_histogram(binwidth = 10)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.dur, colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3 <- p3 + 
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(allyear, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") 


p3



#### REPEAT PLOTS FOR 6 MONTH DURATION REACTIVE VAX ####
# get summary values from proactive vaccination scenarios
reactvax.results.summ <- analyze_infections(scenario.type =  "reactvax_baseline", param.sets = seq(1,32))

# load and join parameter sets
reactvax.param.data <- param.data
reactvax.param.data$param.set <- seq(1, nrow(reactvax.param.data))
reactvax.results.summ <- left_join(reactvax.results.summ, reactvax.param.data, by = "param.set")
reactvax.results.summ$prop.pro.vax_f <- factor(reactvax.results.summ$prop.pro.vax, levels = seq(0.6, 0, -0.2))

reactvax.results.summ <- subset(reactvax.results.summ, reactvax.results.summ$window.dur==T)


#### progressives/mortalities plot ####
p <- ggplot(window, aes(x = total.prog)) + geom_histogram(binwidth = 3)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.pi, colour = "#e41a1c") +
  xlab("Total FeLV mortalities") +
  ylab("Count of occurences")

p <- p + 
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.pi[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.pi[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") 

## (optional): add test regarding median total individuals vaccianted and median vaccines per individuals
# p <- p + geom_text(data = vax.results.window,
#                    mapping = aes(x = Inf, y = Inf, label = med.tv),
#                    vjust = 1.5,
#                    hjust = 1.15
#                    ) +
#           geom_text(data = vax.results.window,
#                     mapping = aes(x = Inf, y = Inf, label = signif(med.vpi, 3)),
#                     vjust = 3.5,
#                     hjust = 1.15,
#                     colour = "purple"
#                     )

p




### duration of epidemics ####

p3 <- ggplot(window, aes(x = dur.time)) + geom_histogram(binwidth = 10)+
  facet_grid(prop.pro.vax_f ~ dist.type + onset, labeller = labeller(prop.pro.vax_f = prop.labs, dist.type = dist.labs, onset = on.labs)) +
  geom_vline(xintercept = no.inter.results$med.dur, colour = "#e41a1c") +
  xlab("Duration of FeLV epidemics (weeks)") +
  ylab("Count of occurences")

p3 <- p3 + 
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.6]), colour = "#377eb8") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.4]), colour = "#377eb8") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2"), aes(xintercept = provax.summ.sub$med.dur[provax.summ.sub$prop.pro.vax==0.2]), colour = "#377eb8") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.6" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.6 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.4" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.4 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0.2" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0.2 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="random" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="random" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="random" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="26"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="26"]), colour = "#4daf4a") +
  geom_vline(data=filter(window, prop.pro.vax_f=="0" & dist.type=="spatial" & onset=="52"), aes(xintercept = reactvax.results.summ$med.dur[reactvax.results.summ$prop.pro.vax==0 & reactvax.results.summ$dist.type=="spatial" & reactvax.results.summ$onset=="52"]), colour = "#4daf4a") 


p3


