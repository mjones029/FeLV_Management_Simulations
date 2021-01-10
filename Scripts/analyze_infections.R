# analyze_infections
#
#========================================================	
# ---
### title: Analyze simulation infections
# author: Marie Gilbertson
# date: "10/11/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function that outputs mean, median, and range results for simulations


analyze_infections <- function(scenario.type = "nointer_baseline",
                               param.sets = seq(1, 1)
){
  
  # empty data.frame for storing results
  full.summary.results <- data.frame(mean.piri = numeric(),
                                     med.piri = numeric(),
                                     ran.piri1 = numeric(),
                                     ran.piri2 = numeric(),
                                     mean.pi = numeric(),
                                     med.pi = numeric(),
                                     ran.pi1 = numeric(),
                                     ran.pi2 = numeric(),
                                     mean.dur = numeric(),
                                     med.dur = numeric(),
                                     ran.dur1 = numeric(),
                                     ran.dur2 = numeric(),
                                     scenario.type = character(),
                                     param.set = numeric(),
                                     stringsAsFactors = F)
  

  for(i in 1:length(param.sets)){
    
    #### READ IN RESULTS FILES ####
    temp.param.num <- param.sets[i]
      
    full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
    full.sims.results <- get(load(full.name))  
    
    
    #### CALCULATE SUMMARY STATISTICS ####
    # mean infectious (prog and regr)
    pi_ri <- full.sims.results$total.prog+full.sims.results$total.lr
    mean.piri <- mean(pi_ri)
    
    # median infectious (prog and regr)
    med.piri <- median(pi_ri)
    
    # range infectious (prog and regr)
    ran.piri <- range(pi_ri)
    
    # mean prog infections (i.e. mortalities)
    pi <- full.sims.results$total.prog
    mean.pi <- mean(pi)
     
    # median prog infections (i.e. mortalities)
    med.pi <- median(pi)
    
    # range prog infections (i.e. mortalities)
    ran.pi <- range(pi)
    
    # mean duration of outbreak
    dur <- full.sims.results$dur.time
    if(any(is.na(full.sims.results$dur.time))){
      print("Warning: NA outbreak duration")
    }
    mean.dur <-  mean(dur)
    
    # median duration 
    med.dur <- median(dur)
    
    # range of durations
    ran.dur <- range(dur)
    
    
    ### temporary results ###
    temp.results <- data.frame(mean.piri = mean.piri,
                               med.piri = med.piri,
                               ran.piri1 = ran.piri[1],
                               ran.piri2 = ran.piri[2],
                               mean.pi = mean.pi,
                               med.pi = med.pi,
                               ran.pi1 = ran.pi[1],
                               ran.pi2 = ran.pi[2],
                               mean.dur = mean.dur,
                               med.dur = med.dur,
                               ran.dur1 = ran.dur[1],
                               ran.dur2 = ran.dur[2],
                               scenario.type = paste(scenario.type),
                               param.set = temp.param.num,
                               stringsAsFactors = F
    )
    
   full.summary.results[i,] <- temp.results
  }
  

  #### RETURN RESULTS ####
  return(full.summary.results)

}
