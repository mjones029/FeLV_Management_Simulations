# analyze_reactive_vaccinations
#
#========================================================	
# ---
### title: Analyze simulation vaccinations (reactive)
# author: Marie Gilbertson
# date: "10/11/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function that outputs mean, median, and range results for simulations; focuses on 
# factors specifically relevant to reactive vaccination scenarios (e.g. number of vaccines used)


analyze_reactive_vaccinations <- function(scenario.type = "reactvax_baseline",
                                          param.sets = seq(1, 32)
){
  
  # empty data.frame for storing results
  full.summary.results <- data.frame(mean.tv = numeric(),
                                     med.tv = numeric(),
                                     ran.tv1 = numeric(),
                                     ran.tv2 = numeric(),
                                     mean.tvu = numeric(),
                                     med.tvu = numeric(),
                                     ran.tvu1 = numeric(),
                                     ran.tvu2 = numeric(),
                                     mean.vpi = numeric(),
                                     med.vpi = numeric(),
                                     ran.vpi1 = numeric(),
                                     ran.vpi2 = numeric(),
                                     scenario.type = character(),
                                     param.set = numeric(),
                                     stringsAsFactors = F)
  
  
  for(i in 1:length(param.sets)){
    
    #### READ IN RESULTS FILES ####
    temp.param.num <- param.sets[i]
    
    full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
    full.sims.results <- get(load(full.name))  
    
    
    #### CALCULATE SUMMARY STATISTICS ####
    
    # mean total.vaxed
    total.vaxed <- full.sims.results$total.vaxed
    mean.tv <- mean(total.vaxed)
    
    # median total.vaxed
    med.tv <- median(total.vaxed)
    
    # range total.vaxed
    ran.tv <- range(total.vaxed)
    
    # mean total.vax.used
    total.vax.used <- full.sims.results$total.vax.used
    mean.tvu <- mean(total.vax.used)
    
    # median total.vax.used
    med.tvu <- median(total.vax.used)
    
    # range total.vax.used
    ran.tvu <- range(total.vax.used)
    
    # mean vax per vaccinated individual
    vpi <- full.sims.results$total.vax.used/full.sims.results$total.vaxed
    
    # check if none were vaccinated
    if(any(is.nan(vpi))){
      nans <- which(is.nan(vpi))
      for(j in 1:length(nans)){
        if(full.sims.results$total.vaxed[nans[j]]==0){
          vpi[nans[j]] <- 0
        }else{
          print("Warning: NaN not resolved")
        }
      }
    }
    
    
    mean.vpi <- mean(vpi)
    
    # median vax per vaccinated individual
    med.vpi <- median(vpi)
    
    # range vax per vaccinated individual
    ran.vpi <- range(vpi)
    
    
    ### temporary results ###
    temp.results <- data.frame(mean.tv = mean.tv,
                               med.tv = med.tv,
                               ran.tv1 = ran.tv[1],
                               ran.tv2 = ran.tv[2],
                               mean.tvu = mean.tvu,
                               med.tvu = med.tvu,
                               ran.tvu1 = ran.tvu[1],
                               ran.tvu2 = ran.tvu[2],
                               mean.vpi = mean.vpi,
                               med.vpi = med.vpi,
                               ran.vpi1 = ran.vpi[1],
                               ran.vpi2 = ran.vpi[2],
                               scenario.type = paste(scenario.type),
                               param.set = temp.param.num,
                               stringsAsFactors = F
    )
    
    full.summary.results[i,] <- temp.results
  }
  
  
  #### RETURN RESULTS ####
  return(full.summary.results)
  
}
