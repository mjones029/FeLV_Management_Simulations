# analyze_reactive_tandr
#
#========================================================	
# ---
### title: Analyze simulation test and removal (reactive)
# author: Marie Gilbertson
# date: "10/30/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function that outputs mean, median, and range results for simulations; specific to results
# of interest/unique to test-and-removal scenarios (e.g. proportion of captures that were "successful")


analyze_reactive_tandr <- function(scenario.type = "reacttandr_baseline",
                                          param.sets = seq(1, 16)
){
  
  # empty data.frame for storing results
  full.summary.results <- data.frame(mean.tr = numeric(),
                                     med.tr = numeric(),
                                     ran.tr1 = numeric(),
                                     ran.tr2 = numeric(),
                                     mean.te = numeric(),
                                     med.te = numeric(),
                                     ran.te1 = numeric(),
                                     ran.te2 = numeric(),
                                     mean.psc = numeric(),
                                     med.psc = numeric(),
                                     ran.psc1 = numeric(),
                                     ran.psc2 = numeric(),
                                     scenario.type = character(),
                                     param.set = numeric(),
                                     stringsAsFactors = F)
  
  
  for(i in 1:length(param.sets)){
    
    #### READ IN RESULTS FILES ####
    temp.param.num <- param.sets[i]
    
    full.name <- paste("Simulation_results/", scenario.type, "/fulldata_", scenario.type, "_paramset", temp.param.num, ".Rdata", sep = "")
    full.sims.results <- get(load(full.name))  
    
    
    #### CALCULATE SUMMARY STATISTICS ####
    
    # mean total removed
    tr <- full.sims.results$total.rems
    mean.tr <- mean(tr)
    
    # median total removed
    med.tr <- median(tr)
    
    # range total removed
    ran.tr <- range(tr)
    
    # mean total euthanized
    te <- full.sims.results$total.euth
    mean.te <- mean(te)
    
    # median total euthanized
    med.te <- median(te)
    
    # range total euthanized
    ran.te <- range(te)
    
    # mean proportion "successful" captures
    psc <- (tr + te)/full.sims.results$total.caps
    
    
    # check if none were successfully captured
    if(any(is.nan(psc))){
      nans <- which(is.nan(psc))
      for(j in 1:length(nans)){
        if(full.sims.results$total.caps[nans[j]]==0){
          psc[nans[j]] <- 0
        }else{
          print("Warning: NaN not resolved")
        }
      }
    }
    
    mean.psc <- mean(psc)
    
    # median proportion "successful" captures
    med.psc <- median(psc)
    
    # range proportion "successful" captures
    ran.psc <- range(psc)
    
        
    
    
    
    ### temporary results ###
    temp.results <- data.frame(mean.tr = mean.tr,
                               med.tr = med.tr,
                               ran.tr1 = ran.tr[1],
                               ran.tr2 = ran.tr[2],
                               mean.te = mean.te,
                               med.te = med.te,
                               ran.te1 = ran.te[1],
                               ran.te2 = ran.te[2],
                               mean.psc = mean.psc,
                               med.psc = med.psc,
                               ran.psc1 = ran.psc[1],
                               ran.psc2 = ran.psc[2],
                               scenario.type = paste(scenario.type),
                               param.set = temp.param.num,
                               stringsAsFactors = F
    )
    
    full.summary.results[i,] <- temp.results
  }
  
  
  #### RETURN RESULTS ####
  return(full.summary.results)
  
}
