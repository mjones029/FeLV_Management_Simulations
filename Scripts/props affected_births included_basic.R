# props affected_births included_basic
#
#========================================================	
# ---
### title: Proportions affected in outbreak with births
# author: Marie Gilbertson
# date: "03/16/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function that processes outbreak results when births/respawnings are included for the purposes of analysis/plotting
# More specifically, calculates proportions in each disease category over time.

props_affected_births_included <- function(m.new = m.new){

  #### prep data ####
  # new dataframe of m.new to work with (and conserve original m.new data)
  p.data <- m.new

  
  # make all not-born-yet individuals class "10" for ease in calculating population sizes
  p.data[is.na(p.data)] <- 10 
  
  # make data frame for population results
  p.results <- data.frame(
    time=numeric(nrow(m.new)),
    pop.size=numeric(nrow(m.new)),
    prop.s=numeric(nrow(m.new)),
    prop.i=numeric(nrow(m.new)),
    prop.lr=numeric(nrow(m.new)),
    prop.r=numeric(nrow(m.new)),
    prop.v=numeric(nrow(m.new)),
    num.d=numeric(nrow(m.new))
  )
  
  p.results$time <- seq(0,nrow(m.new)-1)
  
  
  #### calculate population size per time step ####
  # population size per time, subtracting all dead and not-born-yet individuals
  p.results$pop.size <- apply(p.data!=4 & p.data!=10, 1, sum)
  
  
  
  
  #### calculate propotion of population in each disease category per time step ####
  
  p.results$prop.s <- apply(p.data==0, 1, sum)/p.results$pop.size
  p.results$prop.i <- apply(p.data==1, 1, sum)/p.results$pop.size
  p.results$prop.lr <- apply(p.data==2, 1, sum)/p.results$pop.size
  p.results$prop.r <- apply(p.data==3, 1, sum)/p.results$pop.size
  p.results$num.d <- apply(p.data==4, 1, sum) # just record NUMBER of dead
  p.results$prop.v <- apply(p.data==6, 1, sum)/p.results$pop.size
  
  
  
  #### output data results ####
  
  return(p.results)

}