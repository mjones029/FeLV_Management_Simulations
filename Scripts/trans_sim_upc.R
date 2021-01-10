# trans_sim_upc
#
#========================================================	
# ---
### title: Function for simulating transmission on networks with exponential wait time for duration of infectiousness
# author: Marie Gilbertson
# date: "10/09/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates transmission of FeLV on networks, but with underpass closures factored in
# 2. Adapted from "trans_sim_geodur.R" script from FIV-->FeLV project




trans_sim_upc <- function(g = new.g2,                           # network on which to simulate (should be an igraph object)
                               beta = beta,                     # base transmission rate (applies to all progressive individuals)
                               c.beta = c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                               c.rate = c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                               duration.yrs = duration.yrs,     # number of years to run simulation
                               prop.outcomes = prop.outcomes,   # proportion per each infected state
                               progr.dur = progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                               regr.dur.c = regr.dur.c,         # regressive infection duration multiplier
                               terr.repop = terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                               initiate.upc = 26,               # time point to initiate underpass closures (in weeks); for simplicity, will occur "instantaneously"
                               upc.dur = 4                      # duration of underpass closures (also occurs instantaneously)
                               ) {
  
  # assign pop.size value 
  pop.size <- length(V(g))
  
  
  
  
  
  ############# Network model ########################
  
  ER_full<-get.adjacency(g) #make an adjacency network; convert graph to an adjacency matrix or an edge list
  
  
  netwb <- ER_full
  
  #### assign FeLV model parameters ####
  beta<-beta
  sigma <- c.beta*beta # transmission rate for regressively infected individuals
  regr.dur <- regr.dur.c*progr.dur # multiply constant by the progressive duration distribution to get distribution for regressives
  duration <- 52*duration.yrs # 52 * number of years to run simulation
  
  # track number of failed epidemics
  fails <- c()
  
  
  
  repeat{ # repeat to account for epidemics that "fail"
    
    #### run the disease simulation ####
    
    time<-1
    gp <- ncol(netwb)    # number of columns in network    
    
    # create vector of infection status for each individual
    stat <- V(g)$status
    
    # create vector of vax efficacy status for each individual
    veff.stat <- V(g)$veff
    
    
    
    
    # randomly select a susceptible individual initiating infection in the population
    # for efficiency, only keep instances in which infection initiates in a non-isolate
    # but keep record of how many iterations before a non-isolate is selected
    initiate <- c()
    repeat {
      
      init.inf <- sample(which(stat==0),1)
      initiate <- c(initiate, init.inf)
      
      # exit loop when infection initiates in a non-isolate
      if (any(netwb[init.inf,]==1) & any(netwb[,init.inf]==1)) break
    }
    
    stat[init.inf] <- 1
    # "initiate.dur" gives the number of iterations before a non-isolate was selected
    # in other words, the number of "failed" epidemics + 1
    
    initiate.dur <- length(initiate)
    
    
    m <- matrix(0,duration+1,gp)  # make empty matrix for number of time steps long by number of individuals wide (this will be the progression of infection per time step per indiv.) 52 plus 1 more for time=0
    m[time,] <- stat
    
    
    # create matrix to accommodate duration of infection assignments
    n <- matrix(0,duration*100,gp) # need this to be extra long to accommodate duration of infection past end of monitoring
    
    # assign duration of infection for first infected individual
    recovery.all <- c()
    repeat {
      # repeatedly draw until recovery = 1, recording number of draws until that point
      recovery <- rbinom(1,1,progr.dur)
      recovery.all <- c(recovery.all, recovery)
      
      # exit loop when death = 1
      if (recovery > 0) break
    }
    dur.1 <- length(recovery.all)
    n[c(1:dur.1),which(stat==1)] <- 1
    
    
    write(c(0,stat),file="test_sim.txt",ncolumns=(gp+1),append=F) #writes initial conditions to file
    
    #### proceed with main body of transmission simulations ####
    while ((any(stat==1)|any(stat==2))&time<=duration) {
      #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
      
      # for each (weekly) time step
      statc <- stat  #assign stat to copy of stat so can change 'statc' while still working off of 'stat'
      
      
      #### UNDERPASS CLOSURES/OPENINGS ####
      
      # If time to close underpasses
      if(time==initiate.upc){ 
        # delete edges that cross I75
        # 1 = edge does NOT cross I75; 2 = edge DOES cross I75
        g.upc <- g
        g.upc <- delete.edges(g.upc, which(E(g.upc)$up.weight==2))
        
        # update "netwb" to reflect change in adjacency matrix
        ER_upc<-get.adjacency(g.upc)
        netwb <- ER_upc
        
        print("closed") # can track progress
      }
      
      # if time to reopen underpasses
      if(time==(initiate.upc+upc.dur)){  
        # return to original adjacency matrix
        netwb <- ER_full

        print("re-opened")
      }
      
      
      
      
      
      #### TRANSMISSION TO SUSCEPTIBLES ####
      for (i in (which(stat==0))){     #for every susceptible individual in original stat...
        if(rbinom(1,1,c.rate)==1){    #does that individual have a contact this time step? If so...
          for (j in (1:gp)[-i]){          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
            
            if ((netwb[i,j]>=1) & (((stat[j]==1)&(rbinom(1,1,beta)==1)) | ((stat[j]==2)&(rbinom(1,1,sigma)==1))) ) #is there an edge to an infected individual and does a transmission event take place?  (By looking down column of that indiv in netw for a 1)  
            {
              
              # does the infected individual become progressively, latently, or regressively infected?
              state <- sample(c(1:3), size = 1, prob = prop.outcomes) # 1 = progressive, 2 = regressive, 3 = abortive
              statc[i] <- state
              
              if(state==1){ 
                # assign duration of infection
                recovery.all <- c()
                repeat {
                  # repeatedly draw until recovery = 1, recording number of draws until that point
                  recovery <- rbinom(1,1,progr.dur)
                  recovery.all <- c(recovery.all, recovery)
                  
                  # exit loop when death = 1
                  if (recovery > 0) break
                }
                dur.temp <- length(recovery.all)
                n[c((time+1):(time+dur.temp)),i] <- 1
              }
              if(state==2){
                # assign duration of infection
                recovery.all <- c()
                repeat {
                  # repeatedly draw until recovery = 1, recording number of draws until that point
                  recovery <- rbinom(1,1,regr.dur)
                  recovery.all <- c(recovery.all, recovery)
                  
                  # exit loop when death = 1
                  if (recovery > 0) break
                }
                dur.temp <- length(recovery.all)
                n[c((time+1):(time+dur.temp)),i] <- 2
              }
              if(state==3){ # assign "duration" for abortive infections (for record keeing purposes)
                n[c((time+1):(duration+1)),i] <- 3 # this bit only matters for record keeping
              }
              break # assign outcome of infection to the copy and get out of loop
            }
          }
        }
      }
      
      #### TRANSMISSION TO VACCINATED INDNIVIDUALS ####
      # loop through all vaccinated individuals (same basic loop as for susceptibles, but now accounting for vaccine efficacy)
      for (i in (which(stat==6))){     #for every vaccinated individual in original stat...
        if(rbinom(1,1,c.rate)==1){    #does that individual have a contact this time step? If so...
          for (j in (1:gp)[-i]){          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
            
            if ((netwb[i,j]>=1) & (((stat[j]==1)&(rbinom(1,1,beta)==1)&(rbinom(1,1,veff.stat[i])==0)) | ((stat[j]==2)&(rbinom(1,1,sigma)==1)&(rbinom(1,1,veff.stat[i])==0))) ) #is there an edge to an infected individual and does a transmission event take place and does the vaccine fail to prevent infection? 
            {
              
              # does the infected individual become progressively, latently, or regressively infected?
              state <- sample(c(1:3), size = 1, prob = prop.outcomes) # 1 = progressive, 2 = regressive, 3 = abortive
              statc[i] <- state
              
              if(state==1){ 
                # assign duration of infection
                recovery.all <- c()
                repeat {
                  # repeatedly draw until recovery = 1, recording number of draws until that point
                  recovery <- rbinom(1,1,progr.dur)
                  recovery.all <- c(recovery.all, recovery)
                  
                  # exit loop when death = 1
                  if (recovery > 0) break
                }
                dur.temp <- length(recovery.all)
                n[c((time+1):(time+dur.temp)),i] <- 1
              }
              if(state==2){
                # assign duration of infection
                recovery.all <- c()
                repeat {
                  # repeatedly draw until recovery = 1, recording number of draws until that point
                  recovery <- rbinom(1,1,regr.dur)
                  recovery.all <- c(recovery.all, recovery)
                  
                  # exit loop when death = 1
                  if (recovery > 0) break
                }
                dur.temp <- length(recovery.all)
                n[c((time+1):(time+dur.temp)),i] <- 2
              }
              if(state==3){ # assign "duration" for regressives (for record keeping purposes)
                n[c((time+1):(duration+1)),i] <- 3 # this bit only matters for record keeping
              }
              break # assign outcome of infection to the copy and get out of loop
            }
          }
        }
      }
      
      #### DEATH OF PROGRESSIVES ####
      for (i in (which(stat==1)))     #for every progressively infectious individual in original stat
      {
        if (n[time+1,i]==0) # determine which individuals die from infection
        {
          statc[i] <- 4
          
          # at death, veff.stat drops to 0 so respawns aren't spawned with vaccine efficacy
          veff.stat[i] <- 0
          
          # Assign a "duration of death" to account for new "births" occupying empty territories
          # first, draw binomial trials to determine time until territory is repopulated
          # if/else loop to accommodate no birth/repopulation rate
          if(terr.repop>0){ # if territories are repopulable, assign duration of time until repopulation
            death.all <- c()
            
            repeat {
              # repeatedly draw until death = 1, recording number of draws until that point
              death <- rbinom(1,1,terr.repop)
              death.all <- c(death.all, death)
              
              # exit loop when death = 1
              if (death > 0) break
            }
            death.dur <- length(death.all) # assign time until repopulation of territory 
          }else{ # if territories are not repopulable, just assign out death duration through end of simulation
            death.dur <- duration+1 # will put 4s out too far, but does the job
          }
          
          # then assign that death duration in the n matrix
          n[c((time+1):(time+death.dur)),i] <- 4 
          
        }
      }
      
      #### RECOVERY OF REGRESSIVES ####
      for (i in (which(stat==2)))     #for every latently/regressively infectious individual in original stat
      {
        if (n[time+1,i]==0) # determine which individuals recover from infection
        {
          statc[i] <- 5
          n[c((time+1):(duration+1)),i] <- 5
        }
      }
      
      #### RESPAWNING ####
      # determine which territories are re-occupied by a new, susceptible individual ("births") after prior individual's death
      for (i in (which(stat==4)))     #for every dead individual in original stat
      {
        if (n[time+1,i]==0) # determine which territories are re-occupied
        {
          statc[i] <- 0
          n[c((time+1):(duration+1)),i] <- 0
          # THIS FORMUALTION MEANS THAT A SINGLE COLUMN COULD REPRESENT MORE THAN ONE INDIVIDUAL
          # MUST BE CAPTURED IN "POST-PROCESSING"
          # doing this way because adding new columns for new individuals will also mess things up... (loop functions by looping through original population size)
        }
      }
      
      stat <- statc
      m[time+1,] <- stat
      write(c(time,m[time+1,]),file="test_sim.txt",ncolumns=(gp+1),append=TRUE)
      time<-(time+1)
      # print(time)
      
    }
    
    m <- m[rowSums(m)>0,]
    
    # process outbreak results with external function
    m.new <- process_birth_outbreak(m)
    
    # check if 5 or more individuals were ever infectious; if not, repeat simulation
    #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
    if(sum(apply(m.new, 2, function(r) any(r %in% c(1, 2)))) < 5){
      fails <- c(fails, "f")
      print("fail")
    }else if(sum(apply(m.new, 2, function(r) any(r %in% c(1, 2)))) >= 5) break
  }
  
  
  m.list <- list(m.new, initiate.dur, fails)
  names(m.list) <- c("m.new", "initiate.dur", "fails")
  return(m.list)
  
}