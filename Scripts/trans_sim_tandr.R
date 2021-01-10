# trans_sim_tandr
#
#========================================================	
# ---
### title: Function for simulating transmission on networks; includes test-and-removal (tandr)
# author: Marie Gilbertson
# date: "10/09/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates transmission of FeLV on networks with inclusion of test-and-removal (tandr)
# 2. Adapted from "trans_sim.R" script from FIV-->FeLV project




trans_sim_tandr <- function(g = g,                           # network on which to simulate (should be an igraph object)
                            beta = beta,                     # base transmission rate (applies to all progressive individuals)
                            c.beta = c.beta,                 # constant applied to beta to determine transmission rate of latent/regressive individuals 
                            c.rate = c.rate,                 # "contact rate" which is just a lower probability of transmission occuring
                            duration.yrs = duration.yrs,     # number of years to run simulation
                            prop.outcomes = prop.outcomes,   # proportion per each infected state
                            progr.dur = progr.dur,           # weekly probability of progressive mortality (= 1/progressive infection duration in weeks)
                            regr.dur.c = regr.dur.c,         # regressive infection duration multiplier
                            terr.repop = terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                            initiate.int = initiate.int,     # time point to initiate intervention (in weeks); in this case, to initiate test and removal
                            max.occupancy = max.occupancy,   # maximum number of temporary removals possible at one time
                            dist.type = "random",            # distribution of captures; "random" or "spatial"
                            int.dur = 17                     # duration of capture season each year (17 weeks ~ 4 months)
                            ) {
  
  # assign pop.size value 
  pop.size <- length(V(g))
  
  
  
  
  
  ############# Network model ########################
  
  ER<-get.adjacency(g) #make an adjacency network; convert graph to an adjacency matrix or an edge list
  
  
  netwb <- ER
  
  
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
    
    # set intervention baselines
    int.switch <- "off"
    occupancy <- data.frame(matrix(nrow = 1, ncol = 2)) # track occupancy of "temporary removal" state
    colnames(occupancy) <- c("id", "release.time")
    occupancy <- na.omit(occupancy)
    total.tr <- c() # total test-and-removes
    attempted.tr <- c() # total attempts (i.e. captures)
    attempt.time <- c() # time point of attempts
    release.data <- NULL # tracking info about releases after recovery
    
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
    
    
    # if using spatial captures, class locations of all individuals according to side of freeway of initial infection
    if(dist.type == "spatial"){
      # if using spatial vaccination, set up dataframe of initial vaccination candidates to consider
      locations.df <- data.frame(id = V(g)$vertex.names,
                                 lat = V(g)$latitude,
                                 long = V(g)$longitude,
                                 status = stat,
                                 i75.lat = 26.16)
      
      locations.df$i75_NS <- ifelse(locations.df$lat >= locations.df$i75.lat, "N", "S")
      
      init.lat <- vertex_attr(g, "latitude", index = which(V(g)$vertex.names==init.inf))
      init.side <- ifelse(init.lat>=26.16, "N", "S")
      
      locations.df$candidate <- ifelse(locations.df$i75_NS==init.side, 1, 0)
      locations.df <- locations.df[locations.df$candidate==1,]
    }
    
    
    
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
    
    #### start of main body of transmission simulation ####
    while ((any(stat==1)|any(stat==2))&time<=duration) {
      # for (l in 1:duration) {  
      # print(l)
      #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated; 7 means temporarily removed
      
      # for each (weekly) time step
      statc <- stat  #assign stat to copy of stat so can change 'statc' while still working off of 'stat'
      
      
      #### TEST AND REMOVAL CAPTURES ####
      ## decide if in a capture season ##
      if((time-initiate.int)%%52==0){
        int.switch <- "on"
        print(int.switch) # can monitor status of capture seasons
      }
      if((time-(initiate.int + int.dur))%%52==0){
        int.switch <- "off"
        print(int.switch)
      }
      
      ## if so, decide if any room left for temp removals ##
      if(int.switch=="on" & (nrow(occupancy)<max.occupancy)){
        
        # capture an alive/present individual
        currently.dead <- which(stat==4 | stat==7)
        
        candidates <- which(!V(g)$vertex.names %in% currently.dead)
        
        
        if(dist.type=="random"){  
          # if random captures, capture any alive individual
          cap.ind <- sample(candidates, 1) 
          
        }else if(dist.type=="spatial"){
          # if spatial captures, capture based on side of freeway of initially infected individual
          candidates.df <- locations.df[locations.df$id %in% candidates,]
          cap.ind <- sample(as.numeric(as.character(candidates.df$id)), 1)
        }
        
        # if state = 1, euthanize (convert to state = 4)
        if(stat[cap.ind]==1){
          statc[cap.ind] <- 4
          total.tr_temp <- "e" # track type of "removal"
          
          # at death, veff.stat drops to 0 so respawns aren't spawned with vaccine efficacy
          veff.stat[cap.ind] <- 0
          
          
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
          # (first need to clear out all the "1"s to account for "premature" death)
          n[c((time+1):nrow(n)), cap.ind] <- 0
          n[c((time+1):(time+death.dur)),cap.ind] <- 4 
          
        }else if(stat[cap.ind]==2){
          # if state = 2, add to occupancy pool, convert to special state 7 (treat like a 4, but don't consider a disease mortality)
          statc[cap.ind] <- 7
          unadjusted <- min(which(n[c(time:nrow(n)),cap.ind]==0))
          release.time <- (unadjusted+time)-1 # determine time point of recovery
          temp_occupancy <- data.frame(id = cap.ind,
                                       release.time = release.time)
          occupancy <- rbind(occupancy, temp_occupancy)
          
          total.tr_temp <- "r" # record type of "removal"
          
          # at removal, veff.stat drops to 0 so respawns aren't spawned with vaccine efficacy
          veff.stat[cap.ind] <- 0
          
          
          # Assign a "duration of death" (in this case, un-occupancy) to account for new "births" occupying empty territories
          
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
            death.dur <- duration+1 # will put 7s out too far, but does the job
          }
          
          # then assign that unoccupancy duration in the n matrix
          # (first need to clear all the instances of "2"s)
          n[c((time+1):nrow(n)), cap.ind] <- 0
          n[c((time+1):(time+death.dur)),cap.ind] <- 7 
          
        }else if(stat[cap.ind] %in% c(0, 3, 5, 6)){
          # if state = 0, 3, 5, or 6, do nothing
          total.tr_temp <- stat[cap.ind]
        }else{
          print("Warning: capture state error")
        }
        
        # monitor outcomes of each capture
        total.tr <- c(total.tr, total.tr_temp)
        # record "attempted test and removals" 
        attempted.tr <- c(attempted.tr, cap.ind)
        attempt.time <- c(attempt.time, time)
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
      
      #### TRANSMISSION TO VACCINATES ####
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
              if(state==3){ # assign "duration" for regressives (for record keeing purposes)
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
      
      #### RELEASE OF RECOVERED REMOVALS ####
      # do before respawning steps so no conflict between release and respawning
      if(any(occupancy$release.time<=time) & any(stat==4 | stat==7)){
        release.cands <- occupancy$id[occupancy$release.time<=time] # determine candidates for release based on recovery times
        # only one release per week
        if(length(release.cands)>1){
          old.id <- as.numeric(as.character(sample(release.cands, 1)))
        }else if(length(release.cands)==1){
          old.id <- release.cands
        }else{
          print("Warning: release candidates fail")
        }
        
        # assign individual selected for release to an "unoccupied territory" 
        if(length(which(stat==4 | stat==7))==1){
          new.id <- which(stat==4 | stat==7)
        }else if(length(which(stat==4 | stat==7))>1){
          new.id <- sample(which(stat==4 | stat==7), 1)
        }else{
          print("Warning: new.id assignment error")
        }
        # update state for tracking and record this "release data"
        statc[new.id] <- 5
        n[c((time+1):(duration+1)),new.id] <- 5
        temp.release.data <- data.frame(old.id = old.id,
                                        new.id = new.id,
                                        release.time = time)
        release.data <- rbind(release.data, temp.release.data)

        # update occupancy
        occupancy <- occupancy[occupancy$id != old.id,]
      }
      
      
      #### RESPAWNING OF DECEASED ####
      # determine which territories are re-occupied by a new, susceptible individual ("births") after prior individual's death
      for (i in (which(stat==4)))     #for every dead individual in original stat
      {
        if (n[time+1,i]==0) # determine which territories are re-occupied
        {
          statc[i] <- 0
          n[c((time+1):(duration+1)),i] <- 0
          # THIS FORMUALTION MEANS THAT A SINGLE COLUMN COULD REPRESENT MORE THAN ONE INDIVIDUAL
          # MUST BE CAPTURED IN "POST-PROCESSING"
        }
      }
      
      
      #### RESPAWNING IN REMOVAL TERRITORIES ####
      # determine which territories are re-occupied by a new, susceptible individual ("births") after prior individual's removal
      for (i in (which(stat==7)))     #for every dead individual in original stat
      {
        if (n[time+1,i]==0) # determine which territories are re-occupied
        {
          statc[i] <- 0
          n[c((time+1):(duration+1)),i] <- 0
          # THIS FORMUALTION MEANS THAT A SINGLE COLUMN COULD REPRESENT MORE THAN ONE INDIVIDUAL
          # MUST BE CAPTURED IN "POST-PROCESSING"
        }
      }
      
      
     
      
      
      stat <- statc
      m[time+1,] <- stat
      write(c(time,m[time+1,]),file="test_sim.txt",ncolumns=(gp+1),append=TRUE)
      time<-(time+1)
      # print(time)
    }
    
    m <- m[rowSums(m)>0,]
    
    # process outbreak results
    m.new <- process_birth_outbreak_tandr(m=m,
                                          occupancy = occupancy,
                                          release.data = release.data,
                                          attempted.tr = attempted.tr,
                                          attempt.time = attempt.time,
                                          total.tr = total.tr
                                          )
    # check if 5 or more individuals were ever infectious; if not, repeat simulation
    #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated; 7 means removed
    if(sum(apply(m.new, 2, function(r) any(r %in% c(1, 2)))) < 5){
      fails <- c(fails, "f")
      print("fail")
    }else if(sum(apply(m.new, 2, function(r) any(r %in% c(1, 2)))) >= 5) break
  }
  
  
  m.list <- list(m.new, initiate.dur, fails, attempted.tr, attempt.time, total.tr, release.data, occupancy)
  names(m.list) <- c("m.new", "initiate.dur", "fails", "attempted.tr", "attempt.time", "total.tr", "release.data", "final.occupancy")
  
  return(m.list)
  
}
