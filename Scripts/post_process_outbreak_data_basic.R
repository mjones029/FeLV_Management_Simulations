# process_birth_outbreak_basic
#
#========================================================	
# ---
### title: Process outbreaks with births/respawning
# author: Marie Gilbertson
# date: "10/10/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for processing outbreak data when the outbreak involved "births"
# The point is to produce a new "m" matrix where the individuals that occupied a dead territory now have their own column

process_birth_outbreak <- function(m=m){

  m.orig <- m # matrix output from trarns_sim_forcedur.R function
  colnames(m.orig) <- seq(1, ncol(m.orig))
  
  m.new <- data.frame(m.orig)
  colnames(m.new) <- seq(1,ncol(m.new))
  
  results <- NULL
  
  for(j in 1:ncol(m.orig)){
    
    # print(j)
    
    # pull out a column
    temp.m <- m.orig[,j]
    
    # change vaccine tracking number so that 4-->0 transition is the only one with a different in status number of 4
    temp.m[temp.m==6] <- 12
    
    # determine if that column has a 4-->0 transition
    transitions <- NA
    for(i in 1:(length(temp.m)-1)){
      now <- temp.m[i]
      now.1 <- temp.m[i+1]
      
      transitions <- c(transitions, (now-now.1))
      
    }
    
    if(any(na.omit(transitions)==4)){
    # save the data about which individuals/territories have a 4-->0 transition and the time point(s) of the transition(s)
    temp.results <- data.frame(id = paste(j),
                               trans.time = which(transitions==4)
                                )
  
    results <- rbind(results, temp.results)
    
    # fill in the original column with 4s to represent the "ongoing" death of that individual
    m.new[c(which(transitions==4)[1]:nrow(m.orig)),j] <- 4
    }
  }
  
  if(is.null(results)==F){ # basically, check if any territories were reoccupied. If not, return the original "m" data 
    
    # for all re-occupied "territories", pull the 0 onward data and add to new column, with a column name that reflects the "territory" it took over
    
    m.new.all <- data.frame(matrix(NA, ncol = 0, nrow = nrow(m.orig)))
    ids.to.eval <- unique(results$id)
    
    for(j in 1:length(ids.to.eval)){
      
      temp.id <- ids.to.eval[j]
      temp <- results[results$id==temp.id,]
      
      if(nrow(temp)==1){
        
        m.new_sing <- data.frame(matrix(NA, ncol = 1, nrow = nrow(m.orig)))
        colnames(m.new_sing) <- "new.ind"
        
        temp.m.orig <- m.orig[,paste(temp$id)]
        
        m.new_sing$new.ind[c(temp$trans.time:nrow(m.new_sing))] <- temp.m.orig[c(temp$trans.time:nrow(m.new_sing))]
    
        # rename the column to represent the territory that was taken over
        colnames(m.new_sing) <- paste(temp$id, ".2", sep = "")
        
        }else{
          # need to deal with the possibility that a territory could be re-occupied more than once
          m.new_multi <- data.frame(matrix(NA, ncol = 0, nrow = nrow(m.orig)))
    
        for(k in 1:nrow(temp)){
          m.new_multi_temp <- data.frame(matrix(NA, ncol = 1, nrow = nrow(m.orig)))
          colnames(m.new_multi_temp) <- "new.ind"
          
          temp.m.orig <- m.orig[,unique(paste(temp$id))]
          
          trans.time_multi_start <- temp$trans.time[k]
          trans.time_multi_end <- ifelse(k==nrow(temp), nrow(m.new_multi_temp), (temp$trans.time[k+1]-1))
          
          m.new_multi_temp$new.ind[c(trans.time_multi_start:trans.time_multi_end)] <- temp.m.orig[c(trans.time_multi_start:trans.time_multi_end)]
          
          # for all "not last occupant" individuals, fill in remainder of occupancy with 4s to indicate ongoing death
          if(k!=nrow(temp)){
            m.new_multi_temp$new.ind[c((trans.time_multi_end+1):nrow(m.new_multi_temp))] <- 4
          }
          # rename the column to represent the territory that was taken over
          colnames(m.new_multi_temp) <- paste(unique(temp$id), k+1, sep = ".")
          
          m.new_multi <- cbind(m.new_multi, m.new_multi_temp)
        }
          m.new_sing <- m.new_multi
      }
      
      m.new.all <- cbind(m.new.all, m.new_sing)
    }
    
    
    m.new <- as.matrix(cbind(m.new, m.new.all))
  
  } else {
    m.new <- m.orig
  }

return(m.new)
}