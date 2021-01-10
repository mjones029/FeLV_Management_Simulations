# process_birth_outbreak_tandr
#
#========================================================	
# ---
### title: Process outbreaks with births; specific for test-and-removal
# author: Marie Gilbertson
# date: "10/10/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for processing outbreak data when the outbreak involved "births"; SPECIFIC TO TEST AND REMOVAL SCENARIOS
# The point is to produce a new "m" matrix where the individuals that occupied a dead territory now have their own column
# Test and removal needs unique processes because of the additional way individuals can be removed and replaced within the population

process_birth_outbreak_tandr <- function(m=m,
                                         occupancy = occupancy,
                                         release.data = release.data,
                                         attempted.tr = attempted.tr,
                                         attempt.time = attempt.time,
                                         total.tr = total.tr
                                         ){
  
  m.orig <- m # matrix output from trarns_sim_forcedur.R function
  colnames(m.orig) <- seq(1, ncol(m.orig))
  
  m.new <- data.frame(matrix(ncol = 0, nrow = nrow(m.orig)))
  all.release.records <- data.frame(matrix(ncol = 0, nrow = nrow(m.orig)))
  
  # loop through each column (represents a "territory") and convert to individual histories
  for(i in 1:ncol(m.orig)){
    respawn.state <- "off"
    changes <- c()
    temp.m <- m.orig[,i]
    all.new.temp <- data.frame(temp.m)
    release.records <- data.frame(matrix(nrow = nrow(m.orig), ncol = 0))
    repeat{
      # loop through each time point
      for(j in 2:length(temp.m)){
        # print(j)
        # look for 4-->0 transition (indicates respawning)
        if(temp.m[j-1]==4 & temp.m[j]==0)
          {
          temp.new.ind <- temp.m[c(j:length(temp.m))]
          empty <- rep(20, length(temp.m)-length(temp.new.ind)) # assign a "20" as a filler
          temp.new.ind <- c(empty, temp.new.ind)
          all.new.temp <- cbind(all.new.temp, temp.new.ind)
          temp.m <- temp.new.ind
          respawn.state <- "on"
          changes <- c(changes, j)
          break
        }
        
        
          
        # look for 7-->0 (also indicates respawning)
        if(temp.m[j-1]==7 & temp.m[j]==0)
        {
          temp.new.ind <- temp.m[c(j:length(temp.m))]
          empty <- rep(20, length(temp.m)-length(temp.new.ind)) # assign a "20" as a filler
          temp.new.ind <- c(empty, temp.new.ind)
          all.new.temp <- cbind(all.new.temp, temp.new.ind)
          temp.m <- temp.new.ind
          respawn.state <- "on"
          changes <- c(changes, j)
          break
        }
          
          
        # look for 4-->5 # indicates release of recovered
        if(temp.m[j-1]==4 & temp.m[j]==5)
        {
          temp.new.ind <- temp.m[c(j:length(temp.m))]
          empty <- rep(20, length(temp.m)-length(temp.new.ind)) # assign a "20" as a filler
          temp.new.ind <- c(empty, temp.new.ind)
          release.records <- cbind(release.records, temp.new.ind)
          temp.m <- temp.new.ind
          respawn.state <- "on"
          changes <- c(changes, j)
          break
        }
          
        # look for 7-->5 (also indicates release of recovered)
        if(temp.m[j-1]==7 & temp.m[j]==5)
        {
          temp.new.ind <- temp.m[c(j:length(temp.m))]
          empty <- rep(20, length(temp.m)-length(temp.new.ind)) # assign a "20" as a filler
          temp.new.ind <- c(empty, temp.new.ind)
          release.records <- cbind(release.records, temp.new.ind)
          temp.m <- temp.new.ind
          respawn.state <- "on"
          changes <- c(changes, j)
          break
        }
        
      }

      if(j == nrow(m.orig)) break
    }
    
    # within this column of m.orig (within "i" loop), update column names and extend status of deaths/removals 
    if(respawn.state=="on"){
      
      all.new.temp[c(min(changes):nrow(all.new.temp)),1] <- all.new.temp[(min(changes)-1),1]
      
      if(ncol(all.new.temp)>2){
        for(j in 2:(ncol(all.new.temp)-1)){
          
          all.new.temp[c(changes[j]:nrow(all.new.temp)),j] <- all.new.temp[(changes[j]-1),j]
          
        }
      }
      
    }
    
    for(j in 1:ncol(all.new.temp)){
      colnames(all.new.temp)[j] <- paste(i, j, sep = ".")
    }
    
    if(ncol(release.records)>0){
      for(j in 1:ncol(release.records)){
        colnames(release.records)[j] <- paste(i, j, sep = ".")
      }
    }
    
    
    # store m.new and full release records
    all.new.temp[all.new.temp==20] <- NA
    release.records[release.records==20] <- NA
    
    m.new <- cbind(m.new, all.new.temp)
    all.release.records <- cbind(all.release.records, release.records)
    
    
  } # end of "i" loop through columns
  
  
  # now go back and put "release records" into original territories they were removed from
  
  if(is.null(release.data)==F){

    for(j in 1:nrow(release.data)){
      
      # temporary new.id (represents the release location)
      temp.id <- release.data$new.id[j]
      
      ### start by finding the post-release data for the removed/released individual ###
      # find where new.id match AND where release.time is correct
      right.ind <- data.frame(all.release.records[,which(sapply(strsplit(colnames(all.release.records), "[.]"), `[`, 1)==temp.id)])
      
      right.new <- right.ind[,which(right.ind[release.data$release.time[j]+1,]==5 &
                              is.na(right.ind[release.data$release.time[j],]))]
      # "right.new" is the post-release data matching the j'th observation in "release.data"
      
      
      ### now assign "right.new" to the correct column and time period in m.new (insert post-release data back with the original removed individual) ###
      
      # load capture times to make sure assign release info to correct initial individual
      # original id
      orig.id <- release.data$old.id[j]
      # original capture time for this original id
      old.cands_time <- attempt.time[which(attempted.tr==orig.id & total.tr=="r")]
      
      # which column matches the original ID and the correct capture time
      right.old.index <- which((sapply(strsplit(colnames(m.new), "[.]"), `[`, 1)==orig.id) &
                              m.new[old.cands_time+1,]==7 &  m.new[old.cands_time,]==2)
      
      # for that column, insert the post-release data in the correct location
      m.new[(release.data$release.time[j]+1):nrow(m.new),right.old.index] <- right.new[(release.data$release.time[j]+1):nrow(m.new)]
      
      
      
    }
  
  } # end of processing release records
  

 
  return(m.new)
}