# up_edge_id.R
#
#========================================================	
# ---
### title: Underpass edges identification
# author: Marie Gilbertson
# date: "10/31/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function to identify edges that cross the I-75 freeway, and would therefore be 
# removed with underpass closure


up_edge_id <- function(net = new.g) {
  
  # extract node locations
  locations.df <- data.frame(id = V(net)$vertex.names,
                             lat = V(net)$latitude,
                             long = V(net)$longitude,
                             i75.lat = 26.16)
  
  # label nodes as north vs south of I75
  locations.df$i75_NS <- ifelse(locations.df$lat >= locations.df$i75.lat, "N", "S")
  
  locations_n <- locations.df[locations.df$i75_NS=="N",]
  locations_s <- locations.df[locations.df$i75_NS=="S",]
  
  # create edgelist of network
  el <- data.frame(as_edgelist(net))
  
  # add an attribute per edge in edgelist corresponding to whether or not that edge "crosses" I75
  # 1 = edge does NOT cross I75; 2 = edge DOES cross I75
  el$up.weight <- ifelse(((el$X1 %in% locations_n$id) & (el$X2 %in% locations_s$id)) | ((el$X1 %in% locations_s$id) & (el$X2 %in% locations_n$id)), 2, 1)
  
  # now set underpass "weight" as a new edge attribute in "net"
  net2 <- set_edge_attr(net, "up.weight", index = E(net), value = el$up.weight)
  
  
  return(net2)
}

