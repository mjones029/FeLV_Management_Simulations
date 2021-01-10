# proactive_vax.R
#
#========================================================	
# ---
### title: Proactive vaccination
# author: Marie Gilbertson
# date: "10/11/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function to randomly distribute proactive vaccination in a network


proactive_vax <- function(g = g,
                          prop.vax = 0.1,
                          single.vax.eff = 0.4,
                          boosted.vax.eff = 0.8,
                          single.prob = 0.5
                          ){

  pop.size <- length(V(g))
  
  # sample which individuals to vaccinate
  vaxed <- sample(pop.size, size = prop.vax*pop.size, replace = F)
  
  # randomly assign vaccine efficacy based on relative proportions of boosted vs unboosted
  # randomly select from vax individuals who will receive a single (unboosted) vaccination
  singles <- sample(vaxed, size = ceiling(length(vaxed)*single.prob), replace = F) # take ceiling to be conservative (vax more w/out booster)
  
  vax.data <- data.frame(vaxed = as.character(vaxed),
                         stringsAsFactors = F)
  
  # all others receive two vaccinations (boosted) 
  vax.data$vax.type <- ifelse(vax.data$vaxed %in% singles, single.vax.eff, boosted.vax.eff)
  
  # pad ids with leading zero(s) for ordering purposes
  vax.data$padded.id <- str_pad(vax.data$vaxed, 3, pad = "0")
  
  vax.data2 <- vax.data[order(vax.data$padded.id),]
  
  
  # assign vaccinated status (status = 6) and vaccine efficacy to vaccinated individuals in network
  vaxed <- as.character(vax.data2$vaxed)
  
  g2 <- g
  g2 <- set_vertex_attr(g2, "status", index = V(g)$vertex.names %in% vaxed, value = 6)  
  g2 <- set_vertex_attr(g2, "status", index = !V(g)$vertex.names %in% vaxed, value = 0) 
  # V(g2)$status
  
  g3 <- g2
  g3 <- set_vertex_attr(g3, "veff", index = V(g3)$vertex.names %in% vaxed, value = vax.data2$vax.type)  
  g3 <- set_vertex_attr(g3, "veff", index = !V(g3)$vertex.names %in% vaxed, value = 0)  
  # V(g3)$veff  
  
  return(g3)
}