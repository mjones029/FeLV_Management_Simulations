## simulate_pop_ergm
# 
#========================================================	
# ---
### title: Test simulating from ERGM
# author: Marie Gilbertson
# date: "03/02/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulate population based on ERGM coefficients
# Adapted from simulate_ergm function in FIV-->FeLV project


simulate_pop_ergm <- function(telemetry = tele.ordered,       # telemetry data to use for simulating HR centroids
                              net.den = net.den,           # network density to constrain simulated network to
                              a_prop = a_prop,                # approximate proportion of population that should be categorized as adults (vs juveniles; ignore kittens) 
                              ergm.mod = combo.mod,           # ergm model results to use
                              pop.size = pop.size,                  # size of simulated population
                              net.style = "bestFIV"              # "bestFIV" - determines which ergm to use to generate network; just "bestFIV" right now, but could incorporate other ERGM results
){
  
  # must detach igraph package in order to set vertex attributes later on (also have to make sure tnet is detached or igraph won't detach)
  if( "tnet" %in% (.packages())) {
    detach("package:tnet", unload=TRUE) 
  }
  
  if( "igraph" %in% (.packages())) {
    detach("package:igraph", unload=TRUE) 
  }
  
  #### simulate population representative of contemporary FL panthers ####
  
  #### simulate pairwise distances of simulated population ####
  ### load location data ###
  tele.ordered <- telemetry
  
  
  ### calculate MCP home ranges to calculate HR centroids ###
  # Note: currently doing this every time to generate as there is some randomization in sampling for MCP
  
  cats <- unique(tele.ordered$CATNUMBER)
  
  tele.mcp_30.50 <- NULL # to go into MCP analysis
  tele.mcp_1 <- NULL # save in place of MCP centroid
  
  for(i in 1:length(cats)){
    #print(i)
    temp.cat <- as.character(cats[i])
    temp.tele <- subset(tele.ordered, tele.ordered$CATNUMBER==temp.cat)
    
    if(nrow(temp.tele)>=50){
      
      sample.locs <- sample_n(temp.tele, size = 50, replace = F)
      sample.locs <- sample.locs[order(sample.locs$date.time),]
      tele.mcp_30.50 <- rbind(tele.mcp_30.50, sample.locs)
      
    }else if(nrow(temp.tele)<50 & nrow(temp.tele)>=30){
      
      sample.locs <- temp.tele
      tele.mcp_30.50 <- rbind(tele.mcp_30.50, sample.locs)
      
    }else if(nrow(temp.tele)<30){
      
      sample.locs <- sample_n(temp.tele, size = 1, replace = F)
      tele.mcp_1 <- rbind(tele.mcp_1, sample.locs)
      
    }
  }
  
  
  id <- as.character(tele.mcp_30.50$CATNUMBER)
  xyt <- as.matrix(tele.mcp_30.50[,c(6,7)]) # references UTM83EAST and UTM83NORTH, respectively
  idsp <- data.frame(id)
  coordinates(idsp) <- xyt                 
  class(idsp)
  
  crs.idsp <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  proj4string(idsp) <- crs.idsp
  is.projected(idsp)
  
  
  polys <- mcp(idsp, percent = 95)
  
  trueCentroids = gCentroid(polys,byid=TRUE)
  
  centroids <- data.frame(trueCentroids)
  colnames(centroids) <- c("UTM_East", "UTM_North")
  centroids$tele.id <- rownames(centroids)
  centroids <- centroids[,c(3,1,2)] # just change the order of the columns so tele.id comes first
  
  # add back in the random locations for individuals with fewer than 30 locations
  sing.centroids <- tele.mcp_1[,c("CATNUMBER", "UTM83EAST", "UTM83NORTH")]
  colnames(sing.centroids) <- c("tele.id", "UTM_East", "UTM_North")
  row.names(sing.centroids) <- sing.centroids$tele.id
  
  centroids <- rbind(centroids, sing.centroids)
  
  
  
  
  # convert UTM locations to lat/long
  sputm <- SpatialPoints(cbind(centroids$UTM_East, centroids$UTM_North), proj4string=CRS("+proj=utm +zone=17N +datum=WGS84"))
  spgeo <- data.frame(spTransform(sputm, CRS("+proj=longlat +datum=WGS84")))
  colnames(spgeo) <- c("Longitude", "Latitude")
  
  centroids <- cbind(centroids, spgeo)
  
  
  
  ####### simulate new centroid points using spatial quadrats ######
  
  # make a polygon window based on centroid locations (rather than a rectangular window)
  xyt2 <- as.matrix(centroids[,c("UTM_East","UTM_North")]) # references UTM83EAST and UTM83NORTH, respectively
  idsp2 <- data.frame(matrix("A", nrow = nrow(centroids), ncol = 1)) # assign arbitrary id, as want a single polygon encompassing these points
  coordinates(idsp2) <- xyt2                 
  class(idsp2)
  
  crs.idsp <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  proj4string(idsp2) <- crs.idsp
  is.projected(idsp2)
  
  
  
  poly.window <- mcp(idsp2, percent = 100)
  poly.coords <- data.frame(poly.window@polygons[[1]]@Polygons[[1]]@coords)
  rev.poly.coords_x <- rev(poly.coords$UTM_East) # have to reverse the order of points to work in "owin" function
  rev.poly.coords_y <- rev(poly.coords$UTM_North)
  rev.poly.coords <- data.frame(x = rev.poly.coords_x,
                                y = rev.poly.coords_y)
  
  tele.window <- owin(xrange=c(min(centroids$UTM_East), max(centroids$UTM_East)), yrange = c(min(centroids$UTM_North), max(centroids$UTM_North)), 
                      poly=list(x=rev.poly.coords$x, y=rev.poly.coords$y))
  
  
  
  # convert FeLV period centroids into point pattern object for spatstat
  ss.centroids <- ppp(x = centroids$UTM_East, y = centroids$UTM_North, window = tele.window)
  # plot(ss.centroids)
  
  # calculate quadrats, counts, and proportions
  Q <- quadratcount(ss.centroids, nx=7, ny=10) # started with 4x7; increased to 7x10 to increase clustering (based on comparison to observed centroids)
  Q.prop <- Q/nrow(centroids)
  
  
  quads <- as.tess(Q)
  t.quads <- tiles(quads)
  
  # loop through each tile extracting the polygon coordinates for each tile
  all.sample.coords <- NULL
  pop.size <- pop.size
  
  for(i in 1:length(Q.prop)){
    #print(i)
    
    if(t.quads[[i]]$type=="polygonal"){
      tq.1x <- t.quads[[i]]$bdry[[1]]$x
      tq.1y <- t.quads[[i]]$bdry[[1]]$y
    }else if(t.quads[[i]]$type=="rectangle"){
      tq.1x <- c(rep(t.quads[[i]]$xrange[1], 2), rep(t.quads[[i]]$xrange[2], 2))
      tq.1y <- c(t.quads[[i]]$yrange, t.quads[[i]]$yrange[2], t.quads[[i]]$yrange[1])
    }
    
    
    # convert those points back into polygons for sampling purposes
    x_coord <- tq.1x
    y_coord <- tq.1y
    xym <- cbind(x_coord, y_coord)
    
    # triangles result in a polygon warning; simply append first coordinate of triangle to become the last as well
    if(nrow(xym)==3){
      xym <- rbind(xym, xym[1,])
    }
    
    p = Polygon(xym)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    # plot(sps)
    
    if(Q.prop[i]!=0){
      # sample points based on quadrat counts/proportions from earlier
      temp.sample <- spsample(sps, n = ceiling(Q.prop[i]*pop.size), type = "random", iter = 10)
      # spsample will sometimes error out (likely because of small polygons); increased iterations from default 4 to 10
      
      # extract coordinates from Spatial Points
      temp.coords <- data.frame(coordinates(temp.sample))
      colnames(temp.coords) <- c("x", "y")
      temp.coords$tile <- paste(i)
      
      # store for next round
      all.sample.coords <- rbind(all.sample.coords, temp.coords)
    }
  }
  
  # will generally create too large of a population due to rounding; sample down to proper pop.size
  # using "ceiling" prevents too *small* of a population which would be more difficult to resolve
  if(nrow(all.sample.coords)>pop.size){
    all.sample.coords <- sample_n(all.sample.coords, size = pop.size, replace = F)
  }
  
  # plot(ss.centroids)
  # points(all.sample.coords$x, all.sample.coords$y, col = "blue")
  
  
  ### calculate pairwise distance matrix with quadrat-generated centroids ###
  
  sputm2 <- SpatialPoints(cbind(all.sample.coords$x, all.sample.coords$y), proj4string=CRS("+proj=utm +zone=17N +datum=WGS84"))
  spgeo2 <- data.frame(spTransform(sputm2, CRS("+proj=longlat +datum=WGS84")))
  colnames(spgeo2) <- c("Longitude", "Latitude")
  
  gen.locs <- as.matrix(spgeo2)
  
  # final distance matrices
  dist.mat2 <- distm(gen.locs, fun = distGeo)
  # covert to log(km) as used in original ERGM fitting
  dist.mat2_km <- dist.mat2/1000
  dist.mat2_log <- log(dist.mat2_km)
  dist.mat2_log[dist.mat2_log=="-Inf"] <- 0
  
  
  #### simulate categorical ages of simulated population ####
  
  # based on results from "ideal" parameter sets in FIV/FeLV project
  age.data <- data.frame(matrix(nrow = nrow(dist.mat2_log), ncol = 2))
  colnames(age.data) <- c("id", "age.cat")
  age.data$id <- seq(1, nrow(dist.mat2_log))
  rand.ages <- rbinom(nrow(dist.mat2_log), size = 1, prob = a_prop)
  age.data$age.cat <- ifelse(rand.ages==1, "Adult", "Yearling")
  
  
  
  # assign consistent naming 
  colnames(dist.mat2_log) <- row.names(dist.mat2_log) <- age.data$id
  
  
  
  
  #### simulate network ####
  
  # initialize x to a random undirected network with 80 nodes (min panther count in 2002)
  x <- network(nrow(dist.mat2_log), directed = F, density = net.den)
  
  x <- set.network.attribute(x, "distances", dist.mat2_log)
  x <- set.vertex.attribute(x, "age.cat", age.data$age.cat, v=seq_len(network.size(x)))
  x <- set.vertex.attribute(x, "latitude", spgeo2$Latitude, v=seq_len(network.size(x)))
  x <- set.vertex.attribute(x, "longitude", spgeo2$Longitude, v=seq_len(network.size(x)))
  
  # get combo.mod from ERGM script
  combo.mod <- ergm.mod
  
  
  # based on function input, generate network based on specific ERGM model
  
  if(net.style=="bestFIV"){
    coefs <-  c("edges", "gwesp.fixed.0.7", "altkstar.0.5", "nodefactor.age.cat.Adult", "edgecov.pair.dist") # ergm coefficients to use in simulating new network
    sim.formula <-  x ~ edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base=2) + edgecov("distances")  # ergm formula to use in simulating new network - must be same items and order as coefficients
    constrain <-  ~edges
    
  }else{
    print("WARNING: Undefined ERGM")
  }
  
  sim.net <- simulate(sim.formula,
                      coef = combo.mod$coef[coefs], verbose = F, constraints = constrain)
  # plot(sim.net)
  
  return(sim.net)
  
}