library('move')
library('maps')
library('foreach')
library('lubridate')

#does not work for lon=lat=0 (errors)

rFunction <- function(data, duration=NULL, radius=NULL)
{
  Sys.setenv(tz="UTC")
  
  n.all <- length(timestamps(data))
  data <- data[!duplicated(paste0(round_date(timestamps(data), "5 min"), trackId(data))),]
  logger.info(paste0("For better performance, the data have been thinned to max 5 minute resolution. From the total ",n.all," positions, the algorithm retained ",length(timestamps(data))," positions for calculation."))

  if (is.null(duration) & is.null(radius)) 
  {
    logger.info("You didnt provide any stopover site radius or minimum stopover duration. Please go back and configure them. Here return input data set.") 
    result <- data
  }
  if (is.null(duration) & !is.null(radius)) 
    {
    logger.info(paste0("You have selected a stopover site radius of ",radius,"m, but no minimum stopover duration. We here use 1 day by default. If that is not what you need, please go back and configure the parameters."))
    duration <- 24*3600
    }
  if (!is.null(duration) & is.null(radius))
    {
    logger.info(paste0("You have selected a minimum stopover duration of ",duration,"h, but no radius. We here use 10'000 m = 10 km by default. If that is not what you need, please go back and configure the parameters."))
    radius <- 10000
  }
  if (!is.null(duration) & !is.null(radius))
  {
    logger.info(paste0("You have selected a minimum stopover duration of ",duration," seconds and a radius of ",radius," metres."))
    
    ll2xyz <- function(ll) {
      ## Convert lat/long coordinates (in radians) to 3d Cartesian (on unit sphere)
      c( cos(ll[1])*c(cos(ll[2]),sin(ll[2])) , sin(ll[1]) )
    }
    
    xyz2ll <- function(xyz) {
      ## Convert 3d Cartesian coordinates to lat/long (in radians)
      ## Also computes angles for coordinates not on the unit sphere
      c( atan2(xyz[3], sqrt(sum(xyz[1:2]^2))), atan2(xyz[2],xyz[1]) )
    }
    
    sphereCircle <- function(lat, lon, r, rSphere=6371008.8, n=360) {
      ## Generate a circle centered at (lat,lon) of radius r on the sphere or radius rSphere
      ## Returns matrix of lat/long coordinates for n equally spaced points on circle
      
      ## First compute distance from c to central axis, normalized to unit circle
      r <- sin(r/rSphere)
      
      ## Generate n points equidistant around the North Pole at radius r
      a <- (2*pi/n)*(0:n)
      P <- matrix(c(r*cos(a),r*sin(a),rep(sqrt(1-r^2), n+1)), ncol=3)
      
      ## Compute rotation matrix so North Pole ends up at (lat,lon)
      l <- -pi*(0.5-(lat/180)) # to radian and subtract from 90 degree == pi/2
      Ry <- matrix(c(cos(l), 0, sin(l), 0, 1, 0, -sin(l), 0, cos(l)), nrow=3)
      l <- lon/180*pi
      Rz <- matrix(c(cos(l), sin(l), 0, -sin(l), cos(l), 0, 0, 0, 1), nrow=3)
      R <- Rz %*% Ry
      ## Rotate points to lie around new center
      P <- t(apply(P, 1, function(p) {
        R %*% p
      }))
      
      ## Unproject to lat/long
      P <- t(apply(P, 1, xyz2ll) * (180/pi))
      
      ## Find the break in longitude if any and reorder
      m <- which.max(P[,2])
      if (m != nrow(P)) {
        P <- P[c((m+1):nrow(P),1:m),]
      }
      P
    }
    
    dist.ct <- function(p,q) {
      sqrt(sum((p-q)^2))
    }
    dist <- dist.ct
    
    dist.ll <- function(p,q) {
      ## Great-circle distance between p,q in lat/long coordinates (in radians)
      ## Normalized to unit sphere
      
      # Calculates the geodesic distance using the Haversine formula
      d <- q-p
      a <- sin(d[1]/2)^2 + cos(p[1]) * cos(q[1]) * sin(d[2]/2)^2
      2 * asin(min(1,sqrt(a)))
    }
    
    midpoint.ct <- function(p,q) {
      ## Midpoint of p and q in Cartesian coordinates (arbitrary d)
      (p+q)*0.5
    }
    midpoint <- midpoint.ct
    
    midpoint.ll <- function(p,q) {
      ## Midpoint of shortest great-circle arc between p,q (lat/lon in radians)
      ## Undefined for antipodal points
      
      ## Convert to Cartesian, compute midpoint, convert back
      xyz2ll( midpoint.ct( ll2xyz(p), ll2xyz(q) ) )
    }
    
    circumcircle.ct <- function(p,q,r) {
      ## Compute circumcircle of triangle embedded in higher dimensional space
      ## (p,q,r in Cartesian coordinates, d arbitrary)
      ## Based on https://en.wikipedia.org/wiki/Circumscribed_circle#Higher_dimensions
      l  <- function(x) { sqrt(sum(x^2)) } # shortcut for vector length
      l2 <- function(x) { sum(x^2) } # squared vector length
      
      ## Translate p to origin
      q <- q-p
      r <- r-p
      
      # Compute translated circumcenter
      qxr2 <- l2(q)*l2(r) - sum(q*r)^2 # ||q cross r||^2
      d <- (l2(q)*r - l2(r)*q)
      m <- (sum(d*r)*q - sum(d*q)*r) / (2*qxr2)
      # Compute radius
      r <- l(q) * l(r) * l(q-r) / (2*sqrt(qxr2))
      
      # translate back circumcenter and return
      c(m+p, r)
    }
    circumcircle <- circumcircle.ct
    
    circumcircle.ll <- function(p, q, r) {
      ## Midpoint of smallest circumcircle of p,q,r (lat/long in radians)
      ## Undefined for p,q,r and origin coplanar
      
      ## Convert to Cartesian, compute circumcircle, convert back
      m <- xyz2ll(circumcircle.ct( ll2xyz(p), ll2xyz(q), ll2xyz(r) )[1:3])
      c(m, dist.ll(m,p)) # Recompute radius along great circle
    }
    
    miniDisc <- function(P, rMax = Inf, rSphere = Inf) {
      ## Compute smallest enclosing disc, stop if radius becomes > rMax
      ## P is a point set (without temporal coordinates)
      ## Points can have either cartesian (d=2) or spherical coordinates
      ## All points must be contained in a hemisphere for spherical coordinates
      ## rSphere is the radius of the sphere or Inf for cartesian coordinates
      ## Earth's mean radius is about 6371 km
      
      ## Based on the algorithm described in Section 4.7 of
      ## De Berg, Cheong, Van Kreveld, Overmars: Computational Geometry,
      ## Algorithms and Applications, 3rd ed. Springer, 2008.
      
      if (nrow(P) <= 1) {
        if (nrow(P) == 0) {
          return(NULL)
        }
        disc <- c(P, 0)
      } else {
        miniDiscWithPoints <- function (i, Q) {
          ## Compute SED of P[1:i,] under condition that all points in Q are on boundary
          ## Implemented with this signature to avoid copying P unnecessarily
          
          # Set up initial disc
          init.points <- rbind(Q, P[1:2,])[1:2,] # Get 2 points from either Q or P
          disc <- c(midpoint(init.points[1,], init.points[2,]), 
                    dist(init.points[1,], init.points[2,])/2)
          if (disc[length(disc)] > rMax) {
            return(NULL)
          }
          
          # Iteratively add points
          for (j in (2-nrow(Q)+seq_len(i-2+nrow(Q)))) {
            ## Test if P[j,] is inside disc
            if (dist(P[j,],disc[-length(disc)]) > disc[length(disc)]) {
              ## Recursive call requiring P[j,] on boundary
              if (nrow(Q) >= 2) {
                disc <- circumcircle(Q[1,], Q[2,], P[j,])
                if (disc[length(disc)] > rMax) {
                  return(NULL)
                }
              } else {
                disc <- miniDiscWithPoints(j-1, rbind(Q, P[j,]))
                if (is.null(disc)) { return(disc) }
              }
            }
          }
          disc
        } ## /miniDiscWithPoints()
        
        if (is.finite(rSphere)) {
          # Convert to unit sphere for simplicity of calculations
          rMax <- rMax / rSphere
          # Convert ll coordinates to radians
          P <- P * (pi / 180)
          
          ## switch to dist, midpoint, circumcircle functions for spherical coordinates
          dist <- dist.ll
          midpoint <- midpoint.ll
          circumcircle <- circumcircle.ll
        }
        
        P <- unique(P) # Duplicate points cause trouble; filter them
        # if only one unique location, then same disc as if only one location at the beginning
        if (nrow(P)==1) disc <- c(P, 0) else
        {
          # Randomly permute points
          # Start with subtrajectory endpoints as a heuristic
          perm <- c(1,nrow(P),1+sample(seq_len(nrow(P)-2),size=nrow(P)-2,replace=FALSE))
          P <- P[perm,]
          
          # Compute disc
          disc <- miniDiscWithPoints(i=nrow(P), Q=matrix(double(0), ncol=ncol(P)))
        }
      }
      
      # Clean result and return
      if (!is.null(disc)) {
        if (is.infinite(rSphere)) { 
          names(disc) <- c("x", "y", "r") 
        } else {
          disc[1:2] <- disc[1:2] * (180 / pi) # Convert lat/lon to degrees
          disc[3] <- disc[3] * rSphere # Convert radius from unit sphere
          names(disc) <- c("lat", "long", "r") 
        }
      }
      disc
    }
    
    stopOvers <- function(tr, tMin, rMax) {
      ## Compute stopovers in tr (Move object).
      ## tMin: min stopover duration (s), rMax: max stopover radius (m)
      ## returns matrix with one stopover per row
      ## Columns are 
      ## - The first and last indices of observations contained in the stopover
      ## - Duration in seconds
      ## - Coordinates of disc center
      ## - Radius of the disc (in meters for lat/long coordinates)
      
      if (length(grep("+proj=longlat", tr@proj4string, value = FALSE)) >= 1) {
        ## Lat/long coordinates
        dist <- dist.ll
        rSphere <- 6371008.8 # Mean Earth radius
        rm <- rMax / rSphere # Normalize to unit sphere
      } else {
        ## Projected coordinates, assume same units as rMax
        rSphere <- Inf # Cartesian coordinates
        rm <- rMax # No normalization necessary
      }
      tMin <- as.difftime(tMin, units="secs")
      
      ts <- tr@timestamps
      tc <- tr@coords[,2:1] # Switch to lat/long format
      
      stopOvers <- data.frame(iStart=integer(),
                              iEnd=integer(),
                              duration=integer(),
                              cLat=double(),
                              cLong=double(),
                              radius=double())
      ## Scan over trajectory; test whether tMin-length window exists with endpoints closer than rMax
      start <- end <- 1
      while (end < nrow(tr)) {
        ## Find first point more than tMin time away from start
        while (ts[end] - ts[start] < tMin && end < nrow(tr)) {
          end <- end+1
        }
        
        if (ts[end] - ts[start] >= tMin && dist(tc[start,], tc[end,]) <= rm) {
          ## Potential stopover, compute SED
          disc <- miniDisc(tc[start:end,], rMax, rSphere)
          if (!is.null(disc)) {
            ## Stopover detected
            ## Maximize duration, use exponential and binary search
            while(!is.null(disc) && end < nrow(tc)) { # Exponential search...
              so <- c(start,end, as.double(ts[end]-ts[start], units="secs"), disc)
              end <- min(end + (end-start), nrow(tc)) # ~Double number of points
              disc <- miniDisc(tc[start:end,], rMax, rSphere)
            }
            if (!is.null(disc)) { ## Stopover lasts until end of tr
              so <- c(start, end, as.double(ts[end]-ts[start], units="secs"), disc)
            } else {
              while(end>so[2]) { # ... Binary search
                m <- ceiling((so[2] + end)/2)
                disc <- miniDisc(tc[start:m,], rMax, rSphere)
                if (is.null(disc)) {
                  end <- m-1
                } else {
                  ## Found new lower bound, store longer stopover
                  so <- c(start, m, as.double(ts[m]-ts[start], units="secs"), disc)
                }
              }
            }
            
            stopOvers[nrow(stopOvers)+1,] <- so
            ## Find disjoint stopovers only
            start <- end <- so[2]
          }
        }
        
        start <- start + 1
      }
      
      if (is.infinite(rSphere)) {
        colnames(stopOvers) <- c("iStart", "iEnd", "duration", "cX", "cY", "radius")
      } else {
        colnames(stopOvers) <- c("iStart", "iEnd", "duration", "cLat", "cLong", "radius")
      }
      rownames(stopOvers) <- NULL
      
      stopOvers
    }
    
    data.split <- move::split(data)
    stopover.tab <- data.frame("individual_local_identifier"=character(),"timestamp_arrival"=character(),"timestamps_departure"=character(),"location_long"=numeric(),"location_lat"=numeric(),"duration"=numeric(),"radius"=numeric(),"taxon_canonical_name"=character(),"sensor"=character())
    
    foreach(datai = data.split) %do% {
      res <- stopOvers(datai,duration,radius)
      ARR <- timestamps(datai)[res$iStart]
      DEP <- timestamps(datai)[res$iEnd]
      ID <- rep(namesIndiv(datai),length(res$iStart))
      DUR <- res$duration #in seconds
      
      nr <- length(DUR)
      tax <- rep(idData(datai)$taxon_canonical_name,nr)
      senso <- rep(as.character(sensor(datai)[1]),nr)
      
      resi <- data.frame("individual_local_identifier"=ID,"timestamp_arrival"=ARR,"timestamp_departure"=DEP,"location_long"=res$cLong,"location_lat"=res$cLat,"duration"=DUR,"radius"=res$radius,"taxon_canonical_name"=tax,"sensor"=senso)
      stopover.tab <- rbind(stopover.tab,resi)
    }  
     
    stopover.tab.list <- as.list(data.frame(t(stopover.tab)))
     
    stopover.sites <- foreach(stopoveri = stopover.tab.list) %do% {
      dataj <- data.split[[which(names(data.split)==stopoveri[1])]]
      dataj[timestamps(dataj)>=as.POSIXct(stopoveri[2]) & timestamps(dataj)<=as.POSIXct(stopoveri[3])]
    }
    names(stopover.sites) <- names(stopover.tab.list)      

    stopover.sites.nozero <- stopover.sites[unlist(lapply(stopover.sites, length) > 0)] #remove IDs with no data
    
    if (length(stopover.sites.nozero)==0) 
    {
      logger.info("Your output file contains no positions/stopover sites. No csv saved. Return NULL.")
      result <- NULL
    } else 
    {
      result <- moveStack(stopover.sites.nozero)
      
      names(stopover.tab)[4:7] <- c("mid_longitude","mid_latitude","duration (s)","radius (m)")
      
      write.csv(stopover.tab,file = paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"stopover_sites.csv"),row.names=FALSE) #csv artefakt
      #write.csv(stopover.tab,file = "stopover_sites.csv",row.names=FALSE)
    }
  }
  
  return(result)
}

  
  
  
  
  
  
  
  
  
  
