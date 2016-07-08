SimSpatSocMovtV6 <- function(params.in,
                                 N, N.sites, Timesteps, Site.0, P.move.mat, site.dat,
                                 Site.dist, Site.neigh.ind.vec, Site.random.in, Site.neigh.ind){
  
  New.S <- New.X <- New.Y <- New.T <- occd.site.degree.out <- vector("list", dim(params.in)[1])  
  degree.spattemp <- grstrength.spattemp <- grstrength.spat <- matrix(NA, nrow = N * dim(params.in)[1] * Timesteps, ncol = 3)
  degree.spat <- grstrength.spat <- matrix(NA, nrow = N.sites * dim(params.in)[1] * Timesteps, ncol = 3)
  modularity.spat <- modularity.spattemp <- rep(NA, Timesteps)

  meandist.spat <- meandist.spattemp <- rep(NA, Timesteps)
  transitivity.spat <- transitivity.spattemp <- rep(NA, Timesteps)
  power.spat <- power.spattemp <- poisson.spat <- poisson.spattemp <- rep(NA, Timesteps)
  power.spat.loglik <- power.spattemp.loglik <- poisson.spat.loglik <- poisson.spattemp.loglik <- rep(NA, Timesteps)
  
  for(p in 1:dim(params.in)[1]){
    alpha <- params.in[p, 1]
    gamma <- params.in[p, 2]
    delta <- params.in[p, 3]
    rho <- params.in[p, 4]
    
    New.S[[p]] <- New.X[[p]] <- New.Y[[p]] <- New.T[[p]] <- matrix(NA, nrow = N, ncol = Timesteps)
    P.move <- P.move.mat <- Site.random <- vector("list", Timesteps)
    for(t in 1:(Timesteps - 1)){
      if(t == 1){
        # infect one individual at random 
        Site.random[[t]] <- Site.random.in
        for(i in 1:N){
          if(i %in% group.index.animal){
            P.move[[t]] <- alpha * 1 / Site.dist + delta * Site.random[[t]]
            P.move.mat[[t]] <- matrix(P.move[[t]], nrow = N.sites, byrow = T)
            # NEW 02/21/2016: bump up prob of moving to cell in original neighborhood
            # NEW 03/06/2016: instead of additive shift, make shift multiplicative
            P.move.mat[[t]][Site.0[i], ] <- P.move.mat[[t]][Site.0[i], ] + Site.neigh.ind[Site.0[i], ] * gamma
            New.S[[p]][i, t] <- which(rmultinom(1:max(site.dat$Site.no), 1, 
                                                prob = P.move.mat[[t]][Site.0[i], ] / sum(P.move.mat[[t]][Site.0[i], ])) == 1)
            New.X[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.x
            New.Y[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.y
            New.T[[p]][i, t] <- t
          }
          else{
            P.move[[t]] <- alpha * 1 / Site.dist + delta * Site.random[[t]]
            P.move.mat[[t]] <- matrix(P.move[[t]], nrow = N.sites, byrow = T)
            P.move.mat[[t]][Site.0[i], New.S[[p]][NsIndexAnimal[i], t]] <- P.move.mat[[t]][Site.0[i], New.S[[p]][NsIndexAnimal[i], t]] * rho # bump up prob for site chosen by group index animal
            # NEW 02/21/2016: bump up prob of moving to cell in original neighborhood
            # NEW 03/06/2016: instead of additive shift, make shift multiplicative
            P.move.mat[[t]][Site.0[i], ] <- P.move.mat[[t]][Site.0[i], ] + Site.neigh.ind[Site.0[i], ] * gamma
            
            New.S[[p]][i, t] <- which(rmultinom(1:max(site.dat$Site.no), 1, 
                                                prob = P.move.mat[[t]][Site.0[i], ] / sum(P.move.mat[[t]][Site.0[i], ])) == 1)
            New.X[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.x
            New.Y[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.y
            New.T[[p]][i, t] <- t
          }
        }
        New.S.vec <- as.vector(t(New.S[[p]][ , 1:t]))
        New.S.Time <- rep(seq(1:t), times = N)
        New.S.x <- site.dat[New.S.vec, ]$Site.x
        New.S.y <- site.dat[New.S.vec, ]$Site.y
        New.S.individ <- rep(seq(1:N), each = t)
                
        Full.dat.obsInRows <- as.data.frame(cbind(New.S.Time, New.S.x, New.S.y, New.S.individ, New.S.vec))
        names(Full.dat.obsInRows) <- c("Time", "New.X", "New.Y", "ID", "New.S")
        simmed.affin <- CalcSimmedAffinity(data.in = Full.dat.obsInRows)
        simmed.affin.sp <- CalcSimmedAffinitySpat_Sims(data.in = Full.dat.obsInRows, N.sites = N.sites)
        simgraph <- buildGraph(simmed.affin[[2]], data = simmed.affin[[2]])
        simgraph <- delete_edges(simgraph, which(E(simgraph)$weight == 0))
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1] <- degree(simgraph)
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 2] <- rep(t, N)
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 3] <- rep(p, N)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1] <- graph.strength(simgraph)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 2] <- rep(t, N)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 3] <- rep(p, N)    
        modularity.spattemp[t] <- modularity(simgraph, membership(cluster_walktrap(simgraph)))    
        meandist.spattemp[t] <- mean_distance(simgraph, unconnected = F)
        transitivity.spattemp[t] <- transitivity(simgraph, type = "global")
        k <- try(fit_power_law(degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1]))
        if(class(k) != "try-error"){
           power.spattemp[t] <- k$alpha
           power.spattemp.loglik[t] <- k$logLik
        }
        k2 <- try(fitdistr(degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1], "Poisson"))
        if(class(k2) != "try-error"){
          poisson.spattemp[t] <- k2$estimate
          poisson.spattemp.loglik[t] <- k2$loglik
        }
        
        simgraphsp <- buildGraph(simmed.affin.sp[[2]], data = simmed.affin.sp[[2]])
        simgraphsp <- delete_edges(simgraphsp, which(E(simgraphsp)$weight == 0))
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 1] <- degree(simgraphsp)
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 2] <- rep(t, N.sites)
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 3] <- rep(p, N.sites)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 1] <- graph.strength(simgraphsp)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 2] <- rep(t, N.sites)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 3] <- rep(p, N.sites)    
        modularity.spat[t] <- modularity(simgraphsp, membership(cluster_walktrap(simgraphsp)))   
        meandist.spat[t] <- mean_distance(simgraphsp, unconnected = F)
        transitivity.spat[t] <- transitivity(simgraphsp, type = "global")
        q <- try(fit_power_law(degree.spat[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1]))
        if(class(q) != "try-error"){
          power.spat[t] <- q$alpha
          power.spat.loglik[t] <- q$logLik
        }
        q2 <- try(fitdistr(degree.spat[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1], "Poisson"))
        if(class(q2) != "try-error"){
          poisson.spat[t] <- q2$estimate
          poisson.spat.loglik[t] <- q2$loglik
        }

        print(paste(p, "_", t, sep = ""))
      } 
      else { # timesteps after 1
        Site.random[[t]] <- Site.random.in
        for(i in 1:N){
          if(i %in% group.index.animal){
            P.move[[t]] <- alpha * 1 / Site.dist + delta * Site.random[[t]]
            P.move.mat[[t]] <- matrix(P.move[[t]], nrow = N.sites, byrow = T)
            # NEW 02/21/2016: bump up prob of moving to cell in original neighborhood
            # NEW 03/06/2016: instead of additive shift, make shift multiplicative
            # NEW 03/06/2016: need to replace values in CURRENT row (New.S[[p]][i, t-1], not original row (Site.0[[p]][i, ]))
            P.move.mat[[t]][New.S[[p]][i, t - 1], ] <- P.move.mat[[t]][New.S[[p]][i, t - 1], ] + Site.neigh.ind[Site.0[i], ] * gamma
            
            New.S[[p]][i, t] <- which(rmultinom(1:max(site.dat$Site.no), 1, 
                                                prob = P.move.mat[[t]][New.S[[p]][i, t - 1], ] / sum(P.move.mat[[t]][New.S[[p]][i, t - 1], ])) == 1)
            New.X[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.x
            New.Y[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.y
            New.T[[p]][i, t] <- t
          }
          else{
            P.move[[t]] <- alpha * 1 / Site.dist + delta * Site.random[[t]]
            P.move.mat[[t]] <- matrix(P.move[[t]], nrow = N.sites, byrow = T)
            P.move.mat[[t]][New.S[[p]][i, t - 1], New.S[[p]][NsIndexAnimal[i], t]] <- P.move.mat[[t]][New.S[[p]][i, t - 1], New.S[[p]][NsIndexAnimal[i], t]] * rho # bump up prob for site chosen by group index animal
            # NEW 02/21/2016: bump up prob of moving to cell in original neighborhood
            # NEW 03/06/2016: instead of additive shift, make shift multiplicative
            P.move.mat[[t]][New.S[[p]][i, t - 1], ] <- P.move.mat[[t]][New.S[[p]][i, t - 1], ] + Site.neigh.ind[Site.0[i], ] * gamma
            
            New.S[[p]][i, t] <- which(rmultinom(1:max(site.dat$Site.no), 1, 
                                                prob = P.move.mat[[t]][New.S[[p]][i, t - 1], ] / sum(P.move.mat[[t]][New.S[[p]][i, t - 1], ])) == 1)
            New.X[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.x
            New.Y[[p]][i, t] <- site.dat[New.S[[p]][i, t], ]$Site.y
            New.T[[p]][i, t] <- t
          }
        }
        New.S.vec <- as.vector(t(New.S[[p]][ , 1:t]))
        New.S.Time <- rep(seq(1:t), times = N)
        New.S.x <- site.dat[New.S.vec, ]$Site.x
        New.S.y <- site.dat[New.S.vec, ]$Site.y
        New.S.individ <- rep(seq(1:N), each = t)
                
        Full.dat.obsInRows <- as.data.frame(cbind(New.S.Time, New.S.x, New.S.y, New.S.individ, New.S.vec))
        names(Full.dat.obsInRows) <- c("Time", "New.X", "New.Y", "ID", "New.S")
        simmed.affin <- CalcSimmedAffinity(data.in = Full.dat.obsInRows)
        simmed.affin.sp <- CalcSimmedAffinitySpat_Sims(data.in = Full.dat.obsInRows, N.sites = N.sites)
        simgraph <- buildGraph(simmed.affin[[2]], data = simmed.affin[[2]])
        simgraph <- delete_edges(simgraph, which(E(simgraph)$weight == 0))
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1] <- degree(simgraph)
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 2] <- rep(t, N)
        degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 3] <- rep(p, N)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1] <- graph.strength(simgraph)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 2] <- rep(t, N)
        grstrength.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 3] <- rep(p, N)    
        modularity.spattemp[t] <- modularity(simgraph, membership(cluster_walktrap(simgraph)))  
        meandist.spattemp[t] <- mean_distance(simgraph, unconnected = F)
        transitivity.spattemp[t] <- transitivity(simgraph, type = "global")
        k <- try(fit_power_law(degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1]))
        if(class(k) != "try-error"){
          power.spattemp[t] <- k$alpha
          power.spattemp.loglik[t] <- k$logLik
        }
        k2 <- try(fitdistr(degree.spattemp[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1], "Poisson"))
        if(class(k2) != "try-error"){
          poisson.spattemp[t] <- k2$estimate
          poisson.spattemp.loglik[t] <- k2$loglik
        }

        simgraphsp <- buildGraph(simmed.affin.sp[[2]], data = simmed.affin.sp[[2]])
        simgraphsp <- delete_edges(simgraphsp, which(E(simgraphsp)$weight == 0))
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 1] <- degree(simgraphsp)
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 2] <- rep(t, N.sites)
        degree.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 3] <- rep(p, N.sites)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 1] <- graph.strength(simgraphsp)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 2] <- rep(t, N.sites)
        grstrength.spat[((N.sites * (p - 1) * Timesteps) + N.sites * (t - 1) + 1) : (N.sites * (p - 1) * Timesteps + N.sites * t), 3] <- rep(p, N.sites)    
        modularity.spat[t] <- modularity(simgraphsp, membership(cluster_walktrap(simgraphsp)))     
        meandist.spat[t] <- mean_distance(simgraphsp, unconnected = F)
        transitivity.spat[t] <- transitivity(simgraphsp, type = "global")
        q <- try(fit_power_law(degree.spat[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1]))
        if(class(q) != "try-error"){
          power.spat[t] <- q$alpha
          power.spat.loglik[t] <- q$logLik
        }
        q2 <- try(fitdistr(degree.spat[((N * (p - 1) * Timesteps) + N * (t - 1) + 1) : (N * (p - 1) * Timesteps + N * t), 1], "Poisson"))
        if(class(q2) != "try-error"){
          poisson.spat[t] <- q2$estimate
          poisson.spat.loglik[t] <- q2$loglik
        }

        print(paste(p, "_", t, sep = ""))
      }
      # occupied site degree
      occd.site <- occd.site.degree <- matrix(NA, nrow = N, ncol = Timesteps)
      for(i in 1:N){
        for(t in 2:Timesteps){
          occd.site[i, t] <- ifelse(New.S[[p]][i, t] %in% New.S[[p]][-i, seq(1:(t-1))] == T, New.S[[p]][i, t], NA)
          occd.site.degree[i, t] <- length(unique(na.omit(occd.site[i, 1:t])))
        }
      }    
      occd.site.degree.out[[p]] <- cbind(as.vector(occd.site.degree), rep(1:Timesteps, each = N), rep(p, N * Timesteps))
      
    }
  
    # Get individual core areas and HRs over full timeseries
    data.x <- New.X[[p]]
    data.y <- New.Y[[p]]
    spat.pts <- kernel.in <- vector("list", N)
    kernel.95 <- kernel.50 <- core.to.hr <- rep(NA, N)
    # dimensions are Individual X Timestep
    for(j in 1:N){
      spat.pts[[j]] <- SpatialPoints(cbind(data.x[j, 1:(Timesteps-1)], data.y[j, (Timesteps-1)]))
      kernel.in[[j]] <- kernelUD(spat.pts[[j]], h = "href")
      kernel.95[j] <- kernel.area(kernel.in[[j]], percent = 95)
      kernel.50[j] <- kernel.area(kernel.in[[j]], percent = 50)
      core.to.hr[j] <- kernel.50[j] / kernel.95[j]
     }
    print(p)
  }

  out.list <- list(params.in = params.in,
                   New.S = New.S,
                   degree.spattemp = degree.spattemp,
                   degree.spat = degree.spat,
                   grstrength.spattemp = grstrength.spattemp,
                   grstrength.spat = grstrength.spat,
                   New.X = New.X,
                   New.Y = New.Y,
                   modularity.spattemp = modularity.spattemp,
                   modularity.spat = modularity.spat,
                   transitivity.spattemp = transitivity.spattemp,
                   transitivity.spat = transitivity.spat,
                   meandist.spattemp = meandist.spattemp,
                   meandist.spat = meandist.spat,
                   power.spattemp = power.spattemp,
                   power.spattemp.loglik = power.spattemp.loglik,
                   power.spat = power.spat,
                   power.spat.loglik = power.spat.loglik,
                   poisson.spattemp = poisson.spattemp,
                   poisson.spattemp.loglik = poisson.spattemp.loglik,
                   poisson.spat = poisson.spat,
                   poisson.spat.loglik = poisson.spat.loglik,
                   kernel.95 = kernel.95,
                   kernel.50 = kernel.50,
                   core.to.hr = core.to.hr
                  )
  return(out.list)
}