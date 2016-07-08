# each sim.out object contains one replicate of the FULL parameter space
# store all reps from each parameter combo in their own object

# SPATTEMP DEGREE PLOTTING
alpha.in <- delta.in <- gamma.in <- rho.in <- c(1, 100)
alpha.in <- gamma.in <- .1
delta.in <- c(.01, .1, 1, 10, 100)
gamma.in <- exp(seq(-5, 5, length.out = 20))
rho.in <- exp(seq(-10, 10, length.out = 20))

params.in.full <- expand.grid(alpha.in, gamma.in, delta.in, rho.in)

#-------------------------------#
#-- Read in simulation output --#
#-------------------------------#
delta.low <- which(params.in.full[, 3] <= .01)
params.in.deltlow <- params.in.full[delta.low, ]
delta.high <- which(params.in.full[, 3] >= 1)
params.in.delthigh <- params.in.full[delta.high, ]

require(MASS)
require(igraph)
require(fields)

indices <- seq(1:2000)
power.spat.out <- power.spattemp.out <- poisson.spat.out <- poisson.spattemp.out <- grstrength.spat.out <- grstrength.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
power.spat.ll.out <- power.spattemp.ll.out <- poisson.spat.ll.out <- poisson.spattemp.ll.out <- matrix(NA, nrow = length(indices), ncol = 100)
data.in <- vector("list", length(indices))
modularity.spat.out <- modularity.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
transitivity.spat.out <- transitivity.spattemp.out <- meandist.spat.out <- meandist.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
core.to.hr.out <- kernel.50.out <- kernel.95.out <- matrix(NA, nrow = length(indices), ncol = 50)
groupsize.lambda <- rep(NA, length(indices))
group.size.vec <- vector("list", length(indices))
degree.spat.out <- matrix(NA, nrow = length(indices), ncol = 10000)
degree.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 5000)

for(i in 1:length(indices)){
  k <- try(load(paste("./Data/Simulations_20160629/output_", indices[i], ".RData", sep = "")))
  if(class(k) != "try-error"){
    data.in[[i]] <- sim.out
    degree.spat.out[i, ] <- data.in[[i]]$degree.spat[, 1]
    degree.spattemp.out[i, ] <- data.in[[i]]$degree.spattemp[, 1]
    modularity.spat.out[i, ] <- data.in[[i]]$modularity.spat
    modularity.spattemp.out[i, ] <- data.in[[i]]$modularity.spattemp
#    grstrength.spat.out[i, ] <- data.in[[i]]$grstrength.spat
#    grstrength.spattemp.out[i, ] <- data.in[[i]]$grstrength.spattemp
    power.spattemp.out[i,] <- data.in[[i]]$power.spattemp
    poisson.spattemp.out[i,] <- data.in[[i]]$poisson.spattemp
    power.spat.out[i, ] <- data.in[[i]]$power.spat
    poisson.spat.out[i, ] <- data.in[[i]]$poisson.spat
    transitivity.spat.out[i, ] <- data.in[[i]]$transitivity.spat
    transitivity.spattemp.out[i, ] <- data.in[[i]]$transitivity.spattemp
    meandist.spat.out[i, ] <- data.in[[i]]$meandist.spat
    meandist.spattemp.out[i, ] <- data.in[[i]]$meandist.spattemp
    kernel.50.out[i, ] <- data.in[[i]]$kernel.50
    kernel.95.out[i, ] <- data.in[[i]]$kernel.95
    core.to.hr.out[i, ] <- data.in[[i]]$core.to.hr
    group.size.list <- vector("list", dim(sim.out$New.S[[1]])[1])
    for(j in 1:dim(sim.out$New.S[[1]])[1]){
      group.size.list[[j]] <- table(sim.out$New.S[[1]][, j])
    }
    group.size.vec[[i]] <- do.call("c", group.size.list)
    groupsize.lambda[i] <- fitdistr(group.size.vec[[i]], "Poisson")$estimate
  }
  print(i)
}

# Extract degree, determine number of edges, build Erdos-Renyi graph analog for each graph
# get mean path length and clustering coefficient for each E-R graph match
# get ratio of observed transitivity:ER transivitiy, observed mean dist:ER mean dist

tot.degrees.sp <- tot.edges.sp <- rep(NA, length(indices))
ER.analog.sp <- vector("list", length(indices))
ER.trans.sp <- ER.meandist.sp <- small.world.ness.sp <- rep(NA, length(indices))
tot.degrees.soc <- tot.edges.soc <- rep(NA, length(indices))
ER.analog.soc <- vector("list", length(indices))
ER.trans.soc <- ER.meandist.soc <- small.world.ness.soc <- rep(NA, length(indices))
for(i in 1:length(indices)){
  tot.degrees.sp[i] <- sum(degree.spat.out[i, 9800:9899])
  tot.edges.sp[i] <- tot.degrees.sp[i] / 2
  ER.analog.sp[[i]] <- sample_gnm(100, tot.edges.sp, directed = F, loops = F)
  ER.trans.sp[i] <- transitivity(ER.analog.sp[[i]], type = "global")
  ER.meandist.sp[i] <- average.path.length(ER.analog.sp[[i]])
  small.world.ness.sp[i] <- (transitivity.spat.out[i, 99] / ER.trans.sp[i]) / (meandist.spat.out[i, 99] / ER.meandist.sp[i])

  tot.degrees.soc[i] <- sum(degree.spattemp.out[i, 4900:4949])
  tot.edges.soc[i] <- tot.degrees.soc[i] / 2
  ER.analog.soc[[i]] <- sample_gnm(50, tot.edges.soc, directed = F, loops = F)
  ER.trans.soc[i] <- transitivity(ER.analog.soc[[i]], type = "global")
  ER.meandist.soc[i] <- average.path.length(ER.analog.soc[[i]])
  small.world.ness.soc[i] <- (transitivity.spattemp.out[i, 99] / ER.trans.soc[i]) / (meandist.spattemp.out[i, 99] / ER.meandist.soc[i])
  
  print(i)
}


mean.core.to.hr <- apply(core.to.hr.out, 1, mean)
mean.kernel.50 <- apply(kernel.50.out, 1, mean)
mean.kernel.95 <- apply(kernel.95.out, 1, mean)


# build data frame 
sim.data <- as.data.frame(cbind(params.in.full, power.spat.out[, 99], power.spattemp.out[, 99], 
                                poisson.spat.out[, 99], poisson.spattemp.out[, 99], 
                                modularity.spat.out[, 99], modularity.spattemp.out[, 99], 
                                meandist.spat.out[, 99], meandist.spattemp.out[, 99], 
                                transitivity.spat.out[,99], transitivity.spattemp.out[, 99], groupsize.lambda,
                                mean.core.to.hr, mean.kernel.50, mean.kernel.95, small.world.ness.sp, small.world.ness.soc))

names(sim.data) <- c("alpha", "gamma", "delta", "rho", "power.spat", "power.spattemp", "poisson.spat", "poisson.spattemp", 
                     "modularity.spat", "modularity.spattemp", "meandist.spat", "meandist.spattemp", "transitivity.spat", "transitivity.spattemp",
                     "groupsize.lambda", "mean.core.to.hr", "mean.kernel.50", "mean.kernel.95", "small.world.ness.sp", "small.world.ness.soc")


#----------------------#
#-- Small-world-ness --#
#----------------------#

delta.01 <- subset(sim.data, delta == .01)
delta.1 <- subset(sim.data, delta == .1)
delta1 <- subset(sim.data, delta == 1)
delta10 <- subset(sim.data, delta == 10)
# check.data <- cbind(delta.01$modularity.spat, delta.01$gamma, delta.01$rho)
sw.sp.delt.01.mat <- matrix(delta.01$small.world.ness.sp, nrow = 20, ncol = 20, byrow = F)
sw.sptemp.delt.01.mat <- matrix(delta.1$small.world.ness.soc, nrow = 20, ncol = 20, byrow = F)
sw.sp.delt.1.mat <- matrix(delta.1$small.world.ness.sp, nrow = 20, ncol = 20, byrow = F)
sw.sptemp.delt.1.mat <- matrix(delta.01$small.world.ness.soc, nrow = 20, ncol = 20, byrow = F)
sw.sp.delt1.mat <- matrix(delta1$small.world.ness.sp, nrow = 20, ncol = 20, byrow = F)
sw.sptemp.delt1.mat <- matrix(delta1$small.world.ness.soc, nrow = 20, ncol = 20, byrow = F)
sw.sp.delt10.mat <- matrix(delta10$small.world.ness.sp, nrow = 20, ncol = 20, byrow = F)
sw.sptemp.delt10.mat <- matrix(delta10$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
sw.lims <- range(na.omit(c(sw.sp.delt.01.mat, sw.sptemp.delt.01.mat,
                           sw.sp.delt.1.mat, sw.sptemp.delt.1.mat,
                           sw.sp.delt1.mat, sw.sptemp.delt1.mat,
                           sw.sp.delt10.mat, sw.sptemp.delt10.mat)))
par(mfrow = c(2, 4), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(sw.sptemp.delt.01.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Social small-world-ness, delta = .01")
image.plot(sw.sptemp.delt.1.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(sw.sptemp.delt1.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(sw.sptemp.delt10.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(sw.sp.delt.01.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Spatial small-world-ness, delta = .01")
image.plot(sw.sp.delt.1.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(sw.sp.delt1.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(sw.sp.delt10.mat, zlim = sw.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)
