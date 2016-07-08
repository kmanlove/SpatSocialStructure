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
power.spat.out <- power.spattemp.out <- poisson.spat.out <- poisson.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
power.spat.ll.out <- power.spattemp.ll.out <- poisson.spat.ll.out <- poisson.spattemp.ll.out <- matrix(NA, nrow = length(indices), ncol = 100)
data.in <- vector("list", length(indices))
modularity.spat.out <- modularity.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
transitivity.spat.out <- transitivity.spattemp.out <- meandist.spat.out <- meandist.spattemp.out <- matrix(NA, nrow = length(indices), ncol = 100)
core.to.hr.out <- kernel.50.out <- kernel.95.out <- matrix(NA, nrow = length(indices), ncol = 50)
groupsize.lambda <- rep(NA, length(indices))
group.size.vec <- vector("list", length(indices))

for(i in 1:length(indices)){
  k <- try(load(paste("./Data/Simulations_20160629/output_", indices[i], ".RData", sep = "")))
  if(class(k) != "try-error"){
    data.in[[i]] <- sim.out
    modularity.spat.out[i, ] <- data.in[[i]]$modularity.spat
    modularity.spattemp.out[i, ] <- data.in[[i]]$modularity.spattemp
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

mean.core.to.hr <- apply(core.to.hr.out, 1, mean)
mean.kernel.50 <- apply(kernel.50.out, 1, mean)
mean.kernel.95 <- apply(kernel.95.out, 1, mean)


# build data frame 
sim.data <- as.data.frame(cbind(params.in.full, power.spat.out[, 99], power.spattemp.out[, 99], 
                                poisson.spat.out[, 99], poisson.spattemp.out[, 99], 
                                modularity.spat.out[, 99], modularity.spattemp.out[, 99], 
                                meandist.spat.out[, 99], meandist.spattemp.out[, 99], 
                                transitivity.spat.out[,99], transitivity.spattemp.out[, 99], groupsize.lambda,
                                mean.core.to.hr, mean.kernel.50, mean.kernel.95))
names(sim.data) <- c("alpha", "gamma", "delta", "rho", "power.spat", "power.spattemp", "poisson.spat", "poisson.spattemp", 
                     "modularity.spat", "modularity.spattemp", "meandist.spat", "meandist.spattemp", "transitivity.spat", "transitivity.spattemp",
                     "groupsize.lambda", "mean.core.to.hr", "mean.kernel.50", "mean.kernel.95")

# SPATIAL DEGREE
require(asbio)
pois.spat.lm.fit.sat <- lm(poisson.spat ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spat.lm.fit.nogammadelta <- lm(poisson.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spat.lm.fit.nogammarho <- lm(poisson.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
pois.spat.lm.fit.nodeltarho <- lm(poisson.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
pois.spat.lm.fit.nogamma <- lm(poisson.spat ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spat.lm.fit.nodelta <- lm(poisson.spat ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spat.lm.fit.norho <- lm(poisson.spat ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(pois.spat.lm.fit.nogammadelta, pois.spat.lm.fit.sat)
partial.R2(pois.spat.lm.fit.nogammarho, pois.spat.lm.fit.sat)
partial.R2(pois.spat.lm.fit.nodeltarho, pois.spat.lm.fit.sat)
partial.R2(pois.spat.lm.fit.nogamma, pois.spat.lm.fit.sat)
partial.R2(pois.spat.lm.fit.nodelta, pois.spat.lm.fit.sat)
partial.R2(pois.spat.lm.fit.norho, pois.spat.lm.fit.sat)

# SPATTEMP DEGREE
pois.spattemp.lm.fit.sat <- lm(poisson.spattemp ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spattemp.lm.fit.nogammadelta <- lm(poisson.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spattemp.lm.fit.nogammarho <- lm(poisson.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
pois.spattemp.lm.fit.nodeltarho <- lm(poisson.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
pois.spattemp.lm.fit.nogamma <- lm(poisson.spattemp ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spattemp.lm.fit.nodelta <- lm(poisson.spattemp ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
pois.spattemp.lm.fit.norho <- lm(poisson.spattemp ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(pois.spattemp.lm.fit.nogamma, pois.spattemp.lm.fit.sat)
partial.R2(pois.spattemp.lm.fit.nodelta, pois.spattemp.lm.fit.sat)
partial.R2(pois.spattemp.lm.fit.norho, pois.spattemp.lm.fit.sat)
partial.R2(pois.spattemp.lm.fit.nogammadelta, pois.spattemp.lm.fit.sat)
partial.R2(pois.spattemp.lm.fit.nodeltarho, pois.spattemp.lm.fit.sat)
partial.R2(pois.spattemp.lm.fit.nogammarho, pois.spattemp.lm.fit.sat)

# SPAT MODULARITY
mod.spat.lm.fit.sat <- lm(modularity.spat ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spat.lm.fit.nogammadelta <- lm(modularity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spat.lm.fit.nogammarho <- lm(modularity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
mod.spat.lm.fit.nodeltarho <- lm(modularity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
mod.spat.lm.fit.nogamma <- lm(modularity.spat ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spat.lm.fit.nodelta <- lm(modularity.spat ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spat.lm.fit.norho <- lm(modularity.spat ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(mod.spat.lm.fit.nogamma, mod.spat.lm.fit.sat)
partial.R2(mod.spat.lm.fit.nodelta, mod.spat.lm.fit.sat)
partial.R2(mod.spat.lm.fit.norho, mod.spat.lm.fit.sat)
partial.R2(mod.spat.lm.fit.nogammadelta, mod.spat.lm.fit.sat)
partial.R2(mod.spat.lm.fit.nodeltarho, mod.spat.lm.fit.sat)
partial.R2(mod.spat.lm.fit.nogammarho, mod.spat.lm.fit.sat)

# SPATTEMP MODULARITY
mod.spattemp.lm.fit.sat <- lm(modularity.spattemp ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spattemp.lm.fit.nogammadelta <- lm(modularity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spattemp.lm.fit.nogammarho <- lm(modularity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
mod.spattemp.lm.fit.nodeltarho <- lm(modularity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
mod.spattemp.lm.fit.nogamma <- lm(modularity.spattemp ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spattemp.lm.fit.nodelta <- lm(modularity.spattemp ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
mod.spattemp.lm.fit.norho <- lm(modularity.spattemp ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(mod.spattemp.lm.fit.nogamma, mod.spattemp.lm.fit.sat)
partial.R2(mod.spattemp.lm.fit.nodelta, mod.spattemp.lm.fit.sat)
partial.R2(mod.spattemp.lm.fit.norho, mod.spattemp.lm.fit.sat)
partial.R2(mod.spattemp.lm.fit.nogammadelta, mod.spattemp.lm.fit.sat)
partial.R2(mod.spattemp.lm.fit.nodeltarho, mod.spattemp.lm.fit.sat)
partial.R2(mod.spattemp.lm.fit.nogammarho, mod.spattemp.lm.fit.sat)

# SPAT TRANSITIVITY
trans.spat.lm.fit.sat <- lm(transitivity.spat ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spat.lm.fit.nogammadelta <- lm(transitivity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spat.lm.fit.nogammarho <- lm(transitivity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
trans.spat.lm.fit.nodeltarho <- lm(transitivity.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
trans.spat.lm.fit.nogamma <- lm(transitivity.spat ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spat.lm.fit.nodelta <- lm(transitivity.spat ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spat.lm.fit.norho <- lm(transitivity.spat ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(trans.spat.lm.fit.nogamma, trans.spat.lm.fit.sat)
partial.R2(trans.spat.lm.fit.nodelta, trans.spat.lm.fit.sat)
partial.R2(trans.spat.lm.fit.norho, trans.spat.lm.fit.sat)
partial.R2(trans.spat.lm.fit.nogammadelta, trans.spat.lm.fit.sat)
partial.R2(trans.spat.lm.fit.nodeltarho, trans.spat.lm.fit.sat)
partial.R2(trans.spat.lm.fit.nogammarho, trans.spat.lm.fit.sat)

# SPATTEMP TRANSITIVITY
trans.spattemp.lm.fit.sat <- lm(transitivity.spattemp ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spattemp.lm.fit.nogammadelta <- lm(transitivity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spattemp.lm.fit.nogammarho <- lm(transitivity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
trans.spattemp.lm.fit.nodeltarho <- lm(transitivity.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
trans.spattemp.lm.fit.nogamma <- lm(transitivity.spattemp ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spattemp.lm.fit.nodelta <- lm(transitivity.spattemp ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
trans.spattemp.lm.fit.norho <- lm(transitivity.spattemp ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(trans.spattemp.lm.fit.nogamma, trans.spattemp.lm.fit.sat)
partial.R2(trans.spattemp.lm.fit.nodelta, trans.spattemp.lm.fit.sat)
partial.R2(trans.spattemp.lm.fit.norho, trans.spattemp.lm.fit.sat)
partial.R2(trans.spattemp.lm.fit.nogammadelta, trans.spattemp.lm.fit.sat)
partial.R2(trans.spattemp.lm.fit.nodeltarho, trans.spattemp.lm.fit.sat)
partial.R2(trans.spattemp.lm.fit.nogammarho, trans.spattemp.lm.fit.sat)


# SPAT MEAN DISTANCE
dist.spat.lm.fit.sat <- lm(meandist.spat ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spat.lm.fit.nogammadelta <- lm(meandist.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spat.lm.fit.nogammarho <- lm(meandist.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
dist.spat.lm.fit.nodeltarho <- lm(meandist.spat ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
dist.spat.lm.fit.nogamma <- lm(meandist.spat ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spat.lm.fit.nodelta <- lm(meandist.spat ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spat.lm.fit.norho <- lm(meandist.spat ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(dist.spat.lm.fit.nogamma, dist.spat.lm.fit.sat)
partial.R2(dist.spat.lm.fit.nodelta, dist.spat.lm.fit.sat)
partial.R2(dist.spat.lm.fit.norho, dist.spat.lm.fit.sat)
partial.R2(dist.spat.lm.fit.nogammadelta, dist.spat.lm.fit.sat)
partial.R2(dist.spat.lm.fit.nodeltarho, dist.spat.lm.fit.sat)
partial.R2(dist.spat.lm.fit.nogammarho, dist.spat.lm.fit.sat)

# SPATTEMP MEAN DISTANCE
dist.spattemp.lm.fit.sat <- lm(meandist.spattemp ~ log(gamma) + log(delta) +  log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spattemp.lm.fit.nogammadelta <- lm(meandist.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spattemp.lm.fit.nogammarho <- lm(meandist.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho), data = sim.data)
dist.spattemp.lm.fit.nodeltarho <- lm(meandist.spattemp ~ log(gamma) + log(delta) + log(rho) + log(delta):log(gamma) + log(gamma):log(rho), data = sim.data)
dist.spattemp.lm.fit.nogamma <- lm(meandist.spattemp ~ log(delta) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spattemp.lm.fit.nodelta <- lm(meandist.spattemp ~ log(gamma) + log(rho) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
dist.spattemp.lm.fit.norho <- lm(meandist.spattemp ~ log(gamma) + log(delta) + log(delta):log(gamma) + log(delta):log(rho) + log(gamma):log(rho), data = sim.data)
partial.R2(dist.spattemp.lm.fit.nogamma, dist.spattemp.lm.fit.sat)
partial.R2(dist.spattemp.lm.fit.nodelta, dist.spattemp.lm.fit.sat)
partial.R2(dist.spattemp.lm.fit.norho, dist.spattemp.lm.fit.sat)
partial.R2(dist.spattemp.lm.fit.nogammadelta, dist.spattemp.lm.fit.sat)
partial.R2(dist.spattemp.lm.fit.nodeltarho, dist.spattemp.lm.fit.sat)
partial.R2(dist.spattemp.lm.fit.nogammarho, dist.spattemp.lm.fit.sat)


#---- HOW WELL DO GROUPSIZE/CORE AREA:HR AREA PREDICT METRICS? ----#
#-- SPATIAL DEGREE --#
pois.spat.lm.sat <- lm(poisson.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.nocake <- lm(poisson.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.nogske <- lm(poisson.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.nogsca <- lm(poisson.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.noke <- lm(poisson.spat ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.noca <- lm(poisson.spat ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spat.lm.gs.ca.ker.nogs <- lm(poisson.spat ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(pois.spat.lm.gs.ca.ker.nogs, pois.spat.lm.sat)
partial.R2(pois.spat.lm.gs.ca.ker.noca, pois.spat.lm.sat)
partial.R2(pois.spat.lm.gs.ca.ker.noke, pois.spat.lm.sat)
partial.R2(pois.spat.lm.gs.ca.ker.nogsca, pois.spat.lm.sat)
partial.R2(pois.spat.lm.gs.ca.ker.nogske, pois.spat.lm.sat)
partial.R2(pois.spat.lm.gs.ca.ker.nocake, pois.spat.lm.sat)

# SPATIOTEMPORAL DEGREE
pois.spattemp.lm.sat <- lm(poisson.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.nocake <- lm(poisson.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.nogske <- lm(poisson.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.nogsca <- lm(poisson.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.noke <- lm(poisson.spattemp ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.noca <- lm(poisson.spattemp ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
pois.spattemp.lm.gs.ca.ker.nogs <- lm(poisson.spattemp ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(pois.spattemp.lm.gs.ca.ker.nogs, pois.spattemp.lm.sat)
partial.R2(pois.spattemp.lm.gs.ca.ker.noca, pois.spattemp.lm.sat)
partial.R2(pois.spattemp.lm.gs.ca.ker.noke, pois.spattemp.lm.sat)
partial.R2(pois.spattemp.lm.gs.ca.ker.nogsca, pois.spattemp.lm.sat)
partial.R2(pois.spattemp.lm.gs.ca.ker.nogske, pois.spattemp.lm.sat)
partial.R2(pois.spattemp.lm.gs.ca.ker.nocake, pois.spattemp.lm.sat)

# SPATIAL MODULARITY
mod.spat.lm.sat <- lm(modularity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.nocake <- lm(modularity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.nogske <- lm(modularity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.nogsca <- lm(modularity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.noke <- lm(modularity.spat ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.noca <- lm(modularity.spat ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spat.lm.gs.ca.ker.nogs <- lm(modularity.spat ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(mod.spat.lm.gs.ca.ker.nogs, mod.spat.lm.sat)
partial.R2(mod.spat.lm.gs.ca.ker.noca, mod.spat.lm.sat)
partial.R2(mod.spat.lm.gs.ca.ker.noke, mod.spat.lm.sat)
partial.R2(mod.spat.lm.gs.ca.ker.nogsca, mod.spat.lm.sat)
partial.R2(mod.spat.lm.gs.ca.ker.nogske, mod.spat.lm.sat)
partial.R2(mod.spat.lm.gs.ca.ker.nocake, mod.spat.lm.sat)

# SPATIOTEMPORAL MODULARITY
mod.spattemp.lm.sat <- lm(modularity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.nocake <- lm(modularity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.nogske <- lm(modularity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.nogsca <- lm(modularity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.noke <- lm(modularity.spattemp ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.noca <- lm(modularity.spattemp ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
mod.spattemp.lm.gs.ca.ker.nogs <- lm(modularity.spattemp ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(mod.spattemp.lm.gs.ca.ker.nogs, mod.spattemp.lm.sat)
partial.R2(mod.spattemp.lm.gs.ca.ker.noca, mod.spattemp.lm.sat)
partial.R2(mod.spattemp.lm.gs.ca.ker.noke, mod.spattemp.lm.sat)
partial.R2(mod.spattemp.lm.gs.ca.ker.nogsca, mod.spattemp.lm.sat)
partial.R2(mod.spattemp.lm.gs.ca.ker.nogske, mod.spattemp.lm.sat)
partial.R2(mod.spattemp.lm.gs.ca.ker.nocake, mod.spattemp.lm.sat)

# SPATIAL TRANSITIVITY
trans.spat.lm.sat <- lm(transitivity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.nocake <- lm(transitivity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.nogske <- lm(transitivity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.nogsca <- lm(transitivity.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.noke <- lm(transitivity.spat ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.noca <- lm(transitivity.spat ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spat.lm.gs.ca.ker.nogs <- lm(transitivity.spat ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(trans.spat.lm.gs.ca.ker.nogs, trans.spat.lm.sat)
partial.R2(trans.spat.lm.gs.ca.ker.noca, trans.spat.lm.sat)
partial.R2(trans.spat.lm.gs.ca.ker.noke, trans.spat.lm.sat)
partial.R2(trans.spat.lm.gs.ca.ker.nogsca, trans.spat.lm.sat)
partial.R2(trans.spat.lm.gs.ca.ker.nogske, trans.spat.lm.sat)
partial.R2(trans.spat.lm.gs.ca.ker.nocake, trans.spat.lm.sat)

# SPATIOTEMPORAL TRANSITIVITY
trans.spattemp.lm.sat <- lm(transitivity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.nocake <- lm(transitivity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.nogske <- lm(transitivity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.nogsca <- lm(transitivity.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.noke <- lm(transitivity.spattemp ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.noca <- lm(transitivity.spattemp ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
trans.spattemp.lm.gs.ca.ker.nogs <- lm(transitivity.spattemp ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(trans.spattemp.lm.gs.ca.ker.nogs, trans.spattemp.lm.sat)
partial.R2(trans.spattemp.lm.gs.ca.ker.noca, trans.spattemp.lm.sat)
partial.R2(trans.spattemp.lm.gs.ca.ker.noke, trans.spattemp.lm.sat)
partial.R2(trans.spattemp.lm.gs.ca.ker.nogsca, trans.spattemp.lm.sat)
partial.R2(trans.spattemp.lm.gs.ca.ker.nogske, trans.spattemp.lm.sat)
partial.R2(trans.spattemp.lm.gs.ca.ker.nocake, trans.spattemp.lm.sat)

# SPATIAL MEAN DISTANCE
dist.spat.lm.sat <- lm(meandist.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.nocake <- lm(meandist.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.nogske <- lm(meandist.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.nogsca <- lm(meandist.spat ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.noke <- lm(meandist.spat ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.noca <- lm(meandist.spat ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spat.lm.gs.ca.ker.nogs <- lm(meandist.spat ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(dist.spat.lm.gs.ca.ker.nogs, dist.spat.lm.sat)
partial.R2(dist.spat.lm.gs.ca.ker.noca, dist.spat.lm.sat)
partial.R2(dist.spat.lm.gs.ca.ker.noke, dist.spat.lm.sat)
partial.R2(dist.spat.lm.gs.ca.ker.nogsca, dist.spat.lm.sat)
partial.R2(dist.spat.lm.gs.ca.ker.nogske, dist.spat.lm.sat)
partial.R2(dist.spat.lm.gs.ca.ker.nocake, dist.spat.lm.sat)

# SPATIOTEMPORAL MEAN DISTANCE
dist.spattemp.lm.sat <- lm(meandist.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.nocake <- lm(meandist.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.nogske <- lm(meandist.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.nogsca <- lm(meandist.spattemp ~ groupsize.lambda + mean.core.to.hr + mean.kernel.50 + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.noke <- lm(meandist.spattemp ~ groupsize.lambda + mean.core.to.hr + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.noca <- lm(meandist.spattemp ~ groupsize.lambda + mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
dist.spattemp.lm.gs.ca.ker.nogs <- lm(meandist.spattemp ~ mean.kernel.50 + groupsize.lambda:mean.core.to.hr + groupsize.lambda:mean.kernel.50 + mean.core.to.hr:mean.kernel.50, data = sim.data)
partial.R2(dist.spattemp.lm.gs.ca.ker.nogs, dist.spattemp.lm.sat)
partial.R2(dist.spattemp.lm.gs.ca.ker.noca, dist.spattemp.lm.sat)
partial.R2(dist.spattemp.lm.gs.ca.ker.noke, dist.spattemp.lm.sat)
partial.R2(dist.spattemp.lm.gs.ca.ker.nogsca, dist.spattemp.lm.sat)
partial.R2(dist.spattemp.lm.gs.ca.ker.nogske, dist.spattemp.lm.sat)
partial.R2(dist.spattemp.lm.gs.ca.ker.nocake, dist.spattemp.lm.sat)

## CCA APPROACH 
require(vegan)
