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

for(i in 1:length(indices)){
  k <- try(load(paste("./Data/Simulations_20160629/output_", indices[i], ".RData", sep = "")))
  if(class(k) != "try-error"){
    data.in[[i]] <- sim.out
    modularity.spat.out[i, ] <- data.in[[i]]$modularity.spat
    modularity.spattemp.out[i, ] <- data.in[[i]]$modularity.spattemp
    grstrength.spat.out[i, ] <- data.in[[i]]$grstrength.spat
    grstrength.spattemp.out[i, ] <- data.in[[i]]$grstrength.spattemp
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


#----------------#
#-- Modularity --#
#----------------#

delta.01 <- subset(sim.data, delta == .01)
delta.1 <- subset(sim.data, delta == .1)
delta1 <- subset(sim.data, delta == 1)
delta10 <- subset(sim.data, delta == 10)
# check.data <- cbind(delta.01$modularity.spat, delta.01$gamma, delta.01$rho)
mod.sp.delt.01.mat <- matrix(delta.01$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt.01.mat <- matrix(delta.1$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt.1.mat <- matrix(delta.1$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt.1.mat <- matrix(delta.01$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt1.mat <- matrix(delta1$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt1.mat <- matrix(delta1$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt10.mat <- matrix(delta10$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt10.mat <- matrix(delta10$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.lims <- range(na.omit(c(mod.sp.delt.01.mat, mod.sptemp.delt.01.mat,
                            mod.sp.delt.1.mat, mod.sptemp.delt.1.mat,
                            mod.sp.delt1.mat, mod.sptemp.delt1.mat,
                            mod.sp.delt10.mat, mod.sptemp.delt10.mat)))
par(mfrow = c(4, 2), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(mod.sp.delt.01.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Spatial modularity, delta = .01")
image.plot(mod.sptemp.delt.01.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Social modularity, delta = .01")
image.plot(mod.sp.delt.1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(mod.sptemp.delt.1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(mod.sp.delt1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(mod.sptemp.delt1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(mod.sp.delt10.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(mod.sptemp.delt10.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)

#--------------------#
#-- Poisson degree --#
#--------------------#
pois.sp.delt.01.mat <- matrix(delta.01$poisson.spat, nrow = 20, ncol = 20, byrow = F)
pois.sptemp.delt.01.mat <- matrix(delta.1$poisson.spattemp, nrow = 20, ncol = 20, byrow = F)
pois.sp.delt.1.mat <- matrix(delta.1$poisson.spat, nrow = 20, ncol = 20, byrow = F)
pois.sptemp.delt.1.mat <- matrix(delta.01$poisson.spattemp, nrow = 20, ncol = 20, byrow = F)
pois.sp.delt1.mat <- matrix(delta1$poisson.spat, nrow = 20, ncol = 20, byrow = F)
pois.sptemp.delt1.mat <- matrix(delta1$poisson.spattemp, nrow = 20, ncol = 20, byrow = F)
pois.sp.delt10.mat <- matrix(delta10$poisson.spat, nrow = 20, ncol = 20, byrow = F)
pois.sptemp.delt10.mat <- matrix(delta10$poisson.spattemp, nrow = 20, ncol = 20, byrow = F)
pois.lims <- range(na.omit(c(pois.sp.delt.01.mat, pois.sptemp.delt.01.mat,
                            pois.sp.delt.1.mat, pois.sptemp.delt.1.mat,
                            pois.sp.delt1.mat, pois.sptemp.delt1.mat,
                            pois.sp.delt10.mat, pois.sptemp.delt10.mat)))
par(mfrow = c(4, 2), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(pois.sp.delt.01.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Lambda for spatial degree, delta = .01")
image.plot(pois.sptemp.delt.01.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Lambda for social degree, delta = .01")
image.plot(pois.sp.delt.1.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(pois.sptemp.delt.1.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(pois.sp.delt1.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(pois.sptemp.delt1.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(pois.sp.delt10.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(pois.sptemp.delt10.mat, zlim = pois.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)



#-------------------------#
# -- Spat mod - soc mod --#
#-------------------------#

moddiff.sp.delt.01.mat <- matrix((delta.01$modularity.spat - delta.01$modularity.spattemp), nrow = 20, ncol = 20, byrow = F)
moddiff.sp.delt.1.mat <- matrix((delta.1$modularity.spat - delta.1$modularity.spattemp), nrow = 20, ncol = 20, byrow = F)
moddiff.sp.delt1.mat <- matrix((delta1$modularity.spat - delta1$modularity.spattemp), nrow = 20, ncol = 20, byrow = F)
moddiff.sp.delt10.mat <- matrix((delta10$modularity.spat - delta10$modularity.spattemp), nrow = 20, ncol = 20, byrow = F)
moddiff.lims <- range(na.omit(c(moddiff.sp.delt.01.mat,
                            moddiff.sp.delt.1.mat,
                            moddiff.sp.delt1.mat,
                            moddiff.sp.delt10.mat)))

poisdiff.sp.delt.01.mat <- matrix((delta.01$poisson.spat - delta.01$poisson.spattemp), nrow = 20, ncol = 20, byrow = F)
poisdiff.sp.delt.1.mat <- matrix((delta.1$poisson.spat - delta.1$poisson.spattemp), nrow = 20, ncol = 20, byrow = F)
poisdiff.sp.delt1.mat <- matrix((delta1$poisson.spat - delta1$poisson.spattemp), nrow = 20, ncol = 20, byrow = F)
poisdiff.sp.delt10.mat <- matrix((delta10$poisson.spat - delta10$poisson.spattemp), nrow = 20, ncol = 20, byrow = F)
poisdiff.lims <- range(na.omit(c(poisdiff.sp.delt.01.mat,
                                poisdiff.sp.delt.1.mat,
                                poisdiff.sp.delt1.mat,
                                poisdiff.sp.delt10.mat)))
par(mfrow = c(2,4), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(poisdiff.sp.delt.01.mat, zlim = poisdiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Spat degree - Soc degree, delta = .01")
image.plot(poisdiff.sp.delt.1.mat, zlim = poisdiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(poisdiff.sp.delt1.mat, zlim = poisdiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(poisdiff.sp.delt10.mat, zlim = poisdiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(moddiff.sp.delt.01.mat, zlim = moddiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Spat mod - Soc mod, delta = .01")
image.plot(moddiff.sp.delt.1.mat, zlim = moddiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(moddiff.sp.delt1.mat, zlim = moddiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(moddiff.sp.delt10.mat, zlim = moddiff.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)


#----------------#
#-- Group size --#
#----------------#
gs.delt.01.mat <- matrix(delta.01$groupsize.lambda, nrow = 20, ncol = 20, byrow = F)
gs.delt.1.mat <- matrix(delta.1$groupsize.lambda, nrow = 20, ncol = 20, byrow = F)
gs.delt1.mat <- matrix(delta1$groupsize.lambda, nrow = 20, ncol = 20, byrow = F)
gs.delt10.mat <- matrix(delta10$groupsize.lambda, nrow = 20, ncol = 20, byrow = F)
gs.lims <- range(na.omit(c(gs.delt.01.mat, 
                             gs.delt.1.mat, 
                             gs.delt1.mat, 
                             gs.delt10.mat)))
par(mfrow = c(4, 1), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(gs.delt.01.mat, zlim = gs.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Groupsize lambda, delta = .01")
image.plot(gs.delt.1.mat, zlim = gs.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(gs.delt1.mat, zlim = gs.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(gs.delt10.mat, zlim = gs.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")

#---------------#
#-- Space use --#
#---------------#
core.to.hr.delt.01.mat <- matrix(delta.01$mean.core.to.hr, nrow = 20, ncol = 20, byrow = F)
kernel.95.delt.01.mat <- matrix(delta.1$mean.kernel.95, nrow = 20, ncol = 20, byrow = F)
core.to.hr.delt.1.mat <- matrix(delta.1$mean.core.to.hr, nrow = 20, ncol = 20, byrow = F)
kernel.95.delt.1.mat <- matrix(delta.01$mean.kernel.95, nrow = 20, ncol = 20, byrow = F)
core.to.hr.delt1.mat <- matrix(delta1$mean.core.to.hr, nrow = 20, ncol = 20, byrow = F)
kernel.95.delt1.mat <- matrix(delta1$mean.kernel.95, nrow = 20, ncol = 20, byrow = F)
core.to.hr.delt10.mat <- matrix(delta10$mean.core.to.hr, nrow = 20, ncol = 20, byrow = F)
kernel.95.delt10.mat <- matrix(delta10$mean.kernel.95, nrow = 20, ncol = 20, byrow = F)
c.to.hr.lims <- range(na.omit(c(core.to.hr.delt.01.mat, 
                             core.to.hr.delt.1.mat, 
                             core.to.hr.delt1.mat, 
                             core.to.hr.delt10.mat)))
kernel.95.lims <- range(na.omit(c(kernel.95.delt.01.mat, 
                                  kernel.95.delt.1.mat, 
                                  kernel.95.delt1.mat, 
                                  kernel.95.delt10.mat)))
par(mfrow = c(4, 2), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 2))
image.plot(core.to.hr.delt.01.mat, zlim = c.to.hr.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Core area:HR, delta = .01")
image.plot(kernel.95.delt.01.mat, zlim = kernel.95.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "95% kernel area, delta = .01")
image.plot(core.to.hr.delt.1.mat, zlim = c.to.hr.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(kernel.95.delt.1.mat, zlim = kernel.95.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(core.to.hr.delt1.mat, zlim = c.to.hr.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(kernel.95.delt1.mat, zlim = kernel.95.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(core.to.hr.delt10.mat, zlim = c.to.hr.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(kernel.95.delt10.mat, zlim = kernel.95.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")

mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)



mod.sp.delt.01.mat <- matrix(delta.01$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt.01.mat <- matrix(delta.1$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt.1.mat <- matrix(delta.1$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt.1.mat <- matrix(delta.01$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt1.mat <- matrix(delta1$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt1.mat <- matrix(delta1$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.sp.delt10.mat <- matrix(delta10$modularity.spat, nrow = 20, ncol = 20, byrow = F)
mod.sptemp.delt10.mat <- matrix(delta10$modularity.spattemp, nrow = 20, ncol = 20, byrow = F)
mod.lims <- range(na.omit(c(mod.sp.delt.01.mat, mod.sptemp.delt.01.mat,
                            mod.sp.delt.1.mat, mod.sptemp.delt.1.mat,
                            mod.sp.delt1.mat, mod.sptemp.delt1.mat,
                            mod.sp.delt10.mat, mod.sptemp.delt10.mat)))
par(mfrow = c(4, 2), mar = c(5, 5, 2, 3), oma = c(0, 1, 1, 0))
image.plot(mod.sp.delt.01.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Spatial modularity, delta = .01")
image.plot(mod.sptemp.delt.01.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "Social modularity, delta = .01")
image.plot(mod.sp.delt.1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = .1")
image.plot(mod.sptemp.delt.1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(mod.sp.delt1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 1")
image.plot(mod.sptemp.delt1.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
image.plot(mod.sp.delt10.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity",
           main = "delta = 10")
image.plot(mod.sptemp.delt10.mat, zlim = mod.lims,
           xlab = "Home range affinity", 
           ylab = "Group affinity")
mtext(text = "Habitat heterogeneity INCREASES DOWN rows", side = 2, line = 2, outer = F, las = 0)

#--------------------------#
#-- Degree by core-to-hr --#
#--------------------------#
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in <- c(.5, 1)
ylim.in <- c(0, 100)
ylim.in.spattemp <- c(0, 50)
plot(delta.01$poisson.spat ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$poisson.spat ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$poisson.spat ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$poisson.spat ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$poisson.spattemp ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta.1$poisson.spattemp ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta1$poisson.spattemp ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta10$poisson.spattemp ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)

#-------------------------#
#-- Degree by groupsize --#
#-------------------------#
# SHOULD PROBABLY COLOR THIS ONE BY RHO, NOT GAMMA (CORE-TO-HR COLORED BY GAMMA, NOT RHO)
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in <- c(1, 6)
ylim.in <- c(0, 100)
ylim.in.spattemp <- c(0, 50)
plot(delta.01$poisson.spat ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$poisson.spat ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$poisson.spat ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$poisson.spat ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$poisson.spattemp ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta.1$poisson.spattemp ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta1$poisson.spattemp ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta10$poisson.spattemp ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Social Degree",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)


#------------------------------#
#-- Modularity by core-to-hr --#
#------------------------------#
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in <- c(.5, 1)
ylim.in <- c(-0.008, 1)
ylim.in.spattemp <- c(-0.008, 1)
plot(delta.01$modularity.spat ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$modularity.spat ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$modularity.spat ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$modularity.spat ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$modularity.spattemp ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta.1$modularity.spattemp ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta1$modularity.spattemp ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta10$modularity.spattemp ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)


#------------------------------------#
#-- Modularity by groupsize lambda --#
#------------------------------------#
# AGAIN, CONSIDER RECODING COLOR SO THAT IT RAMPS ON RHO, NOT GAMMA.
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in <- c(1, 6)
ylim.in <- c(-0.008, 1)
ylim.in.spattemp <- c(-0.008, 1)
plot(delta.01$modularity.spat ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$modularity.spat ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$modularity.spat ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$modularity.spat ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$modularity.spattemp ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta.1$modularity.spattemp ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta1$modularity.spattemp ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)
plot(delta10$modularity.spattemp ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Social modularity",
     xlim = xlim.in,
     ylim = ylim.in.spattemp)


#----------------------------------------------#
#-- Degree and Mod by Group size (Spat only) --#
#----------------------------------------------#
# SHOULD PROBABLY COLOR THIS ONE BY RHO, NOT GAMMA (CORE-TO-HR COLORED BY GAMMA, NOT RHO)
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in.mod <- c(1, 6)
ylim.in <- c(0, 100)
ylim.in.mod <- c(-0.008, 1)
plot(delta.01$poisson.spat ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$poisson.spat ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for spatial degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$poisson.spat ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for spatial degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$poisson.spat ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Lambda for spatial degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$modularity.spat ~ delta.01$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.01$rho))], 
     cex = log10(delta.01$gamma * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod)
plot(delta.1$modularity.spat ~ delta.1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta.1$rho))], 
     cex = log10(delta.1$gamma * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod)
plot(delta1$modularity.spat ~ delta1$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta1$rho))], 
     cex = log10(delta1$gamma * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod)
plot(delta10$modularity.spat ~ delta10$groupsize.lambda, 
     col = col.in[as.numeric(factor(delta10$rho))], 
     cex = log10(delta10$gamma * 10^5)/2,
     xlab = "Groupsize lambda", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod)


#----------------------------------#
#-- Degree and Mod by core-to-hr --#
#----------------------------------#
col.in <- rainbow(n = 20)
par(mfrow = c(2, 4), las = 1)
xlim.in <- c(.5, 1)
ylim.in <- c(0, 100)
ylim.in.mod <- c(-0.008, 1)
plot(delta.01$poisson.spat ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .01")
plot(delta.1$poisson.spat ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = .1")
plot(delta1$poisson.spat ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 1")
plot(delta10$poisson.spat ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Lambda for Spatial Degree",
     xlim = xlim.in,
     ylim = ylim.in,
     main = "delta = 10")
plot(delta.01$modularity.spat ~ delta.01$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.01$gamma))], 
     cex = log10(delta.01$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod,
     main = "delta = .01")
plot(delta.1$modularity.spat ~ delta.1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta.1$gamma))], 
     cex = log10(delta.1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod,
     main = "delta = .1")
plot(delta1$modularity.spat ~ delta1$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta1$gamma))], 
     cex = log10(delta1$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod,
     main = "delta = 1")
plot(delta10$modularity.spat ~ delta10$mean.core.to.hr, 
     col = col.in[as.numeric(factor(delta10$gamma))], 
     cex = log10(delta10$rho * 10^5)/2,
     xlab = "Mean core area:HR", 
     ylab = "Spatial modularity",
     xlim = xlim.in,
     ylim = ylim.in.mod,
     main = "delta = 10")
