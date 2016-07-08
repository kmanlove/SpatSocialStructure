#-- specify path to source file --#
root.path <- "~/work/SpatSocNetworks/" 
write.path <- "~/work/SpatSocNetworks/Output/Batch10_20160705/"
index=1
#install.packages("igraph")
require(igraph)
require(adehabitatHR)
require(MASS)
source(paste(root.path, "Source/CalcSimmedAffinity.R", sep = ""))
source(paste(root.path, "Source/CalcSimmedAffinitySpat.R", sep = ""))
source(paste(root.path, "Source/buildGraph.R", sep = ""))
source(paste(root.path, "Source/SimSpatSocMovtV6.R", sep = ""))
# source(paste(root.path, "Source/SimSpatSocEpidemic.R", sep = ""))

x <- rep(1:10, each = 10)
y <- rep(1:10, times = 10)

S <- cbind(x, y)
# S.dist <- dist(S, method = "euclidean", diag = T, upper = T)
S.dist.torus <- matrix(NA, nrow = dim(S)[1], ncol = dim(S)[1])
for(i in 1:dim(S)[1]){
  for(j in 1:dim(S)[1]){
    S.dist.torus[i,j] <- sqrt(min(abs(S[i, 1] - S[j, 1]), max(x) - abs(S[i, 1] - S[j, 1]))^2 +
                                min(abs(S[i, 2] - S[j, 2]), max(y) - abs(S[i, 2] - S[j, 2]))^2)
  }
#  print(i)
}
  
S.dist.mat <- as.matrix(S.dist.torus)
diag(S.dist.mat) <- .1

# build neighborhood indicators
row.blocks <- c(rep(rep(1:5, each = 2), times = 2),
                rep(rep(6:10, each = 2), times = 2),
                rep(rep(11:15, each = 2), times = 2),
                rep(rep(16:20, each = 2), times = 2),
                rep(rep(21:25, each = 2), times = 2)
)
#neigh.ind.mat <- matrix(row.blocks, byrow = F, nrow = length(levels(factor(x))))

# row.blocks <- c(rep(rep(1:5, each = 4), times = 4),
#                 rep(rep(6:10, each = 4), times = 4),
#                 rep(rep(11:15, each = 4), times = 4),
#                 rep(rep(16:20, each = 4), times = 4),
#                 rep(rep(21:25, each = 4), times = 4)
# )
#neigh.ind.mat <- matrix(row.blocks, byrow = F, nrow = 10)

neigh.ind.mat <- matrix(row.blocks, byrow = F, nrow = length(levels(factor(x))))
neigh.ind.vec <- rep(NA, (length(x) * length(y))*(length(levels(factor(x))) * length(levels(factor(y)))))

for(i in 1:max(x)){
  for(j in 1:max(y)){
    z <- length(levels(factor(x)))
    neigh.ind.vec[i + (j - 1) * z] <- ifelse(row.blocks[i] == row.blocks[j], 1, .001)
  }
}

site.dat <- as.data.frame(cbind(seq(1:length(x)), x, y, row.blocks))
names(site.dat) <- c("Site.no", "Site.x", "Site.y", "Site.block")
Site.neigh.ind <- matrix(NA, nrow = dim(site.dat)[1], ncol = dim(site.dat)[1])
for(i in 1:dim(Site.neigh.ind)[1]){
  for(j in 1:dim(Site.neigh.ind)[1]){
    Site.neigh.ind[i, j] <- ifelse(site.dat$Site.block[i] == site.dat$Site.block[j], 1, .001)
  }
#  print(i)
}

Site.neigh.ind.vec <- as.vector(Site.neigh.ind)
Site.dist <- as.vector(S.dist.mat)
Site.random <- round(runif(length(x) * length(x), min = 0, max = 1), 0)

N <- 50
Timesteps <- 100
N.sites <- length(x)

Site.0 <- c(rep(1, 5), rep(30, 5), rep(70, 5), rep(100, 5), rep(15, 5), rep(45, 5), rep(50, 5), rep(85, 5), rep(26, 5), rep(67, 5))
# Site.0 <- c(rep(1, 5), rep(120, 5), rep(280, 5), rep(400, 5), rep(60, 5), rep(180, 5), rep(200, 5), rep(340, 5), rep(104, 5), rep(268, 5))
social.neighbors <- vector("list", N)
for(i in 1:N){
  social.neighbors[[i]] <- which(Site.0 == Site.0[i])
}
group.index.animal <- rep(NA, length(levels(factor(Site.0))))
for(i in 1:(length(x) * length(x))){
  k <- which(Site.0 == i)
  group.index.animal[i] <- ifelse(length(k) == 0, NA, min(k))
}
NsIndexAnimal <- group.index.animal[Site.0]

alpha.in <- delta.in <- gamma.in <- rho.in <- c(1, 100)
alpha.in <- gamma.in <- .1
delta.in <- c(.01, .1, 1, 10, 100)
gamma.in <- exp(seq(-5, 5, length.out = 20))
rho.in <- exp(seq(-10, 10, length.out = 20))
#delta.in <- exp(seq(-5, 5, length.out = 20))
# specify probablity of infection given contact, and recovery
# P.spat.infect <- c(0, .5, 1)
# P.direct.infect <- c(0, .5, 1)
# P.spat.recover <- .1
# P.individ.recover <- .1

params.in.full <- expand.grid(alpha.in, gamma.in, delta.in, rho.in
                              )

sim.out <- SimSpatSocMovtV6(params.in = params.in.full[index, ], N, N.sites, Timesteps, 
                                Site.0, P.move.mat, site.dat,
                                Site.dist, Site.neigh.ind.vec, Site.random, Site.neigh.ind)
save(sim.out, file = paste(write.path, "output_", index, ".RData", sep = ""))
