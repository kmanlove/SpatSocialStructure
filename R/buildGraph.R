buildGraph <- function(social_mat, data){
    # social_mat = output of CalculateSocialAffinities (or other affinity / assoc generating function)
    # data = studysheep (...I think...)
  graphvector <- as.vector(t(as.matrix(social_mat[, c(1, 2)])))
  assoc_graph <- make_graph(graphvector, directed = F)
  minC <- rep(-Inf, vcount(assoc_graph))
  maxC <- rep(Inf, vcount(assoc_graph))
  minC[1] <- maxC[1] <- 0
  E(assoc_graph)$weight <- as.numeric(as.character(social_mat$SocialAffinity))
#  node.size <- rep(NA, length(V(assoc_graph)$name))
#  for(i in 1:length(node.size)){
#    k <- subset(data, DemogGrp == V(assoc_graph)$name[i])
#    node.size[i] <- mean(na.omit(k$Grpsz))
#  }
#  V(assoc_graph)$size <- node.size 
  
  return(assoc_graph)
}