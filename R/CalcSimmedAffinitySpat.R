CalcSimmedAffinitySpat_Sims <- function(data.in, N.sites){
# Builds spatial network: how quickly do all patches become connected? (get connected if by timestep t, 
# one individual has visited both site i and site j)
  ind <- seq(1:N.sites)
  n <-  length(ind) # number of individuals in the dataset
  
  dyads <-  n*(n-1)/2
  TimesTogether <- DyadID <- Ind1 <- Ind2 <- TimesInd1Separate <- rep(NA,length(dyads))
  TimesInd2Separate <- TotalSpInd1 <- TotalSpInd2 <- SimpleAssociation <- SocialAffinity <- cooccur.sp <- rep(NA,length(dyads))
  Avec <- rep(NA, dyads*5)
  dim(Avec) <- c(dyads,5)
  cooccur <- separate <- BothObserved <- vector("list",dyads)
  
  i=1
  j=2  
  for(k in 1:dyads){
    Avec[k,1] <- ind[i]  # individual 1
    Avec[k,2] <- ind[j]  # individual 2
    ind1<-ind[i]
    ind2<-ind[j]
#    temp.a<-subset(data.in,as.character(data.in$ID) ==
#                     as.character(ind[i]))
#    temp.2a<-subset(data.in,as.character(data.in$ID) == as.character(ind[j]))
    temp.a<-subset(data.in,as.character(data.in$New.S) ==
                     as.character(ind[i]))
    temp.2a<-subset(data.in,as.character(data.in$New.S) == as.character(ind[j]))
    
#    cooccur[[k]] <- merge(temp.a, temp.2a, by = c("New.X","New.Y")) #-- HERE'S THE PROBLEM. 
    cooccur[[k]] <- merge(temp.a, temp.2a, by = c("ID")) #-- HERE'S THE PROBLEM. 
#    cooccur.sp[k] <- dim(unique(subset(cooccur[[k]], select = c("New.X", "New.Y"))))[1]
    cooccur.sp[k] <- dim(unique(subset(cooccur[[k]], select = "ID")))[1]
#    BothObserved[[k]] <- levels(factor(temp.a$Time))[which(levels(factor(temp.a$Time)) %in% levels(factor(temp.2a$Time)))]
#     separate[[k]] <- subset(data.in, Time %in% BothObserved[[k]] & !(Time %in% as.character(cooccur[[k]]$Time)) & ID %in% c(as.character(ind[i]), as.character(ind[j])))
#     TimesTogether[k]<-ifelse(dim(cooccur[[k]])[1]==0,0,length(levels(factor(cooccur[[k]]$Time))))
#     TimesInd1Separate[k]<-ifelse(dim(temp.a)[1]==0,0,length(levels(factor(temp.a$Time))))-TimesTogether[k] - length(levels(factor(separate[[k]]$Time)))
#     TimesInd2Separate[k]<-ifelse(dim(temp.2a)[1]==0,0,length(levels(factor(temp.2a$Time))))-TimesTogether[k]- length(levels(factor(separate[[k]]$Time)))
#     TotalInd1[k]<-ifelse(dim(temp.a)[1]==0,0,length(levels(factor(temp.a$Time))))
#     TotalInd2[k]<-ifelse(dim(temp.2a)[1]==0,0,length(levels(factor(temp.2a$Time))))
    TotalSpInd1[k] <- dim(unique(subset(temp.a, select = c("ID"))))[1]
    TotalSpInd2[k] <- dim(unique(subset(temp.2a, select = c("ID"))))[1]
    DyadID[k]<-paste(ind[i],"_",ind[j],sep="")
    Ind1[k]<-as.character(ind1)
    Ind2[k]<-as.character(ind2)
    #-- krm addition 28 Jan 2014 --#
#    SimpleAssociation[k] <- as.numeric(as.character(TimesTogether[k])) / (as.numeric(as.character(TotalSpInd1[k])) + as.numeric(as.character(TotalSpInd2[k])) - as.numeric(as.character(TimesTogether[k])))
#     SocialAffinity[k] <- TimesTogether[k] / 
#       (min((TimesTogether[k] + length(levels(factor(separate[[k]]$Time))) + TimesInd1Separate[k]), 
#            (TimesTogether[k] + length(levels(factor(separate[[k]]$Time))) + TimesInd2Separate[k])))
    SocialAffinity[k] <- cooccur.sp[k] / (TotalSpInd1[k] + TotalSpInd2[k])
    SocialAffinity[k] <- ifelse(is.na(SocialAffinity[k]), 0, SocialAffinity[k])

    
    i_new<-ifelse(j!=n,i,i+1)
    i<-i_new
    j_new<-ifelse(j==n,i_new+1,j+1)
    j<-j_new
#    print(i)
#    print(j)
  }
  out.mat <- as.data.frame(cbind(Ind1,Ind2,DyadID,TotalSpInd1,TotalSpInd2,SocialAffinity))
  
  out.mat$TotalSpInd1 <- as.numeric(as.character(out.mat$TotalSpInd1))
  out.mat$TotalSpInd2 <- as.numeric(as.character(out.mat$TotalSpInd2))
#  out.mat$TimesTogether <- as.numeric(as.character(out.mat$TimesTogether))
  
  out.list<-list(ind,out.mat)
  return(out.list)
}
