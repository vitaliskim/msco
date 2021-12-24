richness.variances <- function(my.files, Boxplot=TRUE){
  myfiles = lapply(my.files, utils::read.csv, header=T)
  richness_cv <- c()
  Archs <- list()
  community_archs <- list()
  richness_cv_archs <- list()
  RNGkind(sample.kind = "Rejection")
  set.seed(39)
  for (j in 1:length(myfiles)) {
    coe <- msco::Jo.eng(myfiles[[j]], nReps = 999)
    Archs[[j]] <- coe$Archetype ### Archetypes
    richness <- colSums(myfiles[[j]])
    richness_cv[j] <- stats::sd(richness)/mean(richness)
    if(Archs[[j]]=="A1"){
      community_archs$A1[j] <- my.files[[j]]
      richness_cv_archs$A1[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A2"){
      community_archs$A2[j] <- my.files[[j]]
      richness_cv_archs$A2[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A3"){
      community_archs$A3[j] <- my.files[[j]]
      richness_cv_archs$A3[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A4"){
      community_archs$A4[j] <- my.files[[j]]
      richness_cv_archs$A4[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A5"){
      community_archs$A5[j] <- my.files[[j]]
      richness_cv_archs$A5[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A6"){
      community_archs$A6[j] <- my.files[[j]]
      richness_cv_archs$A6[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A7"){
      community_archs$A7[j] <- my.files[[j]]
      richness_cv_archs$A7[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A8"){
      community_archs$A8[j] <- my.files[[j]]
      richness_cv_archs$A8[j] <- richness_cv[j]
    }else if(Archs[[j]]=="A9"){
      community_archs$A9[j] <- my.files[[j]]
      richness_cv_archs$A9[j] <- richness_cv[j]
    }else if(Archs[[j]]=="NA"){
      community_archs$`NA`[j] <- my.files[[j]]
      richness_cv_archs$`NA`[j] <- richness_cv[j]
    }
  }

  for (i in 1:length(community_archs)) {
    community_archs[[i]] <- community_archs[[i]][stats::complete.cases(community_archs[[i]])]
    richness_cv_archs[[i]] <- richness_cv_archs[[i]][stats::complete.cases(richness_cv_archs[[i]])]
  }
  richn <- list()
  richn$community_archs <- community_archs
  richn$richness_cv_archs <- richness_cv_archs
  if(Boxplot==TRUE){
  graphics::boxplot(richness_cv_archs[-which(names(richness_cv_archs)=="NA")],
                                       xlab="Archetypes", ylab="Richness coefficient of variation")
  }
return(richn)
}
# st <- Sys.time()
# my.path <- system.file("extdata/myCSVs", package = "msco")
# setwd(my.path)
# my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
# richn.cv <- msco:::richness.variances(my.files)
# et <- Sys.time();et-st

