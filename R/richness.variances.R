st <- Sys.time()
my.path <- system.file("extdata/myCSVs", package = "msco")
setwd(my.path)
my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
myfiles = lapply(my.files, utils::read.csv, header=T)

richness_var <- c()
Archs <- list()
community_archs <- list()
richness_var_archs <- list()
RNGkind(sample.kind = "Rejection")
set.seed(14)
for (j in 1:length(myfiles)) {
  coe <- msco::Jo.eng(myfiles[[j]], nReps = 999)
  Archs[[j]] <- coe$Archetype ### Archetypes
  richness_var[j] <- stats::var(colSums(myfiles[[j]]))
  if(Archs[[j]]=="A1"){
    community_archs$A1[j] <- my.files[[j]]
    richness_var_archs$A1[j] <- richness_var[j]
  }else if(Archs[[j]]=="A2"){
    community_archs$A2[j] <- my.files[[j]]
    richness_var_archs$A2[j] <- richness_var[j]
  }else if(Archs[[j]]=="A3"){
    community_archs$A3[j] <- my.files[[j]]
    richness_var_archs$A3[j] <- richness_var[j]
  }else if(Archs[[j]]=="A4"){
    community_archs$A4[j] <- my.files[[j]]
    richness_var_archs$A4[j] <- richness_var[j]
  }else if(Archs[[j]]=="A5"){
    community_archs$A5[j] <- my.files[[j]]
    richness_var_archs$A5[j] <- richness_var[j]
  }else if(Archs[[j]]=="A6"){
    community_archs$A6[j] <- my.files[[j]]
    richness_var_archs$A6[j] <- richness_var[j]
  }else if(Archs[[j]]=="A7"){
    community_archs$A7[j] <- my.files[[j]]
    richness_var_archs$A7[j] <- richness_var[j]
  }else if(Archs[[j]]=="A8"){
    community_archs$A8[j] <- my.files[[j]]
    richness_var_archs$A8[j] <- richness_var[j]
  }else if(Archs[[j]]=="A9"){
    community_archs$A9[j] <- my.files[[j]]
    richness_var_archs$A9[j] <- richness_var[j]
  }else if(Archs[[j]]=="NA"){
    community_archs$`NA`[j] <- my.files[[j]]
    richness_var_archs$`NA`[j] <- richness_var[j]
  }
}
et <- Sys.time();et-st
richness_var
for (i in 1:length(community_archs)) {
  community_archs[[i]] <- community_archs[[i]][stats::complete.cases(community_archs[[i]])]
  richness_var_arc <- richness_var_archs[[i]][stats::complete.cases(richness_var_archs[[i]])]
  richness_var_archs[[i]] <- richness_var_arc/max(richness_var_arc)
}
community_archs
richness_var_archs
boxplot(richness_var_archs[-which(names(richness_var_archs)=="NA")], xlab="Archetypes", ylab="Richness variance")
# rich <- list()
# richness_var <- richness_var/max(richness_var)
# for (j in 1:length(richness_var_archs)) {
#   rich[[j]] <- richness_var[as.numeric(gsub(".csv", "", community_archs[[j]]))]
# }
# rich <- `names<-`(rich, names(richness_var_archs))
# boxplot(rich)








