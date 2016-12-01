library(doMC)
registerDoMC(8)
set.seed(1)

source("CRM.R")

N_sim <- 100
TR_p1 <- c(0.05, 0.15, 0.30, 0.45, 0.60)  ### scenario 1
TR_p2 <- c(0.05, 0.10, 0.20, 0.30, 0.50)  ### scenario 2
TR_p3 <- c(0.15, 0.30, 0.45, 0.60, 0.85)  ### scenario 3

get.Stats <- function(out,TR_p) {
  # Find % of patients treated at dose level 5
  perc.at.dose <- function(tab, d) {
    hasDose <- as.character(d) %in% names(tab)
    if (!hasDose) 0 else tab[which(names(tab) == as.character(d))] / sum(tab)
  }
  perc.at.each.dose <- function(tab) sapply(1:5,function(d) perc.at.dose(tab,d))
  perc.doses <- apply(sapply(out, function(o) perc.at.each.dose(o$a)),1,mean)
  names(perc.doses) <- 1:5

  # Find % of trials recommending the true MTD (maximum tolerated dose) as the MTD
  chosen.MTD <- sapply(out, function(o) o$c)
  perc.rec.true.MTD <- mean( chosen.MTD == which(TR_p==.3))

  # Find overall % of DLT (Dose-limiting toxicity)
  overall.perc.DLT <- mean(sapply(out, function(o) o$b))

  list(percDoses=perc.doses, percRecTrueMTD=perc.rec.true.MTD, overallPercDLT=overall.perc.DLT)
}

### Simulation:
system.time(out1 <- foreach(i=1:N_sim) %dopar% sim(TR_p1,chrt_size=1))
system.time(out2 <- foreach(i=1:N_sim) %dopar% sim(TR_p2,chrt_size=1))
system.time(out3 <- foreach(i=1:N_sim) %dopar% sim(TR_p3,chrt_size=1))

out <- rbind(unlist(get.Stats(out1,TR_p1)),
             unlist(get.Stats(out2,TR_p2)),
             unlist(get.Stats(out3,TR_p3)))
rownames(out) <- paste0("Scenario",1:3)

print(out)
