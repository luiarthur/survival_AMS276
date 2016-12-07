set.seed(26)

# Constants:
d_CR <- .2; d_TOX <- .05
k <- 4
n_max <- 75; n_sim <- 5; n_mcmc <- 10000
L_CR <- .02; U_TOX <- .8
B <- 20

theta_S <- c(2.037, 6.111, 30.555, 2.037)
theta_E <- theta_S / sum(theta_S) * 4

scenario <- list(list(p_true_CR=.2, p_true_TOX=.4), 
                 list(p_true_CR=.4, p_true_TOX=.2))

p <- sapply(scenario, function(s) c(s$p_true_CR*(1-s$p_true_TOX),     #A1: CR ~TOX
                                    s$p_true_CR*s$p_true_TOX,         #A2: CR  TOX
                                    (1-s$p_true_CR)*(1-s$p_true_TOX), #A3:~CR ~TOX
                                    (1-s$p_true_CR)*s$p_true_TOX))    #A4:~CR  TOX
                                                                      #CR:  A1 U A2
                                                                      #TOX: A2 U A4


oneSim <- function(ps) { # p scenario
  pi.cr <- matrix(0,B,n_sim)
  pi.tox <- matrix(0,B,n_sim)
  fut.stop <- matrix(FALSE,B,n_sim)
  tox.stop <- matrix(FALSE,B,n_sim)

  for (i in 1:n_sim) {
    # Stays constant within the n_sim=5 trials
    eta.s.cr  <- rbeta(n_mcmc, sum(theta_S[1:2]),    sum(theta_S[3:4]))
    eta.s.tox <- rbeta(n_mcmc, sum(theta_S[c(2,4)]), sum(theta_S[c(1,3)]))
    theta.e <- theta_E

    for (j in 1:B) {
      # cohort
      y <- sample(1:k, k, replace=TRUE, prob=ps) #c(rmultinom(1, k, ps))

      # update theta_e
      for (l in 1:k) theta.e[l] <- theta.e[l] + sum(y==l)

      # eta
      eta.e.cr  <- rbeta(n_mcmc, sum(theta.e[1:2]),    sum(theta.e[3:4]))
      eta.e.tox <- rbeta(n_mcmc, sum(theta.e[c(2,4)]), sum(theta.e[c(1,3)]))

      # pi
      pi.cr[j,i]  <- mean(eta.e.cr  > eta.s.cr  + d_CR)
      pi.tox[j,i] <- mean(eta.e.tox > eta.s.tox + d_TOX)

      if (pi.cr[j,i]  < L_CR)  fut.stop[j,i] <- TRUE
      if (pi.tox[j,i] > U_TOX) tox.stop[j,i] <- TRUE
    }

  }
  list(cr=pi.cr, tox=pi.tox, fut.stop=fut.stop, tox.stop=tox.stop)
}

s1 <- oneSim(p[,1])
s2 <- oneSim(p[,2])

pdf("../img/sim1.pdf")
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:n_sim) {
  plot (s1$cr[,i], type='b',ylim=0:1,col='blue',xlab="cohort number")
  lines(s1$tox[,i],type='b',ylim=0:1,col='red')
  
  min.fut.stop <- min(which(s1$fut.stop[,i]),B+1)
  min.tox.stop <- min(which(s1$tox.stop[,i]),B+1)
  col <- ifelse(min.fut.stop < min.tox.stop, 'blue', 'red')
  abline(v=min(min.fut.stop,min.tox.stop),lty=2,col=col)
  abline(h=c(L_CR,U_TOX),lty=2,col='grey')
}
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

pdf("../img/sim2.pdf")
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:n_sim) {
  plot (s2$cr[,i], type='b',ylim=0:1,col='blue',xlab="cohort number")
  lines(s2$tox[,i],type='b',ylim=0:1,col='red')

  min.fut.stop <- min(which(s2$fut.stop[,i]),B+1)
  min.tox.stop <- min(which(s2$tox.stop[,i]),B+1)
  col <- ifelse(min.fut.stop < min.tox.stop, 'blue', 'red')
  abline(v=min(min.fut.stop,min.tox.stop),lty=2,col=col)
  abline(h=c(L_CR,U_TOX),lty=2,col='grey')
}
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()
