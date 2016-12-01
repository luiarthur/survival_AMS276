rm(list=ls(all=TRUE))
#set.seed(5554)
source("fn_CRM.R")

#########################################################################################################
###  I implement CRM with power model  p(d)= (p^0_d)^a where p^0_d: skeleton (initial value of p(d)) and a: parameter
###  I put p(a)=Gamma(a_gam, b_gam)
#########################################################################################################

###  True probability --- p(d)
#TR_p <- c(0.05, 0.15, 0.30, 0.45, 0.60)  ### scenario 1
#TR_p <- c(0.05, 0.10, 0.20, 0.30, 0.50)  ### scenario 2
#TR_p <- c(0.15, 0.30, 0.45, 0.60, 0.85)  ### scenario 3


sim <- function(TR_p, chrt_size=1)
{
  ### design parameters
  design_para <- NULL
  design_para$m <- 5  ### # of dose levels
  design_para$chrt_size <- chrt_size  ## how many subjects for one group (cohort)
  design_para$N_p <- 20    ### max # of subjects in a trial
  design_para$N_chrt <- design_para$N_p/design_para$chrt_size   ### how many cohorts
  design_para$TTL <- 0.3  ## target toxicity level

  ## hyperparameters
  hyper <- NULL
  hyper$d_0 <- c(0.05, 0.15, 0.30, 0.45, 0.60)  ## skeleton
  hyper$a_gam <- 1
  hyper$b_gam <- 1

  ## MCMC parameter
  n_burn <- 3000
  n_sam <- 3000

  #####################################################
  cur_dat <- NULL
  #####################################################
  ###  FIRST cohort
  cur_dat$a <- 1  ### set a=1 for chrt 1
  cur_dat$d_star <- which.min(abs(hyper$d_0 - design_para$TTL))  ## current MTD
  cur_dat$dose <- rep(cur_dat$d_star, design_para$chrt_size)  ## record dose level for chrt 1 (current MTD)
  cur_dat$y <- fn.sam.dat(TR_p[cur_dat$d_star], design_para$chrt_size)  ## simulate y using the true prob  --- 1: DLT, 0: NO DLT


  #####################################################
  ###  SECOND cohort  --- the LAST cohort

  for(i_chrt in 2:design_para$N_chrt)
  {
      ## use the current data and find the posterior dist of a
      ## then find p(d) based on the updated a and find the current MTD
      tmp <- fn.post.a(n_sam, n_burn, cur_dat$a, cur_dat$y, cur_dat$dose, hyper$d_0, hyper$a_gam, hyper$b_gam, design_para$TTL)
      cur_dat$a <- tmp$a  ## we use the post mean as the initial value of a for the next MCMC running
      cur_dat$d_star <- tmp$d_star
      
      ### assign the next cohort to the current MTD and simulate y for the next chrt
      cur_dat$dose <- c(cur_dat$dose, rep(cur_dat$d_star, design_para$chrt_size))  ## record dose level for chrt 1 (current MTD)
      cur_dat$y <- c(cur_dat$y, fn.sam.dat(TR_p[cur_dat$d_star], design_para$chrt_size))  ## simulate y using the true prob  --- 1: DLT, 0: NO DLT
  }


  ## use ALL data and find the posterior dist of a
  ## then find p(d) based on the updated a and find the *FINAL* MTD
  tmp <- fn.post.a(n_sam, n_burn, cur_dat$a, cur_dat$y, cur_dat$dose, hyper$d_0, hyper$a_gam, hyper$b_gam, design_para$TTL)
  cur_dat$a <- tmp$a  ## we use the post mean as the initial value of a for the next MCMC running
  cur_dat$d_star <- tmp$d_star

  out <- NULL

  out$a <- table(cur_dat$dose)  ### % of patients treated at each dose
  out$b <- mean(cur_dat$y)  ### % of patients who experienced DLT
  out$c <- cur_dat$d_star   ### chosen MTD based on all avaialble data

  return(out)
}




