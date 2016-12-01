
####  DO MCMC for burn-in period
fn.MCMC.Burnin <- function(cur_a, y, dose, d_0, a_gam, b_gam, N_sam)
{
    p_cur <- (d_0)^(cur_a)  ## assume power model
    
    for(i_iter in 1:N_sam)
    {
        pro_a <- exp(log(cur_a) + rnorm(1, 0, 0.5))
        p_pro <- (d_0)^(pro_a)
        
        A_cur <- a_gam*log(cur_a) - b_gam*cur_a; ### prior -- gamma(a_gam, b_gam) + jacobian
        A_pro <- a_gam*log(pro_a) - b_gam*pro_a; ### prior--- gamma(a_gam, b_gam) + jacobian
        
        ### likelihood --- independent Ber.
        A_cur <- A_cur +  sum(y*log(p_cur[dose]) + (1-y)*log(1.0-p_cur[dose]))
        A_pro <- A_pro +  sum(y*log(p_pro[dose]) + (1-y)*log(1.0-p_pro[dose]))
        
        ## accept?
        if(log(runif(1)) < (A_pro - A_cur))
        {
            cur_a <- pro_a
            p_cur <- p_pro
        }
    }
    
    return(cur_a)
}

####  DO MCMC and return posterior sample
### cur_a <- cur_dat$a; y <- cur_dat$y; dose <- cur_dat$dose; d_0 <- hyper$d_0; a_gam <- hyper$a_gam; b_gam <- hyper$b_gam; N_sam <- n_sam
fn.MCMC.Save <- function(cur_a, y, dose, d_0, a_gam, b_gam, N_sam)
{
    Save_a <- rep(NA, N_sam)
    p_cur <- (d_0)^(cur_a)
    
    for(i_iter in 1:N_sam)
    {
        pro_a <- exp(log(cur_a) + rnorm(1, 0, 0.5))
        p_pro <- (d_0)^(pro_a)
        
        A_cur <- a_gam*log(cur_a) - b_gam*cur_a; ### prior -- gamma(a_gam, b_gam) + jacobian
        A_pro <- a_gam*log(pro_a) - b_gam*pro_a; ### prior--- gamma(a_gam, b_gam) + jacobian
        
        ### likelihood --- independent Ber.
        A_cur <- A_cur +  sum(y*log(p_cur[dose]) + (1-y)*log(1.0-p_cur[dose]))
        A_pro <- A_pro +  sum(y*log(p_pro[dose]) + (1-y)*log(1.0-p_pro[dose]))
        
        ## accept?
        if(log(runif(1)) < (A_pro - A_cur))
        {
            cur_a <- pro_a
            p_cur <- p_pro
        }
        
        ## save current sample
        Save_a[i_iter] <- cur_a
    }
    
    return(Save_a)
}

###  Update the posteiror mean of a using the current data and find MTD estimate
fn.post.a <- function(N_sam, N_burn, cur_a, y, dose, d_0, a_gam, b_gam, TTL)
{
    ## Run MCMC for burn-in period
    cur_a <- fn.MCMC.Burnin(cur_a, y, dose, d_0, a_gam, b_gam, N_burn)
    
    ## Run MCMC after burn-in period and save MCMC sample
    Save_a <- fn.MCMC.Save(cur_a, y, dose, d_0, a_gam, b_gam, N_sam)
    
    
    ### Let's find the current MTD
    post_m_a <- mean(Save_a)  ### post mean of a
    p_d <- d_0^post_m_a       ### updated p(d)
    d_star <- which.min(abs(p_d - TTL))  ## Update MTD

    
    return(list(a=post_m_a, d_star = d_star))
    
}


### simulate y using the true prob.
fn.sam.dat <- function(tr_prob_d, N)
{
    u <- runif(N)
    
    y <- rep(0, N)
    y[u < tr_prob_d] <- 1
    
    return(y)
}
