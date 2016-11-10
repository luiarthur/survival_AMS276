using RCall

#=
model:
# h(tᵢⱼ| wᵢ,xᵢⱼ) = h₀(tᵢⱼ) wᵢ exp(xᵢⱼ'β)
=# 
R"""
library(survival)
kidney <- read.table("dat/kidney.txt",header=TRUE,skip=5)
kidney$sex <- ifelse(kidney$sex==2, 0, 1) # M=1, F=0
cox_mod <- coxph(Surv(time,nu) ~ sex+age, data=kidney)
frail_mod <- coxph(Surv(time,nu) ~ sex+age+frailty(cluster,theta=1),
                   data=kidney)
"""

