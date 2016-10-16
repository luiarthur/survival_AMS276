---
title: "AMS 276 - Survival Analysis HW 1"
author: Arthur Lui
date: "11 October, 2016"
geometry: margin=1in
fontsize: 12pt

# Uncomment if using natbib:

# bibliography: BIB.bib
# bibliographystyle: plain 

# This is how you use bibtex refs: @nameOfRef
# see: http://www.mdlerch.com/tutorial-for-pandoc-citations-markdown-to-latex.html)

header-includes: 
    - \usepackage{bm}
    - \usepackage{bbm}
    - \pagestyle{empty}
    - \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    - \newcommand{\p}[1]{\left(#1\right)}
    - \newcommand{\bk}[1]{\left[#1\right]}
    - \newcommand{\bc}[1]{ \left\{#1\right\} }
    - \newcommand{\abs}[1]{ \left|#1\right| }
    - \newcommand{\mat}{ \begin{pmatrix} }
    - \newcommand{\tam}{ \end{pmatrix} }
    - \newcommand{\suml}{ \sum_{i=1}^n }
    - \newcommand{\prodl}{ \prod_{i=1}^n }
    - \newcommand{\ds}{ \displaystyle }
    - \newcommand{\df}[2]{ \frac{d#1}{d#2} }
    - \newcommand{\ddf}[2]{ \frac{d^2#1}{d{#2}^2} }
    - \newcommand{\pd}[2]{ \frac{\partial#1}{\partial#2} }
    - \newcommand{\pdd}[2]{ \frac{\partial^2#1}{\partial{#2}^2} }
    - \newcommand{\N}{ \mathcal{N} }
    - \newcommand{\E}{ \text{E} }
    - \def\given{~\bigg|~}
    # Figures in correct place
    - \usepackage{float}
    - \def\beginmyfig{\begin{figure}[H]\center}
    - \def\endmyfig{\end{figure}}
    # 
    - \allowdisplaybreaks
---

\def\tx{\p{\frac{x-\mu}{\sigma}}}

***

##1a 

\begin{align*}
y_i &\sim Exp(\lambda), ~~~ i=1,...,n \\
\\
\mathcal{L} &\propto \prod_{i=1}^k\lambda\exp\p{-\lambda y_i} 
                     \prod_{j=k+1}^n \exp\p{-\lambda y_j} \\
&\propto \lambda^k \exp\p{-\lambda\suml y_i} \\
\\
l &= \log{\mathcal{L}} = k\log(\lambda) -\lambda\suml y_i \\
\\
\df{l}{\lambda} &= \frac{k}{\lambda} - \suml y_i \\
\ddf{l}{\lambda} &= -\frac{k}{\lambda^2} \propto -\frac{1}{\lambda^2} \\
\\
\pi_J(\lambda) &\propto \sqrt{I(\lambda)} \\
               &\propto \sqrt{-\E\p{\ddf{l}{\lambda}}} \\
               &\propto \sqrt{-\E\p{-\frac{1}{\lambda^2}}} \\
               &\propto \sqrt{\frac{1}{\lambda^2}} \\
               \\
\therefore \pi_J(\lambda) &\propto \frac{1}{\lambda} \\
\end{align*}

##1b

\begin{align*}
p(\lambda|y) &\propto \pi_J(\lambda) \mathcal{L}(\lambda|y) \\
             &\propto \lambda^{-1} \lambda^k\exp\p{-\lambda\suml y_i} \\
             &\propto \lambda^{k-1} \exp\p{-\lambda\suml y_i} \\
             \\
\therefore ~~~ \lambda | y &\sim Gamma(k, \suml y_i) 
\end{align*}

##1c

\begin{align*}
p(z | y) &= \int_0^\infty f(z|\lambda) p(\lambda|y) d\lambda \\
&=\int_0^\infty \lambda e^{-\lambda z} \frac{(\suml y_i)^k}{\Gamma(k)}\lambda^{k-1} \exp\p{-\lambda \suml y_i} d\lambda \\
&= \frac{(\suml y_i)^k}{\Gamma(k)}\frac{\Gamma(k+1)}{(z+\suml y_i)^{k+1}} \\
&= \frac{k(\suml y_i)^k}{(z+\suml y_i)^{k+1}}
\end{align*}

##1d
\begin{align*}
y_i | \theta &\sim Exp(rate=\lambda_i=\exp\p{x_i\beta}) = Exp(\exp\p{x_i\log\theta}) \\
\\
\Rightarrow  y_i | \theta &\sim Exp(\theta^{x_i}) \\
\theta &\sim Gamma(\alpha_0,\lambda_0) \\
\\
p(\theta | y) &\propto \theta^{\alpha_0-1}\exp\bc{-\theta\lambda_0}
                       \prod_{i=1}^k\theta^{x_i}\exp\bc{-\theta^{x_i} y_i} 
                       \prod_{i=k+1}^n\exp\bc{-\theta^{x_i} y_i} \\
              &\propto \theta^{\alpha_0 + (\sum_{i=1}^k x_i)-1}
                       \exp\bc{-\theta\lambda_0-\suml\theta^{x_i} y_i}\\
              &\propto \theta^{\alpha_0 + (\sum_{i=1}^k x_i)-1}
                       \exp\bc{-\theta\lambda_0-\sum_{\bc{i:x_i=1}}\theta^{x_i} y_i}\\
              &\propto \theta^{\alpha_0 + (\sum_{i=1}^k x_i)-1}
                       \exp\bc{-\theta\p{\lambda_0+\sum_{\bc{i:x_i=1}} y_i}}\\
\\
\therefore \theta | y,x &\sim Gamma\p{\alpha_0 + \ds\sum_{i=1}^k x_i,
                                   \lambda_0+ \ds\sum_{\bc{i:x_i=1}} y_i}
\end{align*}

***

##2a
Supposing that $\alpha$ is known,
\begin{align*}
\mathcal{L}(\gamma) &\propto \gamma^{\suml\nu_i} \exp\bc{-\gamma y_i^\alpha} \\
\\
p(\gamma|y) &\propto \mathcal{L}(\gamma) \gamma^{-1} \\
            &\propto  \gamma^{\suml\nu_i} \exp\bc{-\gamma \suml y_i^\alpha}\gamma^{-1} \\
            &\propto \gamma^{\suml\nu_i-1} \exp\bc{-\gamma \suml y_i^\alpha} \\
\therefore p(\gamma|y) &\propto Gamma(\suml\nu_i,~ \suml y_i^\alpha) \\
\end{align*}

##2b
\newcommand{\sua}{\suml\nu_i}
\newcommand{\sub}{\suml y_i^\alpha}

Let $a=\sua$ and $b=\sub$. Then,
\begin{align*}
p(z|y) &= \int_0^\infty p(z|\gamma)p(\gamma|y)d\gamma \\
&= \int_0^\infty \alpha\gamma z^{\alpha-1} \exp\p{-\gamma z^\alpha}
\frac{b^a}{\Gamma(a)} \gamma^{a-1}\exp\p{-\gamma b}d\lambda\\
&= \frac{\alpha z^{\alpha-1}b^a}{\Gamma(a)}\int_0^\infty \gamma^a \exp\p{-\gamma(z^\alpha+b)}d\lambda\\
&= \frac{\alpha z^{\alpha-1}b^a\Gamma(a+1)}{\Gamma(a)(z^\alpha+b)^{a+1}} \\
&= \frac{\alpha a z^{\alpha-1}b^a}{(z^\alpha+b)^{a+1}} \\
&= \ds\frac{\alpha \p{\sua} z^{\alpha-1}\p{\sub}^{\sua}}{(z^\alpha+\sub)^{\p{\sua}+1}} \\
\end{align*}

***

## 3a

Let $\nu_i$ be an indicator for whether observation $i$ was right-censored (=1)
or observed (=0). Let $t_i$ (the time variable in the `tongue` data set) be the
time of death or censorship for observation $i$. In addition let $x_i$ be
an indicator for whether observation $i$ belongs to the aneuploid group (=1)
or the diploid group (=0).

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/post_a.pdf){ height=70% }


blablabla

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/post_b.pdf){ height=70% }

**Need to do sensitivity analysis!!!**

## 3b

### Weibull Accelerated Failure Time (AFT) Model for Aneuploid Data

A Weibull AFT model can be constructed by first that if $Y_i|\alpha,\lambda
\sim V(\alpha,\lambda)$, where $V$ denotes the extreme value
distribution, then $\exp(Y_i) \sim Weibull(\alpha,\lambda)$. The p.d.f.s of the
extreme value and Weibull distributions are given as follows:

\begin{align*}
f_T(t) &= \alpha t^{\alpha-1}\exp\bc{\lambda - t^\alpha\exp\p{\lambda}}, &y>0 \\
\\
f_Y(y|\alpha,\lambda) &= \alpha\exp\bc{\lambda + \alpha y - \exp\p{\lambda+\alpha y}}, &y\in\mathbb{R}  \\
\end{align*}

In an AFT model setting, log time $(\log T_i )$ for a given observation 
is modeled by $-x_i^T\beta + \sigma W_i$. So, the likelihood can be expressed
as
$$
\mathcal{L}(\beta,\sigma| D) = \prodl \bk{f_Y\p{y_i\bigg|\frac{1}{\sigma},-\frac{-x_i^T\beta}{\sigma}}^{\nu_i} S_Y\p{y_i\bigg|\frac{1}{\sigma},-\frac{-x_i^T\beta}{\sigma}}^{1-\nu_i}}\\
$$
where $S_Y(y|\alpha,\lambda) = \exp\bc{-\exp(\lambda+\alpha y)}$ is the
survival function under the extreme value distribution.

Note that here I used the knowledge that $-x_i^T\beta$ and $\sigma$ are
respectively the location and scale of the extreme value
distribution using an *alternate* parameterization[^1].
(Note also, that for this problem, $\beta=(\beta_0,\beta_1)$.)

Placing priors on the unknown parameters $\sigma$ and $\beta$ 
completes the model specification. The priors were chosen
to be
$$
\begin{split}
\sigma &\sim \text{Inverse-Gamma}(2,1) \\
\beta &\sim \text{Normal}(0,3I_2) \\
\end{split}
$$

### Prior Specification for Weibull AFT
The prior for $\sigma$ were chosen to reflect the prior variance to be infinite
and having a prior mean of 1. The prior for $\beta$ was chosen to reflect the
prior mean of the covariates is 0, which is reasonable in a regression setting
when there is no prior knowledge related to the experiment. The prior variance
for $\beta$ is $10I_2$ to reflect that the covariates are not related and each
have a variance of 10. This is also reasonable as the covariates are on log
scale. A variance of 10 on the log scale corresponds to a standard deviation of
about 150 on a non-logged scale, which expresses much uncertainty.

### Posterior Distribution for Weibull AFT Parameters and Results
The posterior distribution of all the parameters cannot be analytically
computed jointly. But by writing out the full conditionals for each parameter,
a Gibbs sampler with 3 separate Metropolis steps can be implemented to obtain
samples from the posterior distribution.

Figure 3 summarizes the posterior distributions of all parameters in the
Weibull AFT model. The diagonals display the univariate posterior distributions
of $\sigma$, $\beta_0$, and $\beta_1$. 

> ![Posterior distribution for $\sigma$ and $\beta$ for Weibull AFT model.](../img/aft_weib.pdf){ height=70% }

Note that within each plot, the trace plots are included and seem to indicate
that the Gibbs sampler has converged. Note also that the posterior mean has
been outlined in red, the standard deviation written, and a 95% equal-tailed
credible interval (CI) were included. The upper subplots of Figure 3 display
the bivariate distribution and trace plots. Again, the trace plots seem to
indicate that the chain has converged. Finally, in the lower subplots, the
correlations between each pair of parameters were included. $\beta_0$ and
$\beta_1$ are more-than-moderately negatively correlated (-.765). This is
expected as when the intercept increases, the relative slopes should decrease.

To briefly summarize in words, $\sigma$ has a posterior mean of 1.272, standard
deviation of 0.15, and 95% CI of (1.026, 1.602). $\beta_0$ has a posterior
mean of -2.031, standard deviation of 0.282, and 95% CI of (-2.604, -1.497).
$\beta_1$ has a posterior mean of -0.669, standard deviation of 0.369, and 95%
CI of (-1.421, 0.07).

### Interpretation of Weibull AFT Coefficients
Recall that the relationship between the covariates and the data is
$$
\begin{split}
Y_i = \log T_i &= -\beta_0 -x_i\beta_1 + \sigma W_i \\
T_i | \beta_0, \beta_1, \sigma &\sim Weibull\p{\frac{1}{\sigma},\frac{\beta_0+x_i\beta_1}{\sigma}} \\
\end{split}
$$
Recall also that the hazard function of the Weibull distribution is
$$
\begin{split}
h_T(t|\alpha,\lambda) &=\alpha e^\lambda t^{\alpha-1} \\
\Rightarrow  h_T\p{t\given\frac{1}{\sigma},\frac{\beta_0+x\beta_1}{\sigma}} &=\frac{1}{\sigma} \exp\p{\frac{\beta_0+x\beta_1}{\sigma}}t^{\frac{1}{\sigma}-1} \\
\Rightarrow  \frac{h(t|x=1)}{h(t|x=0)} = \frac{h(t~|\text{ anueploid})}{h(t~|\text{ diploid})} &= \exp\p{\frac{\beta_1}{\sigma}}
\end{split}
$$
A point estimate of $\exp\p{\ds\frac{\beta_1}{\sigma}}$ can be obtained by 
substituting the posterior means for $\beta_1$ and $\sigma$, yielding
$\exp\p{\ds\frac{\beta_1}{\sigma}} = \exp\p{\ds\frac{-0.699}{1.272}} = 0.577$.
We can say that the relative risk of death (ratio of hazards) for
a patient with aneuploid cells compared to one with diploid cells is 
$0.577$. Or the risk of patients with aneuploid cells is lower than patients
with diploid cells.

We can also say that on average, the log survival time of a patient with
diploid cells is 2.031 (=$-\E\bk{\beta_0|y,X}$), and that of a patient with
aneuploid cells is 6.099 (=$-\E\bk{\beta_1|y,X}$) higher.

## 3c

### Prior Specification for Parameters in Log-logistic AFT Model
The prior distributions chosen for parameters in the log-logistic AFT model
were the same as in the Weibull distribution. The reason being that the
interpretation of the parameters remain mostly unchanged. 

### Posterior Distribution of Parameters in Log-logistic AFT Model
The posterior distribution of parameters in the log-logistic AFT model are 
shown in Figure 4.

> ![Posterior distribution for $\sigma$ and $\beta$ for log-logistic AFT model.](../img/aft_loglog.pdf){ height=70% }

To summarize heuristically, the posterior mean for $\sigma$ is $0.987$, with 
standard deviation 0.117. The posterior means for $(\beta_0,\beta_1)$ are
(-1.372,-0.807) with standard deviations (0.327, 0.413). 

The conclusion under this model is that the risk for patients with aneuploid
cells is 0.441 ($=\exp(\hat{\beta_1}/\hat\sigma)$) that of patients with diploid
cells. That is patients with diploid cells are at higher risk.

### Prior Specification for Parameters in Log-normal AFT Model
Again, the prior distributions chosen for parameters in the log-normal AFT model
were the same as in the Weibull distribution. The reason being that the
interpretation of the parameters remain mostly unchanged. 

### Posterior Distribution of Parameters in Log-normal AFT Model
The posterior distribution of parameters in the log-normal AFT model are 
shown in Figure 5.

> ![Posterior distribution for $\sigma$ and $\beta$ for log-normal AFT model.](../img/aft_lognorm.pdf){ height=70% }

To summarize heuristically, the posterior mean for $\sigma$ is $1.721$, with 
standard deviation 0.187. The posterior means for $(\beta_0,\beta_1)$ are
(-1.36,-0.791) with standard deviations (0.332, 0.412). 

The conclusion under this model is that the risk for patients with aneuploid
cells is 0.631 ($=\exp(\hat{\beta_1}/\hat\sigma)$) that of patients with diploid
cells. That is patients with diploid cells are at higher risk.

## 3d Model Comparison with DIC
The three AFT models can be compared using deviance information criterion
(DIC), defined as 
$$ DIC = \bar{D}(\theta) + \widehat{\text{var}}\p{D(\theta)} / 2 $$
where $D(\theta)= -2\log p(\bm y|\theta)$, $\bar{D}(\theta)$ is the average of
$D(\theta)$ over all posterior samples $\theta$, and
$\widehat{\text{var}}(D(\theta))$ is the sample variance of $D(\theta)$. DIC
measures goodness-of-fit and penalizes model complexity; a lower DIC is
preferred.

The following table summarizes the DICs of the three models

| Model | Weibull | Log-Logistic | Log-Normal |
|:-----:|:-------:|:------------:|:----------:|
|DIC    |  253.977|       253.099| **155.778**|




[//]: # (
\beginmyfig
\includegraphics[height=0.5\textwidth]{../img/aft_weib.pdf}
\endmyfig
)


[//]: # (Footnotes:)

[^1]: A common but alternate parameterization used for the extreme value distribution is 
$$
f_Y(y|\mu,\sigma) = \frac{1}{\sigma}\exp\bc{-\tx -\exp\p{-\tx}},~y\in\mathbb{R}
$$
where $\mu$ and $\sigma$ are location and scale parameters respectively.

