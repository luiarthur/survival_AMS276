---
title: "AMS 276 - Survival Analysis HW 2"
author: Arthur Lui
date: "8 November, 2016"
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
    - \usepackage{graphicx}
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
    - \newcommand{\iid}{\overset{iid}{\sim}}
    - \newcommand{\ind}{\overset{ind}{\sim}}
    # 
    - \allowdisplaybreaks
    - \def\M{\mathcal{M}}
---

The following Bayesian models were fit for this assignment: 

- $\M_1$: a parametric proportional hazards model assuming Weibull errors.
- $\M_2$: a piecewise constant hazards (PCH) proportional hazards model
- $\M_3$: a proportional hazards model using a gamma process (PC) to model the cumulative baseline hazard 

# a) Priors, Full Conditionals, and Joint Posteriors

## $\M_1$: Proportional Hazards Model (Weibull hazard)

Let $h(t|\beta,x) = h_0(t)\exp(x'\beta)$. Where $h_0(t)$ is a Weibull (baseline) hazard function of the form $h_0(t) = \alpha \lambda t^{\alpha-1}$. Then $S(t|\beta,x) = \exp(-\lambda t^\alpha)$, and the likelihood is 

\begin{align*}
\mathcal{L}(\beta, \alpha, \lambda | t, x, \nu) &= \prodl \bc{h(t_i|\beta,\alpha, \lambda,x_i)}^{\nu_i} S(t_i|\beta,\alpha,\lambda,x_i),
\end{align*}

where $\nu_i$ is 1 if observation $i$ is observed to be a failure, and 0 if the observation is right-censored. (Here, it is assumed observations belong to one of those classes.)

This model is fully specified after defining prior distributions for $\alpha, \lambda$, and $\beta$. The priors chosen are

\begin{align*}
p(\beta) &\propto 1 \\
\alpha &\sim \text{Gamma}(1/10,1/10) \\
\lambda &\sim \text{Gamma}(1/10,1/10), \\
\end{align*}

where the expected value of a Gamma($a$,$b$) random variable is $a/b$. The parameters are assumed to be independent apriori.

The joint posterior is
$$ 
p(\beta,\alpha,\lambda|t,x,\nu) \propto  \mathcal{L}(\beta, \alpha, \lambda | t, x, \nu) p(\beta,\alpha,\lambda),
$$
which can be sampled from via MCMC.

## $\M_2$: Piecewise Constant Hazard Model

As in the previous model, the likelihood is 

\begin{align*}
\mathcal{L}(\beta, \lambda | t, x, \nu) &= \prodl \bc{h(t_i|\beta,x_i)}^{\nu_i} S(t_i|\beta,x) \\
&= \prodl \bc{h_0(t_i) \exp(x_i'\beta)}^{\nu_i} \exp\bc{-H_0(t_i)\exp(x_i'\beta)}\\
\end{align*}

Below are the definitions for $h_0(t)$ and $H_0(t)$.

### The Hazard $h_0(t)$ and Cumulative Hazards $H_0(t)$

Using the piecewise constant hazard model, flexible (Bayesian) nonparametric models can be constructed via nonparametric hazard functions. Before defining the hazard and cumulative hazard functions, a finite partition on the survival
times axis ($t$) needs to be defined. We will define a grid with $J$ intervals to be $\mathbf{s} = \bc{s_{0}, s_{1}, ..., s_{J}}$, with $(s_0,s_1]$ being the first interval, and $(s_{J-1}, s_J]$ being the $J^{th}$ interval. We also impose the restrictions that $s_0 = 0$, $s_J >max(t_i)$, and $s_i < s_j$ for all $i < j$.
One possible partition is constructed using the empirical quantiles of the observed times. Another possible partition could be the set of sorted, unique 
observed survival times. Now we can define the baseline hazard as

$$h_0(t_i|\lambda) = \sum_{j=1}^J \mathbbm{1}_{\bc{s_{j-1} < t_i \le s_j}} \p{\lambda_j} $$

and the cumulative (baseline) hazard function as

$$H_0(t_i|\lambda) = \lambda_g(t_i-s_{g-1})+ \sum_{\bc{j: s_j < t_i \cap j > 0}} \lambda_j(s_j-s_{j-1}),$$ where $g=\inf\bc{j: s_j > t_i}$.

The model is fully specified after defining the prior distributions for
$\beta$ and $\lambda_j$, for $j \in \bc{1,...,J}$. The priors I chose were

\begin{align*}
p(\beta) &\propto 1 \\
\lambda_j &\ind \text{Gamma}(1/10,1/10)\\
\end{align*}

The joint posterior is
$$ 
p(\beta,\lambda|t,x,\nu) \propto  \mathcal{L}(\beta, \lambda | t, x, \nu) p(\beta,\lambda),
$$
which can be sampled from via MCMC.

Note that the grid chosen for this problem was simply the quantiles of the observed times with 10 intervals (the quantiles being evenly spaced).

## $\M_3$: Proportional hazards model using a Gamma Process for $H_0$
[//]: # ( See slides 8. I still don't fully understand.)
The likelihood can be specified by

\begin{align*}
\mathcal{L}(\beta,h|X,t,\nu) &\propto \prod_{j=1}^J\bc{\exp\bc{-h_j \sum_{k\in R_j-D_j}e^{X_k\beta}} \prod_{l\in D_j}\p{1-\exp\bc{-h_j e^{X_l\beta}}}} \\
&\p{= \prod_{j=1}^J S_j f_j???, \text{ where } f_j = P(l \in D_j)} \\
\end{align*}

with priors for $\beta$ and $h$
\begin{align*}
p(\beta) &\propto 1 \\
h_j &\ind \text{Gamma}(c\eta(s_j^\kappa - s_{j-1}^\kappa),c)
\end{align*}

I chose $(c,\eta,\kappa)$ to be (.001,1,1). A small value for $c$ reflects
great prior uncertainty on the Weibull cumulative hazard baseline
centering distribution in the Gamma Process.

The joint posterior is
$$ 
p(\beta,h|t,x,\nu) \propto  \mathcal{L}(\beta, h | t, x, \nu) p(\beta,h),
$$
which can be sampled from via MCMC.


# b) Posterior Distributions and Comparisons

Tables 1, 2, and 3 summarizes the posterior distributions of each of
the parameters in the 3 respective models. The trace plots for all
parameters were examined to check for MCMC chain
convergence. Note that in each model, the posterior mean of
$\beta$ is negative, and the 95% CI do not contain 0 in $\M_1$ and
$\M_2$. The posterior means differ (-0.484, -0.561, -0.464), but the posterior standard deviations are about 0.2 in each case.

---------------------------------------------------------------
  Parameter       mean       sd    CI 2.5%    CI 97.5%  $\ne0$
-----------   --------    ------  --------   --------- --------
$\beta$        -0.4844    0.2195   -0.8643   -0.0190     *

$\alpha$        0.8968    0.0503    0.7958    0.9943

$\lambda$       0.4407    0.0637    0.3310    0.5627
---------------------------------------------------------------

Table: Posterior summary for parameters in $\M_1$ \label{stats1}

Note that the acceptance rate for $(\beta,\lambda,\alpha)$ was 25%.

----------------------------------------------------------------
     Parameters     mean       sd    CI 2.5%   CI 97.5%  $\ne0$
--------------- --------  --------  -------- ---------- --------
$\beta$          -0.5614    0.2680   -1.1155    -0.0436     *

$\lambda_1$       0.3541    0.2495    0.0530    0.9622
 
$\lambda_2$       0.2241    0.1006    0.0672    0.4807

$\lambda_3$       0.1776    0.0706    0.0719    0.3213 

$\lambda_4$       0.2146    0.0752    0.1082    0.3883 

$\lambda_5$       0.0635    0.0289    0.0247    0.1414 

$\lambda_6$       0.1505    0.0685    0.0329    0.3062 

$\lambda_7$       0.0385    0.0294    0.0046    0.1044 

$\lambda_8$       0.3122    0.1169    0.1302    0.5869 

$\lambda_9$       0.0769    0.0575    0.0117    0.2343 

$\lambda_{10}$    0.0645    0.0397    0.0121    0.1796 
----------------------------------------------------------------

Table: Posterior summary for parameters in $\M_2$ \label{stats2}

Note that the acceptance rate for $\beta$ was 29% and that of $\lambda$ was 32%. 

----------------------------------------------
   mean       std     lower     upper  $\ne0$
-------    ------   -------   ------- --------
-0.4640    0.2642   -0.9970    0.0573
----------------------------------------------

Table: Posterior distribution for $\beta$ in $\M_3$ \label{stats3}

Note that the acceptance rate for $\beta$ is 43% and that of $h$ is 40%. The posterior summary for $h$ was not included as there were over 50 parameters, and they are considered nuisance parameters. They are, however, used to obtain the posterior distribution for the survival functions below.

# c) Survival Function Estimates
\beginmyfig
\includegraphics[height=0.5\textwidth]{../img/survival.pdf}
\caption{Survival function estimates for $\M_1$, $\M_2$, and $\M_3$. The aneuploid group is represented by orange and the
diploid group is represented by blue. The 95\% inner credible
intervals are shaded. The Kaplan-Meier curves (in grey) are
also shown for comparison. The Bayesian results agree with the
frequentist results}
\label{fig:mylabel}
\endmyfig

---



[//]: # ( example image embedding
\beginmyfig
\includegraphics[height=0.5\textwidth]{path/to/img/img.pdf}
\caption{some caption}
\label{fig:mylabel}
% reference by: \ref{fig:mylabel}
\endmyfig
)
[//]: # ( example image embedding
> ![some caption.\label{mylabel}](path/to/img/img.pdf){ height=70% }
)

[//]: # ( example two figs side-by-side
\begin{figure*}
  \begin{minipage}{.45\linewidth}
    \centering \includegraphics[height=1\textwidth]{img1.pdf}
    \caption{some caption}
    \label{fig:myLabel1}
  \end{minipage}\hfill
  \begin{minipage}{.45\linewidth}
    \centering \includegraphics[height=1\textwidth]{img2.pdf}
    \caption{some caption}
    \label{fig:myLabel2}
  \end{minipage}
\end{figure*}
)


[//]: # (Footnotes:)

[^1]: A common but alternate parameterization used for the extreme value distribution is 
$$
f_Y(y|\mu,\sigma) = \frac{1}{\sigma}\exp\bc{-\tx -\exp\p{-\tx}},~y\in\mathbb{R}
$$
where $\mu$ and $\sigma$ are location and scale parameters respectively.

