---
title: "AMS 276 - Survival Analysis Project 1"
author: Arthur Lui
date: "18 November, 2016"
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
    - \renewcommand{\G}{ \text{Gamma} }
    - \newcommand{\E}{ \text{E} }
    - \newcommand{\zero}{ \mathbf{0} }
    - \newcommand{\I}{ \mathbf{I} }
    - \def\given{~\bigg|~}
    # Figures in correct place
    - \usepackage{float}
    - \def\beginmyfig{\begin{figure}[H]\center}
    - \def\endmyfig{\end{figure}}
    - \newcommand{\iid}{\overset{iid}{\sim}}
    - \newcommand{\ind}{\overset{ind}{\sim}}
    # 
    - \allowdisplaybreaks
    - \def\prodj{\prod_{j=1}^{m_i}}
    - \def\sumj{\sum_{j=1}^{m_i}}
    - \def\xij{x_{ij}}
    - \def\tij{t_{ij}}
    - \def\vij{\nu_{ij}}
    - \def\exijb{\exp(x_{ij}'\beta)}
    - \renewcommand{\arraystretch}{1.1}
---


## Full Conditionals

The likelihood is
\begin{align*}
\mathcal{L}(\beta,\gamma,\alpha,w,\eta|t,X,\nu) &= \prodl\prodj h(t_{ij})^{\vij} S(t_{ij}) \\
&= \prodl\prodj \bc{\gamma\alpha t_{ij}^{\alpha-1}w_i\exijb}^{\vij} \exp\bc{-\gamma t_{ij}^\alpha w_i\exijb} \\
\end{align*}

The priors for the parameters are
\begin{align*}
\beta &\sim \N_2(m,S) \\
\gamma &\sim \G(a_\gamma,b_\gamma) \\
\alpha &\sim \G(a_\alpha,b_\alpha) \\
w_i | \eta &\sim \G(\eta,\eta) \\
\eta &\sim \G(a_\eta,b_\eta) \\
\end{align*}
where $\N_p(m,S)$ denotes the $p$-dimensional multivariate Normal distribution 
with mean vector $m$ and covariance matrix $S$, and $\G(a,b)$ denotes the
Gamma distribution with the mean being $a/b$. To reproduce the results from
the study in ICS example 4.3, $m$ and $S$ were set to be
$\zero_2$ and $10^3\I_2$ respectively. In addition, $(a_z, b_z)$ was set to be
$(0.001,0.001)$ for $z\in\bc{\gamma,\alpha,\eta}$. 

The resulting complete conditionals for each of the parameters are
\begin{align*}
p(\beta|\gamma,\alpha,w,\eta,t,x,\nu) &\propto \exp\bc{-\frac{(\beta-m)'S^{-1}(\beta-m)}{2}+\suml\sumj \vij \xij'\beta - \gamma\tij^\alpha w_i\exijb} \\
\gamma | \beta,\alpha,w,\eta,t,x,\nu &\sim \G\p{a_\gamma + \suml\sumj\nu_{ij}, b_\gamma+\suml\sumj\tij^\alpha w_i\exijb}\\
p(\alpha | \beta,\gamma,w,\eta,t,x,\nu) &\propto \alpha^{a_\alpha+(\suml\sumj\nu_{ij}) -1}\times\\
&\exp\bc{-\alpha b_\alpha+\suml\sumj\alpha\vij\log\tij-\gamma\tij^\alpha w_i\exijb}\\
w_i | \beta,\gamma,\alpha,\eta,t,x,\nu &\sim \G\p{\eta+\sumj\vij,\eta+\gamma\sumj\tij^\alpha\exijb}\\
p(\eta | \beta,\gamma,\alpha,w,t,x,\nu) &\propto \p{\frac{\eta^\eta}{\Gamma(\eta)}}^n\p{\prodl w_i}^\eta \exp\p{-\eta(b_\eta+\suml w_i)}\eta^{a_\eta-1} \\
\end{align*}

This suggests that Gibbs-updates can be used to update $\gamma$ and 
$w_i$.


## Posterior Summaries
Table 1 and 2 summarize the posterior distributions of the 
model parameters. Figures \ref{fig:beta} to \ref{fig:eta}
show the posterior distributions along with their traceplots.
Note that the traceplots do not show (strong) evidence that
the MCMC chain has not converged. Also note the acceptance
rates in Table 1 seem to show good mixing. Figure \ref{fig:w} shows the
posterior distribution of the frailties $w$. These results seem 
comparable to the original results presented in the original research.
Note that for the MCMC, 10000 MCMC samples were obtained after a burn-in
period of 1000.

\begin{table}[h]
\begin{center} \begin{tabular}{crrrr}
\hline \hline
Parameter  &   mean   &   std  &  lower   & upper \\
\hline
 $\beta_1$   &     0.0080  & 0.0134 & -0.0187  & 0.0334 \\
 $\beta_2$   &    -1.8112  & 0.5393 & -2.8399  &-0.7397 \\
 $\gamma$   &     0.0173  & 0.0159 &  0.0021  & 0.0626 \\
 $\alpha$   &     1.1961  & 0.1603 &  0.8229  & 1.4754 \\
 $\eta$   &     3.2465  & 4.7624 &  0.7583  &16.0787 \\
 $\kappa$   &     0.5403  & 0.3186 &  0.0621  &1.3197 \\
\hline
 $\beta$ acceptance  &     0.2941  &        &          &   \\
 $\alpha$ acceptance  &     0.3025  &        &          &   \\
 $\eta$ acceptance  &     0.3045  &        &          &  \\
\hline \hline
\end{tabular}\end{center}
\caption {Posterior estimates and summary for Model II}
\end{table}

\beginmyfig
\includegraphics[height=0.5\textwidth]{../img/beta.pdf}
\caption{Posterior distribution for coefficients ($\beta$). Trace plots included in top right corner of each plot.}
\label{fig:beta}
\endmyfig

\begin{figure*}[h]
  \begin{minipage}{.3\linewidth}
    \centering \includegraphics[height=1\textwidth]{../img/alpha.pdf}
    \caption{Posterior distribution of $\alpha$}
    \label{fig:alpha}
  \end{minipage}\hfill
  \begin{minipage}{.3\linewidth}
    \centering \includegraphics[height=1\textwidth]{../img/gamma.pdf}
    \caption{Posterior distribution of $\gamma$}
    \label{fig:gamma}
  \end{minipage}
  \begin{minipage}{.3\linewidth}
    \centering \includegraphics[height=1\textwidth]{../img/eta.pdf}
    \caption{Posterior distribution of $\kappa=1/\eta$}
    \label{fig:eta}
  \end{minipage}
\end{figure*}

\beginmyfig
\includegraphics[height=0.5\textwidth]{../img/w.pdf}
\caption{Posterior distribution for frailties ($w$). 95\% equal-tailed credible intervals are the dashed lines.}
\label{fig:w}
\endmyfig


## Interpretation of Posterior Estimates

- $\beta_{age}$: All else unchanged, patients who are one year older are $\widehat{\exp(\beta_{age})}$=1.008 times more at risk in terms of hazard of dying. Note that this is not very significant.
- $\beta_{sex}$: All else unchanged, female patients are $\widehat{\exp(\beta_{sex})}$=0.189 times more at risk than males in terms of hazard of dying. More simply, **males are 5.28 times more at risk than females**.


\newpage 

## Comparison to Frequentist Frailty Model

Below is the output to a comparable frequentist frailty model.

***

```
coxph(formula = Surv(time, nu) ~ age + sex + frailty(cluster,
    distribution = "gaussian"), data = kidney)

  n= 76, number of events= 58

                          coef      se(coef) se2      Chisq DF    p
age                        0.004713 0.01248  0.008557  0.14  1.00 0.7100
sex                       -1.410822 0.44518  0.315038 10.04  1.00 0.0015
frailty(cluster, distribu                             26.54 14.73 0.0290

Iterations: 6 outer, 39 Newton-Raphson
     Variance of random effect= 0.5691225
Degrees of freedom for terms=  0.5  0.5 14.7
Concordance= 0.82  (se = 0.046 )
Likelihood ratio test= 47.55  on 15.7 df,   p=4.653e-05
```

***

It appears that the estimated coefficients are comparable, as are the standard
errors. The estimate for the variance of the random effect is 0.569, which
is comparable to that of the Bayesian estimate (0.54).

\newpage
\begin{table}
\begin{center}\begin{tabular}{crrrr}
\hline \hline
Parameter  &   mean   &   std  &  lower   & upper \\
\hline
  $w_{1}$      &1.4940    &0.8194    &0.4088    &3.5696  \\
  $w_{2}$      &1.3950    &0.8997    &0.2971    &3.6923  \\
  $w_{3}$      &1.0951    &0.5684    &0.2704    &2.4548  \\
  $w_{4}$      &0.5039    &0.2964    &0.0914    &1.1902  \\
  $w_{5}$      &1.2717    &0.6970    &0.3023    &2.9835  \\
  $w_{6}$      &1.0382    &0.5621    &0.2354    &2.3895  \\
  $w_{7}$      &1.6298    &0.9401    &0.4570    &3.9925  \\
  $w_{8}$      &0.5591    &0.3301    &0.0991    &1.3232  \\
  $w_{9}$      &0.9205    &0.4982    &0.2058    &2.0931  \\
  $w_{10}$     &0.4400    &0.3122    &0.0621    &1.2015  \\
  $w_{11}$     & 0.7988   & 0.4168   & 0.1747   & 1.7890 \\
  $w_{12}$     & 0.9832   & 0.5672   & 0.1817   & 2.3434 \\
  $w_{13}$     & 1.3845   & 0.7058   & 0.3949   & 3.1749 \\
  $w_{14}$     & 0.5747   & 0.4025   & 0.0315   & 1.5126 \\
  $w_{15}$     & 0.5047   & 0.3411   & 0.0560   & 1.3213 \\
  $w_{16}$     & 1.0807   & 0.6429   & 0.2027   & 2.7157 \\
  $w_{17}$     & 0.7684   & 0.4116   & 0.1631   & 1.7272 \\
  $w_{18}$     & 0.7514   & 0.3912   & 0.1720   & 1.6549 \\
  $w_{19}$     & 0.5952   & 0.4283   & 0.0305   & 1.5879 \\
  $w_{20}$     & 1.0435   & 0.6115   & 0.2051   & 2.5577 \\
  $w_{21}$     & 0.1435   & 0.1720   & 0.0101   & 0.7084 \\
  $w_{22}$     & 0.5698   & 0.3500   & 0.0815   & 1.3955 \\
  $w_{23}$     & 1.5296   & 0.8254   & 0.4304   & 3.5982 \\
  $w_{24}$     & 1.1914   & 0.7248   & 0.2403   & 3.0232 \\
  $w_{25}$     & 1.0199   & 0.5374   & 0.2506   & 2.3383 \\
  $w_{26}$     & 0.6493   & 0.3926   & 0.0940   & 1.5733 \\
  $w_{27}$     & 1.0799   & 0.5948   & 0.2560   & 2.5641 \\
  $w_{28}$     & 1.6806   & 0.9493   & 0.4593   & 4.0666 \\
  $w_{29}$     & 1.3085   & 0.6946   & 0.3398   & 3.0083 \\
  $w_{30}$     & 1.2186   & 0.6273   & 0.3264   & 2.7463 \\
  $w_{31}$     & 1.5546   & 0.8539   & 0.4372   & 3.6879 \\
  $w_{32}$     & 1.3231   & 0.8103   & 0.2716   & 3.3657 \\
  $w_{33}$     & 1.1266   & 0.5790   & 0.2834   & 2.5525 \\
  $w_{34}$     & 0.8349   & 0.4784   & 0.1451   & 1.9653 \\
  $w_{35}$     & 1.4338   & 0.7863   & 0.3864   & 3.4240 \\
  $w_{36}$     & 0.8291   & 0.5862   & 0.0575   & 2.2581 \\
  $w_{37}$     & 1.1556   & 0.6788   & 0.2404   & 2.8857 \\
  $w_{38}$     & 0.6065   & 0.3952   & 0.0754   & 1.5552 \\
\hline \hline
\end{tabular}\end{center}
\caption {Posterior estimates and summary for $w$ in Model II}
\end{table}

[//]: # (Footnotes:)

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
