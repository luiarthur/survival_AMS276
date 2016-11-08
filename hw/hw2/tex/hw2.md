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
    # 
    - \allowdisplaybreaks
    - \def\M{\mathcal{M}}
---

The following Bayesian models were fit for this assignment: 

- $\M_1$: a parametric proportional hazards model assuming Weibull errors.
- $\M_2$: a piecewise constant hazards (PCH) proportional hazards model
- $\M_3$: a proportional hazards model using a gamma process (PC) to model the cumulative baseline hazard 

# a) Priors, Full Conditionals, and Joint Posteriors

## $\M_1$: Weibull Model

\begin{align*}
\mathcal{L}(\beta, \alpha, \gamma | t) &\propto \gamma t^\alpha\\
\end{align*}

# b) Posterior Distributions and Comparisons

Tables 1, 2, and 3 summarizes the posterior distributions of 
each of the parameters in the 3 respective models. The trace plots for all parameters were examined to reasonably check
for MCMC chain convergence. Note that in each model, the posterior
mean of $\beta$ is negative, and the 95% CI do not contain
0. The posterior means differ (-0.484, -0.797, -1.66), but the
posterior standard deviations are about 0.22 in each case.

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
$\beta$          -0.7968    0.2399   -1.3291    -0.3665     *

$\lambda_1$       0.4938    0.3036    0.1189     1.1371

$\lambda_2$       0.2790    0.0970    0.1288     0.4827

$\lambda_3$       0.2214    0.0841    0.0940     0.4241

$\lambda_4$       0.2583    0.0808    0.1244     0.4417

$\lambda_5$       0.0726    0.0287    0.0256     0.1361

$\lambda_6$       0.2045    0.0935    0.0692     0.4466

$\lambda_7$       0.0673    0.0495    0.0022     0.1940

$\lambda_8$       0.4171    0.1574    0.1961     0.8030

$\lambda_9$       0.1105    0.0717    0.0284     0.3098

$\lambda_{10}$    0.0908    0.0431    0.0149     0.1947
----------------------------------------------------------------

Table: Posterior summary for parameters in $\M_2$ \label{stats2}

Note that the acceptance rate for $\beta$ was 29% and that of $\lambda$ was 29% also. 

----------------------------------------------------------
 Parameters    mean       std     lower     upper  $\ne0$
----------- -------    ------   -------   ------- --------
 $\beta_1$  -1.6673    0.2239   -2.1009   -1.2251    *

  $h_1$      0.0673    0.0382    0.0185    0.1597

  $h_2$      0.1154    0.0546    0.0420    0.2452

  $h_3$      0.0739    0.0459    0.0168    0.1848

  $h_4$      0.0960    0.0484    0.0261    0.2313

  $h_5$      0.0545    0.0329    0.0138    0.1360

  $h_6$      0.0710    0.0434    0.0183    0.1973

  $h_7$      0.0352    0.0266    0.0056    0.1192

  $h_8$      0.1284    0.0704    0.0383    0.3075

  $h_9$      0.0882    0.0558    0.0257    0.2325

  $h_{10}$   0.0875    0.0454    0.0293    0.1983

  $h_{11}$   0.0623    0.0551    0.0034    0.1837

  $h_{12}$   0.0790    0.0488    0.0119    0.1904

  $h_{13}$   0.1335    0.0751    0.0176    0.2952

  $h_{14}$   0.1411    0.0731    0.0324    0.3066

  $h_{15}$   0.0677    0.0489    0.0104    0.1843

  $h_{16}$   0.1703    0.0847    0.0551    0.3689

  $h_{17}$   0.0687    0.0433    0.0063    0.1728

  $h_{18}$   0.0436    0.0502    0.0057    0.1778

  $h_{19}$   0.1430    0.0827    0.0081    0.3348

  $h_{20}$   0.0678    0.0442    0.0050    0.1714

  $h_{21}$   0.0980    0.0790    0.0097    0.3006

  $h_{22}$   0.0582    0.0505    0.0037    0.1728

  $h_{23}$   0.1344    0.0943    0.0173    0.3696

  $h_{24}$   0.1137    0.0713    0.0175    0.2970

  $h_{25}$   0.1343    0.1021    0.0100    0.3440

  $h_{26}$   0.1235    0.0874    0.0285    0.4053

  $h_{27}$   0.1427    0.1127    0.0087    0.4322

  $h_{28}$   0.1414    0.1029    0.0187    0.3878

  $h_{29}$   0.1772    0.1157    0.0377    0.4619

  $h_{30}$   0.0877    0.0765    0.0078    0.2872

  $h_{31}$   0.0394    0.0645    0.0007    0.2532

  $h_{32}$   0.1821    0.1125    0.0151    0.4852

  $h_{33}$   0.0954    0.0861    0.0043    0.3108

  $h_{34}$   0.0879    0.1026    0.0038    0.3958

  $h_{35}$   0.0628    0.0550    0.0058    0.1901

  $h_{36}$   0.0679    0.0585    0.0044    0.2243

  $h_{37}$   0.1243    0.1057    0.0226    0.4412

  $h_{38}$   0.0653    0.0628    0.0014    0.2197

  $h_{39}$   0.1547    0.1041    0.0365    0.4080

  $h_{40}$   0.1155    0.0888    0.0132    0.3166

  $h_{41}$   0.1187    0.1001    0.0104    0.3785

  $h_{42}$   0.0605    0.0560    0.0082    0.2082

  $h_{43}$   0.1938    0.1147    0.0534    0.4665

  $h_{44}$   0.0715    0.0819    0.0068    0.3287

  $h_{45}$   0.4670    0.1838    0.1456    0.8349

  $h_{46}$   0.1868    0.1902    0.0198    0.7622

  $h_{47}$   0.0571    0.0547    0.0014    0.2100

  $h_{48}$   0.2204    0.1894    0.0233    0.6780

  $h_{49}$   0.1204    0.1099    0.0087    0.3806

  $h_{50}$   0.2358    0.1849    0.0370    0.7326

  $h_{51}$   0.2759    0.2316    0.0361    0.9408

  $h_{52}$   0.1278    0.1740    0.0083    0.7280

  $h_{53}$   0.2469    0.1736    0.0371    0.7041

  $h_{54}$   0.3945    0.2540    0.0768    0.9621

  $h_{55}$   0.1211    0.1964    0.0055    0.7844

  $h_{56}$   0.6103    0.3996    0.1042    1.6281

  $h_{57}$   0.6701    0.5738    0.0912    2.0054

  $h_{58}$   0.5136    0.6741    0.0411    2.2653

  $h_{59}$   0.0000    0.0000    0.0000    0.0000
----------------------------------------------------------

Table: Posterior summary for parameters in $\M_3$ \label{stats3}

Note that the acceptance rate for $\beta$ is 56% and that of $h$ is 26%.

# c) Survival Function Estimates
\beginmyfig
\includegraphics[height=0.5\textwidth]{../img/survival.pdf}
\caption{Survival Function Estimates for $\M_1$, $\M_2$, and $\M_3$. The aneuploid group is represented by orange and the
diploid group is represented by blue. The 95\% inner credible
intervals are shaded. The Kaplan-Meier curves (in grey) are
also shown for comparison.}
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

