---
title: "AMS 276 - Survival Analysis HW 3"
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

# Simulation Results

The tables below summarize the simulation results. For all
scenario , where the true MTD is dose 3, a larger proportion of
patients receive the true MTD level under a cohort size of 3. 
In Scenario 2, the percentage of patients that receive a dose 
level greater than the true MTD is the same for both cohorts.
But in Scenario 3 (where the true MTD is dose 2), the 
percentage of patients that receive a dose level of 3, is greater in the simulation with a cohort of 3. The proportion
of patients that receive a dose level of 4 is much smaller
in the simulation with a cohort of 3, though.

The proportion of instances where the recommended MTD is
the true MTD is the same in both simulations for scenario
1 and 3. But the proportion of instances where recommended MTD is the true MTD is greater in the first simulation in scenario 2
where the true MTD is high (low toxicity scenario). 

The overall percent of DLT is about the same in both simulations.

\begin{table}[ht]
\centering
\caption{Cohort of 1}
\begin{tabular}{rrrrrrcc}
\\
  \hline
 & Dose1 & Dose2 & Dose3 & Dose4 & Dose5 & Recommend True MTD & overall DLT \\ 
  \hline
Scenario1 & 0.04 & 0.15 & 0.46 & 0.29 & 0.06 & 0.58 & 0.32 \\ 
  Scenario2 & 0.02 & 0.08 & 0.30 & 0.43 & 0.18 & 0.60 & 0.28 \\ 
  Scenario3 & 0.16 & 0.39 & 0.31 & 0.12 & 0.01 & 0.58 & 0.35 \\ 
   \hline
\end{tabular}
\end{table}

\begin{table}[ht]
\centering
\caption{Cohort of 3}
\begin{tabular}{rrrrrrcc}
  \\
  \hline
 & Dose1 & Dose2 & Dose3 & Dose4 & Dose5 & Recommend True MTD & overall DLT \\ 
  \hline
Scenario1 & 0.03 & 0.14 & 0.54 & 0.21 & 0.08 & 0.58 & 0.33 \\ 
  Scenario2 & 0.02 & 0.06 & 0.40 & 0.34 & 0.18 & 0.52 & 0.27 \\ 
  Scenario3 & 0.16 & 0.34 & 0.42 & 0.05 & 0.03 & 0.57 & 0.37 \\ 
   \hline
\end{tabular}
\end{table}


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

