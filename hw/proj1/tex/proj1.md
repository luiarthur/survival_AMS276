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
---


## Full Conditionals

The likelihood is
\begin{align*}
\mathcal{L}(\beta,\lambda,\alpha,w,\eta|t,X,\nu) &= \prodl\prodj h(t_{ij})^{\vij} S(t_{ij}) \\
&= \prodl\prodj \bc{\lambda\alpha t_{ij}^{\alpha-1}w_i\exijb}^{\vij} \exp\bc{-\lambda t_{ij}^\alpha \exijb} \\
\end{align*}

The priors for the parameters are
\begin{align*}
\beta &\sim \N_2(m,S) \\
\lambda &\sim \G(a_\lambda,b_\lambda) \\
\alpha &\sim \G(a_\alpha,b_\alpha) \\
w_i | \eta &\sim \G(\eta,\eta) \\
\eta &\sim \G(a_\eta,b_\eta) \\
\end{align*}
where $\N_p(m,S)$ denotes the $p$-dimensional multivariate Normal distribution 
with mean vector $m$ and covariance matrix $S$, and $\G(a,b)$ denotes the
Gamma distribution with the mean being $a/b$. To reproduce the results from
the study in ICS example 4.3, $m$ and $S$ were set to be
$\zero_2$ and $10^3\I_2$ respectively. In addition, $(a_z, b_z)$ was set to be
$(0.001,0.001)$ for $z\in\bc{\lambda,\alpha,\eta}$. 

The resulting complete conditionals for each of the parameters are
\begin{align*}
p(\beta|\lambda,\alpha,w,\eta,t,x,\nu) &\propto \exp\bc{-\frac{(\beta-m)'S^{-1}(\beta-m)}{2}+\suml\sumj \vij \xij'\beta - \lambda\tij^\alpha\exijb} \\
\lambda | \beta,\alpha,w,\eta,t,x,\nu &\sim \G\p{a_\lambda + \suml\sumj\nu_i, b_\lambda+\suml\sumj\tij^\alpha\exijb}\\
p(\alpha | \beta,\lambda,w,\eta,t,x,\nu) &\propto \alpha^{a_\alpha+\suml m_i}\exp\bc{\alpha\p{-b_\alpha+\suml\sumj\vij\log\tij}-\lambda\suml\sumj\tij^\alpha\exijb}\\
w_i | \beta,\lambda,\alpha,\eta,t,x,\nu &\sim \G\p{\eta+\sumj\vij,\eta}\\
p(\eta | \beta,\lambda,\alpha,w,t,x,\nu) &\propto \p{\frac{\eta^\eta}{\Gamma(\eta)}}^N\p{\prodl w_i}^\eta \exp\p{-\eta(b_\eta+\suml w_i)}\eta^{a_\eta-1} \\
\end{align*}


## Posterior Summaries

## Interpretation of Posterior Estimates

## Comparison to Frequentist Frailty Model


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
