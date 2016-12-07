---
title: "AMS 276 - Survival Analysis Project 2"
author: Arthur Lui
date: "6 December, 2016"
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

In the Figures below, the blue lines represent the probability of efficacy
($\pi_n(CR)$) and the red lines represent the probability of toxicity 
($\pi_n(TOX)$). Note that in Scenario 1 (Figure \ref{fig:sim1}), early 
stopping is desired as the true probability of a toxicity event is .4 and
that of an efficacy event is .2. This is indeed the case. Note that the blue
vertical lines represent stops for futility, and the red vertical lines
represent stops for toxicity.

In Scenario 2, (Figure \ref{fig:sim1}), early stopping is not desired as the true
probability of a toxicity event is .2 and that of an efficacy event is .4. 
Again, we can note that only in one of the five simulations was there
an early stop for futility, and there were no other stops. 

\beginmyfig
\includegraphics[height=0.6\textwidth]{../img/sim1.pdf}
\caption{Scenario 1: $p_{\text{TRUE}}(\text{CR}) = .2$; $~p_{\text{TRUE}}(\text{TOX}) = .4$;}
\label{fig:sim1}
\endmyfig

\beginmyfig
\includegraphics[height=0.6\textwidth]{../img/sim2.pdf}
\caption{Scenario 2: $p_{\text{TRUE}}(\text{CR}) = .4$; $~p_{\text{TRUE}}(\text{TOX}) = .2$;}
\label{fig:sim2}
\endmyfig



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

