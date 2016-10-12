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
---

#1a 
$$
\begin{split}
y_i &\sim Exp(\lambda), ~~~ i=1,...,n \\
\\
\mathcal{L} &\propto \prod_{i=1}^k\lambda\exp\p{-\lambda y_i} 
                     \prod_{j=k+1}^n \exp\p{-\lambda y_j} \\
&\propto \lambda^k \exp\p{-\lambda\suml y_i} \\
\\
l &= \log{\mathcal{L}} = k\log(\lambda) -\lambda\suml y_i \\
\\
\df{l}{\lambda} &= \frac{k}{\lambda} - \suml y_i \\
\ddf{l}{\lambda} &= -\frac{k}{\lambda^2} \\
                 &\propto -\frac{1}{\lambda^2} \\
\end{split}
$$

$$
\begin{split}
\pi_J(\lambda) &\propto \sqrt{I(\lambda)} \\
               &\propto \sqrt{-\E\p{\ddf{l}{\lambda}}} \\
               &\propto \sqrt{-\E\p{-\frac{1}{\lambda^2}}} \\
               &\propto \sqrt{\frac{1}{\lambda^2}} \\
               &\propto \frac{1}{\lambda} \\
\therefore \pi_J(\lambda) &\propto \frac{1}{\lambda} \\
\end{split}
$$

#1b
$$
\begin{split}
p(\lambda|y) &\propto \pi_J(\lambda) \mathcal{L}(\lambda|y) \\
             &\propto \lambda^{-1} \lambda^k\exp\p{-\lambda\suml y_i} \\
             &\propto \lambda^{k-1} \exp\p{-\lambda\suml y_i} \\
             \\
\therefore ~~~ \lambda | y &\sim Gamma(k, \suml y_i) 
\end{split}
$$

#1c
$$
\begin{split}
p(z | y) &= \int_0^\infty f(z|\lambda) p(\lambda|y) d\lambda \\
&=\int_0^\infty \lambda e^{-\lambda z} \frac{(\suml y_i)^k}{\Gamma(k)}\lambda^{k-1} \exp\p{-\lambda \suml y_i} d\lambda \\
&= \frac{(\suml y_i)^k}{\Gamma(k)}\frac{\Gamma(k+1)}{(z+\suml y_i)^{k+1}} \\
&= \frac{k(\suml y_i)^k}{(z+\suml y_i)^{k+1}}
\end{split}
$$

#1d

