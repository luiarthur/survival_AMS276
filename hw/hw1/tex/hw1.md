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
    # Figures in correct place
    - \usepackage{float}
    - \def\beginmyfig{\begin{figure}[H]\center}
    - \def\endmyfig{\end{figure}}
    # 
    - \allowdisplaybreaks
---

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

### Weibull Model for Aneuploid Data

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/post_a.pdf){ height=70% }


blablabla

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/post_b.pdf){ height=70% }

## 3b


bla

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/aft_weib.pdf){ height=70% }

bla

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/aft_lognorm.pdf){ height=70% }

> ![Posterior distribution for $\lambda$ and $\alpha$ for aneuploid data.](../img/aft_loglog.pdf){ height=70% }

