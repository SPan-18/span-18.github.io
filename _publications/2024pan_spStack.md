---
title: "spStack: Practical Bayesian Geostatistics Using Predictive Stacking in R"
collection: publications
category: preprints
permalink: /publication/2024-sptStack
excerpt: 'This manuscript introduces the R package spStack that delivers Bayesian inferece for geospatial data using predictive stacking.'
date: 2024-11-13
venue: 'arXiv preprint'
slidesurl: 'https://doi.org/(in-process)'
paperurl: /files/spStack.txt
citation: '<b>Pan, S.</b> & Banerjee, S. (2024). &quot;spStack: Practical Bayesian Geostatistics Using Predictive Stacking in R.&quot; <i>arXiv:stat.CO</i>.'
---

Statistical modeling and analysis for spatially oriented point-referenced outcomes play a crucial role in diverse scientific applications such as earth and environmental sciences, ecology, epidemiology, and economics. With the advent of Markov chain Monte Carlo (MCMC) algorithms, Bayesian hierarchical models have gained massive popularity in analyzing such point-referenced or, geostatistical data. However, these models involve latent spatial processes characterized by spatial process parameters, which besides lacking substantive relevance in scientific contexts, are also weakly identified and hence, impedes convergence of MCMC algorithms. Thus, even for moderately sized datasets, the computation for MCMC becomes too onerous for practical use. In this article, we introduce the R package spStack that delivers fast Bayesian inference for a class of geostatistical models, where we obviate these issues by sampling from analytically available posterior distributions conditional upon candidate values of the spatial process parameters and, subsequently assimilate inference from these individual posterior distributions using Bayesian predictive stacking. Besides delivering competitive predictive performance as compared to fully Bayesian inference using MCMC, our proposed algorithm is executable in parallel, thus drastically improving runtime and elevating the utility of our software to a diverse group of practitioners with limited computational resources at their disposal.

```r
# install the package from CRAN
install.packages("spStack")

# load library
library("spStack")

# see documentations of some core functions
help(spLMstack)          ## Gaussian spatial model
help(spGLMstack)         ## Non-Gaussian spatial models
```

See the standalone script [spStack-code.R](https://span-18.github.io/files/spStack-code.R) in order to reproduce the examples in the [manuscript](https://span-18.github.io/files/spStack-v1.pdf).
