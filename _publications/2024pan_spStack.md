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

Analysing non-Gaussian spatial-temporal data typically requires introducing spatial dependence in generalised linear models through the link function of an exponential family distribution. However, unlike in Gaussian likelihoods, inference is considerably encumbered by the inability to analytically integrate out the random effects and reduce the dimension of the parameter space. Iterative estimation algorithms struggle to converge due to the presence of weakly identified parameters. We devise an approach that obviates these issues by exploiting generalised conjugate multivariate distribution theory for exponential families, which enables exact sampling from analytically available posterior distributions conditional upon some fixed process parameters. More specifically, we expand upon the Diaconis-Ylvisaker family of conjugate priors to achieve analytically tractable posterior inference for spatially-temporally varying regression models conditional on some kernel parameters. Subsequently, we assimilate inference from these individual posterior distributions over a range of values of these parameters using Bayesian predictive stacking. We evaluate inferential performance on simulated data, compare with fully Bayesian inference using Markov chain Monte Carlo and apply our proposed method to analyse spatially-temporally referenced avian count data from the North American Breeding Bird Survey database.

<br/><img src='/images/surfaceplots.png'>
