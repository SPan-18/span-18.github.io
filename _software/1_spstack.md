---
layout: page
title: spStack
description: Bayesian Geostatistics Using Predictive Stacking
img: assets/img/spStack-logo.png
importance: 1
category: software
languages: [R, C++, Fortran]
cran: https://cran.r-project.org/package=spStack
github: https://github.com/SPan-18/spStack-dev/
related_publications: true
---

<div class="spstack-figure-float">
    {% include figure.liquid loading="eager" path="assets/img/spStack-logo.png" title="spStack logo" class="img-fluid" %}
</div>

**[spStack](https://span-18.github.io/spStack-dev/)** fits Bayesian hierarchical spatial process models for point-referenced Gaussian, Poisson, binomial, and binary data using stacking of predictive densities. It is written in C++ with calls to Fortran routines for optimized linear algebra operations.

The package samples from analytically available posterior distributions conditional upon candidate values of the spatial process parameters, then assimilates inference from these individual posterior distributions using Bayesian predictive stacking. This algorithm is highly parallelizable and hence much faster than traditional Markov chain Monte Carlo algorithms, while delivering competitive predictive performance. See {% cite pan2024spstack %} and {% cite pan2025_stacking_ba %} for details.

Core functions include:

- [`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.html): Bayesian spatial linear model using predictive stacking
- [`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.html): Bayesian spatial generalized linear model using predictive stacking
- [`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.html): Bayesian spatially-temporally varying coefficients generalized linear model using predictive stacking
- [`recoverGLMscale()`](https://span-18.github.io/spStack-dev/reference/recoverGLMscale.html): Recover posterior samples of scale parameters of spatial/spatial-temporal generalized linear models
- [`posteriorPredict()`](https://span-18.github.io/spStack-dev/reference/posteriorPredict.html): Prediction of latent process at new spatial or temporal locations

Additional matrix algebra utilities, implemented in C++ and exposed through R wrappers, including functionality not available through standard BLAS routines:

- [`cholUpdateRankOne()`, `cholUpdateDel()`, `cholUpdateDelBlock()`](https://span-18.github.io/spStack-dev/reference/cholUpdate.html): Rank-one update, single row/column deletion update and a block deletion update of a Cholesky factor

The stable version is available on [CRAN](https://cran.r-project.org/package=spStack), with a dev version and full documentation on its [own site](https://span-18.github.io/spStack-dev/).

**Links:** [CRAN](https://cran.r-project.org/package=spStack) &middot; [R-universe](https://span-18.r-universe.dev/spStack) &middot; [GitHub](https://github.com/SPan-18/spStack-dev/) &middot; [Package Website](https://span-18.github.io/spStack-dev/)
