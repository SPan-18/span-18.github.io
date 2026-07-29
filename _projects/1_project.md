---
layout: page
title: spStack
description: Bayesian Geostatistics Using Predictive Stacking
img: assets/img/spStack-logo.png
importance: 1
category: software
related_publications: true
---

**spStack** fits Bayesian hierarchical spatial process models for point-referenced Gaussian, Poisson, binomial, and binary data using stacking of predictive densities. It is written in C++ with calls to Fortran routines for optimized linear algebra operations.

The package samples from analytically available posterior distributions conditional upon candidate values of the spatial process parameters, then assimilates inference from these individual posterior distributions using Bayesian predictive stacking. This algorithm is highly parallelizable and hence much faster than traditional Markov chain Monte Carlo algorithms, while delivering competitive predictive performance. See {% cite pan2024spstack %} and {% cite pan2024sptglmstack %} for details.

Core functions include:

- `spLMexact()`: spatial linear model with fixed values of process parameters
- `spLMstack()`: spatial linear model using Bayesian predictive stacking
- `spGLMexact()`: spatial GLM with fixed values of process parameters
- `spGLMstack()`: spatial GLM using Bayesian predictive stacking

The stable version is available on [CRAN](https://cran.r-project.org/package=spStack), with a dev version and full documentation on its [own site](https://span-18.github.io/spStack-dev/).

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/spStack-logo.png" title="spStack logo" class="img-fluid rounded z-depth-1" %}
    </div>
</div>

**Links:** [CRAN](https://cran.r-project.org/package=spStack) &middot; [R-universe](https://span-18.r-universe.dev/spStack) &middot; [GitHub](https://github.com/SPan-18/spStack-dev/)
