---
title: "spStack: Bayesian Geostatistics Using Predictive Stacking"
permalink: /spStack-dev
excerpt: "Fast Bayesian inference for Gaussian and non-Gaussian geospatial models without using Markov chain Monte Carlo algorithms.
This R package is written in C++ with calls to FORTRAN routines for optimized linear algebra operations.<br/><br/>
Now available on [CRAN](https://cran.r-project.org/package=spStack). Also check out its new [website](https://span-18.github.io/spStack-dev/)!
<br/><br/><img src='/images/spStack-logo.png'>"
collection: software
---

<br/><img src='/images/spStack-logo-small.png'><br/>

Fits Bayesian hierarchical spatial process models for point-referenced Gaussian, Poisson, binomial, and binary data using stacking
of predictive densities. It involves sampling from analytically available posterior distributions conditional upon some candidate
values of the spatial process parameters and, subsequently assimilate inference from these individual posterior distributions using
Bayesian predictive stacking. Our algorithm is highly parallelizable and hence, much faster than traditional Markov chain Monte Carlo
algorithms while delivering competitive predictive performance. See [Zhang, Tang, and Banerjee (2024)](https://www.doi.org/10.48550/arXiv.2304.12414), and, [Pan, Zhang, Bradley, and Banerjee (2024)](https://www.doi.org/10.48550/arXiv.2406.04655) for details.