---
layout: page
title: spStackCOS
description: Bayesian Inference Under Spatial-Temporal Misalignment
img: assets/img/spStackCOS_hex.png
importance: 2
category: software
languages: [R, C++, Fortran]
github: https://github.com/SPan-18/spStackCOS-dev
related_publications: true
---

<div class="spstack-figure-float">
    {% include figure.liquid loading="eager" path="assets/img/spStackCOS_hex.png" title="spStackCOS logo" class="img-fluid" %}
</div>

**spStackCOS** implements Bayesian predictive stacking for hierarchical regression of spatially-temporally misaligned data. It develops a Bayesian hierarchical modeling framework to analyze spatially-temporally misaligned exposure and health outcome data, using predictive stacking to optimally combine multiple spatial-temporal predictive models while avoiding iterative estimation algorithms such as Markov chain Monte Carlo. See {% cite pan2026_envcs %} for details.

**Links:** [GitHub](https://github.com/SPan-18/spStackCOS-dev)

The following functions fit a Bayesian spatial-temporal linear model and generate samples from the posterior predictive distribution at a spatial-temporal resolution different from what the data is observed at:

- `spLMexactCOS()`: spatial-temporal linear model for misaligned data with fixed values of process parameters
- `sptBlockPredictTimeAgg()`: spatial-temporal linear model for temporally aggregated misaligned data
