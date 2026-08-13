---
layout: page
permalink: /publications/
title: publications and preprints
nav_title: publications
description: Publications and preprints in reversed chronological order.
nav: true
nav_order: 1
---

<!-- _pages/publications.md -->

<p class="co-first-legend"><sup><i class="fa-solid fa-asterisk author-marker"></i></sup>Denotes joint/equal first authorship. <sup><i class="fa-regular fa-envelope author-marker"></i></sup>Corresponding author.</p>

<!-- Bibsearch Feature -->

{% include bib_search.liquid %}

<div class="publications">

{% bibliography --query @*[show=true] %}

</div>
