---
layout: page
title: software
permalink: /software/
description: Open-source statistical software.
nav: true
nav_order: 2
horizontal: false
---

<!-- pages/software.md -->
<div class="projects">
{% assign sorted_software = site.software | sort: "importance" %}
<div class="row row-cols-1 row-cols-md-3">
{% for project in sorted_software %}
  {% include projects.liquid %}
{% endfor %}
</div>
</div>
