---
layout: page
permalink: /talks/
title: talks
description: Talks and presentations, in reverse chronological order.
nav: true
nav_order: 4
---

<div class="talks">
{% assign sorted_talks = site.talks | sort: "date" | reverse %}
{% for talk in sorted_talks %}
  <div class="talk-entry" style="margin-bottom: 1.5rem;">
    <h4 style="margin-bottom: 0.2rem;"><a href="{{ talk.permalink | relative_url }}">{{ talk.title }}</a></h4>
    <p style="margin-bottom: 0;">
      <em>{{ talk.type }}</em> &middot; {{ talk.venue }} &middot; {{ talk.location }} &middot; {{ talk.date | date: "%b %-d, %Y" }}
    </p>
    {{ talk.content }}
  </div>
{% endfor %}
</div>
