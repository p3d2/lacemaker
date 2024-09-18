---
title: "Lace patterns"
layout: collection
permalink: /patterns/
collection: patterns
entries_layout: grid
classes: wide
---

Explore the various tessellation patterns available:

{% for pattern in site.pages %}
  {% if pattern.permalink contains '/patterns/' and pattern.title != 'Patterns' %}
  - [{{ pattern.title }}]({{ pattern.url | relative_url }})
  {% endif %}
{% endfor %}