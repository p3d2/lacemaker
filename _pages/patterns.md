---
title: "Lace patterns"
layout: single
permalink: /patterns/
collection: patterns
entries_layout: grid
classes: wide
---

<div class="patterns-grid">

  {% for pattern in site.patterns %}
    <div class="pattern-item">
      <a href="{{ site.baseurl }}{{ pattern.url }}">
        <img src="{{ pattern.mesh_jpg | relative_url }}" alt="{{ pattern.title }}" class="pattern-thumbnail" loading="lazy" onerror="this.closest('.pattern-item').hidden=true">
        <p>{{ pattern.title }}</p>
      </a>
    </div>
  {% endfor %}

</div>
