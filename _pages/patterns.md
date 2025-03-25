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
        <img src="{{ site.baseurl }}{{ pattern.image_path }}" alt="{{ pattern.title }}" class="pattern-thumbnail" loading="lazy">
        <p>{{ pattern.title }}</p>
      </a>
    </div>
  {% endfor %}

</div>
