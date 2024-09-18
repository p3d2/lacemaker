---
title: "Lace patterns"
layout: posts
permalink: /patterns/
---

{% assign sorted_pages = site.pages | where_exp: "page", "page.permalink contains '/patterns/' and page.title != 'Patterns'" | sort: "title" %}

{% for page in sorted_pages %}
- [{{ page.title }}]({{ page.url | relative_url }})
{% endfor %}
