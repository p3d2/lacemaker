---
title: "Create a pattern"
permalink: /create/
layout: single
sidebar:
  nav: "create"
---

# How to Create JSON Files for Yarn Mesh in LAMMPS

This guide will walk you through the process of creating JSON files to represent a mesh of yarns for export to LAMMPS.

## Step 1: Analyze yarn crossings and repetition

In this step, we'll examine the crossings in the yarn pattern and identify how yarns repeat along their path over this pattern. This will help us understand the structure and repetition of yarns in the JSON file.

1. Start by visualizing the yarn pattern, paying close attention to the crossings where the yarns intersect with each other.

2. Examine the path of individual yarns as they weave through the crossings. Note how the yarns navigate over, under, or around other yarns and the specific pattern they follow.

3. Identify any repetition or regularity in the yarn's path as it repeats across the pattern of crossings. This may include the number of times the yarn goes over or under other yarns before returning to a specific position, or any other repeating sequence. Calculate the distances between nodes.

4. Once you have identified the repetition pattern, determine the number of yarns that are involved in this repetition. This will help when defining the structure of the JSON file.

## Step 2: Identify nodes and paths for describing repetitions

In this step, we will identify the nodes necessary for describing the repetitions in the yarn pattern. It's important to note that defining the smallest unit cell is **not needed**.

To make it easier to identify the yarn pattern, it can be helpful to visualize it. You can sketch the pattern on a paper or use a software tool to draw it digitally (you can use the blank pattern [here]({{ '/assets/images/blank_pattern.svg' | relative_url }}) and open it on [Inkscape](https://inkscape.org/)). This will provide a clear understanding of how the yarns are arranged and repeated.

Let's use the following lace example:

{% raw %}
<div>
  <label for="file-select">Select JSON File:</label>
  <select id="file-select"></select>
  <button id="save-button">Save Changes</button>
</div>

<div id="graph-container" style="width:100%; height:800px;"></div>

<script>
  const dataPath = "{{ '/input/json_patterns/' | relative_url }}";
</script>
<script src="https://d3js.org/d3.v7.min.js"></script>
<script src="{{ '/assets/js/app.js' | relative_url }}"></script>