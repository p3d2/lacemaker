---
title: "Simulation Results"
permalink: /results/
layout: single
classes: wide
---

<p>Interactive overview of pore area change across all simulated patterns and yarn contraction combinations. Use the chart to explore the relationship between contraction magnitude, pore opening, and recoverability. Use the table to filter and compare experiments.</p>

<div id="chart3d" style="width:100%;height:520px;"></div>

<hr>

<h2>All experiments</h2>

<div style="display:flex;gap:1rem;flex-wrap:wrap;margin-bottom:1rem;align-items:flex-end;">
  <div>
    <label for="searchInput" style="display:block;font-size:.85em;margin-bottom:.25rem;">Search</label>
    <input id="searchInput" type="text" placeholder="label, name, contract…"
           style="padding:.4rem .6rem;border:1px solid #ccc;border-radius:4px;font-size:.9em;width:220px;">
  </div>
  <div>
    <label for="familyFilter" style="display:block;font-size:.85em;margin-bottom:.25rem;">Family</label>
    <select id="familyFilter" style="padding:.4rem .6rem;border:1px solid #ccc;border-radius:4px;font-size:.9em;">
      <option value="">All</option>
      <option value="Square">Square</option>
      <option value="Diamond">Diamond</option>
      <option value="Hexa">Hexa</option>
      <option value="Rose">Rose</option>
      <option value="Spruce">Spruce</option>
    </select>
  </div>
</div>

<div style="overflow-x:auto;">
<table id="resultsTable" style="width:100%;border-collapse:collapse;font-size:.85em;">
  <thead>
    <tr style="background:#f5f5f5;text-align:left;">
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;cursor:pointer;" onclick="sortTable(0)">Label &#8645;</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;cursor:pointer;" onclick="sortTable(1)">Full name &#8645;</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;">Book ID</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;">Family</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;">Contract</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;cursor:pointer;" onclick="sortTable(5)">&#916;A(t10) &#8645;</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;cursor:pointer;" onclick="sortTable(6)">&#916;A(t20) &#8645;</th>
      <th style="padding:.5rem .75rem;border-bottom:2px solid #ddd;">Best yarn</th>
    </tr>
  </thead>
  <tbody id="tableBody"></tbody>
</table>
</div>

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script>
// Pattern mapping: exp prefix → [label, full_name, book_id, family]
const PATTERN_MAP = {
  "pattern1084":   ["S2",  "Square-2",  "1084",    "Square"],
  "pattern1092":   ["S1",  "Square-1",  "1092",    "Square"],
  "pattern2001":   ["D1",  "Diamond-1", "2001",    "Diamond"],
  "pattern2002":   ["D2",  "Diamond-2", "2002",    "Diamond"],
  "pattern2009":   ["D3",  "Diamond-3", "2009",    "Diamond"],
  "pattern7002":   ["D4",  "Diamond-4", "7002",    "Diamond"],
  "pattern7001":   ["D5",  "Diamond-5", "7001",    "Diamond"],
  "pattern7000":   ["D6",  "Diamond-6", "7000",    "Diamond"],
  "pattern2023":   ["D7",  "Diamond-7", "2023",    "Diamond"],
  "pattern2048":   ["H1",  "Hexa-1",    "2048",    "Hexa"],
  "pattern3023v2": ["R1",  "Rose-1",    "3023v2",  "Rose"],
  "pattern3025":   ["R2",  "Rose-2",    "3025",    "Rose"],
  "pattern3021":   ["R3",  "Rose-3",    "3021",    "Rose"],
  "pattern3024":   ["R4",  "Rose-4",    "3024",    "Rose"],
  "pattern3026":   ["R5",  "Rose-5",    "3026",    "Rose"],
  "pattern3013":   ["R6",  "Rose-6",    "3013",    "Rose"],
  "pattern3018":   ["R7",  "Rose-7",    "3018",    "Rose"],
  "pattern3001":   ["R8",  "Rose-8",    "3001",    "Rose"],
  "pattern3003":   ["R9",  "Rose-9",    "3003",    "Rose"],
  "pattern3034":   ["R10", "Rose-10",   "3034",    "Rose"],
  "pattern3043":   ["R11", "Rose-11",   "3043",    "Rose"],
  "pattern3051_1t":["R12", "Rose-12",   "3051_1t", "Rose"],
  "pattern3051v2": ["R12", "Rose-12",   "3051v2",  "Rose"],
  "pattern3052":   ["R13", "Rose-13",   "3052",    "Rose"],
  "pattern3059":   ["Sp1", "Spruce-1",  "3059",    "Spruce"],
};

const FAMILY_COLORS = {
  Square:  "#FFCC00",
  Diamond: "#0099CC",
  Hexa:    "#66FF66",
  Rose:    "#EE0077",
  Spruce:  "#009900",
};

const LABEL_URLS = {
  S1:"s1",S2:"s2",D1:"d1",D2:"d2",D3:"d3",D4:"d4",D5:"d5",D6:"d6",D7:"d7",
  H1:"h1",R1:"r1",R2:"r2",R3:"r3",R4:"r4",R5:"r5",R6:"r6",R7:"r7",R8:"r8",
  R9:"r9",R10:"r10",R11:"r11",R12:"r12",R13:"r13",Sp1:"sp1"
};

let allRows = [];
let sortCol = -1, sortAsc = true;

function parseCSV(text) {
  const lines = text.trim().split("\n");
  const headers = lines[0].split(",");
  return lines.slice(1).map(line => {
    const vals = line.split(",");
    const obj = {};
    headers.forEach((h, i) => obj[h.trim()] = vals[i] ? vals[i].trim() : "");
    return obj;
  });
}

function enrichRow(row) {
  const info = PATTERN_MAP[row.exp] || ["?", row.sample_full || "?", "?", "?"];
  return {
    label:    info[0],
    fullName: info[1],
    bookId:   info[2],
    family:   info[3],
    contract: row.contract,
    da10:     parseFloat(row.area_change_t10),
    da20:     parseFloat(row.area_change_t20),
    relMax:   parseFloat(row.rel_max_contraction),
    bestYarn: row.best_mol,
    color:    FAMILY_COLORS[info[3]] || "#888",
  };
}

function buildChart(rows) {
  const byFamily = {};
  rows.forEach(r => {
    if (!byFamily[r.family]) byFamily[r.family] = {x:[],y:[],z:[],text:[]};
    byFamily[r.family].x.push(-(r.relMax));
    byFamily[r.family].y.push(r.da10);
    byFamily[r.family].z.push(r.da10 - r.da20);
    byFamily[r.family].text.push(`${r.label} | ${r.fullName}<br>Contract: ${r.contract}<br>&#916;A(t10): ${r.da10.toFixed(3)}<br>&#916;A(t20): ${r.da20.toFixed(3)}`);
  });

  const traces = Object.entries(byFamily).map(([fam, d]) => ({
    type: "scatter3d", mode: "markers",
    name: fam,
    x: d.x, y: d.y, z: d.z,
    text: d.text, hoverinfo: "text",
    marker: { size: 5, color: FAMILY_COLORS[fam] || "#888", opacity: 0.85 }
  }));

  const layout = {
    margin: {l:0,r:0,b:0,t:20},
    scene: {
      xaxis: {title: "Contraction magnitude", range:[-0.05,0.25]},
      yaxis: {title: "&#916;A(t10)",           range:[-0.05,0.75]},
      zaxis: {title: "Recoverability",          range:[-0.05,1.0]},
    },
    legend: {x:0,y:1},
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor:  "rgba(0,0,0,0)",
  };

  Plotly.newPlot("chart3d", traces, layout, {responsive:true});
}

function renderTable(rows) {
  const tbody = document.getElementById("tableBody");
  tbody.innerHTML = "";
  rows.forEach(r => {
    const url = `/lacemaker/patterns/${LABEL_URLS[r.label] || ""}/`;
    const tr = document.createElement("tr");
    tr.style.borderBottom = "1px solid #eee";
    tr.innerHTML = `
      <td style="padding:.4rem .75rem;"><a href="${url}">${r.label}</a></td>
      <td style="padding:.4rem .75rem;">${r.fullName}</td>
      <td style="padding:.4rem .75rem;font-family:monospace;">${r.bookId}</td>
      <td style="padding:.4rem .75rem;">
        <span style="background:${r.color};padding:2px 6px;border-radius:3px;font-size:.8em;color:#000;">${r.family}</span>
      </td>
      <td style="padding:.4rem .75rem;font-family:monospace;">${r.contract}</td>
      <td style="padding:.4rem .75rem;text-align:right;">${r.da10.toFixed(3)}</td>
      <td style="padding:.4rem .75rem;text-align:right;">${r.da20.toFixed(3)}</td>
      <td style="padding:.4rem .75rem;text-align:center;">${r.bestYarn}</td>
    `;
    tbody.appendChild(tr);
  });
}

function filterRows() {
  const query  = document.getElementById("searchInput").value.toLowerCase();
  const family = document.getElementById("familyFilter").value;
  return allRows.filter(r => {
    const matchFamily = !family || r.family === family;
    const matchSearch = !query || [r.label,r.fullName,r.bookId,r.contract,r.family]
      .some(v => v.toLowerCase().includes(query));
    return matchFamily && matchSearch;
  });
}

function applyFilters() { renderTable(filterRows()); }

function sortTable(col) {
  if (sortCol === col) sortAsc = !sortAsc;
  else { sortCol = col; sortAsc = true; }
  const keys = ["label","fullName","bookId","family","contract","da10","da20","bestYarn"];
  allRows.sort((a,b) => {
    const av = a[keys[col]], bv = b[keys[col]];
    if (typeof av === "number") return sortAsc ? av-bv : bv-av;
    return sortAsc ? String(av).localeCompare(String(bv)) : String(bv).localeCompare(String(av));
  });
  applyFilters();
}

const BASE = "{{ site.baseurl }}";
fetch(BASE + "/all_results.csv")
  .then(r => r.text())
  .then(text => {
    allRows = parseCSV(text).map(enrichRow).filter(r => r.label !== "?");
    buildChart(allRows);
    renderTable(allRows);
  })
  .catch(err => {
    document.getElementById("chart3d").innerText = "Could not load results data.";
    console.error(err);
  });

document.getElementById("searchInput").addEventListener("input", applyFilters);
document.getElementById("familyFilter").addEventListener("change", applyFilters);
</script>
