document.addEventListener('DOMContentLoaded', () => {
  const fileSelect = document.getElementById('file-select');
  const graphContainer = document.getElementById('graph-container');
  let graphData;
  let svg;
  let currentView = 'figure1'; // 'figure1' or 'figure2'

  // List of JSON files (update with your actual file names)
  const files = [
    'pattern1024.json',
    'pattern1084.json',
    'pattern1093.json'
    // Add your files here
  ];

  // Populate the select dropdown
  files.forEach(file => {
    const option = document.createElement('option');
    option.value = file;
    option.textContent = file;
    fileSelect.appendChild(option);
  });

  fileSelect.addEventListener('change', () => {
    const fileName = fileSelect.value;
    loadGraph(fileName);
  });

  // Switch view button
  const viewToggleButton = document.getElementById('view-toggle-button');
  viewToggleButton.addEventListener('click', () => {
    if (currentView === 'figure1') {
      currentView = 'figure2';
      viewToggleButton.textContent = 'Pattern View';
    } else {
      currentView = 'figure1';
      viewToggleButton.textContent = 'Mesh View';
    }
    renderGraph();
  });

  // Load the initial graph
  if (files.length > 0) {
    loadGraph(files[0]);
  }

  function loadGraph(fileName) {
    // Clear previous SVG
    d3.select('#graph-container').selectAll('*').remove();

    fetch(`${dataPath}${fileName}`)
      .then(response => response.json())
      .then(data => {
        graphData = data;
        renderGraph();
      })
      .catch(error => {
        console.error('Error loading JSON file:', error);
        alert('Error loading JSON file. Check the console for details.');
      });
  }

  function renderGraph() {
    // Clear previous SVG content
    d3.select('#graph-container').selectAll('*').remove();
    
    const width = graphContainer.clientWidth;
    const height = graphContainer.clientHeight;

    // Set up SVG
    svg = d3.select('#graph-container').append('svg')
      .attr('width', width)
      .attr('height', height);

    if (currentView === 'figure1') {
      renderFigure1();
    } else {
      renderFigure2();
    }
  }

  function renderFigure1() {
    const G = new Map();
    const nodeColors = new Map();
    const labelDict = new Map();
    const posToNodeId = new Map();
    const posBase = new Map();
    const pathColors = ['#ff6666', '#6666ff', '#ffff00', '#ff66ff', '#66ff33', '#ccffb3', '#b300ff', '#33ffff'];

    // Parse nodes
    for (const [nodeId, pos] of Object.entries(graphData.nodes)) {
      const roundedPos = pos.slice(0, 2).map(p => Math.round(p * 100) / 100);
      G.set(nodeId, { pos: roundedPos });
      nodeColors.set(nodeId, 'white');
      labelDict.set(nodeId, nodeId);
      posToNodeId.set(roundedPos.toString(), nodeId);
    }

    for (const [key, value] of posToNodeId.entries()) {
      posBase.set(key, value);
    }

    // Paths
    const paths = {};
    for (let index = 0; index < graphData.paths.length; index++) {
      const pathInfo = graphData.paths[index];
      paths[index] = pathInfo.path;
    }

    // Build the graph
    for (const [yarnId, yarnData] of Object.entries(graphData.unit_yarns)) {
      const pathId = yarnData[0];
      const nodeStartIndex = yarnData[1].toString();
      const pathStartIndex = yarnData[2];
      let pathZ0 = yarnData[3];
      const yarnPath = paths[pathId].map(toInt);

      // Unit repetitions
      const unitRepetition = graphData.unit_repetion[yarnId];
      const rep1 = unitRepetition[0];
      const vector1 = [unitRepetition[1], unitRepetition[2]];
      const vector2 = [unitRepetition[4], unitRepetition[5]];
      const rep2 = unitRepetition[3];

      // Adjust path for starting point within pattern
      const adjustedPath = yarnPath.slice(pathStartIndex).concat(yarnPath.slice(0, pathStartIndex));
      let currentPos = [...G.get(nodeStartIndex).pos];
      let cumulativeShift = [0.0, 0.0];

      const pathZorder = [];
      for (let i = 0; i < adjustedPath.length; i++) {
        pathZorder.push(pathZ0);
        const nodeId = adjustedPath[i].toString();
        const twists = graphData.nodes[nodeId] && graphData.nodes[nodeId][2] ? graphData.nodes[nodeId][2] : 0;
        if (twists % 2 === 0) pathZ0 = -pathZ0;
      }

      const pathSave = [currentPos];
      for (let i = 0; i < adjustedPath.length; i++) {
        const nodeStart = adjustedPath[i % adjustedPath.length].toString();
        const nodeEnd = adjustedPath[(i + 1) % adjustedPath.length].toString();

        const shiftKey = `[${nodeStart}, ${nodeEnd}]`;
        const currentShift = graphData.paths[pathId].shifts[shiftKey] || [0, 0];
        cumulativeShift[0] += currentShift[0];
        cumulativeShift[1] += currentShift[1];

        const endPosBase = G.get(nodeEnd).pos;
        const endPos = [endPosBase[0] + cumulativeShift[0], endPosBase[1] + cumulativeShift[1]];
        const roundedPos = endPos.map(p => Math.round(p * 100) / 100);

        const posKey = roundedPos.toString();
        let newNodeId = posToNodeId.get(posKey);
        if (!newNodeId) {
          newNodeId = `${nodeEnd}_${roundedPos[0]}_${roundedPos[1]}`;
          G.set(newNodeId, { pos: roundedPos });
          nodeColors.set(newNodeId, 'white');
          labelDict.set(newNodeId, nodeEnd);
          posToNodeId.set(posKey, newNodeId);
        }

        // Add edge
        G.get(posToNodeId.get(currentPos.toString())).edges = G.get(posToNodeId.get(currentPos.toString())).edges || [];
        G.get(posToNodeId.get(currentPos.toString())).edges.push({
          target: newNodeId,
          color: pathColors[yarnId % pathColors.length],
        });

        // Change color of the trivial network
        if (posBase.has(currentPos.toString()) && pathZorder[i % adjustedPath.length] > 0) {
          nodeColors.set(nodeStart, pathColors[yarnId % pathColors.length]);
          const nodeTwists = graphData.nodes[nodeStart] && graphData.nodes[nodeStart][2] ? graphData.nodes[nodeStart][2] : 0;
          if (Math.abs(nodeTwists) > 0) {
            nodeColors.set(nodeStart, 'gray');
          }
        }

        pathSave.push(endPos);

        // Update current point
        currentPos = endPos;
      }
    }

    const margin = 20;

    // Scaling
    const allNodes = Array.from(G.keys());
    const nodeData = allNodes.map(nodeId => {
      return {
        id: nodeId,
        x: G.get(nodeId).pos[0],
        y: G.get(nodeId).pos[1],
        color: nodeColors.get(nodeId),
        label: labelDict.get(nodeId),
      };
    });

    const edgeData = [];
    for (const nodeId of allNodes) {
      const node = G.get(nodeId);
      if (node.edges) {
        for (const edge of node.edges) {
          edgeData.push({
            source: nodeId,
            target: edge.target,
            color: edge.color,
          });
        }
      }
    }

    const xExtent = d3.extent(nodeData, d => d.x);
    const yExtent = d3.extent(nodeData, d => d.y);

    const dataWidth = xExtent[1] - xExtent[0] + 4;
    const dataHeight = yExtent[1] - yExtent[0] + 4;

    // Compute scaleFactor to ensure square units
    const scaleFactor = Math.min(
      (graphContainer.clientWidth - 2 * margin) / dataWidth,
      (graphContainer.clientHeight - 2 * margin) / dataHeight
    );

    const xScale = d3.scaleLinear()
      .domain([xExtent[0] - 2, xExtent[1] + 2])
      .range([margin, margin + dataWidth * scaleFactor]);

    const yScale = d3.scaleLinear()
      .domain([yExtent[1] + 2, yExtent[0] - 2])
      .range([margin + dataHeight * scaleFactor, margin]); // Inverted y-axis

    // Draw grid lines first
    svg.append('g')
      .attr('class', 'grid')
      .attr('transform', `translate(0, ${margin + dataHeight * scaleFactor})`)
      .call(d3.axisBottom(xScale)
        .ticks((xExtent[1] - xExtent[0]))
        .tickSize(-dataHeight * scaleFactor)
        .tickFormat(''))
      .selectAll('line')
      .attr('stroke', 'lightgray')
      .attr('stroke-dasharray', '2,2');

    svg.append('g')
      .attr('class', 'grid')
      .attr('transform', `translate(${margin}, 0)`)
      .call(d3.axisLeft(yScale)
        .ticks((yExtent[1] - yExtent[0]))
        .tickSize(-dataWidth * scaleFactor)
        .tickFormat(''))
      .selectAll('line')
      .attr('stroke', 'lightgray')
      .attr('stroke-dasharray', '2,2');

    // Draw base edges
    svg.append('g')
      .selectAll('.base-edge')
      .data(edgeData)
      .enter()
      .append('line')
      .attr('class', 'base-edge')
      .attr('x1', d => xScale(G.get(d.source).pos[0]))
      .attr('y1', d => yScale(G.get(d.source).pos[1]))
      .attr('x2', d => xScale(G.get(d.target).pos[0]))
      .attr('y2', d => yScale(G.get(d.target).pos[1]))
      .attr('stroke', 'black')
      .attr('stroke-width', 12)
      .attr('opacity', 0.5);

    // Draw edges
    svg.append('g')
      .selectAll('.edge')
      .data(edgeData)
      .enter()
      .append('line')
      .attr('class', 'edge')
      .attr('x1', d => xScale(G.get(d.source).pos[0]))
      .attr('y1', d => yScale(G.get(d.source).pos[1]))
      .attr('x2', d => xScale(G.get(d.target).pos[0]))
      .attr('y2', d => yScale(G.get(d.target).pos[1]))
      .attr('stroke', d => d.color)
      .attr('stroke-width', 8);

    // Draw nodes
    svg.append('g')
      .selectAll('.node')
      .data(nodeData)
      .enter()
      .append('circle')
      .attr('class', 'node')
      .attr('cx', d => xScale(d.x))
      .attr('cy', d => yScale(d.y))
      .attr('r', 10)
      .attr('fill', d => d.color)
      .attr('stroke', 'black')
      .attr('stroke-width', 1);

    // Draw labels
    svg.append('g')
      .selectAll('.label')
      .data(nodeData)
      .enter()
      .append('text')
      .attr('class', 'label')
      .attr('x', d => xScale(d.x) - 6)
      .attr('y', d => yScale(d.y) + 4)
      .text(d => d.label)
      .attr('font-size', 12)
      .attr('font-weight', 'bold');

    // Axes
    const xAxis = d3.axisBottom(xScale).ticks((xExtent[1] - xExtent[0]) / 4);
    const yAxis = d3.axisLeft(yScale).ticks((yExtent[1] - yExtent[0]) / 4);

    // Draw x-axis at the bottom
    svg.append('g')
      .attr('transform', `translate(0, ${margin + dataHeight * scaleFactor})`)
      .call(xAxis);

    // Draw y-axis on the left
    svg.append('g')
      .attr('transform', `translate(${margin}, 0)`)
      .call(yAxis);
  }

  function renderFigure2() {
    const pathColors = ['#ff6666', '#6666ff', '#ffff00', '#ff66ff', '#66ff33', '#ccffb3', '#b300ff', '#33ffff'];
    const lace = [];
    const margin = 20;
  
    // Build paths
    for (const [yarnId, yarnData] of Object.entries(graphData.unit_yarns)) {
      const pathId = yarnData[0];
      const nodeStartIndex = yarnData[1].toString();
      const pathStartIndex = yarnData[2];
      const yarnPath = graphData.paths[pathId].path.map(toInt);
  
      // Unit repetitions
      const unitRepetition = graphData.unit_repetion[yarnId];
      const rep1 = unitRepetition[0];
      const vector1 = [unitRepetition[1], unitRepetition[2]];
      const vector2 = [unitRepetition[4], unitRepetition[5]];
      const rep2 = unitRepetition[3];
  
      // Adjust path for starting point within pattern
      const adjustedPath = yarnPath.slice(pathStartIndex).concat(yarnPath.slice(0, pathStartIndex));
      let currentPos = graphData.nodes[nodeStartIndex].slice(0, 2);
      let cumulativeShift = [0.0, 0.0];
  
      const pathSave = [currentPos];
      for (let i = 0; i < adjustedPath.length; i++) {
        const nodeStart = adjustedPath[i % adjustedPath.length].toString();
        const nodeEnd = adjustedPath[(i + 1) % adjustedPath.length].toString();
  
        const shiftKey = `[${nodeStart}, ${nodeEnd}]`;
        const currentShift = graphData.paths[pathId].shifts[shiftKey] || [0, 0];
        cumulativeShift[0] += currentShift[0];
        cumulativeShift[1] += currentShift[1];
  
        const endPosBase = graphData.nodes[nodeEnd].slice(0, 2);
        const endPos = [endPosBase[0] + cumulativeShift[0], endPosBase[1] + cumulativeShift[1]];
  
        pathSave.push(endPos);
  
        // Update current point
        currentPos = endPos;
      }
  
      for (let k2 = 0; k2 < rep2; k2++) {
        for (let k1 = 0; k1 < rep1; k1++) {
          const replicatedPath = pathSave.map(([x, y]) => [
            x + k1 * vector1[0] + k2 * vector2[0],
            y + k1 * vector1[1] + k2 * vector2[1],
          ]);
          lace.push({
            path: replicatedPath,
            color: pathColors[yarnId % pathColors.length],
          });
        }
      }
    }
  
    // Scaling
    const allPaths = lace.map(d => d.path);
    const allPoints = allPaths.flat();
    const xExtent = d3.extent(allPoints, d => d[0]);
    const yExtent = d3.extent(allPoints, d => d[1]);
  
    // Compute dataWidth and dataHeight based on ROI bounds
    const dataWidth = graphData.roi_bounds.x_max - graphData.roi_bounds.x_min;
    const dataHeight = graphData.roi_bounds.y_max - graphData.roi_bounds.y_min;
  
    // Compute scaleFactor to ensure square units
    const scaleFactor = Math.min(
      (graphContainer.clientWidth - 2 * margin) / dataWidth,
      (graphContainer.clientHeight - 2 * margin) / dataHeight
    );
  
    // Adjust scales to ensure square spacing and invert y-axis
    const xScale = d3.scaleLinear()
      .domain([graphData.roi_bounds.x_min, graphData.roi_bounds.x_max])
      .range([margin, margin + dataWidth * scaleFactor]);
  
    const yScale = d3.scaleLinear()
      .domain([graphData.roi_bounds.y_min, graphData.roi_bounds.y_max])
      .range([margin + dataHeight * scaleFactor, margin]); // Inverted y-axis
  
    // Create SVG and group
    const svg = d3.select('#graph-container').append('svg')
      .attr('width', graphContainer.clientWidth)
      .attr('height', graphContainer.clientHeight);
  
    const svgGroup = svg.append('g');
  
    // Define clip path after calculating scales
    svg.append('defs')
      .append('clipPath')
      .attr('id', 'clip-path')
      .append('rect')
      .attr('x', margin)
      .attr('y', margin)
      .attr('width', dataWidth * scaleFactor)
      .attr('height', dataHeight * scaleFactor);
  
    // Apply the clip path to svgGroup
    svgGroup.attr('clip-path', 'url(#clip-path)');
  
    // Draw paths
    const lineGenerator = d3.line()
      .x(d => xScale(d[0]))
      .y(d => yScale(d[1]));
  
    lace.forEach((d) => {
      svgGroup.append('path')
        .attr('d', lineGenerator(d.path))
        .attr('fill', 'none')
        .attr('stroke', 'black')
        .attr('stroke-width', 6)
        .attr('opacity', 0.5);
  
      svgGroup.append('path')
        .attr('d', lineGenerator(d.path))
        .attr('fill', 'none')
        .attr('stroke', d.color)
        .attr('stroke-width', 4);
    });
  
    // Legend
    const legendItems = Object.keys(graphData.unit_yarns).length;
    const legendWidth = 120; // Adjust as needed
    const legendHeight = legendItems * 20 + 20;
  
    // Add background rectangle behind legend
    svg.append('rect')
      .attr('x', graphContainer.clientWidth - margin - legendWidth - 20)
      .attr('y', margin - 20)
      .attr('width', legendWidth + 40)
      .attr('height', legendHeight)
      .attr('fill', 'white')
      .attr('stroke', 'black')
      .attr('stroke-width', 1)
      .attr('opacity', 0.8);
  
    const legend = svg.selectAll('.legend')
      .data(pathColors.slice(0, legendItems))
      .enter()
      .append('g')
      .attr('class', 'legend')
      .attr('transform', (d, i) => `translate(0, ${margin + i * 20})`);
  
    legend.append('rect')
      .attr('x', graphContainer.clientWidth - margin - 18)
      .attr('width', 18)
      .attr('height', 18)
      .style('fill', d => d);
  
    legend.append('text')
      .attr('x', graphContainer.clientWidth - margin - 24)
      .attr('y', margin + 9)
      .attr('dy', '.35em')
      .style('text-anchor', 'end')
      .text((d, i) => `Path ${i + 1}`);
  
    // Bring legend to front
    legend.raise();
  }
  
  function toInt(value) {
    if (typeof value === 'string') {
      return parseInt(value.replace(/[^\d]/g, ''), 10);
    }
    return value;
  }

  // Save button event (if needed)
  document.getElementById('save-button').addEventListener('click', () => {
    saveGraph();
  });

  function saveGraph() {
    // Implement saving functionality if needed
  }
});