document.addEventListener('DOMContentLoaded', () => {
  const fileSelect = document.getElementById('file-select');
  const graphContainer = document.getElementById('graph-container');
  let graphData;
  let simulation;
  let svg;

  // List of JSON files (update with your actual file names)
  const files = [
    'pattern1024.json',
    'pattern1084.json'
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
        // Process data according to your JSON structure
        graphData = convertDataToGraphFormat(data);
        renderGraph();
      })
      .catch(error => {
        console.error('Error loading JSON file:', error);
        alert('Error loading JSON file. Check the console for details.');
      });
  }

  function convertDataToGraphFormat(data) {
    const nodes = [];
    const links = [];
    const nodeMap = {}; // Map from node ID to node object

    // Parse nodes
    for (const [nodeId, pos] of Object.entries(data.nodes)) {
      const node = {
        id: nodeId,
        x: pos[0],
        y: pos[1],
        // Handle optional third element
        extra: pos[2] || null,
      };
      nodes.push(node);
      nodeMap[nodeId] = node;
    }

    // Parse paths and shifts
    data.paths.forEach((pathObj, pathIndex) => {
      const path = pathObj.path;
      const shifts = pathObj.shifts || {};
      let cumulativeShift = [0, 0];

      for (let i = 0; i < path.length - 1; i++) {
        let sourceId = path[i];
        let targetId = path[i + 1];

        // Handle modified node IDs (e.g., "1r", "8l")
        sourceId = cleanNodeId(sourceId);
        targetId = cleanNodeId(targetId);

        // Apply shifts if any
        const shiftKey = `[${path[i]}, ${path[i + 1]}]`;
        if (shifts[shiftKey]) {
          const shift = shifts[shiftKey];
          cumulativeShift[0] += shift[0];
          cumulativeShift[1] += shift[1];
        }

        // Get source and target nodes
        const sourceNode = nodeMap[sourceId];
        const targetNode = nodeMap[targetId];

        // Clone nodes if necessary to account for shifts
        const shiftedTargetNode = getShiftedNode(targetNode, cumulativeShift);

        // Add shifted node to nodes array if not already present
        if (!nodes.includes(shiftedTargetNode)) {
          nodes.push(shiftedTargetNode);
        }

        // Create link
        links.push({
          source: sourceNode,
          target: shiftedTargetNode,
          pathIndex: pathIndex,
        });
      }
    });

    return { nodes, links };
  }

  function cleanNodeId(nodeId) {
    if (typeof nodeId === 'string') {
      // Remove any suffixes like 'r' or 'l'
      return nodeId.replace(/[rl]$/, '');
    }
    return nodeId.toString();
  }

  const shiftedNodeCache = {};

  function getShiftedNode(node, shift) {
    // Create a unique key for the shifted node
    const key = `${node.id}_${shift[0]}_${shift[1]}`;
    if (shiftedNodeCache[key]) {
      return shiftedNodeCache[key];
    } else {
      const shiftedNode = {
        id: key,
        originalId: node.id,
        x: node.x + shift[0],
        y: node.y + shift[1],
        extra: node.extra,
      };
      shiftedNodeCache[key] = shiftedNode;
      return shiftedNode;
    }
  }

  function renderGraph() {
    const width = graphContainer.clientWidth;
    const height = graphContainer.clientHeight;

    svg = d3.select('#graph-container').append('svg')
      .attr('width', width)
      .attr('height', height)
      .call(d3.zoom().on('zoom', (event) => {
        svg.attr('transform', event.transform);
      }))
      .append('g');

    // Define color scale for paths
    const color = d3.scaleOrdinal(d3.schemeCategory10);

    // Add links
    const link = svg.append('g')
      .attr('class', 'links')
      .selectAll('line')
      .data(graphData.links)
      .enter().append('line')
      .attr('class', 'link')
      .attr('stroke', d => color(d.pathIndex));

    // Add nodes
    const node = svg.append('g')
      .attr('class', 'nodes')
      .selectAll('g')
      .data(graphData.nodes)
      .enter().append('g')
      .call(d3.drag()
        .on('start', dragStarted)
        .on('drag', dragged)
        .on('end', dragEnded));

    node.append('circle')
      .attr('r', 10)
      .attr('fill', '#69b3a2');

    node.append('text')
      .attr('dx', 12)
      .attr('dy', '.35em')
      .text(d => d.originalId || d.id);

    // Add node double-click event to edit node properties
    node.on('dblclick', (event, d) => {
      const newId = prompt('Enter new node ID:', d.originalId || d.id);
      if (newId !== null && newId !== d.id) {
        d.originalId = newId;
        d3.select(event.currentTarget).select('text').text(newId);
        // Update nodeMap if necessary
      }
    });

    // Simulation setup
    simulation = d3.forceSimulation(graphData.nodes)
      .force('link', d3.forceLink(graphData.links).id(d => d.id).distance(100))
      .force('charge', d3.forceManyBody().strength(-300))
      .force('center', d3.forceCenter(width / 2, height / 2));

    // Simulation tick
    simulation.on('tick', () => {
      link
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y);

      node.attr('transform', d => `translate(${d.x},${d.y})`);
    });
  }

  // Drag functions
  function dragStarted(event, d) {
    if (!event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
  }

  function dragged(event, d) {
    d.fx = event.x;
    d.fy = event.y;
  }

  function dragEnded(event, d) {
    if (!event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
  }

  // Save button event
  document.getElementById('save-button').addEventListener('click', () => {
    saveGraph();
  });

  function saveGraph() {
    // Prepare the data to be saved
    const nodes = graphData.nodes.map(node => {
      return {
        id: node.originalId || node.id,
        x: node.x,
        y: node.y,
        extra: node.extra,
        // Include other node properties if any
      };
    });

    const links = graphData.links.map(link => {
      return {
        source: link.source.originalId || link.source.id,
        target: link.target.originalId || link.target.id,
        // Include other link properties if any
      };
    });

    const updatedData = {
      // Include the original data structure
      ...graphData.originalData,
      nodes: nodes.reduce((acc, node) => {
        const id = node.id;
        const posArray = [node.x, node.y];
        if (node.extra !== null) {
          posArray.push(node.extra);
        }
        acc[id] = posArray;
        return acc;
      }, {}),
      // Paths and shifts would need to be updated accordingly
      // For simplicity, we're not updating shifts here
    };

    const jsonStr = JSON.stringify(updatedData, null, 2);
    const blob = new Blob([jsonStr], { type: 'application/json' });
    const url = URL.createObjectURL(blob);

    // Create a link and trigger download
    const a = document.createElement('a');
    a.href = url;
    a.download = `updated_${fileSelect.value}`;
    a.click();
    URL.revokeObjectURL(url);
  }
});
