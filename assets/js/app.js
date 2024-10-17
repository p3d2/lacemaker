// assets/js/app.js

document.addEventListener('DOMContentLoaded', () => {
  const fileSelect = document.getElementById('file-select');
  const graphContainer = document.getElementById('graph-container');
  let graphData;
  let simulation;
  let svg;

  // Fetch list of JSON files from the input/json_patterns/ directory
  // Since we can't read the directory contents directly, we'll create a JSON file that lists the available JSON files
  // Alternatively, if the list of files is static or known, you can hardcode it
  const files = [
    'pattern1024.json',
    'pattern1084.json',
    // Add more file names as needed
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
    d3.select('svg').remove();

    fetch(`${dataPath}${fileName}`)
      .then(response => response.json())
      .then(data => {
        graphData = data;
        renderGraph();
      })
      .catch(error => {
        console.error('Error loading JSON file:', error);
      });
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

    // Initialize simulation
    simulation = d3.forceSimulation(graphData.nodes)
      .force('link', d3.forceLink(graphData.links).id(d => d.id).distance(100))
      .force('charge', d3.forceManyBody().strength(-300))
      .force('center', d3.forceCenter(width / 2, height / 2));

    // Add links
    const link = svg.append('g')
      .attr('class', 'links')
      .selectAll('line')
      .data(graphData.links)
      .enter().append('line')
      .attr('class', 'link');

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
      .text(d => d.id);

    // Add node double-click event to edit node properties
    node.on('dblclick', (event, d) => {
      const newId = prompt('Enter new node ID:', d.id);
      if (newId !== null && newId !== d.id) {
        d.id = newId;
        d3.select(event.currentTarget).select('text').text(newId);
        simulation.nodes(graphData.nodes);
        simulation.alpha(0.3).restart();
      }
    });

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
        id: node.id,
        x: node.x,
        y: node.y,
        // Include other node properties if any
      };
    });

    const links = graphData.links.map(link => {
      return {
        source: link.source.id ? link.source.id : link.source,
        target: link.target.id ? link.target.id : link.target,
        // Include other link properties if any
      };
    });

    const updatedData = {
      ...graphData,
      nodes,
      links,
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