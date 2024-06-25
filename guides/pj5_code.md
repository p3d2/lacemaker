let availableLabels = [];
let nextLabel = 0;  // Still needed to handle new labels if the pool is empty
let crossingLabels = {};  // Dictionary to store labels with coordinates as keys
let yarnInputs = [];
let yarnPaths = Array(8).fill([]); // Store parsed paths for each yarn
let grid;
let cols, rows;
let cellSize = 40;
let gridSize = 8;
let canvasSize = 600;
let gridSizeInput, updateGridSizeButton;
let gridType = 'Orthogonal';
let gridTypes = ['Orthogonal', 'Hexagonal', 'Triangular', 'Non-Orthogonal', 'Curvilinear'];
let currentSelection = 0; // Index of the selected grid type in the array
let currentMenu = 0;

function setup() {
    createCanvas(canvasSize, canvasSize + 200);
    cols = rows = gridSize;
    cellSize = canvasSize / cols;
    grid = new Array(cols).fill().map(() => new Array(rows).fill(false));

    // Create input field and button for grid size update
    gridSizeInput = createInput(gridSize.toString());
    gridSizeInput.position(20, canvasSize + 90);
    updateGridSizeButton = createButton('Update Grid Size');
    updateGridSizeButton.position(gridSizeInput.x + gridSizeInput.width + 10, canvasSize + 90);
    updateGridSizeButton.mousePressed(updateGridSize);
  
    if (currentMenu !== 0) {
        gridSizeInput.hide();
        updateGridSizeButton.hide();
    }
  
    // Setup navigation buttons
    button1 = createButton('Grid Selection');
    button1.position(canvasSize + 20, 30);
    button1.mousePressed(() => changeMenu(0));

    button2 = createButton('Define Crossings');
    button2.position(canvasSize + 20, 60);
    button2.mousePressed(() => changeMenu(1));

    // Setup yarn input fields
    for (let i = 0; i < 8; i++) {
        let input2 = createInput('');
        input2.attribute('placeholder', 'Enter path like 0, 1, ...');
        input2.position(canvasSize + 20, 90 + (i+1) * 20);
        yarnInputs.push(input2);
    }

    button3 = createButton('Define Yarns');
    button3.position(canvasSize + 20, 90);
    button3.mousePressed(() => changeMenu(2)); // Switch to the yarn definition menu
}

function updateGridSize() {
    let newGridSize = parseInt(gridSizeInput.value());
    if (!isNaN(newGridSize) && newGridSize > 0) {
        gridSize = newGridSize;
        cols = rows = gridSize;
        cellSize = canvasSize / cols;
        grid = new Array(cols).fill().map(() => new Array(rows).fill(false)); // Reset grid
        redraw(); // Redraw the canvas with the new grid size
    } else {
        console.error('Invalid grid size entered');
    }
}

function draw() {
    background(255);
    drawGrid(gridType);
    if (currentMenu === 0) {
        drawUI();
        gridSizeInput.show();
        updateGridSizeButton.show();
    } else if (currentMenu === 2) {
        drawYarnPaths();
    } else {
        gridSizeInput.hide();
        updateGridSizeButton.hide();
    }
}

function changeMenu(menuId) {
    currentMenu = menuId;
}

function drawUI() {
    // Draw dropdown for grid types
    fill(200);
    rect(0, canvasSize, canvasSize, 40);
    fill(0);
    text("Grid Type: " + gridType, 10, canvasSize + 30);

    // Draw buttons for grid type selection
    gridTypes.forEach((type, index) => {
        fill(index === currentSelection ? 'lightgreen' : 'lightgray');
        rect(10 + 100 * index, canvasSize + 50, 90, 30);
        fill(0);
        text(type, 15 + 100 * index, canvasSize + 70);
    });
}

function drawGrid(type) {
    switch (type) {
        case 'Orthogonal':
            drawOrthogonalGrid();
            break;
        case 'Hexagonal':
            drawHexagonalGrid();
            break;
        case 'Triangular':
            drawTriangularGrid();
            break;
        case 'Non-Orthogonal':
            drawNonOrthogonalGrid();
            break;
        case 'Curvilinear': // Placeholder for future implementation
            break;
    }
}

function toggleCrossing(i, j) {
    let key = `${i}_${j}`;
    if (grid[i][j]) {  // If already a crossing
        grid[i][j] = false;
        availableLabels.push(crossingLabels[key]);  // Add the label back to the pool
        delete crossingLabels[key];  // Remove the label
    } else {
        grid[i][j] = true;
        let label;
        if (availableLabels.length > 0) {
            label = availableLabels.shift();  // Reuse an old label
        } else {
            label = nextLabel++;  // No available labels, use the next new label
        }
        crossingLabels[key] = label;  // Assign the label
    }
}

function drawOrthogonalGrid() {
    textSize(12); // Adjust font size as needed
    textAlign(LEFT, CENTER); // Align text to be centered
    for (let i = 0; i < cols; i++) {
        for (let j = 0; j < rows; j++) {
            stroke(0);
            fill(grid[i][j] ? 'skyblue' : 240); // Ensure good contrast
            let x = i * cellSize;
            let y = j * cellSize;
            rect(x, y, cellSize, cellSize);
            if (grid[i][j]) {
                let labelKey = `${i}_${j}`;
                if (crossingLabels[labelKey] !== undefined) {
                    fill(0); // Text color
                    noStroke();
                    text(crossingLabels[labelKey], x + cellSize / 2, y + cellSize / 2);
                }
            }
        }
    }
}


function drawHexagonalGrid() {
    // Adjusting horizontal spacing and vertical spacing to fit hexagons perfectly
    let horizontalSpacing = cellSize * 3/4;
    let verticalSpacing = sqrt(3) / 2 * cellSize;

    for (let i = 0; i < cols; i++) {
        for (let j = 0; j < rows; j++) {
            let x = i * horizontalSpacing;
            let y = (j+1) * verticalSpacing + (i % 2) * verticalSpacing / 2;
            fill(grid[i][j] ? 'skyblue' : 240);
            stroke(0);
            drawHexagon(x + cellSize / 2, y, cellSize / 2); // Ensure the center is properly offset
        }
    }
}

function drawHexagon(x, y, radius) {
    beginShape();
    for (let angle = 0; angle < TWO_PI; angle += TWO_PI / 6) {
        let sx = x + cos(angle) * radius;
        let sy = y + sin(angle) * radius;
        vertex(sx, sy);
    }
    endShape(CLOSE);
}

function drawYarnPaths() {
    yarnInputs.forEach((input, index) => {
        let path = input.value().split(',').map(Number).filter(num => Object.values(crossingLabels).includes(num));
        drawYarn(path, index);
    });
}

function drawYarn(path, yarnIndex) {
    let colors = ['red', 'green', 'blue', 'yellow', 'purple', 'orange', 'brown', 'pink'];  // Color for each yarn
    push(); 
    beginShape();
    stroke(colors[yarnIndex % colors.length]);
    noFill();  // Ensure no fill is applied to the shapes
    strokeWeight(3);
    path.forEach(num => {
        let key = Object.keys(crossingLabels).find(key => crossingLabels[key] === num);
        if (key) {
            let [i, j] = key.split('_').map(Number);
            let x = i * cellSize + cellSize / 2;
            let y = j * cellSize + cellSize / 2;
            vertex(x, y);
        }
    });
    endShape();
    pop();
}

function mousePressed() {
    if (currentMenu === 0 && mouseY > canvasSize) {
        if (mouseY > canvasSize + 50 && mouseY < canvasSize + 80) {
            let index = Math.floor((mouseX - 10) / 100);
            if (index >= 0 && index < gridTypes.length) {
                currentSelection = index;
                gridType = gridTypes[index];
                redraw();
            }
        }
    } else if (currentMenu === 1) {
        // Handling crossings definition
        let i = Math.floor(mouseX / cellSize);
        let j = Math.floor(mouseY / cellSize);
        if (i >= 0 && i < cols && j >= 0 && j < rows) {
            toggleCrossing(i, j);  // Toggle crossing on or off
            redraw();  // Redraw the canvas to reflect changes
        }
    }
}
