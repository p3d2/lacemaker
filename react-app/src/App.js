// App.js
import React, { useState } from "react";
import Split from "react-split";
import "./App.css";
import PatternAppBar from "./components/PatternAppBar";
import ConfigurationPanel from "./components/ConfigurationPanel";
import VisualizationPanel from "./components/VisualizationPanel";
import patternData from "./data/pattern.json";
import { Box } from "@mui/material";

function App() {
  const [data, setData] = useState(patternData);

  // Handlers for PatternAppBar
  const handleNewPattern = () => {
    setData({
      nodes: {},
      paths: [],
      unit_yarns: {},
      roi_bounds: {
        x_min: 0.0,
        x_max: 128.0,
        y_min: 0.0,
        y_max: 128.0,
        z_min: -4.0,
        z_max: 4.0,
      },
    });
  };

  const handleLoadPattern = (file) => {
    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const json = JSON.parse(e.target.result);
        setData(json);
      } catch (error) {
        console.error("Invalid JSON file:", error);
      }
    };
    reader.readAsText(file);
  };

  const handleSavePattern = () => {
    const blob = new Blob([JSON.stringify(data, null, 2)], {
      type: "application/json",
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "pattern.json";
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        height: "100vh",
      }}
    >
      {/* Top App Bar */}
      <PatternAppBar
        onNewPattern={handleNewPattern}
        onLoadPattern={handleLoadPattern}
        onSavePattern={handleSavePattern}
      />

      {/* Main content area */}
      <Box sx={{ flex: 1, overflow: "hidden" }}>
        <Split
          sizes={[30, 70]} // 30% left, 70% right
          minSize={[400, 400]}
          gutterSize={10}
          direction="horizontal"
          style={{
            display: "flex",
            width: "100%",
            height: "100%", // critical for each panel to have a bounding box
          }}
        >
          {/* Left Panel */}
          <Box
            sx={{
              overflow: "auto", // separate scrollbar for left panel
              borderRight: "1px solid #ccc",
              p: 2,
            }}
          >
            <ConfigurationPanel data={data} setData={setData} />
          </Box>

          {/* Right Panel */}
          <Box
            sx={{
              overflow: "auto", // separate scrollbar for right panel
              borderRight: "5px solid #ccc",
              p: 2,
            }}
          >
            <VisualizationPanel data={data} />
          </Box>
        </Split>
      </Box>
    </Box>
  );
}

export default App;
