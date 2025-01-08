// VisualizationPanel.js
import React from "react";
import { Box, Button } from "@mui/material";
import Generator from "./Generator";
import BuildVisualization from "./BuildVisualization";
import MeshVisualization from "./MeshVisualization";
import FlowOverviewPanel from "./FlowOverviewPanel"; // Import the new 3rd panel

const VisualizationPanel = ({ data }) => {
  const [viewMode, setViewMode] = React.useState("build");

  // Example placeholders for job management

  const [logs, setLogs] = React.useState("");

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        width: "100%",
        height: "100vh",
        flex: 1,
      }}
    >
      {/* Top bar with 3 buttons now */}
      <Box
        display="flex"
        justifyContent="space-between"
        alignItems="center"
        mb={2}
        sx={{ flexShrink: 0 }}
      >
        <Box>
          <Button
            variant={viewMode === "build" ? "contained" : "outlined"}
            color="primary"
            onClick={() => setViewMode("build")}
            sx={{ mr: 1 }}
          >
            Build Mode
          </Button>
          <Button
            variant={viewMode === "mesh" ? "contained" : "outlined"}
            color="primary"
            onClick={() => setViewMode("mesh")}
            sx={{ mr: 1 }}
          >
            Mesh Mode
          </Button>
          <Button
            variant={viewMode === "flow" ? "contained" : "outlined"}
            color="primary"
            onClick={() => setViewMode("flow")}
          >
            Jobs
          </Button>
        </Box>
      </Box>

      {/* Main content area */}

      {/* If in Flow mode, skip generator & show FlowOverviewPanel */}
      {viewMode === "flow" && (
        <FlowOverviewPanel
          logs={logs} // or omit if you don't track logs
          setLogs={setLogs} // optional
        />
      )}

      {viewMode !== "flow" && (
        <Box
          display="flex"
          justifyContent="center"
          alignItems="center"
          sx={{
            flex: 1, // fill leftover vertical space
            border: "1px solid #ccc",
            borderRadius: 1,
            position: "relative", // helps if child uses absolute positioning
          }}
        >
          <Generator config={data}>
            {({ nodes, paths, unit_yarns, roi_bounds }) => (
              <>
                {viewMode === "build" && (
                  <BuildVisualization
                    nodes={nodes}
                    paths={paths}
                    unit_yarns={unit_yarns}
                    roi_bounds={roi_bounds}
                  />
                )}
                {viewMode === "mesh" && (
                  <MeshVisualization
                    nodes={nodes}
                    paths={paths}
                    unit_yarns={unit_yarns}
                    roi_bounds={roi_bounds}
                  />
                )}
              </>
            )}
          </Generator>
        </Box>
      )}
    </Box>
  );
};

export default VisualizationPanel;
