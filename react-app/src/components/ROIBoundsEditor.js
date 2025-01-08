import React, { useEffect, useState } from "react";
import { Box, Card, CardContent, Collapse, TextField } from "@mui/material";

const ROIBoundsEditor = ({ roiBounds, setROIBounds }) => {
  const [expanded, setExpanded] = useState(true);

  // State for each bound value as a string for input control
  const [xMinStr, setXMinStr] = useState(roiBounds.x_min || "");
  const [yMinStr, setYMinStr] = useState(roiBounds.y_min || "");
  const [zMinStr, setZMinStr] = useState(roiBounds.z_min || "");
  const [xMaxStr, setXMaxStr] = useState(roiBounds.x_max || "");
  const [yMaxStr, setYMaxStr] = useState(roiBounds.y_max || "");
  const [zMaxStr, setZMaxStr] = useState(roiBounds.z_max || "");

  // Handle blur event to update the main state
  const handleBlur = () => {
    setROIBounds({
      x_min: parseFloat(xMinStr) || 0,
      y_min: parseFloat(yMinStr) || 0,
      z_min: parseFloat(zMinStr) || 0,
      x_max: parseFloat(xMaxStr) || 0,
      y_max: parseFloat(yMaxStr) || 0,
      z_max: parseFloat(zMaxStr) || 0,
    });
  };

  useEffect(() => {
    setXMinStr(roiBounds.x_min ?? "");
    setYMinStr(roiBounds.y_min ?? "");
    setZMinStr(roiBounds.z_min ?? "");
    setXMaxStr(roiBounds.x_max ?? "");
    setYMaxStr(roiBounds.y_max ?? "");
    setZMaxStr(roiBounds.z_max ?? "");
  }, [roiBounds]);

  return (
    <Card variant="transparent" sx={{ mt: 1 }}>
      <Collapse in={expanded} timeout="auto" unmountOnExit>
        <CardContent sx={{ p: 0.5 }}>
          {/* First row: x_min, y_min, z_min */}
          <Box sx={{ display: "flex", gap: "8px", mb: 0.5 }}>
            <TextField
              size="small"
              label="x_min"
              value={xMinStr}
              onChange={(e) => setXMinStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
            <TextField
              size="small"
              label="y_min"
              value={yMinStr}
              onChange={(e) => setYMinStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
            <TextField
              size="small"
              label="z_min"
              value={zMinStr}
              onChange={(e) => setZMinStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
          </Box>

          {/* Second row: x_max, y_max, z_max */}
          <Box sx={{ display: "flex", gap: "8px" }}>
            <TextField
              size="small"
              label="x_max"
              value={xMaxStr}
              onChange={(e) => setXMaxStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
            <TextField
              size="small"
              label="y_max"
              value={yMaxStr}
              onChange={(e) => setYMaxStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
            <TextField
              size="small"
              label="z_max"
              value={zMaxStr}
              onChange={(e) => setZMaxStr(e.target.value)}
              onBlur={handleBlur}
              sx={{ width: "80px" }}
            />
          </Box>
        </CardContent>
      </Collapse>
    </Card>
  );
};

export default ROIBoundsEditor;
