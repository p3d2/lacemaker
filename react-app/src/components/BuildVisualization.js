// BuildVisualization.js
import React, { useRef, useEffect, useState } from "react";
import { Box } from "@mui/material";
import IconButton from "@mui/material/IconButton";
import ZoomInIcon from "@mui/icons-material/ZoomIn";
import ZoomOutIcon from "@mui/icons-material/ZoomOut";
import YoutubeSearchedForIcon from "@mui/icons-material/YoutubeSearchedFor";

import {
  calcTranslations,
  nodeAngle,
  generateYarns,
  extendPoints,
  toInt,
} from "./Generator";

const YARN_COLORS = [
  "#FF5C5C", // Soft red
  "#5C5CFF", // Soft blue
  "#FFD966", // Warm yellow
  "#FF5CFF", // Soft magenta
  "#5CFF5C", // Soft green
  "#A6FFCC", // Pastel green
  "#995CFF", // Deep purple
  "#5CFFFF", // Soft cyan
];

const LINE_WIDTHS = {
  grid: 0.05,
  yarn: 0.8,
  path: 0.2,
  nodeStroke: 0.2,
  nodeFill: 0.5,
};

const resizeCanvasAndDraw = (canvasRef, drawFn) => {
  if (!canvasRef.current) return;
  const canvas = canvasRef.current;
  const ctx = canvas.getContext("2d");
  if (!ctx) return;

  const ratio = window.devicePixelRatio || 1;
  const styleW = canvas.clientWidth;
  const styleH = canvas.clientHeight;
  canvas.width = styleW * ratio;
  canvas.height = styleH * ratio;
  ctx.scale(ratio, ratio);
  drawFn();
};

const BuildVisualization = ({ nodes, paths, unit_yarns, roi_bounds }) => {
  const canvasRef = useRef(null);
  const [transform, setTransform] = useState({
    scale: 20,
    offsetX: 80,
    offsetY: 80,
  });
  const [isPanning, setIsPanning] = useState(false);
  const [startPan, setStartPan] = useState({ x: 0, y: 0 });

  const [allTrs, setAllTrs] = useState([]);

  useEffect(() => {
    if (!nodes || !paths) return;
    const computedTrs = paths.map((p) => calcTranslations(nodes, p));
    const trsWithAngles = nodeAngle(computedTrs, paths);
    setAllTrs(trsWithAngles);
  }, [nodes, paths]);

  const draw = () => {
    const canvas = canvasRef.current;
    if (!canvas || !nodes || !paths) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const { scale, offsetX, offsetY } = transform;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    ctx.translate(offsetX, offsetY);
    ctx.scale(scale, scale);

    drawGrid(ctx);
    drawROIBounds(ctx);
    drawUnitYarns(ctx);
    drawPaths(ctx);
    drawgrayNodes(ctx);
    drawNodes(ctx);

    ctx.restore();
  };

  const drawGrid = (ctx) => {
    ctx.save();
    const gridSize = 1.0;
    const gridColor = "#ddd";
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = LINE_WIDTHS.grid;

    const width = ctx.canvas.width / transform.scale;
    const height = ctx.canvas.height / transform.scale;

    const startX =
      Math.floor(-transform.offsetX / transform.scale / gridSize) * gridSize;
    const startY =
      Math.floor(-transform.offsetY / transform.scale / gridSize) * gridSize;

    const endX = startX + width + gridSize;
    const endY = startY + height + gridSize;

    for (let x = startX; x < endX; x += gridSize) {
      ctx.beginPath();
      ctx.moveTo(x, startY);
      ctx.lineTo(x, endY);
      ctx.stroke();
    }
    for (let y = startY; y < endY; y += gridSize) {
      ctx.beginPath();
      ctx.moveTo(startX, y);
      ctx.lineTo(endX, y);
      ctx.stroke();
    }

    ctx.restore();
  };

  const drawPaths = (ctx) => {
    if (!paths || !allTrs.length) return;
    ctx.save();
    ctx.strokeStyle = "black";
    ctx.lineWidth = LINE_WIDTHS.path;
    ctx.setLineDash([0.4, 0.4]);

    paths.forEach((pathObj, i) => {
      const trs = allTrs[i];
      if (!trs || trs.length === 0) return;

      const startNodeId = String(pathObj.path[0]);
      const startNode = nodes[startNodeId];
      if (!startNode) return;

      let x = startNode.x;
      let y = startNode.y;

      ctx.beginPath();
      ctx.moveTo(x, y);

      trs.forEach(([dx, dy]) => {
        x += dx;
        y += dy;
        ctx.lineTo(x, y);
      });

      ctx.stroke();
    });

    ctx.setLineDash([]);
    ctx.restore();
  };

  const drawgrayNodes = (ctx) => {
    if (!paths || !allTrs.length) return;
    ctx.save();
    ctx.strokeStyle = "black";
    ctx.lineWidth = LINE_WIDTHS.path;

    paths.forEach((pathObj, i) => {
      const trs = allTrs[i];
      if (!trs || trs.length === 0) return;

      const startNodeId = String(pathObj.path[0]);
      const startNode = nodes[startNodeId];
      if (!startNode) return;

      let x = startNode.x;
      let y = startNode.y;

      for (let j = 1; j < pathObj.path.length; j++) {
        const [dx, dy] = trs[j - 1];
        x += dx;
        y += dy;
        drawPathPoint(ctx, x, y, toInt(pathObj.path[j]));
      }
    });

    ctx.restore();
  };

  const drawROIBounds = (ctx) => {
    if (!roi_bounds) return;
    const { x_min, x_max, y_min, y_max } = roi_bounds;
    if (
      x_min === undefined ||
      x_max === undefined ||
      y_min === undefined ||
      y_max === undefined
    )
      return;
    ctx.save();
    ctx.strokeStyle = "black";
    ctx.lineWidth = 0.2;
    ctx.beginPath();
    ctx.rect(x_min, y_min, x_max - x_min, y_max - y_min);
    ctx.stroke();
    ctx.restore();
  };

  const drawUnitYarns = (ctx) => {
    if (!unit_yarns || !Object.keys(unit_yarns).length) return;
    if (!allTrs || !allTrs.length) return;
    ctx.save();

    Object.entries(unit_yarns).forEach(([yarnKey, yarnObj], index) => {
      const {
        path_id,
        starting_node,
        start_path_index,
        z_sign,
        z_height,
        repetitions,
      } = yarnObj;
      const trs = allTrs[path_id];
      if (!trs) return;

      const rad = 0.5;
      const params = [index, starting_node, start_path_index, z_sign, z_height];
      const { points } = generateYarns(nodes, trs, params, rad);
      if (!points || points.length === 0) return;

      const colorIndex = index % YARN_COLORS.length;
      const yarnColor = YARN_COLORS[colorIndex] || "purple";
      ctx.strokeStyle = yarnColor;
      ctx.lineWidth = LINE_WIDTHS.yarn;

      const rep1 = repetitions[0];
      const ext1 = extendPoints(points, rep1.count, rep1.dx, rep1.dy);

      let previousPoints = ext1;
      const rep2 = repetitions[1];
      for (let i = 0; i < rep2.count; i++) {
        drawPointsAsLine(ctx, previousPoints);
        previousPoints = previousPoints.map(([x, y, z]) => [
          x + rep2.dx,
          y + rep2.dy,
          z,
        ]);
      }
    });

    ctx.restore();
  };

  const drawNodes = (ctx) => {
    if (!nodes) return;
    ctx.save();
    ctx.lineWidth = LINE_WIDTHS.nodeStroke;
    Object.entries(nodes).forEach(([id, nodeData]) => {
      const { x, y, twist = 0 } = nodeData;
      if (twist !== 0) {
        ctx.fillStyle = "black";
        ctx.strokeStyle = "black";
      } else {
        ctx.fillStyle = "white";
        ctx.strokeStyle = "black";
      }

      ctx.beginPath();
      ctx.arc(x, y, LINE_WIDTHS.nodeFill, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();

      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.font = "0.6px Roboto";
      ctx.fillStyle = twist !== 0 ? "white" : "black";
      ctx.fillText(id, x, y);
    });
    ctx.restore();
  };

  const drawPathPoint = (ctx, x, y, number) => {
    ctx.save();
    ctx.fillStyle = "rgba(128, 128, 128, 1.0)";
    ctx.strokeStyle = "black";
    ctx.lineWidth = 0.1;
    const radius = 0.5;

    ctx.beginPath();
    ctx.arc(x, y, radius, 0, 2 * Math.PI);
    ctx.fill();
    ctx.stroke();

    ctx.fillStyle = "white";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.font = "0.6px Roboto";
    ctx.fillText(number, x, y);
    ctx.restore();
  };

  function drawPointsAsLine(ctx, points) {
    if (!points || points.length < 2) return;
    ctx.beginPath();
    ctx.moveTo(points[0][0], points[0][1]);
    for (let i = 1; i < points.length; i++) {
      ctx.lineTo(points[i][0], points[i][1]);
    }
    ctx.stroke();
  }

  const handleMouseDown = (e) => {
    setIsPanning(true);
    setStartPan({ x: e.clientX, y: e.clientY });
  };
  const handleMouseMove = (e) => {
    if (!isPanning) return;
    setTransform((prev) => ({
      ...prev,
      offsetX: prev.offsetX + (e.clientX - startPan.x),
      offsetY: prev.offsetY + (e.clientY - startPan.y),
    }));
    setStartPan({ x: e.clientX, y: e.clientY });
  };
  const handleMouseUp = () => setIsPanning(false);

  const handleResetZoom = () => {
    setTransform({ scale: 20, offsetX: 80, offsetY: 80 });
  };
  resizeCanvasAndDraw(canvasRef, draw);
  const handleZoom = (scaleFactor) => {
    const { scale, offsetX, offsetY } = transform;
    const newScale = scale * scaleFactor;
    setTransform({
      scale: newScale,
      offsetX,
      offsetY,
    });
  };

  useEffect(() => {
    draw();
  }, [transform, nodes, paths, unit_yarns, roi_bounds, allTrs]);

  // Initial setup for canvas resolution
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    resizeCanvasAndDraw(canvasRef, draw);
  }, []);

  // Global mouse events for panning
  useEffect(() => {
    const handleWindowMouseUp = () => setIsPanning(false);
    const handleWindowMouseMove = (e) => {
      if (!isPanning) return;
      setTransform((prev) => ({
        ...prev,
        offsetX: prev.offsetX + (e.clientX - startPan.x),
        offsetY: prev.offsetY + (e.clientY - startPan.y),
      }));
      setStartPan({ x: e.clientX, y: e.clientY });
    };
    window.addEventListener("mouseup", handleWindowMouseUp);
    window.addEventListener("mousemove", handleWindowMouseMove);

    return () => {
      window.removeEventListener("mouseup", handleWindowMouseUp);
      window.removeEventListener("mousemove", handleWindowMouseMove);
    };
  }, [isPanning, startPan]);

  return (
    <Box
      sx={{
        position: "relative",
        width: "100%",
        height: "100%",
      }}
    >
      {/* Zoom Controls - position absolute, top-left corner */}
      <Box
        sx={{
          position: "absolute",
          top: "20px",
          left: "20px",
          width: "40px",
          height: "130px",
          zIndex: 10,
          display: "flex",
          flexDirection: "column",
          gap: "5px",
          backgroundColor: "rgba(255, 255, 255, 0.5)",
          padding: "5px",
          borderRadius: "5px",
          boxShadow: 3,
        }}
      >
        <IconButton onClick={() => handleZoom(1.25)} color="primary">
          <ZoomInIcon />
        </IconButton>
        <IconButton onClick={handleResetZoom} color="primary">
          <YoutubeSearchedForIcon />
        </IconButton>
        <IconButton onClick={() => handleZoom(1 / 1.25)} color="primary">
          <ZoomOutIcon />
        </IconButton>
      </Box>

      {/* Canvas - fill entire parent */}
      <Box
        sx={{
          position: "absolute",
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          overflow: "hidden",
        }}
      >
        <canvas
          ref={canvasRef}
          style={{
            width: "100%",
            height: "100%",
            cursor: isPanning ? "grabbing" : "grab",
            display: "block",
            // remove any border or fixed size
            border: "1px solid #ccc",
            backgroundColor: "#f0f0f0",
          }}
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseUp}
          aria-label="Build Visualization Canvas"
          role="img"
        />
      </Box>
    </Box>
  );
};

export default BuildVisualization;
