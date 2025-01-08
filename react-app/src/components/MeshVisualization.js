// MeshVisualization.js
import React, { useRef, useEffect, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";

/**
 * Splits an array of points into multiple subarrays at the given break indices.
 * E.g. if points=[p0,p1,p2,p3,p4,p5], breakIndices=[2,4],
 * then you get [[p0,p1,p2],[p3,p4],[p5]].
 */
function isWebGLAvailable() {
  try {
    const canvas = document.createElement("canvas");
    return !!(
      window.WebGLRenderingContext &&
      (canvas.getContext("webgl") || canvas.getContext("experimental-webgl"))
    );
  } catch (e) {
    return false;
  }
}

function splitIntoSegments(points, breakIndices) {
  const segments = [];
  let lastIndex = 0;
  breakIndices.forEach((breakAt) => {
    // slice from lastIndex to breakAt (inclusive)
    const segment = points.slice(lastIndex, breakAt + 1);
    segments.push(segment);
    lastIndex = breakAt + 1;
  });
  // The last slice (whatever remains)
  if (lastIndex < points.length) {
    segments.push(points.slice(lastIndex));
  }
  return segments;
}

function autoBreakPoints(points, threshold = 2.0) {
  const breakIndices = [];
  for (let i = 0; i < points.length - 1; i++) {
    const [x1, y1, z1] = points[i];
    const [x2, y2, z2] = points[i + 1];
    const dist = Math.sqrt(
      Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2) + Math.pow(z2 - z1, 2)
    );
    if (dist > threshold) {
      breakIndices.push(i); // break at i
    }
  }
  return breakIndices;
}

class YarnCurve extends THREE.Curve {
  constructor(points) {
    super();
    // Store the points as Three.js Vector3 objects
    this.points3D = points.map(([x, y, z]) => new THREE.Vector3(x, y, z));
    this.segmentCount = this.points3D.length - 1;
  }

  getPoint(t) {
    // t is between 0 and 1
    // Scale it to [0, segmentCount]
    const total = this.segmentCount;
    const scaledT = t * total;

    // Find the two segment endpoints
    const i = Math.floor(scaledT);
    const alpha = scaledT - i;

    // Clamp the index to avoid out-of-bounds
    if (i >= total) return this.points3D[total];
    if (i < 0) return this.points3D[0];

    // Get the start/end points of the current segment
    const start = this.points3D[i];
    const end = this.points3D[i + 1];

    // Linear interpolate
    const x = start.x + alpha * (end.x - start.x);
    const y = start.y + alpha * (end.y - start.y);
    const z = start.z + alpha * (end.z - start.z);

    return new THREE.Vector3(x, y, z);
  }
}

// Import your logic from Generator.js.
// Adjust paths if needed for your setup:
import {
  calcTranslations,
  nodeAngle,
  generateYarns,
  extendPoints,
  smoothYarnPoints,
} from "./Generator";

/**
 * Optional color palette for yarns
 */
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

const MeshVisualization = ({ nodes, paths, unit_yarns, roi_bounds }) => {
  // This mirrors BuildVisualization's approach to storing translations:
  const [allTrs, setAllTrs] = useState([]);
  const mountRef = useRef(null);

  /**
   * Compute allTrs once nodes/paths are available
   * (Same as in BuildVisualization: calcTranslations -> nodeAngle)
   */
  useEffect(() => {
    if (!nodes || !paths) return;

    // 1) Calculate raw translations for each path
    const computedTrs = paths.map((p) => calcTranslations(nodes, p));
    // 2) Assign angles for special 'l'/'r' nodes
    const trsWithAngles = nodeAngle(computedTrs, paths);
    setAllTrs(trsWithAngles);
  }, [nodes, paths]);

  /**
   * Main Three.js setup
   */
  useEffect(() => {
    if (!unit_yarns || !Object.keys(unit_yarns).length) return;
    if (!allTrs || !allTrs.length) return; // same condition as BuildVisualization
    if (!mountRef.current) return;

    if (!isWebGLAvailable()) {
      console.error("WebGL is not supported in this browser or device.");
      // Optionally display an error message to the user
      return;
    }

    // 1) Create scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xffffff);

    // 2) Camera
    const width = mountRef.current.clientWidth;
    const height = mountRef.current.clientHeight;
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 2000);

    // 3) Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    mountRef.current.appendChild(renderer.domElement);

    // 4) Controls
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    // 5) For each yarn, replicate the BuildVisualization logic
    Object.entries(unit_yarns).forEach(([yarnKey, yarnObj], index) => {
      const {
        path_id,
        starting_node,
        start_path_index,
        z_sign,
        z_height,
        repetitions,
      } = yarnObj;

      // Make sure path translations exist
      const trs = allTrs[path_id];
      if (!trs) return;

      // Generate base points
      const rad = 0.5; // same as in BuildVisualization
      const params = [index, starting_node, start_path_index, z_sign, z_height];
      let { points: basePoints } = generateYarns(nodes, trs, params, rad);

      if (!basePoints || basePoints.length === 0) return;
      const colorIndex = index % YARN_COLORS.length;
      const yarnColor = YARN_COLORS[colorIndex] || "#FF00FF";
      let finalPoints = [];

      // Now you can smooth in-place before repeating:
      const smoothed = smoothYarnPoints(basePoints, 0.5, 0.1, 1000);
      basePoints = smoothed;

      // If your data has 2 repetitions:
      const rep1 = repetitions[0];
      if (rep1) {
        // e.g. { count, dx, dy }
        const ext1 = extendPoints(basePoints, rep1.count + 1, rep1.dx, rep1.dy);

        // We'll feed these into rep2 logic
        let previousPoints = ext1;
        const rep2 = repetitions[1];
        if (rep2) {
          // We build up finalPoints by concatenating each iteration
          for (let i = 0; i < rep2.count; i++) {
            // "Drawing" in 2D was drawPointsAsLine, but here we just accumulate them
            finalPoints = finalPoints.concat(previousPoints);

            // Shift them for the next iteration
            previousPoints = previousPoints.map(([x, y, z]) => [
              x + rep2.dx,
              y + rep2.dy,
              z,
            ]);
          }
          // Add the final shift
          finalPoints = finalPoints.concat(previousPoints);
        }
      }

      // (Optional) Filter points outside the ROI, if needed
      if (roi_bounds) {
        const { x_min, x_max, y_min, y_max, z_min, z_max } = roi_bounds;
        finalPoints = finalPoints.filter(([xx, yy, zz]) => {
          return (
            xx >= x_min &&
            xx <= x_max &&
            yy >= y_min &&
            yy <= y_max &&
            zz >= z_min &&
            zz <= z_max
          );
        });
      }

      const threshold = 20.0; // e.g. if distance > 2 means a gap
      const breakIndices = autoBreakPoints(finalPoints, threshold);
      const segments = splitIntoSegments(finalPoints, breakIndices);
      // If finalPoints is empty, skip
      if (!finalPoints.length) return;

      segments.forEach((segmentPoints) => {
        // If the segment is too short (e.g. only 1 point), skip
        if (segmentPoints.length < 2) return;
        //console.log("Segment points:", finalPoints);

        // Create a YarnCurve from just this segment
        const curve = new YarnCurve(segmentPoints);
        const tubeMaterial = new THREE.MeshBasicMaterial({
          color: yarnColor,
          side: THREE.DoubleSide,
        });
        // Build a TubeGeometry
        const tubeGeometry = new THREE.TubeGeometry(
          curve,
          finalPoints.length, // tubularSegments
          0.5, // radius
          8, // radialSegments
          false // not closed
        );
        const tubeMesh = new THREE.Mesh(tubeGeometry, tubeMaterial);
        scene.add(tubeMesh);
      });
    });

    // 8) If you want to visualize the ROI bounds as a box or lines:
    if (roi_bounds) {
      const { x_min, x_max, y_min, y_max, z_min = 0, z_max = 0 } = roi_bounds;
      const boxWidth = x_max - x_min;
      const boxHeight = y_max - y_min;
      const boxDepth = Math.max(z_max - z_min, 0.1); // at least some depth
      const boxGeom = new THREE.BoxGeometry(boxWidth, boxHeight, boxDepth);
      const boxEdges = new THREE.EdgesGeometry(boxGeom);
      const boxMat = new THREE.LineBasicMaterial({ color: 0x000000 });
      const box = new THREE.LineSegments(boxEdges, boxMat);
      // Center it
      const cx = (x_min + x_max) / 2;
      const cy = (y_min + y_max) / 2;
      const cz = (z_min + z_max) / 2;
      box.position.set(cx, cy, cz);
      scene.add(box);

      // Position the camera slightly “above” that center (pick any offset you like)
      camera.position.set(cx, cy, cz + 300);

      // Aim the camera at the center
      camera.lookAt(cx, cy, cz);

      // Also update orbit controls so they rotate around the ROI center
      controls.target.set(cx, cy, cz);
      controls.update();
    }

    // 9) Basic lights (optional)
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.7);
    scene.add(ambientLight);
    const dirLight = new THREE.DirectionalLight(0xffffff, 0.4);
    dirLight.position.set(0, 0, 100);
    scene.add(dirLight);

    // 10) Animation loop
    const animate = () => {
      requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    };
    animate();

    // 11) Handle resize
    const handleResize = () => {
      if (!mountRef.current) return;
      const newW = mountRef.current.clientWidth;
      const newH = mountRef.current.clientHeight;
      renderer.setSize(newW, newH);
      camera.aspect = newW / newH;
      camera.updateProjectionMatrix();
    };
    window.addEventListener("resize", handleResize);

    // Cleanup
    return () => {
      window.removeEventListener("resize", handleResize);
      if (mountRef.current && renderer.domElement.parentNode) {
        mountRef.current.removeChild(renderer.domElement);
      }
      renderer.dispose();
    };
  }, [allTrs, unit_yarns, roi_bounds, nodes]);

  return (
    <div
      ref={mountRef}
      style={{
        width: "100%",
        height: "100%",
        border: "1px solid #ccc",
      }}
    />
  );
};

export default MeshVisualization;
