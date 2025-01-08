// Generator.js
// Provides logic to process node and path data (calcTranslations, nodeAngle),
// compute angles, generate yarn points (generateYarns), and extend them (extendPoints).

/**
 * Converts a string like "12" or "12r" into an integer (12).
 * If it ends with 'l' or 'r', we remove that char and parse the prefix.
 */

import React from "react";
import CubicSpline from "cubic-spline";

const Generator = ({ config, children }) => {
  // Processing logic
  const { nodes, paths, unit_yarns, roi_bounds } = config;

  return children({ nodes, paths, unit_yarns, roi_bounds });
};

export default Generator;

export function toInt(value) {
  const val = String(value);
  if (val.endsWith("l") || val.endsWith("r")) {
    return parseInt(val.slice(0, -1), 10);
  }
  return parseInt(val, 10);
}

/**
 * Returns +1 if the string ends with 'r', -1 if it ends with 'l', else 0.
 */
export function lr(value) {
  const val = String(value);
  const lastChar = val.charAt(val.length - 1);
  return { l: -1, r: 1 }[lastChar] || 0;
}

/**
 * Generate the local 'twist' shape points for an integer twist value 'twVal',
 * rotated by 'ang'.
 */
export function twistPoints(twVal, ang) {
  const points = [];
  const twist = Math.abs(twVal);

  // Shift angle by -π/4
  const theta = ang - Math.PI / 4;

  let x = -twist;
  let y = -twist;

  for (let k = 0; k < 2 * twist + 1; k++) {
    points.push([x, y]);

    // This matches logic from Python script: int(((k - sign) % 4) / 2) < 1
    const sign = Math.sign(twVal); // +1 or -1 for positive/negative twist
    const modValue = (((k - sign) % 4) + 4) % 4; // positive modulo
    const condition = Math.floor(modValue / 2) < 1;

    if (condition) {
      y += 2;
    } else {
      x += 2;
    }
  }

  // Rotate all generated points by theta
  const rotatedPoints = points.map(([px, py]) => {
    const rx = px * Math.cos(theta) - py * Math.sin(theta);
    const ry = px * Math.sin(theta) + py * Math.cos(theta);
    return [rx, ry];
  });

  return rotatedPoints;
}

/**
 * Given your 'nodes' object and a single path (with possible special shifts),
 * compute the list of translations [dx, dy, twistValue, angle] for each step.
 *
 * @param {Object} nodes  e.g. { "0": {x:0, y:0, twist:0}, "1": {...}, ... }
 * @param {Object} pathObj e.g. { path: [...], shifts: [{ from, to, dx, dy }... ] }
 */
export function calcTranslations(nodes, pathObj) {
  const { path, shifts } = pathObj;

  // Map of "from,to" => [dx, dy], for special user-defined shifts
  const shiftMap = new Map();
  shifts.forEach((shift) => {
    const key = `${shift.from},${shift.to}`;
    shiftMap.set(key, [shift.dx, shift.dy]);
  });

  const spatialShifts = [];

  for (let i = 0; i < path.length - 1; i++) {
    const currentCrossing = toInt(path[i]);
    const nextCrossing = toInt(path[i + 1]);

    // Skip if node is missing
    if (!nodes[String(nextCrossing)] || !nodes[String(currentCrossing)]) {
      continue;
    }

    // Base shift = difference in x,y
    const shiftX =
      nodes[String(nextCrossing)].x - nodes[String(currentCrossing)].x;
    const shiftY =
      nodes[String(nextCrossing)].y - nodes[String(currentCrossing)].y;

    let finalShiftX = shiftX;
    let finalShiftY = shiftY;

    // If there's a custom shift for (current -> next), add it
    const pairShiftKey = `${currentCrossing},${nextCrossing}`;
    if (shiftMap.has(pairShiftKey)) {
      const [dx, dy] = shiftMap.get(pairShiftKey);
      finalShiftX += dx;
      finalShiftY += dy;
    }

    // Twist value depends on left/right marker
    const twVal = lr(path[i + 1]) * (nodes[String(nextCrossing)].twist || 0);
    // We store [dx, dy, twist, angle=0], angle might be assigned later in nodeAngle
    spatialShifts.push([finalShiftX, finalShiftY, twVal, 0]);
  }

  return spatialShifts;
}

/**
 * Processes "special" nodes (ending in 'l' or 'r') to assign angles in the
 * path translations (the 4th element in each translation array).
 *
 * If a node has 'r' or 'l' in more than one place, or not exactly a pair,
 * we skip it with a warning.
 */
export function nodeAngle(trs, paths) {
  const specialNodeIndices = {};
  const seenNodes = new Set();

  // 1) Collect indices where 'r' or 'l' appear
  paths.forEach((pathObj, pathIndex) => {
    const path = pathObj.path;
    path.forEach((node, nodeIndex) => {
      if (
        typeof node === "string" &&
        (node.endsWith("l") || node.endsWith("r"))
      ) {
        if (seenNodes.has(node)) {
          console.warn(
            `Duplicate special node '${node}' found in multiple paths. Skipping duplicate.`
          );
          return; // Optionally, skip further processing for this node
        }
        seenNodes.add(node);

        const numericKey = node.slice(0, -1); // remove last char
        if (!specialNodeIndices[numericKey]) {
          specialNodeIndices[numericKey] = [];
        }
        specialNodeIndices[numericKey].push([pathIndex, nodeIndex]);
      }
    });
  });

  // 2) Check that each special node index has exactly 2 references
  Object.entries(specialNodeIndices).forEach(([key, indices]) => {
    if (indices.length !== 2) {
      console.warn(
        `Skipping special node '${key}' because it doesn't have exactly one 'l' and one 'r' entry.`
      );
    }
  });

  // 3) Assign angles to these special pairs
  Object.entries(specialNodeIndices).forEach(([key, indices]) => {
    if (!indices || indices.length < 2) {
      return; // skip if not exactly two
    }
    const [[l1, u1], [l2, u2]] = indices; // e.g. [ [pathIndex, nodeIndex], [pathIndex, nodeIndex] ]
    const size1 = trs[l1].length;
    const size2 = trs[l2].length;
    if (u1 - 1 < 0 || u2 - 1 < 0) {
      return; // skip if referencing out of range
    }

    // Calculate angle from the sum of translations in adjacent steps
    const sumY1 = trs[l1][u1 - 1][1] + trs[l2][u2 - 1][1];
    const sumX1 = trs[l1][u1 - 1][0] + trs[l2][u2 - 1][0];
    const angle1 = Math.atan2(sumY1, sumX1);

    // Use modulo to handle wrap-around
    const u1_mod = u1 % size1;
    const u2_mod = u2 % size2;

    const sumY2 = trs[l1][u1_mod][1] + trs[l2][u2_mod][1];
    const sumX2 = trs[l1][u1_mod][0] + trs[l2][u2_mod][0];
    const angle2 = Math.atan2(sumY2, sumX2);

    // average angle
    const averageAngle = (angle1 + angle2) / 2;

    // Assign angle to the relevant translation steps
    trs[l1][u1 - 1] = [
      trs[l1][u1 - 1][0],
      trs[l1][u1 - 1][1],
      trs[l1][u1 - 1][2],
      averageAngle,
    ];
    trs[l2][u2 - 1] = [
      trs[l2][u2 - 1][0],
      trs[l2][u2 - 1][1],
      trs[l2][u2 - 1][2],
      averageAngle,
    ];
  });

  return trs;
}

/**
 * generateYarns: produce the 3D (x,y,z) points for a single yarn, given:
 *  - `nodes`: object of node data
 *  - `trs`: the array of [dx, dy, twist, angle] for each step
 *  - `yarnParams`: an array [yarnIndex, startNode, startL, crossingVal, z0]
 *  - `rad`: radius for twisting arcs
 *
 * returns { points, translations }
 */
export function generateYarns(nodes, trs, yarnParams, rad) {
  // yarnParams e.g. [0, startN=0, startL=0, crossingVal=1, z0=0.5]
  const [, startN, startL, crossingVal, z0] = yarnParams;

  let ptX = nodes[String(startN)].x;
  let ptY = nodes[String(startN)].y;
  let ptZ = crossingVal * z0;
  let crossing = -crossingVal; // sign flipping logic

  const points = [];
  const translations = [];
  //points.push([ptX, ptY, ptZ]);

  for (let i = 0; i < trs.length; i++) {
    // Loop in path order, applying twist logic
    let currentL = startL + i;
    if (currentL > trs.length - 1) {
      currentL = currentL % trs.length;
    }

    const [dx, dy, twists, angle = 0] = trs[currentL];

    // If first node is special and has odd twist, invert crossing again
    if (
      i === 0 &&
      Math.abs(nodes[String(startN)].twist || 0) > 0 &&
      (nodes[String(startN)].twist || 0) % 2 !== 0
    ) {
      crossing = -crossing;
    }

    ptX += dx;
    ptY += dy;
    ptZ = crossing * z0;

    // If twist is non-zero
    if (twists !== 0) {
      const auxPts = twistPoints(twists, angle);
      auxPts.forEach((ap, k) => {
        const px = ptX + ap[0] * rad;
        const py = ptY + ap[1] * rad;
        // slightly vary z based on index k, if you want 3D effect
        const pz = crossing * z0 * Math.cos((k * Math.PI) / 2);
        points.push([px, py, pz]);
      });
    } else {
      // no twist
      points.push([ptX, ptY, ptZ]);
    }

    if (twists % 2 === 0) crossing = -crossing;
    translations.push([dx, dy]);
  }

  return { points, translations };
}

/**
 * extendPoints: duplicates the given `points` array 'n' times,
 * each time offset by (vx, vy).
 * e.g. if n=2, you get:
 *   1st batch: original points
 *   2nd batch: each [x+vx, y+vy, z]
 */
export function extendPoints(points, n, vx, vy) {
  const extended = [];
  // if n=2, you get 2 repeats total: one at (0,0), one at (vx, vy)
  for (let k = -1; k < n - 1; k++) {
    points.forEach(([x, y, z], i) => {
      // Skip the first point for repetitions after the first
      extended.push([x + k * vx, y + k * vy, z]);
    });
  }
  return extended;
}

// If you had a Generator React component or other code here, you can keep or remove it.
// The important part is: no references to "unit_repetition" remain.

/**
 * Smooths the yarn points using cubic spline interpolation.
 * @param {Array} points - Array of [x, y, z] points.
 * @param {number} arcLength - Minimum distance between points.
 * @param {number} smoothness - Smoothing factor.
 * @param {number} numPoints - Number of points after smoothing.
 * @returns {Array} - Smoothed array of [x, y, z] points.
 */
export function smoothYarnPoints(
  points,
  arcLength = 1.0,
  smoothness = 0.1,
  numPoints = 10000
) {
  if (points.length < 2) return points;

  // Extract x, y, z components
  const x = points.map((p) => p[0]);
  const y = points.map((p) => p[1]);
  const z = points.map((p) => p[2]);

  // Parameterize the points
  const t = [0];
  for (let i = 1; i < points.length; i++) {
    const dx = x[i] - x[i - 1];
    const dy = y[i] - y[i - 1];
    const dz = z[i] - z[i - 1];
    t.push(t[i - 1] + Math.sqrt(dx * dx + dy * dy + dz * dz));
  }

  // Normalize t
  const tMin = t[0];
  const tMax = t[t.length - 1];
  const tNormalized = t.map((ti) => (ti - tMin) / (tMax - tMin));

  // Create splines for each dimension
  const splineX = new CubicSpline(tNormalized, x);
  const splineY = new CubicSpline(tNormalized, y);
  const splineZ = new CubicSpline(tNormalized, z);

  // Generate smoothed points
  const smoothed = [];
  const step =
    (tNormalized[tNormalized.length - 1] - tNormalized[0]) / numPoints;
  let currentT = tNormalized[0];
  let lastPoint = null;

  for (let i = 0; i < numPoints; i++) {
    const xi = splineX.at(currentT);
    const yi = splineY.at(currentT);
    const zi = splineZ.at(currentT);
    const point = [xi, yi, zi];

    if (!lastPoint || distance(lastPoint, point) >= arcLength) {
      smoothed.push(point);
      lastPoint = point;
    }

    currentT += step;
    if (currentT > tNormalized[tNormalized.length - 1]) break;
  }

  return smoothed;
}

/**
 * Calculates Euclidean distance between two points.
 * @param {Array} p1 - [x, y, z]
 * @param {Array} p2 - [x, y, z]
 * @returns {number} - Distance
 */
function distance(p1, p2) {
  const dx = p2[0] - p1[0];
  const dy = p2[1] - p1[1];
  const dz = p2[2] - p1[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}
