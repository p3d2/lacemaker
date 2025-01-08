import React, { useState, useEffect } from "react";
import {
  Box,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Typography,
  IconButton,
  Button,
  TextField,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import DeleteIcon from "@mui/icons-material/DeleteOutline";
import AddCircleOutlineIcon from "@mui/icons-material/AddCircleOutline";

const PathsEditor = ({ paths, setPaths }) => {
  const [expanded, setExpanded] = useState(false);

  // State to hold raw input strings for each path
  const [tempPaths, setTempPaths] = useState({});

  // Initialize tempPaths based on the initial paths prop
  useEffect(() => {
    const initialTempPaths = {};
    paths.forEach((path, index) => {
      // Join the path array into a comma-separated string
      initialTempPaths[index] = path.path.join(", ");
    });
    setTempPaths(initialTempPaths);
  }, [paths]);

  // Handle accordion expansion
  const handleAccordionChange = (panel) => (event, isExpanded) => {
    setExpanded(isExpanded ? panel : false);
  };

  // Handle adding a new path
  const handleAddPath = () => {
    const newPath = { path: [], shifts: [] };
    setPaths([...paths, newPath]);
    // Initialize the corresponding tempPath as an empty string
    setTempPaths((prev) => ({ ...prev, [paths.length]: "" }));
  };

  // Handle deleting a path
  const handleDeletePath = (index) => {
    const updatedPaths = paths.filter((_, i) => i !== index);
    setPaths(updatedPaths);

    // Remove the corresponding tempPath
    setTempPaths((prev) => {
      const updatedTemp = { ...prev };
      delete updatedTemp[index];
      // Reindex the remaining tempPaths
      const newTemp = {};
      Object.keys(updatedTemp).forEach((key) => {
        const newKey = key > index ? key - 1 : key;
        newTemp[newKey] = updatedTemp[key];
      });
      return newTemp;
    });
  };

  // Handle changing the nodes in a path (raw input)
  const handlePathChange = (index, newValue) => {
    setTempPaths((prev) => ({ ...prev, [index]: newValue }));
  };

  // Handle blur event to parse and update the path
  const handlePathBlur = (index, value) => {
    const input = typeof value === "string" ? value : "";
    if (input.trim() === "") {
      // If input is empty, set path to an empty array
      updatePath(index, []);
      return;
    }

    // Split the input by commas and trim each item
    const items = input.split(",").map((item) => item.trim());

    // Process each item
    const cleaned = items.map((trimmed) => {
      if (trimmed.endsWith("r") || trimmed.endsWith("l")) {
        return trimmed; // Keep as string (e.g., "1r")
      } else {
        const parsed = parseInt(trimmed, 10);
        return isNaN(parsed) ? trimmed : parsed; // Convert to number or keep as string if NaN
      }
    });

    // Update the path in the paths state
    updatePath(index, cleaned);
  };

  // Helper function to update a specific path
  const updatePath = (index, newPathArray) => {
    const updatedPaths = [...paths];
    updatedPaths[index] = { ...updatedPaths[index], path: newPathArray };
    setPaths(updatedPaths);
  };

  // Handle adding a shift
  const handleAddShift = (pathIndex) => {
    const newShift = { from: 0, to: 0, dx: 0, dy: 0 };
    const updatedPaths = [...paths];
    updatedPaths[pathIndex].shifts.push(newShift);
    setPaths(updatedPaths);
  };

  // Handle changing a shift field
  const handleShiftChange = (pathIndex, shiftIndex, field, value) => {
    const updatedPaths = [...paths];
    updatedPaths[pathIndex].shifts[shiftIndex][field] = value;
    setPaths(updatedPaths);
  };

  const handleShiftBlur = (pathIndex, shiftIndex, field, value) => {
    const updatedPaths = [...paths];
    let updatedValue = value;

    if (field === "from" || field === "to") {
      updatedValue = parseInt(value, 10);
      updatedPaths[pathIndex].shifts[shiftIndex][field] = isNaN(updatedValue)
        ? 0
        : updatedValue;
    } else if (field === "dx" || field === "dy") {
      updatedValue = parseFloat(value.replace(",", "."));
      updatedPaths[pathIndex].shifts[shiftIndex][field] = isNaN(updatedValue)
        ? 0
        : updatedValue;
    }
    setPaths(updatedPaths);
  };

  // Handle deleting a shift
  const handleDeleteShift = (pathIndex, shiftIndex) => {
    const updatedPaths = [...paths];
    updatedPaths[pathIndex].shifts.splice(shiftIndex, 1);
    setPaths(updatedPaths);
  };

  // Handle saving all paths and shifts
  const handleSave = () => {
    const updatedPaths = paths.map((p) => {
      const newShifts = p.shifts.map((s) => ({
        ...s,
        dx: parseFloat(String(s.dx)) || 0,
        dy: parseFloat(String(s.dy)) || 0,
        from: parseInt(String(s.from), 10) || 0,
        to: parseInt(String(s.to), 10) || 0,
      }));
      return { ...p, shifts: newShifts };
    });
    setPaths(updatedPaths);
    // then serialize to JSON, etc.
  };

  return (
    <Box>
      {/* Paths Accordions */}
      {paths.map((path, pathIndex) => (
        <Accordion
          key={pathIndex}
          expanded={expanded === `panel${pathIndex}`}
          onChange={handleAccordionChange(`panel${pathIndex}`)}
        >
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography>Path {pathIndex}</Typography>
            <IconButton
              color="error"
              onClick={(e) => {
                e.stopPropagation(); // Prevent accordion toggle
                handleDeletePath(pathIndex);
              }}
              sx={{ marginLeft: "auto" }}
            >
              <DeleteIcon />
            </IconButton>
          </AccordionSummary>

          <AccordionDetails>
            {/* Path Nodes */}
            <Box mb={2}>
              <TextField
                label="Nodes (comma-separated)"
                value={tempPaths[pathIndex] ?? path.path.join(", ")}
                onChange={(e) => handlePathChange(pathIndex, e.target.value)}
                onBlur={(e) => handlePathBlur(pathIndex, e.target.value)}
                fullWidth
                multiline
                rows={4}
                variant="outlined"
              />
            </Box>
            <TableContainer component={Paper}>
              <Table>
                <TableHead>
                  <TableRow>
                    <TableCell sx={{ padding: "4px 4px" }}>From</TableCell>
                    <TableCell sx={{ padding: "4px 4px" }}>To</TableCell>
                    <TableCell sx={{ padding: "4px 4px" }}>dx</TableCell>
                    <TableCell sx={{ padding: "4px 4px" }}>dy</TableCell>
                    <TableCell
                      align="center"
                      sx={{ padding: "8px 12px" }}
                    ></TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {path.shifts.map((shift, shiftIndex) => (
                    <TableRow key={shiftIndex}>
                      <TableCell
                        sx={{
                          width: "40px",
                          padding: "2px 2px",
                          whiteSpace: "nowrap",
                        }}
                      >
                        <TextField
                          value={shift["from"]}
                          onChange={(e) =>
                            handleShiftChange(
                              pathIndex,
                              shiftIndex,
                              "from",
                              e.target.value
                            )
                          }
                          onBlur={(e) =>
                            handleShiftBlur(
                              pathIndex,
                              shiftIndex,
                              "from",
                              e.target.value
                            )
                          }
                          variant="outlined"
                          size="small"
                          fullWidth
                        />
                      </TableCell>
                      <TableCell
                        sx={{
                          width: "40px",
                          padding: "2px 2px",
                          whiteSpace: "nowrap",
                        }}
                      >
                        <TextField
                          value={shift["to"]}
                          onChange={(e) =>
                            handleShiftChange(
                              pathIndex,
                              shiftIndex,
                              "to",
                              e.target.value
                            )
                          }
                          onBlur={(e) =>
                            handleShiftBlur(
                              pathIndex,
                              shiftIndex,
                              "to",
                              e.target.value
                            )
                          }
                          variant="outlined"
                          size="small"
                          fullWidth
                        />
                      </TableCell>
                      <TableCell
                        sx={{
                          width: "40px",
                          padding: "2px 2px",
                          whiteSpace: "nowrap",
                        }}
                      >
                        <TextField
                          value={shift["dx"]}
                          onChange={(e) =>
                            handleShiftChange(
                              pathIndex,
                              shiftIndex,
                              "dx",
                              e.target.value
                            )
                          }
                          onBlur={(e) =>
                            handleShiftBlur(
                              pathIndex,
                              shiftIndex,
                              "dx",
                              e.target.value
                            )
                          }
                          variant="outlined"
                          size="small"
                          fullWidth
                        />
                      </TableCell>
                      <TableCell
                        sx={{
                          width: "40px",
                          padding: "2px 2px",
                          whiteSpace: "nowrap",
                        }}
                      >
                        <TextField
                          value={shift["dy"]}
                          onChange={(e) =>
                            handleShiftChange(
                              pathIndex,
                              shiftIndex,
                              "dy",
                              e.target.value
                            )
                          }
                          onBlur={(e) =>
                            handleShiftBlur(
                              pathIndex,
                              shiftIndex,
                              "dy",
                              e.target.value
                            )
                          }
                          variant="outlined"
                          size="small"
                          fullWidth
                        />
                      </TableCell>
                      <TableCell
                        align="center"
                        sx={{ width: "20px", padding: "4px 4px" }}
                      >
                        <IconButton
                          color="error"
                          onClick={(e) => {
                            e.stopPropagation(); // Prevent row selection
                            handleDeleteShift(pathIndex, shiftIndex);
                          }}
                        >
                          <DeleteIcon />
                        </IconButton>
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>

            {/* Add Shift Button */}
            <Box
              display="flex"
              justifyContent="flex-end"
              sx={{ padding: "8px 0px" }}
            >
              <Button
                variant="outlined"
                color="primary"
                startIcon={<AddCircleOutlineIcon />}
                onClick={() => handleAddShift(pathIndex)}
              >
                Add Shift
              </Button>
            </Box>
          </AccordionDetails>
        </Accordion>
      ))}
      <Box display="flex" justifyContent="flex-end" sx={{ padding: "8px 0px" }}>
        <Button
          variant="outlined"
          color="primary"
          startIcon={<AddCircleOutlineIcon />}
          onClick={() => handleAddPath()}
        >
          Add Path
        </Button>
      </Box>
    </Box>
  );
};

export default PathsEditor;
