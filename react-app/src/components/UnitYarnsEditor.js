// UnitYarnsEditor.js
import React, { useState } from "react";
import {
  Box,
  Typography,
  TextField,
  IconButton,
  Button,
  Accordion,
  AccordionSummary,
  AccordionDetails,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import DeleteIcon from "@mui/icons-material/DeleteOutline";
import AddCircleOutlineIcon from "@mui/icons-material/AddCircleOutline";

const UnitYarnsEditor = ({ unitYarns, setUnitYarns }) => {
  // Track which accordion is expanded: null or the Yarn ID string
  const [expanded, setExpanded] = useState(null);

  // Define which fields are numeric for consistent data types
  const numericFields = [
    "path_id",
    "starting_node",
    "start_path_index",
    "z_sign",
    "z_height",
    "name",
    "count",
    "dx",
    "dy",
  ];

  const isNumericField = (field) => numericFields.includes(field.split(".")[0]);

  // Handle field changes directly in unitYarns
  const handleChange = (yarnId, fieldPath, value) => {
    // Create a shallow copy of unitYarns
    const updatedYarns = { ...unitYarns };

    // Access the specific yarn to update
    const yarn = { ...updatedYarns[yarnId] };

    // Split the fieldPath to handle nested fields (e.g., "repetitions.0.count")
    const fields = fieldPath.split(".");

    if (fields.length === 1) {
      // Top-level field
      yarn[fields[0]] = value;
    } else if (fields.length === 3 && fields[0] === "repetitions") {
      // Nested field within repetitions array
      const repIndex = parseInt(fields[1], 10);
      const repField = fields[2];

      // Ensure the repetitions array exists
      if (!Array.isArray(yarn.repetitions)) {
        yarn.repetitions = [];
      }

      // Ensure the repetition at repIndex exists
      while (yarn.repetitions.length <= repIndex) {
        yarn.repetitions.push({ count: 2, dx: 0, dy: 0 });
      }

      // Update the specific field within the repetition
      yarn.repetitions[repIndex] = {
        ...yarn.repetitions[repIndex],
        [repField]: value,
      };
    } else {
      console.warn(`Unhandled fieldPath: ${fieldPath}`);
      return;
    }

    // Update the yarn in the updatedYarns object
    updatedYarns[yarnId] = yarn;

    // Update the state
    setUnitYarns(updatedYarns);
  };

  // Handle field blur to parse and clean values
  const handleBlur = (yarnId, fieldPath, rawValue) => {
    const cleaned = rawValue.replace(",", ".").replace(/[^\d.\-+]/g, "");
    const parsed = parseFloat(cleaned) || 0;
    handleChange(yarnId, fieldPath, parsed);
  };

  // Delete a yarn by ID
  const handleDeleteYarn = (yarnId) => {
    // Create a shallow copy of unitYarns
    const updatedYarns = { ...unitYarns };
    delete updatedYarns[yarnId];
    setUnitYarns(updatedYarns);

    // If the deleted yarn was expanded, collapse the accordion
    if (expanded === yarnId) {
      setExpanded(null);
    }
  };

  // Add a new yarn (exactly as specified)
  const handleAddYarn = () => {
    const existingIds = Object.keys(unitYarns).map(Number);
    const newIdNum = existingIds.length ? Math.max(...existingIds) + 1 : 0;
    const newId = String(newIdNum); // Ensure newId is a string

    setUnitYarns({
      ...unitYarns,
      [newId]: {
        path_id: 0,
        starting_node: 0,
        start_path_index: 0,
        z_sign: 1,
        z_height: 0.5,
        name: parseInt(newId, 10) + 1,
        repetitions: [
          { count: 4, dx: 0, dy: 0 },
          { count: 4, dx: 0, dy: 0 },
        ],
      },
    });
  };

  // Convert unitYarns object to array for iteration
  const yarnEntries = Object.entries(unitYarns);

  return (
    <Box>
      {yarnEntries.map(([yarnId, data]) => {
        const {
          path_id,
          starting_node,
          start_path_index,
          z_sign,
          z_height,
          name,
          repetitions,
        } = data;

        const isExpanded = expanded === yarnId;

        return (
          <Accordion
            key={yarnId}
            expanded={isExpanded}
            onChange={(e, newState) => setExpanded(newState ? yarnId : null)}
            sx={{ mb: 2 }}
          >
            {/* Header: Yarn #X + Delete */}
            <AccordionSummary expandIcon={<ExpandMoreIcon />}>
              <Box
                display="flex"
                width="100%"
                justifyContent="space-between"
                alignItems="center"
              >
                <Typography variant="subtitle1" sx={{ mr: 1 }}>
                  Yarn {yarnId}
                </Typography>
                <IconButton
                  color="error"
                  onClick={(e) => {
                    e.stopPropagation(); // Prevent accordion toggle
                    handleDeleteYarn(yarnId);
                  }}
                >
                  <DeleteIcon />
                </IconButton>
              </Box>
            </AccordionSummary>

            <AccordionDetails>
              {/* 5 Top-Level Fields */}
              <Box mb={2}>
                <Typography variant="body2">Path ID</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={path_id}
                  onChange={(e) =>
                    handleChange(yarnId, "path_id", e.target.value)
                  }
                  onBlur={(e) => handleBlur(yarnId, "path_id", e.target.value)}
                />
              </Box>

              <Box mb={2}>
                <Typography variant="body2">Starting Node</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={starting_node}
                  onChange={(e) =>
                    handleChange(yarnId, "starting_node", e.target.value)
                  }
                  onBlur={(e) =>
                    handleBlur(yarnId, "starting_node", e.target.value)
                  }
                />
              </Box>

              <Box mb={2}>
                <Typography variant="body2">Start Path Index</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={start_path_index}
                  onChange={(e) =>
                    handleChange(yarnId, "start_path_index", e.target.value)
                  }
                  onBlur={(e) =>
                    handleBlur(yarnId, "start_path_index", e.target.value)
                  }
                />
              </Box>

              <Box mb={2}>
                <Typography variant="body2">Z Sign</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={z_sign}
                  onChange={(e) =>
                    handleChange(yarnId, "z_sign", e.target.value)
                  }
                  onBlur={(e) => handleBlur(yarnId, "z_sign", e.target.value)}
                />
              </Box>

              <Box mb={2}>
                <Typography variant="body2">Z Height</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={z_height}
                  onChange={(e) =>
                    handleChange(yarnId, "z_height", e.target.value)
                  }
                  onBlur={(e) => handleBlur(yarnId, "z_height", e.target.value)}
                />
              </Box>

              <Box mb={2}>
                <Typography variant="body2">Name</Typography>
                <TextField
                  variant="outlined"
                  size="small"
                  fullWidth
                  value={data.name ?? ""}
                  onChange={(e) => handleChange(yarnId, "name", e.target.value)}
                />
              </Box>

              {/* Repetition Rendering */}
              {repetitions.map((rep, index) => (
                <React.Fragment key={index}>
                  <Typography variant="subtitle2" mb={1}>
                    Repetition #{index + 1}
                  </Typography>
                  <Box display="flex" gap={2} mb={2}>
                    <Box flex={1}>
                      <TextField
                        label="Count"
                        variant="outlined"
                        size="small"
                        fullWidth
                        value={rep.count}
                        onChange={(e) =>
                          handleChange(
                            yarnId,
                            `repetitions.${index}.count`,
                            e.target.value
                          )
                        }
                        onBlur={(e) =>
                          handleBlur(
                            yarnId,
                            `repetitions.${index}.count`,
                            e.target.value
                          )
                        }
                      />
                    </Box>
                    <Box flex={1}>
                      <TextField
                        label="dx"
                        variant="outlined"
                        size="small"
                        fullWidth
                        value={rep.dx}
                        onChange={(e) =>
                          handleChange(
                            yarnId,
                            `repetitions.${index}.dx`,
                            e.target.value
                          )
                        }
                        onBlur={(e) =>
                          handleBlur(
                            yarnId,
                            `repetitions.${index}.dx`,
                            e.target.value
                          )
                        }
                      />
                    </Box>
                    <Box flex={1}>
                      <TextField
                        label="dy"
                        variant="outlined"
                        size="small"
                        fullWidth
                        value={rep.dy}
                        onChange={(e) =>
                          handleChange(
                            yarnId,
                            `repetitions.${index}.dy`,
                            e.target.value
                          )
                        }
                        onBlur={(e) =>
                          handleBlur(
                            yarnId,
                            `repetitions.${index}.dy`,
                            e.target.value
                          )
                        }
                      />
                    </Box>
                  </Box>
                </React.Fragment>
              ))}
            </AccordionDetails>
          </Accordion>
        );
      })}

      <Box display="flex" justifyContent="flex-end">
        <Button
          variant="outlined"
          color="primary"
          startIcon={<AddCircleOutlineIcon />}
          onClick={handleAddYarn}
        >
          Add Yarn
        </Button>
      </Box>
    </Box>
  );
};

export default UnitYarnsEditor;
