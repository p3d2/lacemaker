import React, { useState } from "react";
import {
  Box,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextField,
  IconButton,
  Button,
} from "@mui/material";
import DeleteIcon from "@mui/icons-material/DeleteOutline";
import AddCircleOutlineIcon from "@mui/icons-material/AddCircleOutline";

const NodesEditor = ({ nodes, setNodes }) => {
  const [height, setHeight] = useState(300); // Default fixed height

  // Convert nodes object to an array of rows
  const rows = Object.entries(nodes).map(([id, { x, y, twist }]) => ({
    id,
    x,
    y,
    twist,
  }));

  // Handle field changes
  const handleFieldChange = (id, field, value) => {
    setNodes({
      ...nodes,
      [id]: { ...nodes[id], [field]: value },
    });
  };

  const handleFieldBlur = (id, field, value) => {
    const normalizedValue = parseFloat(value.replace(",", "."));
    setNodes({
      ...nodes,
      [id]: {
        ...nodes[id],
        [field]: isNaN(normalizedValue) ? 0 : normalizedValue,
      },
    });
  };

  // Handle delete node
  const handleDelete = (id) => {
    const updatedNodes = { ...nodes };
    delete updatedNodes[id];
    setNodes(updatedNodes);
  };

  // Handle add new node
  const handleAddRow = () => {
    const newId = getNextId();
    setNodes({ ...nodes, [newId]: { x: 0, y: 0, twist: 0 } });
  };

  // Get next available ID
  const getNextId = () => {
    const ids = Object.keys(nodes).map(Number);
    return ids.length > 0 ? Math.max(...ids) + 1 : 0;
  };

  // Handle slider change
  const handleHeightChange = (event, newValue) => {
    setHeight(newValue);
  };

  return (
    <Box>
      <TableContainer
        component={Paper}
        sx={{ height: `${height}px`, overflowY: "auto", mb: 2 }}
      >
        <Table size="small" stickyHeader>
          <TableHead>
            <TableRow>
              <TableCell align="center" sx={{ padding: "4px 8px" }}>
                ID
              </TableCell>
              <TableCell align="center" sx={{ padding: "4px 8px" }}>
                x
              </TableCell>
              <TableCell align="center" sx={{ padding: "4px 8px" }}>
                y
              </TableCell>
              <TableCell align="center" sx={{ padding: "4px 8px" }}>
                Twist
              </TableCell>
              <TableCell align="center" sx={{ padding: "4px 0px" }}></TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {rows.map((row) => (
              <TableRow key={row.id}>
                <TableCell align="center" sx={{ padding: "2px 2px" }}>
                  {row.id}
                </TableCell>
                <TableCell sx={{ padding: "2px 2px" }}>
                  <TextField
                    value={row.x}
                    onChange={(e) =>
                      handleFieldChange(row.id, "x", e.target.value)
                    }
                    onBlur={(e) => handleFieldBlur(row.id, "x", e.target.value)}
                    variant="outlined"
                    size="small"
                    fullWidth
                  />
                </TableCell>
                <TableCell sx={{ padding: "2px 2px" }}>
                  <TextField
                    value={row.y}
                    onChange={(e) =>
                      handleFieldChange(row.id, "y", e.target.value)
                    }
                    onBlur={(e) => handleFieldBlur(row.id, "y", e.target.value)}
                    variant="outlined"
                    size="small"
                    fullWidth
                  />
                </TableCell>
                <TableCell sx={{ padding: "2px 2px" }}>
                  <TextField
                    value={row.twist}
                    onChange={(e) =>
                      handleFieldChange(row.id, "twist", e.target.value)
                    }
                    onBlur={(e) =>
                      handleFieldBlur(row.id, "twist", e.target.value)
                    }
                    variant="outlined"
                    size="small"
                    fullWidth
                  />
                </TableCell>
                <TableCell align="center" sx={{ padding: "2px 2px" }}>
                  <IconButton
                    color="error"
                    onClick={() => handleDelete(row.id)}
                  >
                    <DeleteIcon />
                  </IconButton>
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      <Box display="flex" justifyContent="flex-end">
        <Button
          variant="contained"
          color="primary"
          startIcon={<AddCircleOutlineIcon />}
          onClick={handleAddRow}
        >
          Add Node
        </Button>
      </Box>
    </Box>
  );
};

export default NodesEditor;
