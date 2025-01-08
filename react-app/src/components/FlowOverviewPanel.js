import React, { useState, useEffect, useCallback } from "react";
import { Box, Button, Typography } from "@mui/material";
import Papa from "papaparse";
// Import DataGrid and the grouping types
import {
  DataGrid,
  useGridApiRef,
  GridToolbarContainer,
} from "@mui/x-data-grid";

// Example: All columns from your CSV
//   pattern_id, status, dist_particles, units, mass, threshold, ks, kb,
//   bond0_1, bond0_2, bond0_3, bond0_4, bond0_5, bond0_6, bond0_7, bond0_8,
//   visc, tstep, thermodump, nframes, xmin, xmax, ymin, ymax, etc.
const FlowOverviewPanel = ({ logs = "", setLogs = () => {} }) => {
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [selectedRows, setSelectedRows] = useState([]);

  const CustomToolbar = ({ selectedRows, handleAddRow, handleDeleteRows }) => {
    return (
      <GridToolbarContainer>
        <Button variant="contained" color="primary" onClick={handleAddRow}>
          Add Row
        </Button>
        <Button
          variant="contained"
          color="error"
          onClick={handleDeleteRows}
          disabled={selectedRows.length === 0}
        >
          Delete Selected
        </Button>
        <Button variant="contained" onClick={handleSaveCSV}>
          Save CSV
        </Button>
      </GridToolbarContainer>
    );
  };

  /**
   * Handle adding a new row to the Data Grid.
   * Appends a new job with default values at the end of the table.
   */
  const handleAddRow = useCallback(() => {
    // Define default values for a new job
    const newJob = {
      pattern_id: "ChangeMe",
      status: false, // Default status
      dist_particles: 0.25, // Default or empty value
      units: 1.0, // Default or empty value
      mass: 1.0, // Default or empty value
      threshold: 0.0, // Default or empty value
      ks: 30.0, // Default or empty value
      kb: 0, // Default or empty value
      bond0_1: 0,
      bond0_2: 0,
      bond0_3: 0,
      bond0_4: 0,
      bond0_5: 0,
      bond0_6: 0,
      bond0_7: 0,
      bond0_8: 0,
      visc: 0,
      tstep: 0,
      thermodump: 0,
      nframes: 0,
      xmin: 0,
      xmax: 0,
      ymin: 0,
      ymax: 0,
      visc2: 0,
      tstep2: 0,
      thermodump2: 0,
      nframes2: 0,
      exp: "",
      contract: 0,
      "Obs.": "",
    };
    setJobs((prevJobs) => [...prevJobs, newJob]);

    // Optionally, log the addition
    setLogs(
      (prevLogs) =>
        prevLogs + `Added new job with Pattern ID ${newJob.pattern_id}\n`
    );
  }, [setJobs, setLogs]);

  const handleDeleteRows = useCallback(() => {
    if (selectedRows.length === 0) return;

    // Optional: Add confirmation dialog
    const confirmDelete = window.confirm(
      `Are you sure you want to delete ${selectedRows.length} selected row(s)?`
    );
    if (!confirmDelete) return;

    setJobs((prevJobs) =>
      prevJobs.filter((job) => !selectedRows.includes(String(job.pattern_id)))
    );

    // Clear the selection after deletion
    setSelectedRows([]);

    // Optionally, log the deletion
    setLogs(
      (prevLogs) =>
        prevLogs +
        `Deleted ${
          selectedRows.length
        } row(s) at ${new Date().toLocaleString()}\n`
    );
  }, [selectedRows, setJobs, setSelectedRows, setLogs]);

  /**
   * Load the CSV from /data/jobs.csv once on mount.
   */
  useEffect(() => {
    Papa.parse("/data/jobs.csv", {
      download: true,
      header: true, // Parse first row as column headers
      dynamicTyping: true, // Convert numeric fields to numbers
      complete: (results) => {
        //console.log("CSV Parsing Complete:", results);
        // Filter out any empty rows (optional)
        const parsedJobs = results.data.filter(
          (job) => Object.keys(job).length > 1
        );
        setJobs(parsedJobs);
        setLoading(false);
      },
      error: (err) => {
        console.error("Error loading jobs.csv:", err);
        setError(err);
        setLoading(false);
      },
    });
  }, []);

  /**
   * 2) Define columns for the Data Grid.
   *    Mark 'editable: true' for fields the user can change.
   *    For your boolean status, we can use 'type: "boolean"'.
   */
  const columns = [
    {
      field: "pattern_id",
      headerName: "Pattern ID",
      width: 100,
      editable: true,
      disableColumnMenu: true,
    },
    {
      field: "status",
      headerName: "Status",
      type: "boolean", // Renders checkbox for boolean
      width: 75,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "dist_particles",
      headerName: "d_Particles",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "units",
      headerName: "units",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "mass",
      headerName: "mass",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "threshold",
      headerName: "thresh",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "ks",
      headerName: "ks",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "kb",
      headerName: "kb",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    // Bond columns
    {
      field: "bond0_1",
      headerName: "b_1",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_2",
      headerName: "b_2",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_3",
      headerName: "b_3",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_4",
      headerName: "b_4",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_5",
      headerName: "b_5",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_6",
      headerName: "b_6",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_7",
      headerName: "b_7",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "bond0_8",
      headerName: "b_8",
      type: "number",
      width: 50,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "visc",
      headerName: "Visc",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "tstep",
      headerName: "Tstep",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "thermodump",
      headerName: "Dump",
      type: "number",
      width: 80,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "nframes",
      headerName: "NFrames",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "xmin",
      headerName: "xmin",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "xmax",
      headerName: "xmax",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "ymin",
      headerName: "ymin",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "ymax",
      headerName: "ymax",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "visc2",
      headerName: "Visc2",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "tstep2",
      headerName: "Tstep2",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "thermodump2",
      headerName: "Dump2",
      type: "number",
      width: 80,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "nframes2",
      headerName: "NFrames2",
      type: "number",
      width: 60,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "exp",
      headerName: "Exp",
      width: 150,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "contract",
      headerName: "Contract",
      type: "number",
      width: 100,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
    {
      field: "Obs.",
      headerName: "Obs.",
      width: 200,
      editable: true,
      disableColumnMenu: true,
      sortable: false,
    },
  ];

  /**
   * Column Grouping Model
   * Groups related columns under a common header for better organization.
   * Adjust group names and children as per your CSV structure.
   */
  const columnGroupingModel = [
    {
      groupId: "basic",
      headerName: "Basic Info",
      children: [
        { field: "pattern_id" },
        { field: "status" },
        { field: "dist_particles" },
        { field: "units" },
        { field: "mass" },
        { field: "threshold" },
        { field: "ks" },
        { field: "kb" },
      ],
    },
    {
      groupId: "bonds",
      headerName: "Bonds",
      children: Array.from({ length: 8 }, (_, i) => ({
        field: `bond0_${i + 1}`,
      })),
    },
    {
      groupId: "sim0",
      headerName: "Simulation 0",
      children: [
        { field: "visc" },
        { field: "tstep" },
        { field: "thermodump" },
        { field: "nframes" },
        { field: "xmin" },
        { field: "xmax" },
        { field: "ymin" },
        { field: "ymax" },
      ],
    },
    {
      groupId: "sim1",
      headerName: "Simulation 1",
      children: [
        { field: "visc2" },
        { field: "tstep2" },
        { field: "thermodump2" },
        { field: "nframes2" },
        { field: "exp" },
        { field: "contract" },
      ],
    },
    {
      groupId: "obs",
      headerName: "Observations",
      children: [{ field: "Obs." }],
    },
  ];

  /**
   * Handle cell edits to update the 'jobs' state accordingly.
   */
  const handleCellEditCommit = useCallback(
    (params) => {
      const { id, field, value } = params;
      setJobs((prevJobs) =>
        prevJobs.map((row) => {
          if (String(row.pattern_id) !== String(id)) return row;
          // Handle boolean conversion for 'status' field
          if (field === "status") {
            return { ...row, [field]: Boolean(value) };
          }
          return { ...row, [field]: value };
        })
      );
    },
    [setJobs]
  );

  /**
   * Handle saving the edited jobs back to a CSV file.
   */
  const handleSaveCSV = useCallback(async () => {
    // Convert the jobs array back to CSV format
    const csv = Papa.unparse(jobs);

    try {
      const response = await fetch("/save-jobs", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ csv }),
      });

      if (!response.ok) {
        throw new Error("Failed to save the file.");
      }

      // Optionally, log the save action
      setLogs(
        (prevLogs) =>
          prevLogs +
          `File saved successfully at ${new Date().toLocaleString()}\n`
      );
    } catch (error) {
      console.error("Error saving file:", error);
      setLogs(
        (prevLogs) =>
          prevLogs +
          `Error saving file at ${new Date().toLocaleString()}: ${
            error.message
          }\n`
      );
    }
  }, [jobs, setLogs]);

  /**
   * Define the unique row identifier.
   */
  const getRowId = (row) => String(row.pattern_id);

  /**
   * Render loading or error states.
   */
  if (loading) {
    return (
      <Box sx={{ p: 2 }}>
        <Typography variant="h6">Loading jobs...</Typography>
      </Box>
    );
  }

  if (error) {
    return (
      <Box sx={{ p: 2 }}>
        <Typography variant="h6" color="error">
          Error loading jobs.csv: {error.message}
        </Typography>
      </Box>
    );
  }

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        height: "100%",
      }}
    >
      {/* Header */}
      <Box sx={{ mb: 1 }}>
        <Typography variant="h5">Jobs Editor</Typography>
        <Typography variant="body2">
          Edit job parameters below. After making changes, click &quot;Save
          CSV&quot; to export the updated file. Commit the changes to your
          repository and run your HPC jobs as usual.
        </Typography>
      </Box>

      {/* Data Grid */}
      <Box
        sx={{
          flex: 1,
          mb: 2,
          width: "100%",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <DataGrid
          rows={jobs}
          columns={columns}
          getRowId={getRowId}
          onCellEditCommit={handleCellEditCommit}
          disableSelectionOnClick={false} // Allows row selection by clicking anywhere on the row
          rowSelectionModel={selectedRows}
          onRowSelectionModelChange={(newSelection) => {
            setSelectedRows(newSelection);
          }}
          slots={{
            toolbar: CustomToolbar,
          }}
          slotProps={{
            toolbar: {
              selectedRows, // Pass only relevant props
              handleAddRow,
              handleDeleteRows,
            },
          }}
          pageSizeOptions={[10, 25, 100]}
          sx={{
            "& .MuiDataGrid-columnHeaders": {
              backgroundColor: "#f0f0f0",
            },
            // Enable horizontal scrolling
            "& .MuiDataGrid-root": {
              overflowX: "auto",
            },
          }}
        />
      </Box>
    </Box>
  );
};

export default FlowOverviewPanel;
