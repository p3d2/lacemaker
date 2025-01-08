import React, { useRef } from "react";
import {
  AppBar,
  Toolbar,
  Typography,
  Button,
  IconButton,
  Tooltip,
} from "@mui/material";
import AddIcon from "@mui/icons-material/Add";
import UploadIcon from "@mui/icons-material/Upload";
import SaveIcon from "@mui/icons-material/Save";

const PatternAppBar = ({ onNewPattern, onLoadPattern, onSavePattern }) => {
  const fileInputRef = useRef(null);

  const handleLoadClick = () => {
    fileInputRef.current.click();
  };

  const handleFileSelected = (event) => {
    const file = event.target.files[0];
    if (file) {
      onLoadPattern(file);
    }
  };

  return (
    <AppBar position="static">
      <Toolbar>
        {/* Title */}
        <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
          Pattern Editor
        </Typography>

        {/* New Pattern Button */}
        <Tooltip title="New Pattern">
          <Button
            color="inherit"
            startIcon={<AddIcon />}
            onClick={onNewPattern}
          >
            New Pattern
          </Button>
        </Tooltip>

        {/* Load Pattern Button */}
        <Tooltip title="Load Pattern">
          <Button
            color="inherit"
            startIcon={<UploadIcon />}
            onClick={handleLoadClick}
          >
            Load Pattern
          </Button>
        </Tooltip>

        {/* Hidden File Input for Load */}
        <input
          type="file"
          ref={fileInputRef}
          style={{ display: "none" }}
          onChange={handleFileSelected}
        />

        {/* Save Pattern Button */}
        <Tooltip title="Save Pattern">
          <Button
            color="inherit"
            startIcon={<SaveIcon />}
            onClick={onSavePattern}
          >
            Save Pattern
          </Button>
        </Tooltip>
      </Toolbar>
    </AppBar>
  );
};

export default PatternAppBar;
