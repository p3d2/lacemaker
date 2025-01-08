import React, { useState } from "react";
import {
  Box,
  Typography,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Button,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import NodesEditor from "./NodesEditor";
import PathsEditor from "./PathsEditor";
import UnitYarnsEditor from "./UnitYarnsEditor";
import ROIBoundsEditor from "./ROIBoundsEditor";

const ConfigurationPanel = ({ data, setData }) => {
  const [expanded, setExpanded] = useState({
    nodes: false,
    paths: false,
    unitYarns: false,
    roiBounds: false,
  });

  // Toggle all menus
  const toggleAll = (expand) => {
    setExpanded({
      nodes: expand,
      paths: expand,
      unitYarns: expand,
      roiBounds: expand,
    });
  };

  return (
    <Box
      sx={{
        borderRight: "1px solid #ccc",
        overflow: "auto",
        p: 2,
        background: "#fff",
      }}
    >
      <Box
        display="flex"
        alignItems="center"
        justifyContent="space-between"
        mb={1}
      >
        <Typography variant="h5" gutterBottom>
          Configuration
        </Typography>
        <Button
          variant="outlined"
          size="small"
          onClick={() => toggleAll(!expanded.nodes)}
        >
          {expanded.nodes ? "Collapse All" : "Expand All"}
        </Button>
      </Box>

      {/* Nodes Accordion */}
      <Accordion
        expanded={expanded.nodes}
        onChange={() => setExpanded({ ...expanded, nodes: !expanded.nodes })}
      >
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Nodes</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <NodesEditor
            nodes={data.nodes}
            setNodes={(nodes) => setData({ ...data, nodes })}
          />
        </AccordionDetails>
      </Accordion>

      {/* Paths Accordion */}
      <Accordion
        expanded={expanded.paths}
        onChange={() => setExpanded({ ...expanded, paths: !expanded.paths })}
      >
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Paths</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <PathsEditor
            paths={data.paths}
            setPaths={(paths) => setData({ ...data, paths })}
          />
        </AccordionDetails>
      </Accordion>

      {/* Unit Yarns Accordion */}
      <Accordion
        expanded={expanded.unitYarns}
        onChange={() =>
          setExpanded({ ...expanded, unitYarns: !expanded.unitYarns })
        }
      >
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Unit Yarns</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <UnitYarnsEditor
            unitYarns={data.unit_yarns}
            setUnitYarns={(unit_yarns) => setData({ ...data, unit_yarns })}
          />
        </AccordionDetails>
      </Accordion>

      {/* ROI Bounds Accordion */}
      <Accordion
        expanded={expanded.roiBounds}
        onChange={() =>
          setExpanded({ ...expanded, roiBounds: !expanded.roiBounds })
        }
      >
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>ROI Bounds</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <ROIBoundsEditor
            roiBounds={data.roi_bounds}
            setROIBounds={(roiBounds) =>
              setData({ ...data, roi_bounds: roiBounds })
            }
          />
        </AccordionDetails>
      </Accordion>
    </Box>
  );
};

export default ConfigurationPanel;
