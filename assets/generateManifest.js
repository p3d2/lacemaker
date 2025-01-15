// assets/generateManifest.js
const fs = require("fs");
const path = require("path");

const folderPath = path.join(__dirname, "..", "input", "json_patterns");
const manifest = [];

// Adjust these values as needed.
const username = "p3d2";
const repository = "lacemaker";
const branch = "main"; // or your default branch

fs.readdirSync(folderPath).forEach((file) => {
  if (file.endsWith(".json") && file !== "index.json") {
    manifest.push({
      name: path.basename(file, ".json"),
      url: `https://raw.githubusercontent.com/${username}/${repository}/${branch}/input/json_patterns/${file}`
    });
  }
});

fs.writeFileSync(
  path.join(folderPath, "index.json"),
  JSON.stringify(manifest, null, 2)
);

console.log("Manifest updated!");
