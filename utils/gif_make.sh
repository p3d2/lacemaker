#!/bin/bash

# Variables
input_folder="$WRKDIR/lacemaker/temp"
output_file="$WRKDIR/lacemaker/temp/out.gif"
temp_folder="$input_folder/temp"

# Create temporary folder
mkdir -p "$temp_folder"

# Create an array to hold file names
declare -a files

for file in "$input_folder"/*.png; do

  # Step 1: Make non-transparent pixels black
  convert "$file" -colorspace sRGB -fill black +opaque none "$temp_folder/$(basename "$file" .png)_black.png"

  # Step 2: Expand black pixels by 2px with -draw
  convert "$temp_folder/$(basename "$file" .png)_black.png" -morphology Dilate Disk:10 "$temp_folder/$(basename "$file" .png)_expanded.png"

  # Step 3: Composite expanded image with original image
  convert "$temp_folder/$(basename "$file" .png)_expanded.png" -colorspace sRGB "$file" -compose Over -composite "$temp_folder/$(basename "$file" .png)_final.png"

  # Add the file name to array
  files+=("$temp_folder/$(basename "$file" .png)_final.png")
  
  # Remove intermediate files
  rm "$temp_folder/$(basename "$file" .png)_black.png"
  rm "$temp_folder/$(basename "$file" .png)_expanded.png"
done

# Convert images to GIF
convert -dispose Background -delay 10 -loop 0 "${files[@]}" -delay 100 "${files[-1]}" "$output_file"

# Clean up temporary folder
rm -rf "$temp_folder"