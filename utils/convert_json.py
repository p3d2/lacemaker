# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 16:11:51 2025

@author: silvap1
"""

import os
import json
import argparse

def convert_old_to_new(old_data):
    """
    Convert the old JSON format into the new JSON format.
    Returns a new dictionary in the new format.
    """

    # -------------------------------------------------
    # 1) Convert old "nodes" from something like:
    #    {
    #      "0": [x, y],
    #      "1": [x, y, twist],
    #      ...
    #    }
    #    to new style:
    #    {
    #      "0": {"x": ..., "y": ..., "twist": ...},
    #      ...
    #    }
    # -------------------------------------------------
    old_nodes = old_data["nodes"]  # dict: { "0": [x, y], "1": [x, y, twist], ... }
    new_nodes = {}
    for node_id_str, node_list in old_nodes.items():
        x = node_list[0]
        y = node_list[1]
        # If there's a third element, assume it's twist; else default to 0
        if len(node_list) >= 3:
            twist = node_list[2]
        else:
            twist = 0
        new_nodes[node_id_str] = {
            "x": x,
            "y": y,
            "twist": twist
        }

    # -------------------------------------------------
    # 2) Convert each path's "shifts" from a dict with keys like
    #    "[7, 8]": [dx, dy]
    #    into a list of objects:
    #    [
    #       {"from": 7, "to": 8, "dx": 0.0, "dy": 16.0},
    #       ...
    #    ]
    # -------------------------------------------------
    new_paths = []
    for path_obj in old_data["paths"]:
        old_shifts_dict = path_obj.get("shifts", {})
        shift_list = []
        for key_str, shift_val in old_shifts_dict.items():
            # key_str might look like "[7, 8]"
            clean_key = key_str.strip("[]")
            parts = [x.strip() for x in clean_key.split(",")]
            from_node = int(parts[0])
            to_node   = int(parts[1])
            dx, dy = shift_val
            shift_list.append({
                "from": from_node,
                "to":   to_node,
                "dx":   dx,
                "dy":   dy
            })

        # Keep "path" the same, but replace "shifts" with the new list form
        new_paths.append({
            "path": path_obj["path"],
            "shifts": shift_list
        })

    # -------------------------------------------------
    # 3) Convert old "unit_yarns" from arrays to dicts
    #    e.g. old style:
    #      "0": [path_id, starting_node, start_path_index, z_sign, z_height]
    #    into new style:
    #      "0": {
    #        "path_id": <int>,
    #        "starting_node": <int>,
    #        "start_path_index": <int>,
    #        "z_sign": <int>,
    #        "z_height": <float>,
    #        "name": <int or str>,
    #        "repetitions": [
    #           {"count": ..., "dx": ..., "dy": ...},
    #           {"count": ..., "dx": ..., "dy": ...}
    #        ]
    #      }
    # -------------------------------------------------
    old_unit_yarns = old_data.get("unit_yarns", {})
    old_unit_rep   = old_data.get("unit_repetion", {})

    new_unit_yarns = {}
    for yarn_id_str, yarn_array in old_unit_yarns.items():
        # Suppose your old yarn_array is exactly 5 elements:
        # [path_id, starting_node, start_path_index, z_sign, z_height]
        path_id          = yarn_array[0]
        starting_node    = yarn_array[1]
        start_path_index = yarn_array[2]
        z_sign           = yarn_array[3]
        z_height         = yarn_array[4]
        
        # Start building the new-style object
        yarn_info = {
            "path_id":          path_id,
            "starting_node":    starting_node,
            "start_path_index": start_path_index,
            "z_sign":           z_sign,
            "z_height":         z_height,
            "name": 0,                # We'll fill from old_unit_rep if possible
            "repetitions": []        # We'll fill from old_unit_rep if possible
        }

        # Check if we have an entry in old_unit_rep for this yarn
        if yarn_id_str in old_unit_rep:
            rep_list = old_unit_rep[yarn_id_str]
            # If we expect exactly 7 elements, for example:
            #   [count1, dx1, dy1, count2, dx2, dy2, name]
            if len(rep_list) == 7:
                yarn_info["repetitions"].append({
                    "count": rep_list[0],
                    "dx":    rep_list[1],
                    "dy":    rep_list[2]
                })
                yarn_info["repetitions"].append({
                    "count": rep_list[3],
                    "dx":    rep_list[4],
                    "dy":    rep_list[5]
                })
                yarn_info["name"] = rep_list[6]  # last entry is "name"
            # If the user had 1 repetition plus a name or any other variation,
            # adapt as needed:
            elif len(rep_list) == 4:
                # example: [count, dx, dy, name]
                yarn_info["repetitions"].append({
                    "count": rep_list[0],
                    "dx":    rep_list[1],
                    "dy":    rep_list[2]
                })
                yarn_info["name"] = rep_list[3]
            # etc. (You can add more conditions if your data varies.)

        new_unit_yarns[yarn_id_str] = yarn_info

    # -------------------------------------------------
    # 4) The rest can remain unchanged: roi_bounds, etc.
    # -------------------------------------------------
    new_data = {
        "nodes":       new_nodes,
        "paths":       new_paths,
        "unit_yarns":  new_unit_yarns,
        "roi_bounds":  old_data["roi_bounds"]
    }

    return new_data

def main():
    parser = argparse.ArgumentParser(
        description="Convert old-format JSON files to new-format JSON files."
    )
    parser.add_argument(
        "input_folder",
        type=str,
        help="Path to folder containing old-format JSON files."
    )
    parser.add_argument(
        "--output_folder",
        type=str,
        default="converted_json",
        help="Folder to save new-format JSON files."
    )
    args = parser.parse_args()

    input_folder = args.input_folder
    output_folder = args.output_folder
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all .json files in input_folder
    for filename in os.listdir(input_folder):
        if filename.lower().endswith(".json"):
            old_json_path = os.path.join(input_folder, filename)

            # Load the old JSON
            with open(old_json_path, 'r') as f:
                old_data = json.load(f)

            # Convert to the new format
            new_data = convert_old_to_new(old_data)

            # Build the new filename
            new_filename = os.path.splitext(filename)[0] + "_new.json"
            new_json_path = os.path.join(output_folder, new_filename)

            # Save the new JSON
            with open(new_json_path, 'w') as f:
                json.dump(new_data, f, indent=2)

            print(f"Converted '{filename}' -> '{new_filename}'")

if __name__ == "__main__":
    main()
