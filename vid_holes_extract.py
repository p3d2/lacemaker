import cv2
import numpy as np
import json
import os
import argparse
import subprocess
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import glob

# Create the colormap
cmap = plt.cm.turbo

def process_video(video_path, output_json, vname, frame_skip=1):
    # Initialize video capture
    cap = cv2.VideoCapture(video_path)
    frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fps = cap.get(cv2.CAP_PROP_FPS)
    
    areas_data = []
    
    frame_idx = 0
    frame_c = 0
    frames = []
    while cap.isOpened():
        # Read frame
        ret, frame = cap.read()
        if not ret:
            break

         # Skip initial frames
        if frame_idx < skip_ini:
            frame_idx += 1
            continue

        # Skip every `frame_skip` frames after initial skip
        if frame_idx % frame_skip != 0:
            frame_idx += 1
            continue

        # Step 1: Crop the frame to remove borders
        px = 1
        if frame_c == 0:
            boundary = crop_frame(frame)
            [x, y, w, h] = cv2.boundingRect(boundary)
            px2mm = 50 / np.max((w, h))
            boundary = boundary - np.array([x, y]) + np.array([xalign, yalign])
            #norm_area = ((w + h)/(128*2)) ** 2 # convert vid area to sim area

        # Crop the frame to this rectangle
        cropped_frame = frame[y-px:y+h+px, x-px:x+w+px]

        # Step 2: Convert to grayscale using the value channel
        gray_frame = convert_to_grayscale(cropped_frame)
        
        # Step 3: Apply thresholding
        binary_frame = apply_threshold(gray_frame)

        # Step 4: Perform erosion and dilation
        processed_frame = remove_noise(binary_frame)
        
        # Step 5: Extract contours (holes) and compute areas
        if frame_c == 0: mask = bound_find(processed_frame, boundary)
        (areas, frame_n) = extract_areas(cropped_frame, processed_frame, mask, px2mm, frame_c)
        frames.append(frame_n)

        # Save areas for this frame
        frame_data = {
            'time': frame_c * frame_skip / fps,
            'areas': areas
        }
        areas_data.append(frame_data)
        frame_idx += 1
        frame_c += 1

    # Generate GIF from frames
    gif_filename = os.path.join(subfolder_path, 'plots', vname + '_area_dist.gif')
    with imageio.get_writer(gif_filename, mode='I', fps=10, loop=0) as writer:
        for k in range(len(frames)):
            image = imageio.imread(frames[k])
            writer.append_data(image)

    num_images = 5
    step = (len(frames) - 1) // (num_images - 1)
    png_filename = os.path.join(subfolder_path, 'plots', vname + '_area_dist.png')
    images = []
    for k in range(0, len(frames), step):
        image = imageio.imread(frames[k])
        cv2.putText(image, f"{k * frame_skip / fps:.1f} s", (10, 440), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (0, 0, 0), 15)
        cv2.putText(image, f"{k * frame_skip / fps:.1f} s", (10, 440), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (255, 255, 255), 5)
        images.append(image)
    # Combine images side by side (horizontally)
    thumbnail = np.hstack(images)
    # Save the combined image
    thumbnail_rgb = cv2.cvtColor(thumbnail, cv2.COLOR_BGR2RGB)
    jpeg_filename = png_filename.replace('.png', '.jpg')
    cv2.imwrite(jpeg_filename, thumbnail_rgb, [cv2.IMWRITE_JPEG_QUALITY, 80]) 
    imageio.imwrite(png_filename, cv2.cvtColor(thumbnail_rgb, cv2.COLOR_RGB2BGR))

    # Clean up temporary frame files
    for filename in frames:
        os.remove(filename)
    # Release video capture
    cap.release()

    # Save areas data to JSON
    with open(output_json, 'w') as f:
        json.dump(areas_data, f, indent=2)

    print(f"Processing complete. Areas data saved to {output_json}")

def crop_frame(frame):
    # Convert the frame from BGR to HSV color space
    hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
    
    # Extract the Saturation channel (S channel)
    
    channel = hsv[:, :, 2]
    
    _, binary = cv2.threshold(channel, 120, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    kernel = np.ones((3, 3), np.uint8)
    dilated = cv2.dilate(binary, kernel, iterations=11)
    eroded = cv2.erode(dilated, kernel, iterations=91)
    dilated2 = cv2.dilate(eroded, kernel, iterations=80)
    
    # Perform Canny edge detection on the saturation channel
    edges = cv2.Canny(dilated2, 50, 150)
    
    # Find contours from edges
    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Assume the largest contour is the border
    if contours:
        largest_contour = max(contours, key=cv2.contourArea)
        # Get bounding rectangle of the largest contour
        return largest_contour
    else:
        # If no contours found, return the original frame
        return frame

def convert_to_grayscale(frame):
    if True:
        # Convert to HSV color space
        hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
        h, s, v = cv2.split(hsv)
        
        # Normalize S and V channels to [0,1] range
        s_norm = s.astype(np.float32) / 255.0
        v_norm = v.astype(np.float32) / 255.0

        # Nonlinear combination: Multiply S and V channels
        # Multiplying emphasizes strongly colored regions without losing whites
        sv_combined = np.sqrt(v_norm * (0.5 * v_norm + 0.5 * s_norm))
        
        # Scale to [0,255] and convert to uint8
        sv_combined = cv2.normalize(sv_combined, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        
        # Apply CLAHE for maximum local contrast
        clahe = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(16, 16))
        contrast_enhanced = clahe.apply(sv_combined)
    else:
        hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
        h, s, v = cv2.split(hsv)
        
        # Step 1: Extract red yarn using HSV color masking
        # Red typically wraps around hue=0 (thus two masks needed)
        lower_red1 = np.array([0, 5, 0])
        upper_red1 = np.array([25, 255, 255])
        lower_red2 = np.array([155, 5, 0])
        upper_red2 = np.array([180, 255, 255])
        
        mask_red1 = cv2.inRange(hsv, lower_red1, upper_red1)
        mask_red2 = cv2.inRange(hsv, lower_red2, upper_red2)
        red_mask = cv2.bitwise_or(mask_red1, mask_red2)
        
        # Optional: Slight blur to reduce noise
        red_mask = cv2.GaussianBlur(red_mask, (5, 5), 0)
        
        # Step 2: Enhance red yarn visibility using masked CLAHE
        # Combine S and V to keep both color and brightness information
        sv_combined = cv2.normalize((s.astype(np.float32) / 255.0) * (v.astype(np.float32) / 255.0), 
                                    None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        
        # Apply CLAHE to masked red region
        clahe = cv2.createCLAHE(clipLimit=4.0, tileGridSize=(16,16))
        clahe_red = clahe.apply(sv_combined)
        
        # Mask application: Only enhance red yarn regions
        enhanced_red = cv2.bitwise_and(clahe_red, clahe_red, mask=red_mask)
        
        # Step 3: Enhance non-red (white) yarn separately for clarity
        # White yarn (low saturation, high value):
        white_mask = cv2.inRange(hsv, np.array([0, 0, 120]), np.array([180, 255, 255]))
        clahe_white = clahe.apply(v)  # Enhance brightness channel directly
        enhanced_white = cv2.bitwise_and(clahe_white, clahe_white, mask=white_mask)
        
        # Step 4: Combine red and white enhanced images
        contrast_enhanced = cv2.max(enhanced_red, enhanced_white)

    return contrast_enhanced

def apply_threshold(gray_frame):
    # Apply adaptive thresholding or Otsu's method
    # Here, we use Otsu's method
    _, binary = cv2.threshold(
        gray_frame, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
    )
    # Invert the image so that fibers are white and background is black
    binary = cv2.bitwise_not(binary)
    return binary

def remove_noise(binary_frame):
    # Define the kernel size for erosion and dilation
    kernel = np.ones((3, 3), np.uint8)
    # Perform erosion
    eroded = cv2.erode(binary_frame, kernel, iterations=1)
    # Perform dilation
    dilated = cv2.dilate(eroded, kernel, iterations=1)
    return dilated

def bound_find(processed_frame, boundary):
    
    hull = cv2.convexHull(boundary)
    M = cv2.moments(hull)
    cx = int(M["m10"] / M["m00"]) if M["m00"] != 0 else 0
    cy = int(M["m01"] / M["m00"]) if M["m00"] != 0 else 0
    centroid = np.array([cx, cy])
    
    shrink_factor = border  # Adjust for more/less shrinkage
    
    # Shrink the hull by moving each point towards the centroid
    shrunken_hull = np.array([
        [[int((p[0][0] - cx) * shrink_factor + cx), int((p[0][1] - cy) * shrink_factor + cy)]]
        for p in hull
    ], dtype=np.int32)

    bound_mask = np.zeros_like(processed_frame)
    cv2.drawContours(bound_mask, [shrunken_hull], -1, 255, thickness=cv2.FILLED)
    
    return cv2.bitwise_not(bound_mask)
    
def extract_areas(original, processed, mask, div, time):
    # Find contours in the binary image
    contours, _ = cv2.findContours(processed, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
        
    valid_areas = []
    fig, ax = plt.subplots(figsize=(6 / 2.54, 6 / 2.54))
    ax.set_facecolor('white')
    
    filled_contours = np.zeros((processed.shape[0], processed.shape[1], 4), dtype=np.uint8) + 255
    filled_contours[:,:,:-1] = original
    max_area = 200
    for cnt in contours:
        # Create a mask for the current contour
        mask_current = np.zeros_like(processed)
        cv2.drawContours(mask_current, [cnt], -1, 255, thickness=cv2.FILLED)

        # Check for intersection with the shrunken contour mask
        intersection = cv2.bitwise_and(mask_current, mask)
        if not np.any(intersection):
            # No intersection; include this contour's area
            area = cv2.contourArea(cnt)
            if area > 0:
                valid_areas.append(area * div**2)
                color = cmap(np.log10(area * div**2 + 1) / np.log10(max_area))
                color_rgba = [int(c * 255) for c in color]
                color_bgra = color_rgba[2::-1] + [int(255*1)]  # B, G, R, A
                cv2.drawContours(filled_contours, [cnt], -1, color_bgra, thickness=cv2.FILLED)
                
    ax.imshow(cv2.cvtColor(filled_contours, cv2.COLOR_BGRA2RGBA))
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    frame_filename = os.path.join(subfolder_path, 'plots', f'frame_{time}.png')
    fig.savefig(frame_filename, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    return (valid_areas, frame_filename)

def process_mov_files(mov_files, subfolder_path, overwrite_flag, frame_skip):
    """
    Processes .MOV files in the given subfolder by running video processing and holes_analysis.py.
    
    Parameters:
        mov_files (list): List of .MOV file names.
        subfolder_path (str): Path to the subfolder containing the .MOV files.
        overwrite_flag (int): If 1, overwrite existing output; otherwise, skip existing outputs.
        frame_skip (int): Frame skipping parameter for processing.
    """
    for mov_file in mov_files:
        vname = os.path.splitext(mov_file)[0]

        # Build paths
        video_file = os.path.join(subfolder_path, mov_file)
        plots_folder = os.path.join(subfolder_path, 'plots')
        output_json_file = os.path.join(plots_folder, f'{vname}_holes.json')

        # Create plots directory if it doesn't exist
        os.makedirs(plots_folder, exist_ok=True)

        # Decide whether to process or skip
        if overwrite_flag == 1 or not os.path.exists(output_json_file):
            print(f'Processing {video_file}{" (overwriting existing outputs)" if overwrite_flag == 1 else ""}.')
            process_video(video_file, output_json_file, vname, frame_skip)

            # Run holes_analysis.py using srun
            #command = ['srun', '--mem=8G', 'python', 'holes_analysis.py', subfolder_path, vname]
            #print(f'Running holes_analysis.py for {vname}')
            #subprocess.run(command)
        else:
            print(f'Skipping {video_file}, output already exists.')

if __name__ == "__main__":

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Process videos in subfolders and extract holes areas.')
    parser.add_argument('folder_location', help='The folder that contains multiple subfolders. Each subfolder contains a MOV file')
    parser.add_argument('--overwrite', type=int, default=1, help='Set to 1 to overwrite existing outputs.')
    parser.add_argument('--frame_skip', type=int, default=100, help='Number of frames to skip.')
    parser.add_argument('--border', type=float, default=0.98, help='Remove border regions.')
    parser.add_argument('--skip', type=float, default=0, help='Number of initial frames to skip.')
    parser.add_argument('--xalign', type=int, default=0, help='Align border x.')
    parser.add_argument('--yalign', type=int, default=0, help='Align border y.')
    args = parser.parse_args()

    folder_location = args.folder_location
    overwrite_flag = args.overwrite
    frame_skip = args.frame_skip
    border = args.border
    skip_ini = args.skip
    xalign = args.xalign
    yalign = args.yalign

    # Find all subfolders
    subfolders = [f for f in os.listdir(folder_location) if os.path.isdir(os.path.join(folder_location, f))]
    filtered_subfolders = [f for f in subfolders if "plots" not in os.path.basename(f)]
    
    # If there are subfolders, iterate over them
    if filtered_subfolders:
        for pattern_folder in filtered_subfolders:
            subfolder_path = os.path.join(folder_location, pattern_folder)
            mov_files = [f for f in os.listdir(subfolder_path) if f.lower().endswith('.mov')]
            print("Processing:", subfolder_path, "| MOV files:", mov_files)
            process_mov_files(mov_files, subfolder_path, overwrite_flag, frame_skip)

    else:
        # No subfolders → Treat folder_location as the subfolder
        subfolder_path = folder_location
        mov_files = [f for f in os.listdir(subfolder_path) if f.lower().endswith('.mov')]
        process_mov_files(mov_files, subfolder_path, overwrite_flag, frame_skip)
        print("Processing:", folder_location, "| MOV files:", mov_files)