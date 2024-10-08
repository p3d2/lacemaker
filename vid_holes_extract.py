import cv2
import numpy as np
import json
import os
import matplotlib.pyplot as plt
import imageio.v2 as imageio

# Create the colormap
cmap = plt.cm.turbo

def process_video(video_path, output_json, frame_skip=1):
    # Initialize video capture
    cap = cv2.VideoCapture(video_path)
    frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    
    areas_data = []
    
    frame_idx = 0
    frame_c = 0
    frames = []
    while cap.isOpened():
        # Read frame
        ret, frame = cap.read()
        if not ret:
            break

        # Skip frames if necessary
        if frame_idx % frame_skip != 0:
            frame_idx += 1
            continue

        # Step 1: Crop the frame to remove borders
        px = 3
        if frame_c == 0:
            [x, y, w, h] = crop_frame(frame)
            norm_area = ((w + h)/(128*2)) ** 2 # convert vid area to sim area

        # Crop the frame to this rectangle
        cropped_frame = frame[y-px:y+h+px, x-px:x+w+px]

        # Step 2: Convert to grayscale using the value channel
        gray_frame = convert_to_grayscale(cropped_frame)
        
        # Step 3: Apply thresholding
        binary_frame = apply_threshold(gray_frame)

        # Step 4: Perform erosion and dilation
        processed_frame = remove_noise(binary_frame)
        
        # Step 5: Extract contours (holes) and compute areas
        
        if frame_c == 0: mask = bound_find(processed_frame)
        (areas, frame_n) = extract_areas(cropped_frame, processed_frame, mask, norm_area, frame_c)
        frames.append(frame_n)

        # Save areas for this frame
        frame_data = {
            'time': frame_c,
            'areas': areas
        }
        areas_data.append(frame_data)
        frame_idx += 1
        frame_c += 1

    # Generate GIF from frames
    gif_filename = os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_area_dist.gif')
    with imageio.get_writer(gif_filename, mode='I', fps=10, loop=0) as writer:
        for filename in frames:
            image = imageio.imread(filename)
            writer.append_data(image)

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
    
    saturation = hsv[:, :, 1]
    
    _, binary = cv2.threshold(saturation, 150, 255, cv2.THRESH_BINARY)
    kernel = np.ones((3, 3), np.uint8)
    dilated = cv2.dilate(binary, kernel, iterations=17)
    eroded = cv2.erode(dilated, kernel, iterations=15)
    
    # Perform Canny edge detection on the saturation channel
    edges = cv2.Canny(eroded, 50, 150)
    
    # Find contours from edges
    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Assume the largest contour is the border
    if contours:
        largest_contour = max(contours, key=cv2.contourArea)
        # Get bounding rectangle of the largest contour
        return cv2.boundingRect(largest_contour)
    else:
        # If no contours found, return the original frame
        return frame

def convert_to_grayscale(frame):
    # Convert to HSV color space
    hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
    # Extract the V channel (value)
    _, _, v = cv2.split(hsv)
    return v

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

def bound_find(processed_frame):
    
    contours, _ = cv2.findContours(processed_frame, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    boundary = max(contours, key=cv2.contourArea)
    hull = cv2.convexHull(boundary)
    bound_mask = np.zeros_like(processed_frame)
    cv2.drawContours(bound_mask, [hull], -1, 255, thickness=cv2.FILLED)
    kernel_erode =np.ones((21, 21), np.uint8)
    
    return cv2.bitwise_not(cv2.erode(bound_mask, kernel_erode, iterations=1))
    
def extract_areas(original, processed, mask, div, time):
    # Find contours in the binary image
    contours, _ = cv2.findContours(processed, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
        
    valid_areas = []
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_facecolor('white')
    
    filled_contours = np.zeros((processed.shape[0], processed.shape[1], 4), dtype=np.uint8) + 255
    filled_contours[:,:,:-1] = original
    max_area = 1000
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
                valid_areas.append(area/div)
                color = cmap(np.log10(area/div + 1) / np.log10(max_area))
                color_rgba = [int(c * 255) for c in color]
                color_bgra = color_rgba[2::-1] + [int(255*1)]  # B, G, R, A
                cv2.drawContours(filled_contours, [cnt], -1, color_bgra, thickness=cv2.FILLED)
                
    ax.imshow(cv2.cvtColor(filled_contours, cv2.COLOR_BGRA2RGBA))
    frame_filename = os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + f'frame_{time}.png')
    fig.savefig(frame_filename, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    return (valid_areas, frame_filename)

if __name__ == "__main__":
    # Example usage
    pattern_folder = 'Pattern_3025'
    dump_file = 'contract_'
    video_file = os.path.join('output', 'simulations', pattern_folder, '3025-01 20240905', 'DSC_0181.MOV')
    output_json_file = os.path.join('output', 'simulations', pattern_folder, 'plots', 'DSC_0181_holes.json')
    frame_skip = 100 

    process_video(video_file, output_json_file, frame_skip)