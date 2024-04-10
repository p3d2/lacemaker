import os
import numpy as np
from PIL import Image, ImageChops, ImageFilter, ImageSequence

# Variables
input_folder = os.environ['WRKDIR'] + "/lacemaker/temp"
output_file = os.environ['WRKDIR'] + "/lacemaker/temp/out.gif"

# Process images
frames = []
file_list = sorted([f for f in os.listdir(input_folder) if f.endswith('.png')])

for filename in file_list:

    img_path = os.path.join(input_folder, filename)
    img = Image.open(img_path).convert('RGBA')

    # Convert the image into numpy array
    pixels = np.array(img)

    # Set non transparent pixels to 204 (0.2 * 255) gray and transparent ones black
    pixels[..., :3][pixels[..., 3] > 0] = [51, 51, 51]
    result = Image.fromarray(pixels, 'RGBA')

    # Expand non transparent pixels
    result = result.filter(ImageFilter.MaxFilter(3))
    final_img = Image.alpha_composite(result, img)

    # For the GIF
    final_img = final_img.convert('RGB').convert('P', palette=Image.ADAPTIVE)
    frames.append(final_img)

# Create GIF
durations = [10] * (len(frames) - 1) + [1000]  # Last frame has duration of 1000.
frames[0].save(output_file, save_all=True, append_images=frames[1:], transparency=255, disposal=2, loop=0, duration=durations)