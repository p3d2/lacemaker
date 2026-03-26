import os, glob
import numpy as np
import argparse
import matplotlib.cm as cm
from ovito.io import import_file
from ovito.pipeline import Pipeline
from ovito.modifiers import ColorCodingModifier
from ovito.vis import Viewport, TachyonRenderer
from wand.image import Image as Image
from wand.color import Color
from wand.drawing import Drawing

def process_lammps_data(folder, pattern, file_pattern, w, h, g, border, e_min, e_max, particles_r, gradient):
    output_dir = os.path.join('output', 'simulations', folder)
    temp_dir = os.path.join('temp')
    gifs_dir = os.path.join(output_dir, 'gifs')

    os.makedirs(gifs_dir, exist_ok=True)
    path_data = os.path.join('output', 'lammps_data', pattern + '.data')
    path_trj = glob.glob(os.path.join(output_dir, file_pattern))

    pipeline = import_file(path_data, atom_style='molecular')
    pipeline.add_to_scene()

    # Not currently used since trajectory is loaded to inputs and loses bond info
    particles_vis = pipeline.source.data.particles.vis
    particles_vis.radius = 0.5  # Set the default radius for particles
    bonds_vis = pipeline.source.data.particles.bonds.vis
    bonds_vis.width = 1.0

    color_mod = ColorCodingModifier(
        property='c_peratompe',
        gradient=ColorCodingModifier.Gradient(gradient),
        start_value=e_min,
        end_value=e_max
    )
    pipeline.modifiers.append(color_mod)

    for traj_file in path_trj:
        gif_path = os.path.join(gifs_dir, 'pe_' + os.path.splitext(os.path.basename(traj_file))[0] + '.gif')
        if os.path.exists(gif_path):
            continue

        #traj_mod = LoadTrajectoryModifier() <- implementation that uses trj files as trajectory modifier and not Load (preserves bonds)
        #traj_mod.source.load(path_trj[k])
        #pipeline.modifiers.append(traj_mod)
        pipeline.source.load(traj_file)
        pipeline.compute()

        particles_vis = pipeline.source.data.particles.vis
        particles_vis.radius = particles_r  # Set the default radius for particles

        viewport = Viewport()
        viewport.type = Viewport.Type.TOP
        viewport.zoom_all(size=(w, h))

        frames = []
        for frame in range(pipeline.source.num_frames):
            image_file_path = os.path.join(temp_dir, f'frame_{frame:04d}.png')
            rend_engine = TachyonRenderer(ambient_occlusion_brightness=0.5, direct_light_intensity=2.0)
            viewport.render_image(filename=image_file_path, size=(w, h), alpha=True, renderer=rend_engine, frame=frame)
            frames.append(image_file_path)

        prepare_and_save_gif(frames, gif_path, g, border)

def process_frame(img_path, g, border):
    with Image(filename=img_path) as original:
        # Clone the original image to create the grayscale effect
        with Image(width=original.width, height=original.height, background=Color(f'rgb({int(g * 255)}, {int(g * 255)}, {int(g * 255)})')) as gray:

            # Copy the alpha channel from the original image to the gray image
            gray.composite(original, operator='copy_alpha', left=0, top=0)

            # Applying a maximum filter to simulate the border effect using morphology
            gray.morphology(method='dilate', kernel='disk:'+str(border))

            # Composite the gray image under the original image
            # The 'dst_over' operator places the original on top of the gray image
            original.composite(gray, left=0, top=0, operator='dst_over')
        
        # Ensure the image is in the proper format before returning
        original.format = 'png'
        
        return original.make_blob()

def prepare_and_save_gif(frame_paths, path, g, border):

    # Load all frames and coalesce
    with Image() as gif:
        gif.coalesce()
        background_color = Color('transparent')
        with Image(width=800, height=800, background=Color('transparent')) as frame:
            alpha_channel = 'activate'  # Ensure it's transparent
            frame.delay = 1
            frame.alpha_channel = 'activate'
            frame.background_color = Color('transparent')
            frame.transparent_color = Color('transparent')
            frame.dispose = 'background'  # Ensuring each frame starts with a clean slate if needed
            gif.sequence.append(frame.clone())  # Add as the first frame
        for img_path in frame_paths:
            blob = process_frame(img_path, g, border)
            with Image(blob=blob) as frame:
                frame.delay = 10
                frame.alpha_channel = 'activate'
                frame.background_color = Color('transparent')
                frame.transparent_color = Color('transparent')
                frame.dispose = 'previous'  # Ensuring each frame starts with a clean slate if needed
                gif.sequence.append(frame.clone())
        
        gif.sequence[-1].delay = 100  # Extended delay for the last frame
        # Coalesce frames to ensure each frame is a full snapshot

        gif.loop = 0
        gif.type = 'palette'
        gif.transparent_color = Color('transparent')  # Set the transparent color

        del gif.sequence[0]
        # Save the GIF
        gif.save(filename=path)

def main():
    parser = argparse.ArgumentParser(description='Process LAMMPS data files.')
    parser.add_argument('--folder', required=True, help='Folder name for the simulation output')
    parser.add_argument('--pattern', required=True, help='Pattern for the data file naming')
    parser.add_argument('--file_pattern', required=True, help='Pattern for the trajectory files')
    parser.add_argument('--width', type=int, default=800, help='Viewport width')
    parser.add_argument('--height', type=int, default=800, help='Viewport height')
    parser.add_argument('--gray', type=float, default=0.1, help='Gray scale value for transparency')
    parser.add_argument('--border', type=int, default=3, help='Border size for non-transparent pixels')
    parser.add_argument('--emin', type=float, default=0.0, help='Minimum of the colorbar scale')
    parser.add_argument('--emax', type=float, default=0.2, help='Maximum of the colorbar scale')
    parser.add_argument('--radius', type=float, default=0.5, help='Radius of the fibers')
    parser.add_argument('--cmap', type=str, default='rainbow', help='Colormap name from matplotlib')

    args = parser.parse_args()
    cmap = cm.get_cmap(args.cmap)
    gradient = [cmap(x)[:3] for x in np.linspace(0, 1, cmap.N)]

    process_lammps_data(
        folder=args.folder,
        pattern=args.pattern,
        file_pattern=args.file_pattern,
        w=args.width,
        h=args.height,
        g=args.gray,
        border=args.border,
        e_min=args.emin,
        e_max=args.emax,
        particles_r=args.radius,
        gradient=gradient
    )

if __name__ == "__main__":
    main()