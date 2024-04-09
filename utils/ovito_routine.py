import os
from ovito.io import import_file
from ovito.pipeline import Pipeline, FileSource
from ovito.modifiers import LoadTrajectoryModifier, ColorCodingModifier
from ovito.vis import Viewport, TachyonRenderer

path_data = os.path.join('/scratch','work','silvap1','lacemaker','output','lammps_data','pattern3023_0.25_1.0_2.5_30.0_50.0.data')
path_trj = os.path.join('/scratch','work','silvap1','lacemaker','output','simulations','Pattern_3023','pattern3023_0.25_1.0_2.5_30.0_50.0_contract_fix_1_2_3_4_5_6_7_8.lammpstrj')
path_temp = os.path.join('/scratch','work','silvap1','lacemaker','temp')
w, h = 800, 800

# Create a new ovito.pipeline.FileSource to load data from a LAMMPS dump file
pipeline = import_file(path_data, atom_style='molecular')
pipeline.add_to_scene()

# Access and modify the particles' visual appearance
particles_vis = pipeline.source.data.particles.vis
particles_vis.radius = 0.5  # Set the default radius for particles
bonds_vis = pipeline.source.data.particles.bonds.vis
bonds_vis.width = 1.0

# Load a lammps trajectory file
traj_mod = LoadTrajectoryModifier()
traj_mod.source.load(path_trj)
pipeline.modifiers.append(traj_mod)

# ***Force an initial pipeline computation***
pipeline.compute(0)  # This ensures the pipeline is updated with the loaded trajectory
# Create a viewport
viewport = Viewport()
viewport.type = Viewport.Type.ORTHO
viewport.type = Viewport.Type.TOP
viewport.zoom_all()  # Ensure we are viewing all particles

# Add AssignColorModifier to the pipeline that colors particles based on potential energy
modifier = ColorCodingModifier(
    property='c_peratompe',
    gradient = ColorCodingModifier.Hot()
)
pipeline.modifiers.append(modifier)

# Now safely access the number of frames
n_frames = traj_mod.source.num_frames
for frame in range(n_frames):
    # Compute the pipeline for the current frame
    data = pipeline.compute(frame)

    # Specify the output image file path for this frame
    image_file_path = os.path.join(path_temp, f'frame_{frame:04d}.png')

    # Render and save the image
    viewport.render_image(filename=image_file_path, size=(w, h), alpha=True, renderer=TachyonRenderer())

print("Rendering of all frames completed.")