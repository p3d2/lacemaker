import os
from lammps import PyLammps

# Initialize PyLammps
L = PyLammps()
L.units("lj")
L.atom_style("molecular")
L.boundary("m m m")

# Set bond and angle styles
L.bond_style("harmonic")
L.angle_style("harmonic")

# Read data from file
folderdata = "output/lammps_data/pattern1024_0.5_1.0_0.0_30.0_5.0.data"  # Define your folder path here
L.read_data(f"{folderdata}")

# Define force field settings
L.pair_style("soft_exclude", 1.0, 16)
L.pair_coeff("* *", 200.0)
L.comm_modify("cutoff", 3.5)

# Output settings
folderdump = "test/"  # Define your output path here
filename = "pattern1024_0.5_1.0_0.0_30.0_5.0"  # Define your filename here
thermodump = 2500  # Define thermo output frequency
L.variable("filedump string", f"{folderdump}{filename}.lammpstrj")
L.variable("logdump string", f"{folderdump}{filename}.log")
L.compute("perAtomPE all pe/atom")
L.dump("myDump all custom", thermodump, "${filedump} id mol type x y z fx fy fz c_perAtomPE")
L.log("${logdump}")

# Define regions and groups

xmin, xmax, ymin, ymax = 16.5, 143.5, -143.5, -13.5  # Define dimensions as per your simulation box
L.region("left_block block INF", xmin, "INF INF INF INF units box")
L.region("right_block block", xmax, "INF INF INF INF INF units box")
L.region("bottom_block block INF INF INF", ymin, "INF INF units box")
L.region("top_block block INF INF", ymax, "INF INF INF units box")

L.group("left_wall region left_block")
L.group("right_wall region right_block")
L.group("bottom_wall region bottom_block")
L.group("top_wall region top_block")
L.group("boundary union left_wall right_wall bottom_wall top_wall")

# Apply fixes
L.fix("freezef boundary setforce 0.0 0.0 0.0")
visc = 0.2  # Define viscosity parameter
L.fix("dragU all viscous", visc)
L.fix("nve_all all nve")

# Minimization settings
L.minimize("1.0e-6 1.0e-7 10000 100000")
L.neigh_modify("every 100")

# Integrator settings
tstep = 0.01  # Define timestep
L.timestep(tstep)
L.thermo(thermodump)
L.thermo_style("custom step temp pe ke etotal")

# Run simulation
nframes = 10  # Define number of frames
L.variable("iter equal", f"{thermodump}*{nframes}")
filerestart = f"{folderdump}sim1.restart"
L.run("${iter}")
L.write_restart(f"{filerestart}")
