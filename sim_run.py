import os
from lammps import PyLammps

os.environ['OMP_NUM_THREADS'] = '4'

# Input

folderdata = "output/lammps_data/pattern1024_0.5_1.0_0.0_30.0_5.0.data"  # Define your folder path here
folderdump = "test/"  # Define your output path here
filename = "pattern1024_0.5_1.0_0.0_30.0_5.0"  # Define your filename here
thermodump = 10000  # Define thermo output frequency
xmin, xmax, ymin, ymax = 1.5, 63.0, 1.5, 63.0  # Define dimensions as per your simulation box
visc1 = 0.1  # Define viscosity parameter for simulation 1
visc2 = 0.05  # Define viscosity parameter for simulation 2
tstep = 0.005  # Define timestep
nframes1 = 1  # Define number of frames for simulation 1
nframes2 = 19  # Define number of frames for simulation 2
kb = 30.0
r0 = 0.5  # Initial bond length
dr0 = 0.01  # Change in bond length per frame

# Initialize PyLammps
L = PyLammps()
L.units("lj")
L.atom_style("molecular")
L.boundary("m m m")

# Set bond and angle styles
L.bond_style("harmonic")
L.angle_style("harmonic")

# Read data from file
L.read_data(f"{folderdata}")

#for k in range(6):
#    L.bond_coeff(f"{k+1} {kb} 0.4")

# Define force field settings
L.pair_style("soft_exclude", 1.0, 8)
L.pair_coeff("* *", 100.0)
L.comm_modify("cutoff", 3.5)

# Output settings
L.variable("filedump string", f"{folderdump}{filename}.lammpstrj")
L.variable("logdump string", f"{folderdump}{filename}.log")
L.compute("perAtomPE all pe/atom")
L.dump("myDump all custom", thermodump, "${filedump} id mol type x y z")
L.log("${logdump}")

# Define regions and groups
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
L.fix("dragU all viscous", visc1)
L.fix("nve_all all nve")

# Minimization settings
L.minimize("1.0e-6 1.0e-7 10000 100000")
L.neigh_modify("every 100")

# Integrator settings
L.timestep(tstep)
L.thermo(thermodump)
L.thermo_style("custom step temp pe ke etotal")

# Run 1st simulation -> Stabilization
L.variable("iter equal", f"{thermodump*nframes1}")
filerestart = f"{folderdump}sim1.restart"
L.run("${iter}")

# Run 2nd simulation -> Contraction
L.unfix("dragU")
L.fix("dragU all viscous", visc2)

for i in range(nframes2):
    r0diff = r0 - (i+1) * dr0
    L.command(f"bond_coeff {1} {kb} {r0diff}")
    L.command(f"bond_coeff {2} {kb} {r0diff}")
    L.run(thermodump)