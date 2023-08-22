import sys
import os

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), "src"))

from imports import *
from tekmc import TEKMC
from msd import mean_squared_displacement


# Creating the TEKMC object.

# Parameters
topology_file = 'bulk_spce.gro'
trajectory_file = 'bulk_spce.dcd'
atom_name = 'O'
timestep = 1e-6 * 1000

# Optional parameters
stride = 10
md_file = None
symmetrization_type = "average"
threshold = 0.0
n_cpus = cpu_count()
cmap = "terrain"

if None in [topology_file, trajectory_file, atom_name, timestep]:
    print(
        "Required: Topology, trajectory as supported by \033[1mmdtraj\033[0m",
        "          Atom name for the species of interest",
        "          Time interval between trajectory frames",
        "  Units:  trajectory - A",
        "          timestep   - ns",
        "Also modify other parameters to your specific case. Look up https://github.com/Thanush1/TEKMC.git for details about all the methods.",
        sep="\n",
        flush=True,
    )
    sys.exit()

# Creating TEKMC object.
tekmc = TEKMC(
    topology_file,
    trajectory_file,
    atom_name,
    timestep,
    stride=stride,
    symmetrization_type=symmetrization_type,
    md_file=md_file,
    threshold=threshold,
    n_cpus=n_cpus,
    cmap=cmap,
)

# Various methods of TEKMC object.
# Modify the parameters as required. Look up https://github.com/Thanush1/TEKMC.git for more details.

spacings_list = [[0.25],[0.26],[0.27]]
save_dir = None

for spacings in spacings_list:
    n_walks = 2000
    n_steps = 3000
    tekmc.TEKMC(
        spacings,
        n_walks,
        n_steps=n_steps,
        save_prob=True,
        create_prob=False,
        traj_steps=1000,
        traj_n=3,
        traj_colors=None,
        save_dir=None,
    )
    n_steps = 1000
    tekmc.plot_trajectory(
        spacings,
        n_steps,
        n_traj=3,
        _colors=None,
        save_dir=None,
        unwrapped=True,
        trajectories=None,
        tekmc_file=None,
    )
    tekmc.plot_Hist(spacings, n_bins=100, save_dir=None)
    n_components = 3
    tekmc.plot_connected_components(
        spacings,
        n_components,
        _colors=None,
        save_dir=None,
        Threshold=5e-2,
    )
    tekmc.plot_all_voxels(spacings, save_dir=None, color_map="terrain")

tekmc.find_D(spacings_list=None, save_dir=None)
tekmc.plot_MSDs(spacings_list=None, _colors=None, save_dir=None)
tekmc.plot_Diffusion_coefficient(spacings_list=None, _colors=None, save_dir=None)
tekmc.normalized_MSE(save_dir=None)
tekmc.plot_prob_dist_velocity(spacings_list=None, save_dir=None, _colors=None)

