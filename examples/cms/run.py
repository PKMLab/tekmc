from importlib.util import find_spec
import os
import sys

required_packages = [
    "mdtraj",
    "scipy",
    "sklearn",
    "matplotlib",
    "numpy",
    "p_tqdm",
]

missing_packages = []
for package in required_packages:
    spec = find_spec(package)
    if spec is None:
        missing_packages.append(package)

if len(missing_packages):
    print("These packages are required:")
    print(4 * " ", *missing_packages)
    print("Try installing with:")
    for i in missing_packages:
        print(
            f"    python3 -m pip install {i}\n"
            f"    conda install -c conda-forge {i}\n",
            flush=True,
        )
    sys.exit()

from datetime import datetime
start_date = datetime.now()
print(f"Started on {start_date.ctime()}")

from inputs import *
sys.path.append(path_to_src)

from imports import *
from tekmc import TEKMC
from msd import mean_squared_displacement

if (
    ".pdb" in trajectory_file
    or ".h5" in trajectory_file
    or ".lh5" in trajectory_file
):
    topology_file = -1

if None in [
    topology_file,
    trajectory_file,
    atom_name,
    timestep,
    spacings_list,
    n_walks,
    n_components,
]:
    print(
        "Required: Topology, trajectory as supported by \033[1mmdtraj\033[0m",
        "          Atom name for the species of interest",
        "          Time interval between trajectory frames",
        "          List of spacings to perform TEKMC",
        "          Number of walks during TEKMC for a given grid size",
        "  Units:  trajectory - nm",
        "          timestep   - ns",
        "          grid_size  - nm",
        "Also modify other parameters to your specific case. Look up https://github.com/Thanush1/TEKMC.git for details about all the methods.",
        sep="\n",
        flush=True,
    )
    sys.exit()

def get_parameters(line_no):
    file = open(os.path.join(path_to_src,'defaults.py'),'r')
    all_lines = file.readlines()
    value = all_lines[line_no - 1].split()[2]
    value = eval(value)
    return value

# Creating TEKMC object.
tekmc = TEKMC(
    topology_file,
    trajectory_file,
    atom_name,
    timestep,
    stride=stride,
    md_file=md_file,
    n_cpus=n_cpus,
    cmap=cmap,
    save_dir=save_dir,
    threshold=get_parameters(3),
    symmetrization_type=get_parameters(4),
    hopping_distance_analysis=hopping_distance_analysis,
)

# Various methods of TEKMC object.
# Modify the parameters as required. Look up https://github.com/Thanush1/TEKMC.git for more details.

spacings_list = sorted(spacings_list, key = lambda x: x[0] if type(x) == list else x, reverse=True)

for spacings in spacings_list:
    tekmc.TEKMC(
        spacings,
        n_walks,
        n_steps=n_steps,
        save_prob=get_parameters(12),
        create_prob=get_parameters(13),
        traj_steps=get_parameters(14),
        traj_n=get_parameters(15),
        traj_colors=get_parameters(16),
        save_random_walks=save_random_walks,
    )

    tekmc.plot_histogram(spacings, n_bins=get_parameters(37))

    tekmc.plot_trajectory(
        spacings,
        plot_steps=get_parameters(24),
        n_traj=get_parameters(25),
        _colors=get_parameters(26),
        unwrapped=get_parameters(27),
        trajectories=get_parameters(28),
        tekmc_file=get_parameters(29),
    )

    tekmc.plot_pathways(
        spacings,
        n_components,
        _colors=get_parameters(45),
        Threshold=get_parameters(46),
    )
    tekmc.plot_voxels_probability(spacings, color_map=get_parameters(54))


tekmc.compute_diffusion(spacings_list=get_parameters(62))
tekmc.plot_msd(spacings_list=get_parameters(70), _colors=get_parameters(71), y_scale = 'log', x_scale = 'log')
tekmc.plot_diffusion(spacings_list=get_parameters(79), _colors=get_parameters(80), y_scale = 'log', x_scale = 'log')
tekmc.normalized_MSE(spacings_list=get_parameters(97))

#tekmc.plot_displacement_probability(spacings_list=get_parameters(88), _colors=get_parameters(89))

end_date = datetime.now()

print(f"Finished on {end_date.ctime()}")
print(f"Total runtime {(end_date - start_date)}")
