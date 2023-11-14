from importlib.util import find_spec
import os
import sys

required_packages = [
    "mdtraj",
    "scipy",
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
    print(4*" ", *missing_packages)
    print("Try installing with:")
    for i in missing_packages:
        print(
            f"    python3 -m pip install {i}\n"
            f"    conda install -c conda-forge {i}\n",
            flush=True,
        )
    sys.exit()

from p_tqdm import p_map
from matplotlib import cm
from datetime import datetime
from collections import defaultdict
from scipy.optimize import curve_fit
from matplotlib.colors import Normalize
from scipy.sparse import csr_matrix, save_npz, load_npz, vstack

import numpy as np
import mdtraj as md
import matplotlib as mpl
import matplotlib.pyplot as plt

from msd import mean_squared_displacement, section, msd_fft

class TEKMC:
    def __init__(
        self,
        topology_file,
        trajectory_file,
        atom_name,
        timestep,
        stride=1,
        md_file=None,
        n_cpus=None,
        cmap="terrain",
        save_dir=None,
        symmetrization_type="average",
        threshold=0.0,
        hopping_distance_analysis=False,
    ):
        # key will be used various save filenames
        key = trajectory_file.split(os.sep)[-1].split(".")[0]
        self.key = key
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.atom_name = atom_name
        self.stride = stride
        self.timestep = timestep  # in ns
        self.symmetrization_type = symmetrization_type
        self.threshold = threshold
        self.md_file = md_file
        self.cmap = cmap
        self.nm_to_A = 10
        self.ns_to_ps = 1e3

        if n_cpus is None:
            n_cpus = os.cpu_count()

        if save_dir is None:
            save_dir = os.path.join("outputs")

        save_dir = os.path.join(save_dir, key)

        os.makedirs(save_dir, exist_ok=True)
        self.n_cpus = n_cpus
        self.save_dir = save_dir
        self.hopping_distance_analysis = hopping_distance_analysis

        self.load_trajectory()

        if self.md_file is None:

            self.find_MSD()  # Estimate MSD from MD trajectory

            self.md_file = os.path.join(
                save_dir, f"{self.key}_{self.atom_name}_msd.dat"
            )  # MSD of MD trajectory is saved in this file

    def load_trajectory(self):
        # Load and analyse the MD trajectory
        print(f"Topology:")
        print(f"{self.topology_file.split(os.sep)[-1]}")
        print(f"Trajectory:")
        print(f"{self.trajectory_file.split(os.sep)[-1]}", flush=True)

        if (
            ".pdb" in self.trajectory_file
            or ".h5" in self.trajectory_file
            or ".lh5" in self.trajectory_file
        ):
            traj = md.load(self.trajectory_file)
        else:
            traj = md.load(self.trajectory_file, top=self.topology_file)

        top = traj.topology

        atom_name_set = frozenset(atom.name for atom in top.atoms)

        if self.atom_name not in atom_name_set:
            print("The available atoms are:")
            print(f"\n".join(map(str, atom_name_set)))
            raise ValueError(f"Topology does not contain {atom_name}.")

        n_atoms = traj.n_atoms
        n_frames = traj.n_frames
        unit_cell = traj.unitcell_lengths[1]
        total_time = self.timestep * n_frames
        cell_volume = unit_cell[0] * unit_cell[1] * unit_cell[2]

        section("System Info")

        print(
            f"{'Box Length =':>15s} {' '.join(map(str, np.round(self.nm_to_A * unit_cell, 2)))} A",
            f"{'Volume =':>15s} {self.nm_to_A**3 * cell_volume:.2f} A^3",
            f"{'Atoms =':>15s} {n_atoms}",
            f"{'Frames =':>15s} {n_frames}",
            f"{'Timestep =':>15s} {self.timestep:.4g} ns",
            f"{'Total time =':>15s} {total_time:.4g} ns",
            sep="\n",
            flush=True,
        )

        def geometry(group):
            # Estimate the Center of Mass of a given group of atoms

            xlo = 999999
            xhi = -xlo
            ylo = xlo
            yhi = xhi
            zlo = xlo
            zhi = xhi

            atoms = len(group)

            for i in range(atoms):

                tempx = group[i][0]
                tempy = group[i][1]
                tempz = group[i][2]

                xlo = min(tempx, xlo)
                ylo = min(tempy, ylo)
                zlo = min(tempz, zlo)

                xhi = max(tempx, xhi)
                yhi = max(tempy, yhi)
                zhi = max(tempz, zhi)

            geoCenter = [(xlo + xhi) / 2, (ylo + yhi) / 2, (zlo + zhi) / 2]

            return atoms, xlo, xhi, ylo, yhi, zlo, zhi, np.asarray(geoCenter)

        particles = [
            traj.xyz[0, atom.index] for atom in top.atoms_by_name(self.atom_name)
        ]
        atoms, xlo, xhi, ylo, yhi, zlo, zhi, center = geometry(particles)

        # Center of mass of the particles in the first frame is considered the origin of the MD box
        origin = center
        self.xlo = origin[0]
        self.xhi = self.xlo + unit_cell[0]
        self.ylo = origin[1]
        self.yhi = self.ylo + unit_cell[1]
        self.zlo = origin[2]
        self.zhi = self.zlo + unit_cell[2]
        self.traj = traj
        self.top = top

    def get_ns(self, x_spacing, y_spacing, z_spacing):
        # Estimate and return the number of voxels along each dimension for given spacings

        xlo, xhi, ylo, yhi, zlo, zhi = (
            self.xlo,
            self.xhi,
            self.ylo,
            self.yhi,
            self.zlo,
            self.zhi,
        )

        n_x = (xhi - xlo) / x_spacing
        n_y = (yhi - ylo) / y_spacing
        n_z = (zhi - zlo) / z_spacing

        if n_x != int(n_x):
            n_x = int(n_x) + 1
        if n_y != int(n_y):
            n_y = int(n_y) + 1
        if n_z != int(n_z):
            n_z = int(n_z) + 1

        n_x = int(n_x)
        n_y = int(n_y)
        n_z = int(n_z)

        return n_x, n_y, n_z

    def format(self, value):
        # Format a float

        return f"{value:.4f}"

    def generate_colors(self, n):
        # Generate n distinct colors from the given color map

        vals = np.linspace(0.1, 0.9, n)
        cmap = plt.get_cmap(self.cmap)
        colors = [cmap(x) for x in vals]
        return colors

    def powerLaw(self, x, a, nu):
        # Function for finding the exponent of x dependence

        return a * x**nu

    def linearMSD(self, x, D, c):
        # Function for linear fit

        return D * x + c

    def find_MSD(self):
        # Estimate MSD from topology and trajectory files of MD
        topology_file = self.topology_file
        trajectory_file = [self.trajectory_file]
        timestep = self.timestep  # ns
        atom_name = self.atom_name
        mean_squared_displacement(
            topology_file,
            trajectory_file,
            atom_name,
            timestep,
            start=0,
            end=-1,
            stride=1,
            save_dir=self.save_dir,
            n_cpus=self.n_cpus,
            data_type=np.float16,
            plot=True,
            system_info=False,
        )
        print("\n")

    def compute_diffusion(self, spacings_list=None):
        # Estimate and write the diffusion coefficients to a file for the given list of spacings

        save_dir = self.save_dir

        n = len([x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")])
        if spacings_list:
            n = len(spacings_list)

        if spacings_list is None:
            msd_files = [
                x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")
            ]
        else:
            msd_files = []
            for spac in spacings_list:
                if type(spac) != list:
                    spac = [spac]
                if len(spac) == 1:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[0], spac[0]
                    msd_files.append(f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat")
                else:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[1], spac[2]
                    msd_files.append(
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    )

        msd_files = sorted(msd_files)

        directory = "msd"
        with open(
            os.path.join(save_dir, "Diffusion_coefficients.dat"), "w"
        ) as diff_coeff:
            print("Estimating Diffusion Coefficients ")

            for file in msd_files:

                with open(os.path.join(save_dir, directory, file), "r") as f:
                    msd = []
                    time = []
                    comment = f.readline()
                    while True:
                        line = list(f.readline().split())
                        if line == []:
                            break
                        msd.append(float(line[1]))
                        time.append(float(line[0]))

                msd = np.array(msd)
                time = np.array(time)
                [D, c], pcov = curve_fit(
                    self.linearMSD, time[len(msd) // 4 :], msd[len(msd) // 4 :]
                )  # linear fit to estimate the Diffusion coefficient
                D /= 6
                D *= 10 ** (-7)  # converting to cm^2/s
                D_SD = np.sqrt(np.diag(pcov)[0])
                D_SD *= 10 ** (-7)  # converting to cm^2/s

                [a, nu], pcov = curve_fit(
                    self.powerLaw,
                    time[len(msd) // 2 :],
                    msd[len(msd) // 2 :],
                    p0=[1e2, 1],
                )
                nu_SD = np.sqrt(np.diag(pcov)[1])

                n = f"{self.key}_TEKMC_MSD.dat".count("_")
                spacing = "_".join([x[3:] for x in file.split("_")[n + 1 :]])[:-4]

                diff_coeff.write(f"spacing: {spacing}    ")
                diff_coeff.write(f"Diffusion (D) = {D:.4e} +/- {D_SD:.2e} cm^2/s    ")
                diff_coeff.write(f"Exponent (nu) = {nu:.2f} +/- {nu_SD:.2e}\n")

    def TEKMC(
        self,
        spacings,
        n_walks,
        n_steps=None,
        save_prob=True,
        create_prob=False,
        traj_steps=1000,
        traj_n=3,
        traj_colors=None,
        save_random_walks=False,
    ):
        """Main method that performs TEKMC on the loaded trajectories for the given spacings,
        number of steps, number of random walks and creates MSD file, plot of few longest random walks, etc"""

        section(f"TEKMC: Grid Size {spacings} nm")
        start_date = datetime.now()
        print(f"Started on {start_date.ctime()}", flush=True)

        if type(spacings) != list:
            spacings = [spacings]

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        if n_walks < 10:
            print(f"*** WARNING***    {n_walks = } is too low. Defaulting to 10")
            n_walks = 10

        equal_spacings = len(spacings) == 1

        if n_steps is None:
            n_steps = self.traj.n_frames * 3  # default time for TEKMC is 3 * MD time
        
        if n_steps < 100:
            print(f"*** WARNING***    {n_steps = } is too low. Defaulting to 100")
            n_steps = 100

        if n_steps // self.stride < 100 and save_random_walks:
            print("*** WARNING ***    Trajectory to be saved is very small")
            print("*** WARNING ***    Consider saving at least 100 timesteps")

        save_dir = self.save_dir

        print("Generating voxels", flush=True)

        voxels = {}  # map of the voxel number to the center of the voxel

        n_x, n_y, n_z = self.get_ns(x_spacing, y_spacing, z_spacing)

        ind = 0
        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    voxels[(i, j, k)] = ind
                    ind += 1

        N = (n_x) * (n_y) * (n_z)

        os.makedirs(os.path.join(save_dir, "probability_matrices"), exist_ok=True)
        try:
            if create_prob:  # user request for not loading from saved files
                raise FileNotFoundError("On request...")
            if equal_spacings:
                prob = load_npz(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"probSymm_{self.key}_rgx{x_spacing}.npz",
                    )
                )
            else:
                prob = load_npz(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"probSymm_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.npz",
                    )
                )
            print("Loading saved probability matrix")

        except FileNotFoundError:
            print("Creating probability matrix")
            visits = [defaultdict(int) for i in range(N)]

            def map_to_voxel(coord):  # map from coordinate to voxel
                x = coord[0] - self.xlo
                y = coord[1] - self.ylo
                z = coord[2] - self.zlo
                i = x % unit_cell[0]
                j = y % unit_cell[1]
                k = z % unit_cell[2]

                i = i // x_spacing
                j = j // y_spacing
                k = k // z_spacing

                return voxels[(i, j, k)]

            n_frames = self.traj.n_frames
            unit_cell = self.traj.unitcell_lengths[1]

            for atom in self.top.atoms_by_name(self.atom_name):
                ind = atom.index
                for frame in range(0, n_frames - 1):

                    coord_current = np.asarray(self.traj.xyz[frame, ind])
                    coord_next = np.asarray(self.traj.xyz[frame + 1, ind])
                    voxel_num1 = map_to_voxel(coord_current)
                    voxel_num2 = map_to_voxel(coord_next)
                    visits[voxel_num1][voxel_num2] += 1

            print("Artificial symmetrization", flush=True)

            # Artificial symmetrization of rate of transitions to ensure detailed balance
            for i in range(N):
                visited_voxels = list(visits[i].keys())
                for j in visited_voxels:
                    if self.symmetrization_type == "max":
                        value = max(visits[i][j], visits[j][i])
                    elif self.symmetrization_type == "min":
                        value = min(visits[i][j], visits[j][i])
                    else:  # default is average
                        value = sum([visits[i][j], visits[j][i]]) // 2
                    visits[i][j] = value
                    visits[j][i] = value

            for i in range(N):
                value = sum(visits[i].values())
                if value != 0:
                    visited_voxels = list(visits[i].keys())
                    for j in visited_voxels:
                        visits[i][j] /= value

            # Removing entries with probability < self.threshold
            if self.threshold != 0:
                print(
                    "Removing all entries with value < threshold in the probability matrix"
                )
                for i in range(N):
                    visited_voxels = list(visits[i].keys())
                    for j in visited_voxels:
                        if visits[i][j] < self.threshold:
                            visits[i][j] = 0.0

                # Rescaling the remanining probabilites so that sum(prob[i]) remains 1
                for i in range(N):
                    value = sum(visits[i].values())
                    if value != 0:
                        visited_voxels = list(visits[i].keys())
                        for j in visited_voxels:
                            visits[i][j] /= value

            # Saving the probability matrix
            if save_prob:
                prob = None
                for i in range(N):
                    val = np.zeros(shape=(1, N))
                    for j in visits[i].keys():
                        val[0][j] = visits[i][j]
                    val = csr_matrix(val)
                    if prob is None:
                        prob = val
                    else:
                        prob = vstack([prob, val])

                prob = csr_matrix(prob)
                if equal_spacings:
                    save_npz(
                        os.path.join(
                            save_dir,
                            "probability_matrices",
                            f"probSymm_{self.key}_rgx{x_spacing}",
                        ),
                        prob,
                    )
                else:
                    save_npz(
                        os.path.join(
                            save_dir,
                            "probability_matrices",
                            f"probSymm_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}",
                        ),
                        prob,
                    )
                print("Probability matrix saved")

        # Random walks
        indices = {}  # a map from indices to voxels.
        print("Mapping volxels to indices", flush=True)
        for ele in voxels:
            indices[voxels[ele]] = ele

        try:
            if create_prob:  # User request for not loading from saved files
                raise FileNotFoundError("On request...")
            if equal_spacings:
                visited = np.loadtxt(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"visited_{self.key}_rgx{x_spacing}.dat",
                    ),
                    dtype=np.uint32,
                )
            else:
                visited = np.loadtxt(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"visited_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat",
                    ),
                    dtype=np.uint32,
                )

        except OSError:
            print("Indexing visited voxels", flush=True)
            visited = [
                i for i in range(0, N) if round(sum(prob[i].toarray()[0]), 1) == 1.0
            ]

            if save_prob:
                if equal_spacings:
                    np.savetxt(
                        os.path.join(
                            save_dir,
                            "probability_matrices",
                            f"visited_{self.key}_rgx{x_spacing}.dat",
                        ),
                        visited,
                        fmt="%d",
                    )
                else:
                    np.savetxt(
                        os.path.join(
                            save_dir,
                            "probability_matrices",
                            f"visited_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat",
                        ),
                        visited,
                        fmt="%d",
                    )
                print("Visited voxels saved")

        print("Visited voxels:", f"{np.count_nonzero(visited)/N*100:.2f} %", flush=True)

        values = []  # indices of possible voxels
        prob_values = []  # probabilities for corresponding indices in values
        ns = [n_x, n_y, n_z]

        def squared_distance(pos1, pos2, boxes):
            # Estimate squared distance between given positions
            # Assumes pos2 is always in the 0th periodic image

            pos1 = list(indices[pos1])
            spacings = [x_spacing, y_spacing, z_spacing]

            # Boxes denotes the periodic image to which pos1 belongs
            pos1[0] += boxes[0] * ns[0]
            pos1[1] += boxes[1] * ns[1]
            pos1[2] += boxes[2] * ns[2]

            dis = ((pos1[0] - pos2[0]) * spacings[0]) ** 2
            dis += ((pos1[1] - pos2[1]) * spacings[1]) ** 2
            dis += ((pos1[2] - pos2[2]) * spacings[2]) ** 2

            return dis

        for i in range(N):
            values.append(prob[i].nonzero()[1])
            prob_values.append(prob[i].data)

        # Map of lengths of travel between in each step of random walk
        final_dist = defaultdict(int)

        def random_walk(n_walks):

            all_trajectories = []

            if self.hopping_distance_analysis:
                trajectory = []
                final_dist = defaultdict(int)

            for walker in range(1, n_walks + 1):

                temp = []

                # Starting voxel is randomly chosen from the set of visited voxels
                start = np.random.choice(visited)
                pos = start
                prev = start
                boxes = [0, 0, 0]

                previ = indices[start]
                temp.append(
                    np.array(previ, dtype=np.float16)
                    * np.array([x_spacing, y_spacing, z_spacing], dtype=np.float16)
                )

                for i in range(n_steps):

                    pos = np.random.choice(values[pos], p=prob_values[pos])
                    position = indices[pos]
                    previous = indices[prev]

                    # Crossing > half of the box dimension implies
                    # moving to the previous periodic image
                    for m in range(3):
                        if position[m] > previous[m]:
                            if position[m] - previous[m] > (ns[m] // 2):
                                boxes[m] -= 1
                        else:
                            if (previous[m] - position[m]) > (ns[m] // 2):
                                boxes[m] += 1

                    if self.hopping_distance_analysis:
                        # Total distance travelled in A in the current step
                        curr = squared_distance(pos, previ, boxes)
                        final_dist[(curr ** (0.5)) * self.nm_to_A] += 1
 
                    previ = list(position)

                    previ[0] += boxes[0] * ns[0]
                    previ[1] += boxes[1] * ns[1]
                    previ[2] += boxes[2] * ns[2]

                    prev = pos

                    temp.append(
                        np.array(previ, dtype=np.float16)
                        * np.array([x_spacing, y_spacing, z_spacing], dtype=np.float16)
                    )

                temp = np.array(temp, dtype=np.float16)
                all_trajectories.append(temp)

                if self.hopping_distance_analysis:

                    # Keep the top traj_n trajectories that were highly spread out for plots
                    dis1 = np.dot(temp[0] - temp[-1], temp[0] - temp[-1])
                    trajectory.append([dis1, temp])
                    trajectory.sort(key=lambda x: x[0], reverse=True)
                    trajectory = trajectory[:traj_n]

            if self.hopping_distance_analysis:
                return all_trajectories, final_dist, trajectory

            else:
                return [all_trajectories]

        print(f"Performing {n_steps} steps by {n_walks} walkers", flush=True)
        if n_walks < np.count_nonzero(visited) * 0.75:
            print(
                f"*** WARNING ***     Significantly more voxels than ghost walkers"
            )
            print(
                f"*** WARNING ***     Consider keeping {(np.count_nonzero(visited)*0.76//50)*50 :.0f} or more walkers",
                flush=True,
            )

        walks_list = [n_walks // (self.n_cpus) for i in range(self.n_cpus)]

        for i in range(n_walks % (self.n_cpus)):
            walks_list[i] += 1

        all_result = p_map(random_walk, walks_list, num_cpus=self.n_cpus)

        print("\nGathering random walks")

        all_trajectories = []
        if self.hopping_distance_analysis:
            final_dist_list = []
            trajectory = []

        for result in all_result:
            all_trajectories += result[0]
            if self.hopping_distance_analysis:
                final_dist_list.append(result[1])
                trajectory += result[2]

        all_trajectories = np.array(all_trajectories, dtype=np.float16) * self.nm_to_A
        print(f"Memory usage for all trajectories: {all_trajectories.nbytes/1e6} MB")

        print(f"Computing MSD ", end = "", flush=True)

        msd_list = [msd_fft(all_trajectories[i]) for i in range(n_walks)]

        if self.hopping_distance_analysis:
            for dicti in final_dist_list:
                for ele in dicti:
                    final_dist[ele] += dicti[ele]

            trajectory.sort(key=lambda x: x[0], reverse=True)
            trajectory = np.array([x[1] for x in trajectory[:traj_n]])

            if equal_spacings:
                spacings = [x_spacing]
                tekmc_file = f"traj_{self.key}_rgx{x_spacing}.svg"
            else:
                spacings = [x_spacing, y_spacing, z_spacing]
                tekmc_file = (
                    f"traj_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.svg"
                )

            if traj_colors is None:
                traj_colors = self.generate_colors(traj_n)

        msd = np.average(msd_list, axis=0)  # [1:]
        time = np.arange(0, self.timestep * len(msd), self.timestep)

        [D, c], pcov = curve_fit(
            self.linearMSD, time[len(msd) // 2 :], msd[len(msd) // 2 :]
        )  # Linear fit to estimate the Diffusion coefficient
        D /= 6
        D *= 10 ** (-7)  # Converting to cm^2/s
        D_SD = np.sqrt(np.diag(pcov)[0])
        D_SD *= 10 ** (-7)  # Converting to cm^2/s

        [a, nu], pcov = curve_fit(
            self.powerLaw, time[len(msd) // 2 :], msd[len(msd) // 2 :], p0=[1e2, 1]
        )  # Finding the exponent of t for the MSD from TEKMC
        nu_SD = np.sqrt(np.diag(pcov)[1])

        print()
        print(f"{'Diffusion (D) =':>25s} {D:.3g} +/- {D_SD:.3g} cm^2/s", flush=True)
        print(f"{'Exponent (nu) =':>25s} {nu:.2g} +/- {nu_SD:.2g} ", flush=True)

        os.makedirs(os.path.join(save_dir, "msd"), exist_ok=True)
        t = 0
        if equal_spacings:
            msd_filename = os.path.join(
                save_dir, "msd", f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat"
            )
        else:
            msd_filename = os.path.join(
                save_dir,
                "msd",
                f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat",
            )

        # Saving MSD data
        header = f" t (ns) MSD (A^2)"
        header += f" D = {D:.4e} {D_SD:.2e} cm^2/s"
        header += f" nu = {nu:.4e} {nu_SD:.2e}"
        np.savetxt(
            msd_filename,
            np.c_[time, msd],
            fmt="%.6f",
            delimiter=4 * " ",
            header=header,
        )

        # Plot repeatedly during runtime for inspection
        self.plot_msd(spacings_list=None, _colors=None)
        self.plot_diffusion(spacings_list=None, _colors=None)

        if self.hopping_distance_analysis:
            # Plot traj_n trajectories (most spread out)
            print("Plotting few longest TEKMC trajectories", flush=True)
            self.plot_trajectory(
                spacings,
                plot_steps=traj_steps,
                n_traj=traj_n,
                _colors=traj_colors,
                trajectories=trajectory,
                tekmc_file=tekmc_file,
                unwrapped=True,
            )

            # Write all displacements between steps of random walks to a file
            os.makedirs(os.path.join(save_dir, "hopping_distances"), exist_ok=True)
            if equal_spacings:
                fname = f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat"
            else:
                fname = f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
            with open(os.path.join(save_dir, "hopping_distances", fname), "w") as f:
                for ele in final_dist:
                    # displacement(A), n_occurences
                    f.write(
                        str(self.format(ele)) + "    " + str(final_dist[ele]) + "\n"
                    )

        # Saving trajectories
        if save_random_walks:
            print(
                f"Saving topology, trajectory of ghost walkers, every {self.stride} timestep",
                flush=True,
            )
            os.makedirs(os.path.join(save_dir, "random_walks"), exist_ok=True)
            if equal_spacings:
                xyz_filename = os.path.join(
                    save_dir,
                    "random_walks",
                    f"random_walks_{self.key}_rgx{x_spacing}.xyz",
                )

                gro_filename = os.path.join(
                    save_dir,
                    "random_walks",
                    f"random_walks_{self.key}_rgx{x_spacing}.gro",
                )

            else:
                xyz_filename = os.path.join(
                    save_dir,
                    "random_walks",
                    f"random_walks_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.xyz",
                )
                gro_filename = os.path.join(
                    save_dir,
                    "random_walks",
                    f"random_walks_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.gro",
                )

            with open(xyz_filename, "w") as xyz_file:
                for frame in range(0, n_steps + 1, self.stride):
                    xyz_file.write(f"{n_walks}\n")
                    xyz_file.write(f"id x y z # {frame * self.timestep:.3f} ns\n")
                    for walker in range(n_walks):
                        xyz_file.write(
                            f"G {' '.join(map(lambda x: str(f'{x: .2f}'), all_trajectories[walker, frame]))}\n"
                        )

            with open(gro_filename, "w") as gro_file:
                gro_file.write("Ghost walkers\n")
                gro_file.write(f"{n_walks}\n")
                frame = 1
                for walker in range(n_walks):
                    gro_file.write(
                        f"{walker+1:5g}{'G  ':>10s}{walker+1:5g}{''.join(map(lambda x: str(f'{x:8.3f}'), all_trajectories[walker, frame]/self.nm_to_A))}\n"
                    )

                gro_file.write(
                    f"{self.traj.unitcell_lengths[1,0]:>10.5f}"
                    + f"{self.traj.unitcell_lengths[1,1]:>10.5f}"
                    + f"{self.traj.unitcell_lengths[1,2]:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"{0:>10.5f}"
                    + f"\n"
                )

        else:
            print("Trajectories of ghost walkers not saved", flush=True)

        # Clean up memory
        del all_trajectories

        end_date = datetime.now()
        print(f"Finished on {end_date.ctime()}", end = " || ")
        print(f"Runtime: {(end_date - start_date)}")

    def plot_msd(
        self, spacings_list=None, _colors=None, y_scale="linear", x_scale="linear"
    ):
        # Plots MSD vs time for given list of spacings
        print("Plotting MSD vs time ", flush=True)

        save_dir = self.save_dir
        md_file = self.md_file

        n = len([x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith("dat")])
        if spacings_list:
            n = len(spacings_list)

        if _colors is None:
            _colors = self.generate_colors(n)
        elif len(_colors) < n:
            print("Insufficient colors passed, Automatic colors generated")
            _colors = self.generate_colors(n)

        # Fetch MSD data of the MD simulation.
        with open(md_file, "r") as f:
            data = np.loadtxt(f, comments="#")
            md_time = data[:, 0]
            md_msd = data[:, 1]

        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.figure(figsize=[4, 4], dpi=72)

        if spacings_list is None:
            msd_files = [
                x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")
            ]
        else:
            msd_files = []
            for spac in spacings_list:
                if type(spac) != list:
                    spac = [spac]
                if len(spac) == 1:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[0], spac[0]
                    msd_files.append(f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat")
                else:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[1], spac[2]
                    msd_files.append(
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    )
        ind = 0
        msd_files = sorted(msd_files)
        for filename in msd_files:
            n = f"{self.key}_TEKMC_MSD.dat".count("_")
            spacing = "_".join([x[3:] for x in filename.split("_")[n + 1 :]])[:-4]
            with open(os.path.join(save_dir, "msd", filename), "r") as f:
                data = np.loadtxt(f, comments="#")
                time = data[:, 0]
                msd = data[:, 1]

            plt.plot(
                time[1 : len(msd)],
                msd[1 : len(msd)],
                label=f"{spacing}",
                color=_colors[ind],
                linewidth=2,
                alpha=0.8,
            )
            ind += 1

        # For MD points, color is black(convention).
        plt.scatter(
            md_time[1:],
            md_msd[1:],
            label="MD",
            color="black",
            s=16.0,
            marker="o",
            alpha=0.8,
            zorder=100,
        )


        plt.yscale(y_scale)
        plt.xscale(x_scale)

        plt.tick_params(axis="both", direction="in", which="major", size=8, labelsize="14")
        plt.tick_params(axis="both", direction="in", which="minor", size=5, labelsize="14")

        plt.xlabel("t (ns)", fontsize=14)
        plt.ylabel(r"MSD $(\AA ^2)$", fontsize=14)
        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=8, loc="best", ncol=3)
        plt.tight_layout()
        plt.savefig(
            os.path.join(save_dir, "msd", f"MSD_vs_t_{self.key}.svg"),
            dpi=72,
            transparent=True,
        )
        plt.close()

    def plot_diffusion(
        self, spacings_list=None, _colors=None, y_scale="linear", x_scale="linear"
    ):
        # Plots Diffusion Coefficient vs time for given list of spacings.
        print("Plotting D vs time", flush=True)

        save_dir = self.save_dir
        md_file = self.md_file

        n = len([x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")])
        if spacings_list:
            n = len(spacings_list)

        if _colors is None:
            _colors = self.generate_colors(n)
        elif len(_colors) < n:
            print("Insufficient colors passed, Automatic colors generated")
            _colors = self.generate_colors(n)

        # Fetch MSD data of the MD simulation.
        with open(md_file, "r") as f:
            data = np.loadtxt(f, comments="#")
            md_time = data[1:, 0]
            md_Ds = data[1:, 1] / 6 / md_time * 1e-7

        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.figure(figsize=[4, 4], dpi=72)

        if spacings_list is None:
            msd_files = [
                x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")
            ]
        else:
            msd_files = []
            for spac in spacings_list:
                if type(spac) != list:
                    spac = [spac]
                if len(spac) == 1:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[0], spac[0]
                    msd_files.append(f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat")
                else:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[1], spac[2]
                    msd_files.append(
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    )

        ind = 0
        msd_files = sorted(msd_files)
        for filename in msd_files:
            n = f"{self.key}_TEKMC_MSD.dat".count("_")
            spacing = "_".join([x[3:] for x in filename.split("_")[n + 1 :]])[:-4]
            with open(os.path.join(save_dir, "msd", filename), "r") as f:
                data = np.loadtxt(f, comments="#")
                time = data[1:, 0]
                Ds = data[1:, 1] / 6 / time * 1e-7

            plt.plot(
                time,
                Ds,
                label=f"{spacing}",
                color=_colors[ind],
                linewidth=2,
                alpha=0.8,
            )
            ind += 1

        # for MD 'black' is default color (convention).
        plt.scatter(md_time, md_Ds, label="MD", color="black", s=16.0, alpha=0.8, zorder=100)

        plt.yscale(y_scale)
        plt.xscale(x_scale)

        plt.ylim([0,Ds[-1]*10])
        plt.xlim([md_time[len(md_time)//5], time[-1]])

        plt.tick_params(axis="both", direction="in", which="major", size=8, labelsize="14")
        plt.tick_params(axis="both", direction="in", which="minor", size=5, labelsize="14")

        plt.xlabel("t (ns)", fontsize=14)
        plt.ylabel(r"D $(cm^2/s)$", fontsize=14)
        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=8, loc="best", ncol=3)
        plt.tight_layout()
        plt.savefig(
            os.path.join(save_dir, "msd", f"D_vs_t_{self.key}.svg"),
            dpi=72,
            transparent=True,
        )
        plt.close()

    def plot_trajectory(
        self,
        spacings,
        plot_steps=1000,
        n_traj=3,
        _colors=None,
        unwrapped=True,
        trajectories=None,
        tekmc_file=None,
    ):
        # Plot of trajectories of few random walks (wrapped or unwrapped) for the given spacings.
        # If trajectories are given, they are plotted and saved to tekmc_file.
        # If trajectories are not provided, a few TEKMC moves are first made.

        save_dir = self.save_dir

        if _colors is None:
            _colors = self.generate_colors(n_traj)

        if type(spacings) != list:
            spacings = [spacings]

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        equal_spacings = len(spacings) == 1

        n_x, n_y, n_z = self.get_ns(x_spacing, y_spacing, z_spacing)

        if trajectories is None:
            print("Plotting few random walks", flush=True)

            if equal_spacings:
                prob = load_npz(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"probSymm_{self.key}_rgx{x_spacing}.npz",
                    )
                )
            else:
                prob = load_npz(
                    os.path.join(
                        save_dir,
                        "probability_matrices",
                        f"probSymm_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.npz",
                    )
                )

            N = n_x * n_y * n_z
            ns = [n_x, n_y, n_z]
            visited = [i for i in range(0, N) if sum(prob[i].toarray()[0]) > 0.99]
            values = []  # indices of possible voxels
            prob_values = []  # probabilities for corresponding indices in values.

            for i in range(N):
                values.append(prob[i].nonzero()[1])
                prob_values.append(prob[i].data)

            indices = {}
            ind = 0
            for i in range(n_x):
                for j in range(n_y):
                    for k in range(n_z):
                        indices[ind] = (i, j, k)
                        ind += 1
            coords = []

            def coordinates_unwrapped(pos, boxes, ns):
                # function to estimate unwrapped coordinates for position pos.
                i, j, k = pos
                i += boxes[0] * ns[0]
                j += boxes[1] * ns[1]
                k += boxes[2] * ns[2]
                return [i * x_spacing, j * y_spacing, k * z_spacing]

            def coordinates_wrapped(pos, *args, **kwargs):
                # function to estimate wrapped coordinates of pos.
                i, j, k = pos
                return [i * x_spacing, j * y_spacing, k * z_spacing]

            def get_coordinates():
                if unwrapped:
                    return coordinates_unwrapped
                return coordinates_wrapped

            coordinates = get_coordinates()

            for i_traj in range(n_traj):
                # starting voxel is randomly selected from visited voxels.
                start = np.random.choice(visited)
                pos = start
                prev = start
                boxes = [0, 0, 0]
                new_traj = [coordinates(indices[pos], boxes, ns)]
                for i in range(plot_steps):
                    pos = np.random.choice(values[pos], p=prob_values[pos])
                    position = indices[pos]
                    previous = indices[prev]
                    for m in range(3):
                        if position[m] > previous[m]:
                            if position[m] - previous[m] > (ns[m] // 2):
                                boxes[m] -= 1
                        else:
                            if (previous[m] - position[m]) > (ns[m] // 2):
                                boxes[m] += 1
                    new_traj.append(coordinates(indices[pos], boxes, ns))
                    prev = pos
                coords.append(new_traj)
        else:
            coords = [x[:plot_steps] for x in trajectories]

        mpl.rcParams.update(mpl.rcParamsDefault)
        fig = plt.figure(figsize=[4, 4], dpi=72)
        ax = fig.add_subplot(111, projection="3d")
        # x,y,z = [0,(n_x * x_spacing)],[0,(n_y * y_spacing)],[0,(n_z * z_spacing)]
        x, y, z = (
            [0, round(self.xhi - self.xlo, 4)],
            [0, round(self.yhi - self.ylo, 4)],
            [0, round(self.zhi - self.zlo, 4)],
        )
        corners = [np.array([a, b, c]) for a in x for b in y for c in z]
        for i in range(8):  # Going over all  corners for a cuboid.
            for j in range(i + 1, 8):
                val1, val2 = corners[i], corners[j]
                if sum(val1 == val2) == 2:  # Only one of the coordinates is different
                    xs, ys, zs = (
                        [val1[0], val2[0]],
                        [val1[1], val2[1]],
                        [val1[2], val2[2]],
                    )
                    ax.plot3D(xs, ys, zs, color="black", alpha=0.8)

        ax.grid(False)
        ax.set_xlabel("x (nm)", fontsize=12)
        ax.set_ylabel("y (nm)", fontsize=12)
        ax.set_zlabel("z (nm)", fontsize=12)

        coords = np.array(coords)

        box_lim = [
            min(np.min(coords), x[0], y[0], z[0]),
            max(np.max(coords), x[1], y[1], z[1]),
        ]

        ax.text(x[1] / 2, -0.5 * y[1], z[0], "a")
        ax.text(1.2 * x[1], 0.25 * y[1], z[0], "b")
        ax.text(1.6 * x[1], y[0], z[1], "c")

        if _colors is None or _colors is True:
            _colors = self.generate_colors(n_traj)

        for i_traj in range(n_traj):
            ax.plot3D(
                coords[i_traj][:, 0],
                coords[i_traj][:, 1],
                coords[i_traj][:, 2],
                color=_colors[i_traj],
                linewidth=1,
                alpha=0.8,
            )

        ax.set_xlim(box_lim)
        ax.set_ylim(box_lim)
        ax.set_zlim(box_lim)
        ax.axis("off")
        ax.view_init(20, -60)
        # ax.set_title("Trajectories")

        plt.tight_layout()
        if tekmc_file:  # file to save the plot when trajectories is not None.
            os.makedirs(os.path.join(save_dir, "trajectories"), exist_ok=True)

            plt.savefig(
                os.path.join(save_dir, "trajectories", tekmc_file),
                dpi=72,
                transparent=True,
            )

        else:
            os.makedirs(os.path.join(save_dir, "trajectories"), exist_ok=True)
            if equal_spacings:
                plt.savefig(
                    os.path.join(
                        save_dir,
                        "trajectories",
                        f"traj_{self.key}_rgx{x_spacing}.svg",
                    ),
                    dpi=72,
                    transparent=True,
                )

            else:
                plt.savefig(
                    os.path.join(
                        save_dir,
                        "trajectories",
                        f"traj_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.svg",
                    ),
                    dpi=72,
                    transparent=True,
                )

        plt.close()

    def normalized_MSE(self, spacings_list=None):
        # Normalized MSE for all spacings in ‘msd_files’ folder w.r.t MSD of MD is estimated.
        print("Estimating normalized MSE ")

        save_dir = self.save_dir
        md_file = self.md_file

        key = self.key

        # we are using normalized MSE as an evaluation metric.
        if spacings_list is None:
            msd_files = [
                x for x in os.listdir(os.path.join(save_dir, "msd")) if x.endswith(".dat")
            ]
        else:
            msd_files = []
            for spac in spacings_list:
                if type(spac) != list:
                    spac = [spac]
                if len(spac) == 1:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[0], spac[0]
                    msd_files.append(f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat")
                else:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[1], spac[2]
                    msd_files.append(
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    )

        filenames = msd_files
        ind = f"{key}_TEKMC_MSD_.dat".count("_")
        spacings = ["_".join((x.split("_")[ind:]))[:-4] for x in filenames]
        spacings = [list(float(x[3:]) for x in spac.split("_")) for spac in spacings]
        spacings.sort()
        fMSE = open(os.path.join(save_dir, "normalized_MSE.dat"), "w")
        fMSE.write("# Normalized Mean squared error:\n\n")
        mse_list = []
        for spac in spacings:
            if len(spac) == 1:
                fname = os.path.join(
                    save_dir, "msd", f"{key}_TEKMC_MSD_rgx{spac[0]}.dat"
                )
            else:
                fname = os.path.join(
                    save_dir,
                    "msd",
                    f"{key}_TEKMC_MSD_rgx{spac[0]}_rgy{spac[1]}_rgz{spac[2]}.dat",
                )
            arr1 = np.loadtxt(fname, skiprows=1) # TEKMC MSD
            fname = md_file
            arr2 = np.loadtxt(fname, skiprows=1) # MD MSD

            dicti1 = {}
            for ele in arr2:
                dicti1[round(ele[0], 2)] = ele[1]

            msd1 = []
            msd2 = []
            for ele in arr1:
                try:
                    value = dicti1[round(ele[0], 2)]
                    msd2.append(value)
                    msd1.append(ele[1])
                except KeyError:
                    # One MSD can be longer than the other
                    continue

            end = len(msd2)
            start = 5*end//10

            msd1 = msd1[start:end]
            msd2 = msd2[start:end]

            minima = min(msd2)
            maxima = max(msd2)
            msd1 = (msd1 - minima) / (maxima - minima)

            minima = min(msd2)
            maxima = max(msd2)
            msd2 = (msd2 - minima) / (maxima - minima)

            mse = self.format((np.square(np.subtract(msd1, msd2))).mean())

            fMSE.write(f"Spacing: {spac}         MSE: {mse}\n")
            print(f"Spacing: {spac}         MSE: {mse}", flush=True)
            mse_list.append(float(mse))

        print(
            f"Minimum MSE obtained for spacing {spacings[np.argmin(mse_list)]}",
            flush=True,
        )
        fMSE.write(f"\n# {spacings[np.argmin(mse_list)]} shows minimum MSE")
        fMSE.close()

    def plot_histogram(self, spacings, n_bins=100):

        if not self.hopping_distance_analysis:
            print("Hopping distances not computed")
            return None

        # Plots histogram of displacements between timesteps for the given spacings.
        print("Plotting histograms of displacement", flush=True)

        save_dir = self.save_dir

        if type(spacings) != list:
            spacings = [spacings]

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        equal_spacings = len(spacings) == 1

        directory = "hopping_distances"
        if equal_spacings:
            files = [f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat"]
        else:
            files = [
                f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
            ]

        os.makedirs(os.path.join(save_dir, "hist_displacement"), exist_ok=True)

        for fname in files:
            file = open(os.path.join(save_dir, directory, fname), "r")
            arr = []
            while True:
                value = list(file.readline().split())
                if value == []:
                    break
                val1, val2 = float(value[0]), float(value[1])
                arr.append([val1, val2])
            arr = np.array(sorted(arr))
            weights = arr[:, 1] / (sum(arr[:, 1]))

            for key in ["displacement"]:

                xlabel = "d (\305)"
                ylabel = "Histogram displacement"
                value = arr[:, 0]
                mpl.rcParams.update(mpl.rcParamsDefault)
                plt.figure(figsize=[3, 3], dpi=72)
                plt.hist(value[value > 0], weights=weights[value > 0], bins=n_bins, histtype="bar", alpha=0.8)
                plt.xlabel(xlabel, fontsize=14)
                plt.ylabel(ylabel, fontsize=14)
                plt.tick_params(axis="both", direction="in", which="major", size=8, labelsize="14")
                plt.tick_params(axis="both", direction="in", which="minor", size=5, labelsize="14")
                plt.minorticks_on()
                plt.tight_layout()
                plt.xlim(left=0)

                if equal_spacings:
                    plt.savefig(
                        os.path.join(
                            save_dir,
                            f"hist_{key}",
                            f"hist_{key}_{self.key}_rgx{x_spacing}.svg",
                        ),
                        dpi=72,
                        transparent=True,
                    )

                else:
                    plt.savefig(
                        os.path.join(
                            save_dir,
                            f"hist_{key}",
                            f"hist_{key}_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.svg",
                        ),
                        dpi=72,
                        transparent=True,
                    )

                plt.close()

    def plot_displacement_probability(self, spacings_list=None, _colors=None):

        if not self.hopping_distance_analysis:
            print("Hopping distances not computed")
            return None

        # Scatter plot of probability of displacements for given list of spacings.
        print("Scatter plot of probability of displacements")

        save_dir = self.save_dir

        if spacings_list is None:
            dist_files = [
                x
                for x in os.listdir(os.path.join(save_dir, "hopping_distances"))
                if x.endswith(".dat")
            ]
        else:
            dist_files = []
            for spac in spacings_list:
                if type(spac) != list:
                    spac = [spac]
                if len(spac) == 1:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[0], spac[0]
                    dist_files.append(f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat")
                else:
                    x_spacing, y_spacing, z_spacing = spac[0], spac[1], spac[2]
                    dist_files.append(
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    )

        dist_files = sorted(dist_files)

        N = len(dist_files)

        if _colors is None:
            _colors = self.generate_colors(N)

        for key in ["displacement"]:
            xlabel = "d (\305)"
            ylabel = "P (d)"

            mpl.rcParams.update(mpl.rcParamsDefault)
            plt.figure(figsize=[4, 4], dpi=72)

            plt.xlabel(xlabel, fontsize=14)
            plt.ylabel(ylabel, fontsize=14)
            plt.minorticks_on()

            for i in range(N):
                filename = dist_files[i]
                n = f"{self.key}_TEKMC_MSD.dat".count("_")
                spacing = "_".join(filename.split("_")[n + 1 :])[:-4]
                spacings = [float(x[3:]) for x in spacing.split("_")]
                if len(spacings) == 1:
                    x_spacing, y_spacing, z_spacing = (
                        spacings[0],
                        spacings[0],
                        spacings[0],
                    )
                elif len(spacings) == 3:
                    x_spacing, y_spacing, z_spacing = (
                        spacings[0],
                        spacings[1],
                        spacings[2],
                    )
                equal_spacings = len(spacings) == 1
                directory = "hopping_distances"
                if equal_spacings:
                    files = [f"{self.key}_TEKMC_MSD_rgx{x_spacing}.dat"]
                else:
                    files = [
                        f"{self.key}_TEKMC_MSD_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.dat"
                    ]
                for fname in files:
                    file = open(os.path.join(save_dir, directory, fname), "r")
                    arr = []
                    while True:
                        value = list(file.readline().split())
                        if value == []:
                            break
                        val1, val2 = float(value[0]), float(value[1])
                        arr.append([val1, val2])
                    arr = np.array(sorted(arr))
                    weights = arr[:, 1] / (sum(arr[:, 1]))
                    value = arr[:, 0]
                    plt.scatter(
                        value[1:],
                        weights[1:],
                        color=_colors[i],
                        alpha=0.6,
                    )

                    plt.plot(
                        value[1:],
                        weights[1:],
                        color=_colors[i],
                        label="_".join(np.array(spacings, dtype=str)),
                        alpha=0.5,
                    )

            plt.xlim(left=0)
            plt.tick_params(axis="both", direction="in", which="major", size=8, labelsize="14")
            plt.tick_params(axis="both", direction="in", which="minor", size=5, labelsize="14")

            plt.legend(fontsize=10)
            plt.tight_layout()

            plt.savefig(
                os.path.join(save_dir, "hopping_distances", f"hop_prob_{key}.svg"),
                dpi=72,
                transparent=True,
            )
            plt.close()

    def pathways(self, spacings, Threshold=1e-6):
        # Estimating all pathways of 3D grid of voxels assuming each voxel corresponds to a node in a graph.

        save_dir = self.save_dir

        class Graph:
            def __init__(self, N, prob, Threshold):
                self.N = N
                self.prob = prob
                # 2 nodes are connected only if transition probability (in both directions) between them is > Threshold.
                self.Threshold = Threshold

            def DFS_iterative(self, start):
                stack = [start]
                self.visited[start] = True
                component = []
                while stack:
                    i = stack.pop()
                    component.append(i)
                    for j in list(prob[i].nonzero()[1]):
                        if (
                            self.visited[j] == False
                            and self.prob[i, j] > self.Threshold
                            and self.prob[j, i] > self.Threshold
                        ):
                            self.visited[j] = True
                            stack.append(j)
                return component

            def pathways(self):
                self.visited = [False for i in range(self.N)]
                components = []
                for i in range(self.N):
                    if self.visited[i] == False:
                        component = self.DFS_iterative(i)
                        components.append(component)
                return components

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        equal_spacings = len(spacings) == 1

        if equal_spacings:
            prob = load_npz(
                os.path.join(
                    save_dir,
                    "probability_matrices",
                    f"probSymm_{self.key}_rgx{x_spacing}.npz",
                )
            )
        else:
            prob = load_npz(
                os.path.join(
                    save_dir,
                    "probability_matrices",
                    f"probSymm_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.npz",
                )
            )

        n_x, n_y, n_z = self.get_ns(x_spacing, y_spacing, z_spacing)
        N = n_x * n_y * n_z

        # a node will be part of the connected component if transition prob in both directions > threshold
        graph = Graph(N, prob, Threshold)
        components = graph.pathways()
        components = sorted(components, key=lambda x: len(x), reverse=True)

        return components

    def plot_pathways(
        self,
        spacings,
        n_components,
        _colors=None,
        Threshold=1e-6,
    ):
        # Plot of few longest pathways of the 3D grid of voxels.
        print(f"Plotting {n_components} longest pathways ", flush=True)

        save_dir = self.save_dir

        if _colors is None:
            _colors = self.generate_colors(n_components)

        if type(spacings) != list:
            spacings = [spacings]

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        equal_spacings = len(spacings) == 1

        n_x, n_y, n_z = self.get_ns(x_spacing, y_spacing, z_spacing)
        mpl.rcParams.update(mpl.rcParamsDefault)
        fig = plt.figure(figsize=[4, 4], dpi=72)
        ax = fig.add_subplot(111, projection="3d")

        ax.grid(False)
        origin = [0, 0, 0]
        data = np.empty((n_x, n_y, n_z), dtype=object)
        colors = np.empty(data.shape, dtype=object)
        indices = {}
        ind = 0
        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    indices[ind] = (i, j, k)
                    ind += 1
        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    data[i][j][k] = False
                    colors[i][j][k] = "#00000000"

        # 2 nodes are connected only if the probability of transition (in both directions) > Threshold
        components = self.pathways(spacings, Threshold=Threshold)

        vol_voxel = x_spacing * y_spacing * z_spacing
        components_vol = [f"{len(x) * vol_voxel:.4f}" for x in components]
        components_vol = sorted(components_vol, reverse=True)

        n_components = min(n_components, len(components))

        components_vol = components_vol[:n_components]
        print(
            f"Volume of {n_components} longest pathways (nm^3):",
            *components_vol,
            end="\n",
        )

        for z in range(n_components):
            ele = components[z]
            for ind in ele:
                i, j, k = indices[ind]
                colors[i][j][k] = _colors[z]
                data[i][j][k] = True

        x, y, z = np.indices((n_x + 1, n_y + 1, n_z + 1))
        ax.voxels(
            x * x_spacing,
            y * y_spacing,
            z * z_spacing,
            data,
            facecolors=colors,
            edgecolor="#00000011",
            alpha=0.8,
        )
        ax.set_xlabel("x (nm)", fontsize=11)
        ax.set_ylabel("y (nm)", fontsize=11)
        ax.set_zlabel("z (nm)", fontsize=11)

        os.makedirs(os.path.join(save_dir, "pathways"), exist_ok=True)
        plt.tight_layout()
        if equal_spacings:
            plt.savefig(
                os.path.join(
                    save_dir,
                    "pathways",
                    f"connected_{self.key}_rgx{x_spacing}.svg",
                ),
                dpi=72,
                transparent=True,
            )

        else:
            plt.savefig(
                os.path.join(
                    save_dir,
                    "pathways",
                    f"connected_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.svg",
                ),
                dpi=72,
                transparent=True,
            )
        plt.close()

    def plot_voxels_probability(self, spacings, color_map="terrain"):
        # Plot of all voxels with color gradient based on probability of transition to itself.
        print(
            "Plotting probability of residing in each voxel",
            flush=True,
        )

        save_dir = self.save_dir

        if type(spacings) != list:
            spacings = [spacings]

        if len(spacings) == 1:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[0], spacings[0]
        elif len(spacings) == 3:
            x_spacing, y_spacing, z_spacing = spacings[0], spacings[1], spacings[2]

        equal_spacings = len(spacings) == 1

        if equal_spacings:
            prob = load_npz(
                os.path.join(
                    save_dir,
                    "probability_matrices",
                    f"probSymm_{self.key}_rgx{x_spacing}.npz",
                )
            )
        else:
            prob = load_npz(
                os.path.join(
                    save_dir,
                    "probability_matrices",
                    f"probSymm_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.npz",
                )
            )

        n_x, n_y, n_z = self.get_ns(x_spacing, y_spacing, z_spacing)
        N = n_x * n_y * n_z
        voxels = {}
        ind = 0
        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    voxels[(i, j, k)] = ind
                    ind += 1

        mpl.rcParams.update(mpl.rcParamsDefault)
        fig = plt.figure(figsize=[6, 4], dpi=72)
        ax = fig.add_subplot(111, projection="3d")

        ax.grid(False)
        origin = [0, 0, 0]
        data = np.empty((n_x, n_y, n_z), dtype=object)
        colors = np.empty(data.shape, dtype=object)
        cmap = plt.get_cmap(color_map)
        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    ind = voxels[(i, j, k)]
                    if prob[ind, ind] != 0:
                        data[i][j][k] = True
                        colors[i][j][k] = cmap(prob[ind, ind])

        x, y, z = np.indices((n_x + 1, n_y + 1, n_z + 1))
        ax.voxels(
            x * x_spacing,
            y * y_spacing,
            z * z_spacing,
            data,
            facecolors=colors,
            edgecolor="#00000011",
            alpha=0.6,
        )
        ax.set_xlabel("x (nm)", fontsize=12)
        ax.set_ylabel("y (nm)", fontsize=12)
        ax.set_zlabel("z (nm)", fontsize=12)
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        ax.tick_params(axis="z", labelsize=12)

        diagonal_values = [prob[i, i] for i in range(N)]
        norm = Normalize(vmin=min(diagonal_values), vmax=max(diagonal_values))
        m = cm.ScalarMappable(cmap=cmap, norm=norm)
        plt.colorbar(m, pad=0.25, shrink=0.6)

        os.makedirs(os.path.join(save_dir, "voxels_probability"), exist_ok=True)

        plt.tight_layout()

        if equal_spacings:
            plt.savefig(
                os.path.join(
                    save_dir,
                    "voxels_probability",
                    f"voxels_prob_{self.key}_rgx{x_spacing}.svg",
                ),
                dpi=72,
                transparent=True,
            )
        else:
            plt.savefig(
                os.path.join(
                    save_dir,
                    "voxels_probability",
                    f"voxels_prob_{self.key}_rgx{x_spacing}_rgy{y_spacing}_rgz{z_spacing}.svg",
                ),
                dpi=72,
                transparent=True,
            )

        plt.close()


# ---------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":

    # Required Parameters

    topology_file = None
    trajectory_file = None
    atom_name = None
    timestep = None
    spacings_list = None
    n_walks = None
    n_components = None

    # Optional parameters

    stride = 1
    md_file = None
    n_cpus = None
    cmap = "terrain"
    save_dir = None
    n_steps = None
    save_random_walks = False
    hopping_distance_analysis = False

    # To modify other parameters, modify parameters passed in the methods below.
    # Look up https://github.com/Thanush1/TEKMC.git for more details on the parameters.

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
            "Also modify other parameters to your specific case.",
            "Look up https://github.com/Thanush1/TEKMC.git for details about all the methods.",
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
        md_file=md_file,
        n_cpus=n_cpus,
        cmap=cmap,
        save_dir=save_dir,
        threshold=0.0,
        symmetrization_type="average",
    )

    # Various methods of TEKMC object.
    # Modify the parameters as required. Look up https://github.com/Thanush1/TEKMC.git for more details.

    spacings_list = sorted(
        spacings_list, key=lambda x: x[0] if type(x) == list else x, reverse=True
    )

    for spacings in spacings_list:
        tekmc.TEKMC(
            spacings,
            n_walks,
            n_steps=n_steps,
            save_prob=True,
            create_prob=False,
            traj_steps=1000,
            traj_n=3,
            traj_colors=None,
            save_random_walks=save_random_walks,
            hopping_distance_analysis=hopping_distance_analysis,
        )

        tekmc.plot_trajectory(
            spacings,
            plot_steps=1000,
            n_traj=3,
            _colors=None,
            unwrapped=True,
            trajectories=None,
            tekmc_file=None,
        )

        tekmc.plot_histogram(spacings, n_bins=100)

        tekmc.plot_pathways(
            spacings,
            n_components,
            _colors=None,
            Threshold=1e-6,
        )
        tekmc.plot_voxels_probability(spacings, color_map="terrain")

    tekmc.compute_diffusion(spacings_list=None)
    tekmc.plot_diffusion(spacings_list=None, _colors=None)
    tekmc.plot_msd(spacings_list=None, _colors=None)
    tekmc.plot_displacement_probability(spacings_list=None, _colors=None)
    tekmc.normalized_MSE(spacings_list=None)
