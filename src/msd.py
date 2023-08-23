#!/bin/env python3
# ---------------------------------------------------------------------------- #
# Computing the mean squared displacement of a species
# obtained from molecular dynamics simulation, and then
# finding its exponent and diffusion constant
# ---------------------------------------------------------------------------- #
# Subhadeep Dasgupta, Arun kumar S., Prabal K. Maiti
# Dept. of Physics, IISc Bangalore, India
# ---------------------------------------------------------------------------- #
# Version | Comments
# ---------------------------------------------------------------------------- #
# v0.7    | Option to only perform essential calculations       | Sep 30, 2022
#         | Timestep is read in ns and converted to ps          |
#         | Fixed RuntimeWarning whle computing nu              |
# v0.6    | Support for multiple trajectory files               | May 19, 2022
# v0.5    | Option to analyze a portion of trajectory           | Mar  8, 2022
# v0.4    | Package management and error handling               | Feb 26, 2022
# v0.3    | Time and space optimization at sake of our sanity   | Feb 24, 2022
# v0.2    | Parallel implementation implementation              | Feb 23, 2022
# v0.1    | Basic functionality with parallel implementation    | Feb 23, 2022
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Imports
# ---------------------------------------------------------------------------- #

import importlib
import os
import sys

required_packages = ["mdtraj", "scipy.optimize", "joblib", "matplotlib", "numpy"]

missing_packages = []
for package in required_packages:
    try:
        importlib.import_module(package)
    except ModuleNotFoundError:
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

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from joblib import Parallel, delayed, cpu_count

# ---------------------------------------------------------------------------- #
# Module
# ---------------------------------------------------------------------------- #


def section(heading: str = "-", sep: str = "-"):
    print("".center(60, "-"))
    print(heading.center(60, " "))
    print("".center(60, "-"), flush=True)


def readme():
    section("Usage")
    print(
        "Expects: Topology, trajectory as supported by \033[1mmdtraj\033[0m",
        "         Atom name for the species of interest",
        "         Time interval between trajectory frames",
        "  Units: trajectory - A",
        "         timestep   - ns",
        "         MSD        - A^2/ns",
        "         D          - cm^2/s",
        "   Note: Multiple trajectory files are space separated",
        "         For e.g. <path_to_file_1> <path_to_file_2>",
        sep="\n",
        flush=True,
    )


def msd_fft(trajectory: np.ndarray) -> np.ndarray:
    def autocorr_fft(x: np.ndarray) -> np.ndarray:
        N = len(x)
        F = np.fft.fft(x, n=2 * N)
        PSD = F * F.conjugate()
        res = np.fft.ifft(PSD)
        n = N * np.ones(N, dtype=np.uint32) - np.arange(0, N, dtype=np.uint32)
        return res[:N].real / n

    frames = len(trajectory)
    D = np.square(trajectory, dtype=np.float32).sum(axis=1)
    D = np.array(np.append(D, 0), dtype=np.float32)
    S2 = sum([autocorr_fft(trajectory[:, i]) for i in range(trajectory.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(frames, dtype=np.float32)

    for m in range(frames):
        Q = Q - D[m - 1] - D[frames - m]
        S1[m] = Q / (frames - m)

    msd = np.array(S1[: frames // 2] - 2 * S2[: frames // 2], dtype=np.float32)
    return msd


def mean_squared_displacement(
    topology: str,
    trajectory: str,
    atom_name: str,
    timestep: str,
    start: float = 0,
    end: float = -1,
    stride: int = 1,
    save_dir=None,
    n_cpus: int = cpu_count() // 2,
    data_type=np.float32,
    plot: bool = True,
    system_info: bool = True,
    analyze: bool = True,
) -> None:

    # -------------------------------------- #
    # Guard clauses
    # -------------------------------------- #

    if not os.path.isfile(topology):
        readme()
        raise FileNotFoundError("Topology not found.")

    for filename in trajectory:
        if not os.path.isfile(filename):
            readme()
            raise FileNotFoundError(f"Trajectory not found: {filename}")

    if data_type(timestep) < 0:
        readme()
        raise ValueError(f"Invalid timestep: {timestep}")

    if not isinstance(atom_name, str):
        readme()
        raise ValueError(f"Invalid atom name: {atom_name}")

    if not isinstance(stride, int):
        readme()
        raise ValueError(f"Invalid stride: {stride}")

    if save_dir is None:
        save_dir = os.getcwd()

    if n_cpus < 1:
        readme()
        raise ValueError(f"Invalid number of CPUs: {n_cpus}")

    # -------------------------------------- #
    # Initialize
    # -------------------------------------- #

    nm_to_A = 1e1
    A_to_cm = 1e-8
    ps_to_ns = 1e-3
    ns_to_s = 1e-9

    key = (topology.split(os.sep)[-1])[:-4]
    timestep = stride * data_type(timestep) / ps_to_ns  # converted to ps

    if save_dir:
        data_file_name = os.path.join(save_dir, f"{key}_{atom_name}_msd.dat")
        plot_file_name = os.path.join(save_dir, f"{key}_{atom_name}_msd.pdf")

    # -------------------------------------- #
    # Functions
    # -------------------------------------- #

    def done():
        print("done", flush=True)

    def power_law(x: float, a: float, nu: float) -> float:
        return a * x**nu

    def linear_msd(x: float, D: float, c: float) -> float:
        return D * x + c

    def group_msd(group: np.ndarray, start: int = 0, end: int = -1) -> np.ndarray:
        msd_per_group = None
        for index in group:
            msd_per_atom = (
                msd_fft(traj.xyz[start:end, index] - center_of_mass[start:end])
                / n_track_atoms
            )
            if msd_per_group is None:
                msd_per_group = msd_per_atom
            else:
                msd_per_group += msd_per_atom
        return msd_per_group

    # -------------------------------------- #
    # Load system
    # -------------------------------------- #

    if system_info:
        print(f"Loading trajectory in topology... ", flush=True)
        print(f"Topology:\n{topology.split(os.sep)[-1]}")
        print(f"Trajectory path: {os.path.dirname(trajectory[0])}", flush=True)
        for filename in trajectory:
            print(f"{filename.split(os.sep)[-1]}", flush=True)

    traj = md.load(trajectory, top=topology, stride=stride)
    top = traj.topology

    atom_name_set = frozenset(atom.name for atom in top.atoms)

    if atom_name not in atom_name_set:
        print("The available atoms are:")
        print(f"\n".join(map(str, atom_name_set)))
        raise ValueError(f"Topology does not contain {atom_name}.")

    n_atoms = traj.n_atoms
    n_frames = traj.n_frames
    total_time = timestep * n_frames
    dimensions = traj.xyz[0, 1].shape[0]

    if system_info:
        section("System Info")
        print(
            f"     Atoms = {n_atoms}",
            f"    Frames = {n_frames}",
            f"  timestep = {timestep:.4g} ns",
            f"Total time = {total_time*ps_to_ns:.4g} ns",
            sep="\n",
            flush=True,
        )

    # -------------------------------------- #
    # Atom selection
    # -------------------------------------- #

    track_list = [atom.index for atom in top.atoms_by_name(atom_name)]
    n_track_atoms = len(track_list)

    if n_track_atoms < 72:
        n_cpus = 1

    group_list = np.array_split(track_list, n_cpus)

    """
    Create a binary row matrix corresponding to indices of
    tracked atoms. Its product with the trajectory provides
    the center of mass coordinates for each frame.
    """
    survival_list = np.zeros(n_atoms, dtype=bool)
    survival_list[track_list] = True

    center_of_mass = np.matmul(survival_list, traj.xyz) / n_track_atoms

    # -------------------------------------- #
    # MSD for half the trajectory time
    # -------------------------------------- #

    start_frame = int(start / timestep)
    end_frame = int(end / timestep) if end > 0 else n_frames

    if system_info:
        print(
            f"  Analysis = {start*ps_to_ns:.2f} -> {end_frame*timestep*ps_to_ns:.2f} ns",
            flush=True,
        )

        print(
            f"Tracking {n_track_atoms} compounds",
            f"of {atom_name} with {n_cpus} cores",
            f"in {len(group_list)} groups... ",
            end="",
            flush=True,
        )

    msd_group = Parallel(n_jobs=n_cpus)(
        delayed(group_msd)(group, start=start_frame, end=end_frame)
        for group in group_list
    )

    msd_avg = np.sum(msd_group, axis=0) * nm_to_A**2
    duration = (
        timestep * np.arange(start_frame, (start_frame + end_frame) // 2, 1) * ps_to_ns
    )
    if system_info:
        done()

    header = f" t (ns) MSD (A^2)"

    if analyze:
        # -------------------------------------- #
        # Diffusion coefficient and MSD exponent
        # -------------------------------------- #

        [D, c], pcov = curve_fit(linear_msd, duration, msd_avg)

        D *= A_to_cm**2 / ns_to_s / (2 * dimensions)
        D_SD = np.sqrt(np.diag(pcov)[0])
        D_SD *= A_to_cm**2 / ns_to_s / (2 * dimensions)

        [a, nu], pcov = curve_fit(
            power_law, duration[2:] - duration[1], msd_avg[2:], p0=[1e2, 1]
        )
        fit_function = [power_law(t - duration[1], a, nu) for t in duration[2:]]
        nu_SD = np.sqrt(np.diag(pcov)[1])

        # -------------------------------------- #
        # MSD output
        # -------------------------------------- #

        section("MSD")
        print(f"{'Diffusion (D) =':>25s} {D:.3g} +/- {D_SD:.3g} cm^2/s", flush=True)
        print(f"{'Exponent (nu) =':>25s} {nu:.2g} +/- {nu_SD:.2g} ", flush=True)

        header += f" D = {D:.4e} {D_SD:.2e} cm^2/s"
        header += f" nu = {nu:.4e} {nu_SD:.2e}"

    else:
        if system_info:
            print("MSD not analysed")

    if save_dir:
        np.savetxt(
            data_file_name,
            np.c_[duration, msd_avg],
            fmt="%.6f",
            delimiter=4 * " ",
            header=header,
        )

    if system_info and save_dir:
        print(f"MSD file: {data_file_name}", flush=True)

    # -------------------------------------- #
    # Plotting
    # -------------------------------------- #

    if plot and save_dir:
        from matplotlib.ticker import (
            MultipleLocator,
            FormatStrFormatter,
            AutoMinorLocator,
        )

        if system_info:
            print(f"Plotting... ", end="", flush=True)
        fig = plt.figure(figsize=(4, 4), dpi=72)
        ax = fig.gca()

        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.tick_params(axis="both", which="major", direction = "in", length = 8)
        ax.tick_params(axis="both", which="minor", direction = "in", length = 5)

        # for X axis
        min_x_limit = np.ceil(min(duration) / 1) * 1
        max_x_limit = np.ceil(max(duration) / 1) * 1
        majorLocator = MultipleLocator((max_x_limit - min_x_limit) / 5)
        majorFormatter = FormatStrFormatter("%.1f")
        minorLocator = MultipleLocator(max_x_limit / 20)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        ax.set_xlabel("t (ns)", fontsize=20)

        # for Y axis
        max_y_limit = np.ceil(max(msd_avg) / 1000) * 1000
        majorLocator = MultipleLocator(max_y_limit / 5)
        majorFormatter = FormatStrFormatter("%.0f")
        minorLocator = MultipleLocator(max_y_limit / 20)

        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.set_ylabel(r"MSD $(\AA ^2)$", fontsize=14)

        ax.plot(duration, msd_avg, linewidth=2.5, color="black")

        ax.plot(duration[2:], fit_function, linewidth=2.5, linestyle="--", color="red")
        if analyze:
            ax.text(
                0.3 * max(duration),
                0.1 * max(msd_avg),
                f"{D = :.2e} $cm^2s^{{-1}}$\n" + rf"$\nu = {nu:.2f}$",
                fontsize=9,
            )

        plt.tight_layout()
        plt.savefig(plot_file_name)

        if system_info:
            done()

    if analyze and system_info:
        if nu - nu_SD < 0.85:
            print(
                f'\n{"*** CAUTION ***":>25s}',
                f"The obtained exponent of MSD is small,",
                f"and the system is probably sub diffusive.",
                sep="\n",
                flush=True,
            )

    if analyze:
        return msd_avg, D, D_SD, nu, nu_SD
    else:
        return msd_avg


# ---------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":
    if (len(sys.argv)) < 5:
        print(
            "  Usage:",
            "<topology>",
            "<trajectory>",
            "<atom name>",
            "<time interval>",
            sep=4 * " ",
            flush=True,
        )
        readme()

    topology_file = sys.argv[1]
    trajectory_file = sys.argv[2:-2]
    atom_name = sys.argv[-2]
    timestep = sys.argv[-1]

    mean_squared_displacement(
        topology_file,
        trajectory_file,
        atom_name,
        timestep,
        start=0,
        end=-1,
        data_type=np.float32,
    )
