# TEKMC - Trajectory-Extending Kinetic Monte Carlo

Please read the docs for [Installation Instructions](docs/installation_instructions.md) and description of other methods and variables.

# About

The kinetic Monte Carlo algorithm is used to sample a short simulation trajectory of small particles and extend it to longer timescales with low computation cost.
The simulation box is sub-divided into small voxels having equal grid sizes.
The input trajectory is analysed to construct a transition probability matrix, representative of the system and all its interactions.
A particle is then assumed to perform Markovian random walks inside the box, with its displacements biased with the probability matrix.
The timestep between two successive random walks is assumed to be constant.
This matrix depends on the chosen grid sizes, which affects the dynamics of the random walkers.
The grid size is tuned until the averaged mean squared displacement of the random walkers correspond to that of the input trajectory.
The advantage of Trajectory Extending Kinetic Monte Carlo technique lies in utilisation of the final probability matrix, to perform significantly long simulations that represent the true dynamics of the input system.
This enables studying long-time dynamics of slowly diffusing systems that would otherwise take significant computation time.

# Article
Details regarding the main work: [Journal](https://pubs.acs.org/doi/10.1021/acs.jpcb.3c05661), [arXiv](https://doi.org/10.48550/arXiv.2311.02878) \
Please cite us if you have used any part of this code.

# Usage
To run, modify [run.py](tools/run.py) for your simulation trajectories and run using python.
Object of TEKMC class requires the following inputs.

**Parameters**
<pre>
trajectory_file                       :  str 
					 name of the trajectory file (.dcd, etc)
topology_file	                      :  str 
					 name of the topology file (.gro, etc)
timestep	                      :  float
                                         timestep between MD frames in ns
stride		                      :  int
                                         stride to add between the frames of the trajectory
atom_name 	                      :  str
                                         atom to track during TEKMC
</pre>

**Optional parameters**
<pre>
md_filename	                      :  str, default: None
                                         MSD file of MD trajectories (should contain time and msd in ns and A^2)
		       	                 This file will be created if not provided
symmetrization_type                   :  {‘min’, ‘max’, ‘average’}, default: ‘average’
                                         type of artificial symmetrization imposed during TEKMC
threshold	                      :  float, default: 0.0
			                 Only entries with prob > threshold will be retained in the probability matrix
n_cpus		                      :  int , default: number of cores in the system
			                 number of CPU cores to utilize during TEKMC run
cmap		                      :  str, default: ‘terrain’
			                 Color map used to generate colors whenever required
</pre>


## References:
[1] Dasgupta, S., KS, A., Ayappa, K. G., & Maiti, P. K. (2023). Trajectory-Extending Kinetic Monte Carlo Simulations to Evaluate Pure and Gas Mixture Diffusivities through a Dense Polymeric Membrane. The Journal of Physical Chemistry B. \
[2] Neyertz, S., & Brown, D. (2010). A trajectory-extending kinetic Monte Carlo (TEKMC) method for estimating penetrant diffusion coefficients in molecular dynamics simulations of glassy polymers. Macromolecules, 43(21), 9210-9214.
