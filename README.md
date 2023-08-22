# TEKMC - Trajectory Extending Kinetic Monte Carlo

## Authors

Subhadeep Dasgupta \
Arun K. S \
Prabal K. Maiti \
*Affiliation - Department of Physics, Indian Institute of Science, Bangalore 560012, India*

## [Installation Instructions](docs/installation_instructions.md)

In tools folder, modify [run_tekmc.py](tools/run_tekmc.py) for your case. \
Object of TEKMC class requires the following parameters.

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
**TEKMC object has the following methods**

#### [tekmc.TEKMC(spacings, n_walks, n_steps = None, save_prob = True, create_prob = False, traj_steps = 1000, traj_n = 3, traj_colors = traj_colors, save_dir = None)](docs/method_1.md)
