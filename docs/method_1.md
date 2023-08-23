## tekmc.TEKMC(spacings, n_walks, n_steps = None, save_prob = True, create_prob = False, traj_steps = 1000, traj_n = 3, traj_colors = traj_colors, save_dir = None)

### Performs TEKMC with the given spacings and creates MSD file, plot of few longest random walks, etc.

**Parameters**
<pre>
spacings	: list of floats
                  The list can contain 3 floats corresponding to grid size along x, y and z direction (in nm). 
		  If list has one float, grid size is assumed equal along all directions. Ex: [0.23,0.24,0.24], [0.24], etc.

n_walks		: int
		  Number of random walks during TEKMC.
</pre>

**Optional parameters**
<pre>
n_steps		: int
                  Number of steps in each random walk. Default n_steps is 3 times the number of steps of MD trajectory.

save_prob	: boolean, default: 'True'
                  Whether the probability matrix is to be saved or not. Probability matrices are saved to ‘saved’ folder.

create_prob	: boolean, default: 'False'
                  Whether to create the probability matrix or load saved matrix.
		  If True, probability matrix is created irrespective of it is present or absent in ‘saved’ folder.

traj_steps	: int, default: 1000
		  Number of steps of random walk (of TEKMC) plotted.

traj_n		: int, default: 3
                  Number of random walks plotted from TEKMC. Longest traj_n trajectories will be plotted.

traj_colors	: list of colors, default: 'None'
                  Colors for the traj_n TEKMC trajectories. If None, colors are generated using the cmap of tekmc object.

save_dir	: str, default: None
                  Directory where plot of TEKMC trajectories will be saved. If None, save_dir = ‘visualizations’.

hopping_distance_analysis : boolean, default: 'False'
			    Analysis of hopping distances of individual random walks.

</pre>

**Outputs**
<pre>
Probability matrix corresponding to given spacings is saved in 'saved' folder. 
MSD file (containing t(ns) and MSD(Å²)) is saved in 'All_MSDs' folder.
Displacement_file (containing displacement between timesteps and number of occurence of the displacement) is saved in 'Distances' folder.
Plot of few (traj_n) longest trajectories of TEKMC is saved in save_dir.
</pre>
