## tekmc.plot_trajectory( spacings, n_steps, n_traj = 4, save_dir = None, unwrapped = True, \_colors = None, trajectories = None, tekmc_file = None)

### Plots the trajectories of few random walks of TEKMC for the given spacings. 

**Parameters**
<pre>
spacings	   : list of floats
		     The list can contain 3 floats corresponding to grid size along x, y and z direction (in nm). 
		     If the list has one float, grid size is assumed equal along all directions. Ex: [0.23,0.24,0.24], [0.24], etc.
		    
n_steps		   : int
		     Number of steps in each random walk.
		     
n_traj		   : int
		     Number of random walks.
</pre>

**Optional parameters**
<pre>
save_dir	   : str, default: None
  		     Directory where plot will be saved. If None, save_dir = ‘visualizations’.
		     
unwrapped	   : boolean, default: True
		     Whether to plot the unwrapped or wrapped trajectories. 
		     
_colors		   : list of colors, default: None
		     List of colors for the plot. If None, colors are generated using the cmap of tekmc object.
		     
trajectories       : list of np.ndarray, default: None
		     list of trajectories of random walks. If None, random walks will be generated from the probability matrix corresponding to spacings. 
	 	     This is used by TEKMC method to save TEKMC trajectories.
		     
tekmc_file	   : str, default: None
		     Plot will be saved in this file in directory save_dir if trajectories is provided. 
		     This is used by TEKMC method to save TEKMC trajectories.   
</pre>
**Outputs**
<pre>
Plot of few random walks of TEKMC corresponding to the given spacings is saved in save_dir directory.
</pre>
