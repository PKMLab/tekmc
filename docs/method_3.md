## tekmc.plot_Hist(spacings, n_bins=100, save_dir=save_dir)

### Plots histogram of displacements and velocities between timesteps for the given spacings.

**Parameters**
<pre>
spacings	 : list of floats
                   The list can contain 3 floats corresponding to grid size along x, y and z direction (in nm). 
                   If list has one float, grid size is assumed equal along all directions. Ex: [0.23,0.24,0.24], [0.24], etc.
</pre>
**Optional parameters**
<pre>
n_bins		 : int, default: 100
		   Bins of the Histogram.
		   
save_dir	 : str, default: None
  		   Directory where plot will be saved. If None, save_dir = ‘visualizations’.
</pre>
**Outputs**
<pre>
Histogram plots of displacement and velocity between timesteps is saved in save_dir directory.
</pre>
