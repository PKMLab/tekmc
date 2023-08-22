## tekmc.plot_all_voxels(spacings, save_dir=save_dir, color_map="terrain")

### Plot of all voxels with color gradient based on probability of transition to the same voxel for given spacings.

**Parameters**
<pre>
spacings	  : list of floats
                    The list can contain 3 floats corresponding to grid size along x, y and z direction (in nm). 
                    If the list has one float, grid size is assumed equal along all directions. Ex: [0.23,0.24,0.24], [0.24], etc.
</pre>
**Optional parameters**
<pre>
save_dir	  : str, default: None
  		    Directory where plot will be saved. If None, save_dir = ‘visualizations’.

color_map	  : str, default: ‘terrain’
		    Color map used in the plot.
</pre>
**Outputs**
<pre>
Plot is saved in save_dir directory.
</pre>
