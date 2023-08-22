## tekmc.plot_connected_components( spacings, n_components, \_colors = None, save_dir = None, Threshold = 1e-6)

### Plot of few longest connected components of the 3D grid of voxels.

**Parameters**
<pre>
spacings	 : list of floats
		   The list can contain 3 floats corresponding to grid size along x, y and z direction (in nm). 
		   If the list has one float, grid size is assumed equal along all directions. Ex: [0.23,0.24,0.24], [0.24], etc.
		   
n_components 	 : int 
		   number of connected components to plot.  
</pre>
**Optional parameters**
<pre>
_colors		 : list of colors, default: None
		   colors for the connected components. 
		   If None, colors are generated using the cmap of tekmc object.
		   
Threshold	 : float, default:1e-6
		   2 nodes will be in the same connected component only if probability of transition between them is greater than Threshold.
		   
save_dir	 : str, default: None
  		   Directory where plot will be saved. If None, save_dir = ‘visualizations’.
</pre>
**Outputs**
<pre>
Plot of few longest connected components is saved in save_dir directory.
</pre>

