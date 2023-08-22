## tekmc.plot_Diffusion_coefficient(spacings_list=None, \_colors=None, save_dir=None)

### Plots Diffusion Coefficient vs time for given list of spacings.

**Optional parameters**
<pre>
spacings_list         : list of lists, default: None
                        List of spacings. Ex: [[0.23],[0.22,0.23,0.24]], [[0.22],[0.13],[0.5]]. 
                        If None, plot will be generated for all spacings in ‘All_MSDs’ folder.
                        
_colors	              : list of colors, default: None
                        List of colors for the plot. If None, colors are generated using the cmap of tekmc object.
                        
save_dir              : str, default: None
                        Directory where plot will be saved. If None, save_dir = ‘visualizations’.
</pre>
**Outputs**
<pre>
D vs t plot is saved as 'D_vs_t.pdf' in save_dir directory.
</pre>
