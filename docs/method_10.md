## tekmc.plot_prob_dist_velocity(spacings_list=None, save_dir=save_dir, \_colors=None)

### Scatter plot of probability of displacements and velocities for given list of spacings.

**Optional parameters**
<pre>
spacings_list           : list of lists, default: None
                          List of spacings. Ex: [[0.23],[0.22,0.23,0.24]], [[0.22],[0.13],[0.5]]. 
                          If None, plot will be generated for all spacings in ‘All_MSDs’ folder.
                          
_colors                 : list of colors, default: None
                          List of colors for the plot. If None, colors are generated using the cmap of tekmc object.
                          
save_dir                : str, default: None
                          Directory where plot will be saved. If None, save_dir = ‘visualizations’.
</pre>
**Outputs**
<pre>
Scatter plot for displacements is saved as 'prob_hopp_displacement.pdf' in save_dir directory.
Scatter plot for velocities is saved as 'prob_hopp_velocity.pdf' in save_dir directory.
</pre>
