## tekmc.find_D(save_dir=None, spacings_list=None)

### Estimate and write the diffusion coefficients to a file for the given list of spacings.

**Optional parameters**
<pre>
spacings_list		: list of lists, default: None
			  List of spacings. Ex: [[0.23],[0.22,0.23,0.24]], [[0.22],[0.13],[0.5]]
                  	  If None, D is estimated for all grid sizes present in 'All_MSDs' folder.
                  
save_dir		: str, default: None
  			  Diffusion coefficients are written to 'Diffusion_coefficients.dat' in this directory. 
			  If None, save_dir = ‘visualizations’.
</pre>
**Outputs**
<pre>
Diffusion coefficients are written to 'Diffusion_coefficients.dat' in save_dir directory.
</pre>
