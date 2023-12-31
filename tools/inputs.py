import os

# ****************************** Important parameters *************************************

# Required Parameters

topology_file       = 'topo.gro'
trajectory_file     = 'traj.dcd'
atom_name           = 'A'
timestep            = 0.01
spacings_list       = [0.25, 0.26, 0.27]
n_walks             = 5000
n_components        = 3
path_to_src         = "/home/subhadeepd/tekmc/src"

# Some optional parameters

stride              = 1
md_file             = None
n_cpus              = 8
cmap                = "tab20"
save_dir            = 'production_5us'
n_steps             = 500000
save_random_walks   = False
hopping_distance_analysis = True

# ****************************** Important parameters *************************************

# To modify other default parameters, modify defaults.py in src folder.
# Look up https://github.com/PKMLab/tekmc for more details on the parameters.
