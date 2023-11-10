import os

# ****************************** Important parameters *************************************

# Required Parameters

topology_file       = 'system/50mer_s_CO2_308K_2e6Pa_CO2.gro'
trajectory_file     = 'system/50mer_s_CO2_308K_2e6Pa_CO2.dcd'
atom_name           = 'CO2'
timestep            = 0.01
spacings_list       = [0.18,0.19,0.2]
n_walks             = 500
n_components        = 3
path_to_src         = "/home/subhadeepd/tekmc/src"

# Some optional parameters

stride              = 1
md_file             = None
n_cpus              = 24
cmap                = "tab20"
save_dir            = 'production_5us'
n_steps             = 100
save_random_walks   = False
hopping_distance_analysis = False

# ****************************** Important parameters *************************************

# To modify other default parameters, modify defaults.py in src folder.
# Look up https://github.com/PKMLab/tekmc for more details on the parameters.
