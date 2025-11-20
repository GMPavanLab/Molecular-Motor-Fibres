Figure 4 a and b. 
from molecular_muscle/mkfiber_cluster/run/run & molecular_muscle/mkfiber_cluster/system

To render bonds, use cg_bonds. Have all .itp files in location specified in .top
To rendor the trajectory, use .pdb + .xtc, remove W and CA from .top file. 
To render the full system, use .gro and whole .top file
cg_bonds -top 5_layers_topol_sol.top 