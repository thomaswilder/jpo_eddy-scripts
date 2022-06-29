# jpo_eddy-scripts

The following files are intended to allow anyone to replicate the data used in https://doi.org/10.1175/JPO-D-22-0044.1 entitled 'The Response of a Baroclinic Anticyclonic Eddy to Relative Wind Stress Forcing'. The scripts are mostly written in MATLAB, except the vertical grid script written by Stewart and detailed in Stewart et al (2017). 

## Initialisation scripts

- init_mitgcm.m - Initialisation data for baroclinic eddy. Requires ocean_vertical_grid.nc and geostrophic_vels.m

- geostrophic_vels.m - Calculates horizontal velocities in thermal wind balance using sea surface height and temperature.

- build_vertical_grid_kds.py - Generates vertical grid using Python including parameters chosen for this experiment. Written by Stewart https://github.com/kialstewart/vertical_grid_for_ocean_models

## Energy budget scripts

The following scripts use a 16 day time-mean, which is somewhat inherently built in, though modifications to this time can be made but require some fiddling about. The choice of time-mean needs to ensure a full rotation of the wind vector has taken place, e.g. 16 days = 5 wind rotations. Each script also includes the option to choose a domain integral or horizontal profile. Depending on horizontal resolution chosen, the scripts can either be run locally or on a hpc.

- eddy_energy_mean.m - Calculates mean eddy energy terms.

- eddy_energy_turb.m - Calculates turbulent eddy energy terms. 
