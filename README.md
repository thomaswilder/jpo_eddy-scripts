# jpo_eddy-scripts

The following files are intended to allow anyone to replicate the data used in https://doi.org/10.1175/JPO-D-22-0044.1 entitled 'The Response of a Baroclinic Anticyclonic Eddy to Relative Wind Stress Forcing'. The scripts are mostly written in MATLAB, except the vertical grid script written by Stewart and detailed in Stewart et al (2017). The energy scripts assumes the model data is on an Arakawa C-grid such as the MITgcm. 

If anything isn't working, do get in touch. I am happy to help.

## Initialisation scripts/functions

- init_mitgcm.m - Initialisation data script for baroclinic eddy. Requires ocean_vertical_grid.nc and geostrophic_vels.m

- geostrophic_vels.m - Function that calculates horizontal velocities in thermal wind balance using sea surface height and temperature.

- build_vertical_grid_kds.py - Function that generates vertical grid using Python including parameters chosen for this experiment. Written by Stewart https://github.com/kialstewart/vertical_grid_for_ocean_models

## Energy budget scripts

The following scripts use a 16 day time-mean, which is somewhat inherently built in, though modifications to this time can be made but require some fiddling about. The choice of time-mean needs to ensure a full rotation of the wind vector has taken place, e.g. 16 days = 5 wind rotations. Each script also includes the option to choose a domain integral or horizontal profile. Depending on horizontal resolution chosen, the scripts can either be run locally or on a hpc.

- eddy_energy_mean.m - Script that calculates mean eddy energy terms based on MITgcm input data.

- eddy_energy_turb.m - Script that calculates turbulent eddy energy terms. Can be split into separate jobs e.g. day 31 to day 150, then day 151 to 300 and so on. 


## General use functions

- geostrophic_uv.m - Similar to geostrophic_vels.m, but puts values at corner of grid cells since MITgcm output tracers at grid box centre.

- dvald.m - Calculates first x, y, and z derivative.


## MITgcm data files
Folder 'MITgcm_datfiles' contains the namelist and CPP header files that were used for this experiment. File 'job.slurm' is an example slurm script. 


