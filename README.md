# ArtemisTools

This is a rough package that manipulates data from the ARTEMIS simulation. I've included some functionality that might not be necessary 
(e.g. loading IllustrisTNG data in the LoadData.py script), so just comment it out if you don't need it. 

The key packages you'll need are:
* numpy
* scipy
* matplotlib
* pandas
* h5py
* eagle_IO

The first packages are easily installable, for eagle_IO refer to https://github.com/christopherlovell/eagle_IO. Alternatively, replace
the relevant sections of code with your chosen method of loading EAGLE-like data.

# Artemis.py

This script does the heavy lifting. For a specified ARTEMIS run and TAG, you have the option to load in star, gas and dark matter particle data.
The data is reduced to be that related to the main galaxy in the simulation, and centres the particle coordinates on the galaxy's centre of potential. 
Following this, you have the option to rotate particle coordinates such that they are in the reference frame of the best-fitting ellipsoid of the 
stellar particle distribution, defined by the iterative reduced inertia tensor. You can also compute the NFW profile, retrieve the dimensions of the 
best-fitting ellipsoid, and save the particle positions to a txt file.

# LoadData.py

This script loads in the simulation data. File locations are relevant to the Astrophysics Research Institute. Each class is related to a different simulation,
and the methods therein load in different aspects of the simulation data, such as partdata, header information, subhalo and FOF group information. 

# tens_3d.py

This script handles the computation of the interative reduced inertia tensor. See e.g. https://ui.adsabs.harvard.edu/abs/2021MNRAS.505...65H/abstract for 
more information.
