# VLEO Aerodynamics CORE
This repository contains functions to calculate the aerodynamic forces and moments acting on a satellite in a VLEO environment. The code is based on the Aerodynamic Database for Satellites [ADBSat](https://github.com/nhcrisp/ADBSat) and the corresponding references [1] and [2].

This repository allows for a simple rotation of geometries without providing new geometry files. Instead, rotation directions and hinge points are provided for each body. Using these, a rotation can be performed by specifying the rotation angle. This allows for a usage of this aerodynamics model within a simulator where the orientation of aerodynamic control surfaces is controlled by the attitude control system and needs to be updated continuously.

The functions in this repository expect the space_math_utilities class to be available. It is provided by the following repository: [Space Math Utilities](https://git.ifr.uni-stuttgart.de/ifr-space/space_math_utilities).

## References
[1]: D. Mostaza-Prieto, “CHARACTERISATION AND APPLICATIONS OF AERODYNAMIC TORQUES ON SATELLITES,” University of Manchester, 2017.

[2]: L. A. Sinpetru, N. H. Crisp, D. Mostaza-Prieto, S. Livadiotti, and P. C. E. Roberts, “ADBSat: Methodology of a novel panel method tool for aerodynamic analysis of satellites,” Computer Physics Communications, vol. 275, p. 108326, 2022, doi: https://doi.org/10.1016/j.cpc.2022.108326.

