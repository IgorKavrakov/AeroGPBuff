# AeroGPBuff
AeroGPBuff: Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Gaussian Processes

Cite:
To cite the paper, please use:
```
@article{Kavrakov2024AeroGPBuff,
title = {Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Gaussian Processes},
journal = {Journal of Wind Engineering and Industrial Aerodynamics},
pages = {103534},
year = {2023},
issn = {0266-8920},
doi = {https://doi.org/10.1016/j.probengmech.2023.103534},
author = {Igor Kavrakov, Guido Morgenthal and Allan McRobie}
}
```

AeroGPBuff is a Matlab-based computer code that verifies an aerodynamic force model based on Gaussian Processes (GP) with aerodynamic priors for aerodynamic analyses of linelike structures in turbulent wind conditions.
The scripts `Example1c_FlatPlateAerodynamicPrior_Forced.m`, `Example1d_FlatPlateCoupledMean_Flutter.m` and `Example1e_FlatPlateCoupledMean_Buffeting.m` contain the verification of the model for a Flat Plate.
The folder GP includes an implementation of Gaussian Processes, while the folder Example1_FlatPlateAnalytical contains the details for the analytical flat plate analysis and plots.
The details on of the utilized methods are given in the abovementioned article.

The accepted manuscript can be found on arXiv:
https://arxiv.org/abs/2406.15603 

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 

AeroGPBuff is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AeroGPBuff is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with AeroGPBuff.  If not, see <https://www.gnu.org/licenses/>.

Copyright (c) Igor Kavrakov, Guido Morgenthal, Allan McRobie 2024

