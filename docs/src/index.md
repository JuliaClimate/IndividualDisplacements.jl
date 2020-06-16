# IndividualDisplacements.jl

**IndividualDisplacements.jl** computes elementary point displacements over a gridded Earth domain (e.g. a climate model `C-grid`). A typical application is the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons).

Inter-operability with common climate model grids and their representation in [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) is a central element. The package can also read and plot trajectory simulation output from e.g. the [MITgcm](https://mitgcm.readthedocs.io/en/latest/?badge=latest). It was originally developed using [ECCOv4](https://eccov4.readthedocs.io/en/latest/) and [CBIOMES](https://cbiomes.readthedocs.io/en/latest/) ocean model simulations ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

[![simulated particle movie (300m)](https://user-images.githubusercontent.com/20276764/84767001-b89a4400-af9f-11ea-956f-2e207f892c4f.png)](https://youtu.be/M6vAUtIsIIY)
