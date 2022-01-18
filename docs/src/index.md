# IndividualDisplacements.jl

[IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) computes elementary point displacements over a gridded Earth domain (e.g. a climate model `C-grid`). A typical application is the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons).

Inter-operability with common climate model grids and their representation in [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) is a central element ([docs](https://juliaclimate.github.io/MeshArrays.jl/dev/)). [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) can also read and plot trajectories using [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) output or [Ocean Drifting Buoy](https://doi.org/10.1002/2016JC011716) data ([movie](https://youtu.be/82HPnYBtoVo)). It was originally developed using [ECCOv4](https://eccov4.readthedocs.io/en/latest/) and [CBIOMES](https://cbiomes.readthedocs.io/en/latest/) ocean model simulations ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015), [Forget, 2018](http://doi.org/10.5281/zenodo.1343303)).

[![simulated particle movie (95m)](https://user-images.githubusercontent.com/20276764/90925145-ca90cc80-e3be-11ea-8eed-559307dcb925.png)](https://youtu.be/tsdf4fmYt1k)

[![Global Drifter Program data](https://user-images.githubusercontent.com/20276764/90924860-41799580-e3be-11ea-96bd-9a5784d00ecc.png)](https://youtu.be/82HPnYBtoVo)

