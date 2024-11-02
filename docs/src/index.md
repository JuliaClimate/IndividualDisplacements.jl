# IndividualDisplacements.jl

[IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) computes elementary point displacements over gridded domains. Inter-operability with global climate model output (on [Arakawa](https://en.wikipedia.org/wiki/Arakawa_grids) `C-grids` in general) that can be represented via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) (see [these docs](https://juliaclimate.github.io/MeshArrays.jl/dev/)) is a central goal of [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl). 

The initial example suite relies on gridded ocean transports estimated in [OCCA](https://doi.org/10.7910/DVN/RNXA2A) ([Forget 2010](http://dx.doi.org/10.1175/2009JPO4043.1)), [ECCOv4](https://eccov4.readthedocs.io/en/latest/) ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)), and [CBIOMES](https://cbiomes.readthedocs.io/en/latest/) ([Forget, 2018](http://doi.org/10.5281/zenodo.1343303)). For the Atmosphere, an additional example based on an idealized model simulation is provided in [MITgcm.jl](https://gaelforget.github.io/MITgcm.jl/dev/). 

Typical applications include the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons). An illustration of the type of observations that [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl) can simulate is provided below. These data collected by [Ocean Drifting Buoy](https://doi.org/10.1002/2016JC011716) can be accessed via [OceanRobots.jl](https://gaelforget.github.io/OceanRobots.jl/dev/).

[![Global Drifter Program data](https://user-images.githubusercontent.com/20276764/90924860-41799580-e3be-11ea-96bd-9a5784d00ecc.png)](https://youtu.be/82HPnYBtoVo)

