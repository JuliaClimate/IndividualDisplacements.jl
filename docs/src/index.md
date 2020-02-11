# IndividualDisplacements.jl

**IndividualDisplacements.jl** computes elementary point displacements over a gridded Earth domain (e.g. a climate model C-grid). A typical application is the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons).

Inter-operability with popular climate model grids and their representation in [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) is a central element. The package can also read and plot trajectory simulation output from e.g. the [MITgcm](https://mitgcm.readthedocs.io/en/latest/?badge=latest). It was originally developed using [ECCOv4](https://eccov4.readthedocs.io/en/latest/) and [CBIOMES](https://cbiomes.readthedocs.io/en/latest/) ocean model simulations ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).

The `VelComp!` and `VelComp` functions compute the velocity of tracked points. `tests/runtests.jl` uses solid body rotation as a benchmark (see below).

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/SolidBodyRotation.png)

## List Of Examples

Solid body rotation is used for unit testing:

```
test/runtests.jl
examples/SolidBodyRotation.jl
```

The two examples documented under `API` + extenstions to `example2`:

```
src/examples.jl
examples/example2fleet.jl
examples/example2more.jl
```

An intermediate example: `examples/PeriodicChannel_fleet.jl`

Two examples use `VelComp!` and `update_locations.jl`:

```
examples/PeriodicDomainRandom_fleet.jl
examples/GlobalDomain_fleet.jl
```

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/LatLonCap300mDepth.png)

## API Guide

```@index
```

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type,:function]
```

