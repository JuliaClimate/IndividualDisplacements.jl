# IndividualDisplacements.jl

**IndividualDisplacements.jl** computes elementary point displacements over a gridded Earth domain (e.g. a climate model `C-grid`). A typical application is the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons).

Inter-operability with common climate model grids and their representation in [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) is a central element. The package can also read and plot trajectory simulation output from e.g. the [MITgcm](https://mitgcm.readthedocs.io/en/latest/?badge=latest). It was originally developed using [ECCOv4](https://eccov4.readthedocs.io/en/latest/) and [CBIOMES](https://cbiomes.readthedocs.io/en/latest/) ocean model simulations ([Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).

The `⬡` and `⬡!` functions compute the tracked point / individual / agent velocities. 

## List Of Examples

A solid-body-rotation example is used for unit testing:

```
test/runtests.jl
examples/SolidBodyRotation.jl
```

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

Two examples using `⬡!` and `update_locations.jl`:

```
examples/RandomFlow_fleet.jl
examples/GlobalDomain_fleet.jl
```

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

![LatLonCap300mDepth](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/LatLonCap300mDepth.png)

Plotting recipes using three popular plotting packages:

```
examples/recipes_plots.jl
examples/recipes_makie.jl
examples/recipes_pyplot.jl	
```

Three other examples using `⬡`  are documented in the `API Guide`:

```
examples/examples123.jl
examples/example2fleet.jl
examples/example2more.jl
```

## API Guide

```@index
```

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type,:function]
```

