## `random_flow_field.jl`

A random flow field is generated on a doubly periodic grid, and used to advect a cloud of points. This illustrates defining a grid from scracth, and then simulating many trajectories at once.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## `global_ocean_circulation.jl`

A simulation over the global ocean based on a data-constrained,realistic, model:

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## `solid_body_rotation.jl`

A very simple example, based on solid body rotation, which is also used for unit testing (via `test/runtests.jl`):

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## `detailed_look.jl`, `particle_cloud.jl`

Two examples which reproduce trajectories computed by an earlier implementation in Fortran ([MITgcm/pkg/flt](https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#)). `detailed_look.jl` illustrates package features in more detail. `particle_cloud.jl` illustrates a computation of many trajectories at once.

## Plotting recipes 

Examples using three popular plotting packages:

```
examples/recipes_plots.jl
examples/recipes_makie.jl
examples/recipes_pyplot.jl	
```

