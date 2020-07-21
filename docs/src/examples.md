## `random_flow_field.jl`

A random flow field is generated on a doubly periodic grid, and used to advect a cloud of points. This illustrates `⬡!` and `update_locations.jl` to compute many trajectories at once.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## `global_ocean_circulation.jl`

A real ocean simulation:

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## `solid_body_rotation.jl`

A very simple solid-body-rotation example which is used for unit testing (via `test/runtests.jl`):

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## `detailed_look.jl`, `particle_cloud.jl`

Two examples which use `⬡` and reproduce trajectories computed by a different (fortran) tool. `detailed_look.jl` illustrates package features in more detail. `particle_cloud.jl` illustrates a computation of many trajectories at once.

## Plotting recipes 

Using three popular plotting packages:

```
examples/recipes_plots.jl
examples/recipes_makie.jl
examples/recipes_pyplot.jl	
```

