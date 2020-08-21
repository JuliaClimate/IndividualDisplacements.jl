## `solid_body_rotation.jl`

An idealized flow example, based on solid body rotation, also useful for unit testing.

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## `random_flow_field.jl`

A random flow field is generated on a doubly periodic grid, and used to advect a cloud of points. This illustrates defining a grid from scracth, and then simulating many trajectories at once.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## `global_ocean_circulation.jl`

A simulation over the global ocean based on a data-constrained,realistic, model:

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## More Examples

Two examples which reproduce trajectories computed by an earlier implementation in Fortran ([MITgcm/pkg/flt](https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#)): `detailed_look.jl` illustrates package features in more detail; `particle_cloud.jl` illustrates a computation of many trajectories at once. Also: `example123.jl`, `helper_functions.jl`, `example_CyclicArray.jl`

## Plotting Recipes 

Examples using three popular plotting packages (see `recipes_plots.jl`, `recipes_makie.jl`, `recipes_pyplot.jl`) are provided in the various notebooks.
