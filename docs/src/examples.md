In summary ... 

## `solid_body_rotation.jl`

An idealized flow example, based on solid body rotation, also useful for unit testing.

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## `random_flow_field.jl`

A random flow field is generated on a doubly periodic grid, and used to advect a cloud of points. This illustrates defining a grid from scracth, and then simulating many trajectories at once.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## `global_ocean_circulation.jl`

A simulation over the global ocean based on a data-constrained, realistic, model:

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## `three_dimensional_ocean.jl`

A simulation over the global ocean based on a data-constrained, realistic, model:

[![simulated particle movie (3D)](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)](https://youtu.be/twAAE_WUs_g)

## Tools And More


- Plotting: `recipes_plots.jl`, `recipes_makie.jl`, `recipes_pyplot.jl` demonstrate three popular plotting packages (see examples).

- Tools (see `compute.jl`, `data_wrangling.jl`):
	- Velocity interpolaton functions for several array / grid types.
	- Preprocessing and postprocessing methods.
	- I/O routines to read / write results from / to file.

- Two examples that reproduce trajectories computed earlier in Fortran ([MITgcm/pkg/flt](https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#)): `detailed_look.jl` illustrates package features in more detail; `particle_cloud.jl` illustrates a computation of many trajectories at once. 

- For more see also: `example_CyclicArray.jl`, `example123.jl`, `helper_functions.jl` in the `examples/` folder.

## How-To Run The Examples

#### 1. Download the examples folder:

```
git clone https://github.com/JuliaClimate/IndividualDisplacements.jl
julia --project=IndividualDisplacements.jl/docs/
```

#### 2. Get all needed packages and `IndividualDisplacements.jl`:

```
using Pkg
Pkg.activate("IndividualDisplacements.jl/docs/")
Pkg.instantiate()
Pkg.add("IndividualDisplacements")
```

#### 3. Generate jupyter notebook using `Literate.jl`:

```
using Literate
Literate.notebook("IndividualDisplacements.jl/examples/basics/solid_body_rotation.jl", ".", execute = true, documenter = false)
Literate.notebook("IndividualDisplacements.jl/examples/worldwide/three_dimensional_ocean.jl", ".", execute = true, documenter = false)
```

And so on and so forth -- see documentation for a list of examples.
