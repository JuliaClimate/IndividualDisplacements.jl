## Summary

The four examples outlined below form a tutorial, and thus complement the rest of the package documentation. Afterwards, this section provides a listing of the other examples, plotting recipes, and tools which are included in the package. 

Users aiming to set up `IndividualDisplacements` for a different problem might also find examples of data wrangling codes in `helper_functions.jl` useful. These define grids and ingest velocity fields for the examples below. 

For plotting result, the examples use `recipes_plots.jl`, `recipes_makie.jl`, and `recipes_pyplot.jl` which demo three popular plotting packages.


## Single Particle Example

This example starts with a three-dimensional flow field `u,v,w`, initializes a single particle / individual position `📌`, and wraps everything up as data structure `𝐼`.

It then displaces the individual(s) position in `𝐼` by integrating its instantaneous velocity, [moving through space with the flow](https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field), over time `𝑇`. 

This is generally the main computation done in this package -- interpolating `u,v,w` to individual positions `𝐼.📌` on the fly, using `𝐼.🚄`, and integrating through time, using `𝐼.∫`.

Here, the idealized flow field consists of [rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), plus a convergent term, plus a sinking term. This generates a downward, converging spiral -- a relevant case in the Ocean.

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## Particle Set Example

A random flow field is generated on a doubly periodic grid, and used to advect a cloud of points. This illustrates defining a grid from scracth, and then simulating many trajectories at once.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## Global Ocean Circulation

A simulation over the global ocean based on a data-constrained, realistic, model:

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## Three Dimensional Paths

A simulation over the global ocean based on a data-constrained, realistic, model:

[![simulated particle movie (3D)](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)](https://youtu.be/twAAE_WUs_g)

## Tool Box, Etc.

- Tools included in `src/compute.jl` and `data_wrangling.jl`:
	- Velocity interpolaton functions for several array / grid types.
	- Preprocessing and postprocessing methods.
	- I/O routines to read / write results from / to file.

-  Examples reproducing trajectories that had been computed earlier in Fortran ([MITgcm/pkg/flt](https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#)) are `detailed_look.jl` and `particle_cloud.jl`. 

- For more see also: `example_CyclicArray.jl`, `example123.jl`, `helper_functions.jl`; and for plotting : `recipes_plots.jl`, `recipes_makie.jl`, `recipes_pyplot.jl` also in the `examples/` folder.

## Running The Examples

Running the examples requires `julia` and its relevant packages. Inputs get downloaded as needed upon running the examples. The first three steps below do this, and generate jupyter notebook versions of the examples (easy to rerun afterwards). Once everything is setup then user can just call examples directly in `julia` (i.e., skip to last step below).

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

#### 4. or alternatively, in the julia REPL:

```
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(p*"/../examples/worldwide/global_ocean_circulation.jl")
```
