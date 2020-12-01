## Summary

The four examples outlined below form a tutorial, and thus complement the rest of the package documentation. Afterwards, we provide a listing of other examples, plotting recipes, and tools included in the package. 

A user seeking to configure `IndividualDisplacements.jl` for a new problem might also find, hopefully useful, examples of data wrangling codes in `helper_functions.jl`. They define grids and ingest velocity fields for the examples below. 

Output is in [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format which readily provides powerful and convenient analysis methods. Plotting trajectories in space and time, for example, can be done as in `recipes_plots.jl`, `recipes_makie.jl`, and `recipes_pyplot.jl` (see examples).

## Single Particle Example

Here we start with a three-dimensional flow field `u,v,w`, initialize a single particle / individual position `📌`, and wrap everything up within a custom data structure `𝐼`.

The individual(s) in `𝐼` is then displaced by integrating its instantaneous velocity, [moving along through space](https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field), over time `𝑇`. 
This is generally the main computation done in this package -- interpolating `u,v,w` to individual positions `𝐼.📌` on the fly, using `𝐼.🚄`, and integrating through time, using `𝐼.∫`.

The flow field consists of [rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), plus a convergent term, plus a sinking term in the third direction. This generates a downward, converging spiral -- a idealized version of a relevant case in the Ocean.

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## Particle Set Example

Here we illustrate (1) the simulation of an ensemble of particles and (2) how one simply goes from a velocity array to solving for trajectories using `IndividualDisplacements.jl`. 

The included convenience function (`setup_point_cloud`) defines a grid based on input array dimensions, adds the initial condition and time range, and returns the `Individuals` data structure `𝐼`. All that is left to do at this stage is call `∫!(𝐼)`.

The prototype for this example is based on a randomly generated flow field in a doubly periodic gridded domain. Exercises include the non-periodic domain case, statistics made easy via `DataFrames.jl`, and replacing the flow field with your own.

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
