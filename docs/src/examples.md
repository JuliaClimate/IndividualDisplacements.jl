## Summary

The four examples outlined below form a tutorial of sort, and thus complement the rest of the package documentation. The corresponing data wrangling examples in `flow_fields.jl` define grids and ingest velocity fields. They hopefully provide a useful jumping off point in order to configure `IndividualDisplacements.jl` for new problems.

Output is in [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format which comes with powerful and convenient analysis methods. Plotting results in space and time can be done as in `recipes_plots.jl`, `recipes_makie.jl`, and `recipes_pyplot.jl` -- see the examples.

## Particle Set Example

Simulation of an ensemble of displacements (and trajectories) in a simple 2D configuration. 

The `convert_to_FlowFields` convenience function defines a simplified gridded domain that matches the velocity array size, adds a time range, and returns a `FlowFields`
data structure `𝐹`. 
All that is left to do at this stage is to define initial conditions, put them together with `𝐹` within the `Individuals` data structure `𝐼`, and call `∫!(𝐼)`.

The prototype for this example is based on a randomly generated flow field in a doubly periodic gridded domain. 
Exercises include the non-periodic domain case, statistics made easy via `DataFrames.jl`, and replacing the flow field with your own.

![RandomFlow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

## Single Particle Example

Setup a three-dimensional flow field `u,v,w`, initialize a single particle / individual position `📌`, and wrap everything up within an `Individuals` data structure `𝐼`.

`𝐼` is displaced by integrating the individual velocity, [moving along through space](https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field), over time `𝑇`.  This is the main computation done in this package -- interpolating `u,v,w` to individual positions `𝐼.📌` on the fly, using `𝐼.🚄`, and integrating through time, using `𝐼.∫`.

The flow field consists of [rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), plus a convergent term, plus a sinking term in the vertical direction. This flow field generates a downward, converging spiral -- a idealized version of a relevant case in the Ocean.

![SolidBodyRotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

## Global Ocean Circulation

A simulation of floating particles over the Global Ocean which illustrates (1) using time variable velocity fields, (2) global connections, (3) particle re-seeding, and (4) output statistics. 

The flow field is based on a data-constrained ocean model solution. The problem is configured in a way to mimic, albeit very crudely, the near-surface tranport of plastics or planktons.

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

## Three Dimensional Pathways

A simulation of particles that follow the three-dimensional ocean circulation. This example illustrates (1) the 3D case in a relatistic configuration, (2) tracking the advent or origin of a water patch, and (3) multifacted visualizations in 3D.

The flow field is based on a data-constrained, realistic, ocean model. The problem configuration mimics, albeit very approximately, ocean tracers / coumpounds transported by water masses .

[![simulated particle movie (3D)](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)](https://youtu.be/twAAE_WUs_g)

## Running The Examples

Running the examples requires `julia` and its relevant packages. Inputs get downloaded as needed upon running the examples. The following commands ought to run on either your local computer or e.g [here in the cloud](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/master).

Maybe the simplest option is to generate jupyter notebooks for all examples as part of the documentation build (see below). The generated notebooks should be easy to re-run afterwards. To generate and run the notebooks, start `julia` in the terminal or a jupyter notebook at the `IndividualDisplacements.jl/` folder. Then try:

```
using Pkg
Pkg.activate("docs/")
Pkg.instantiate()
Pkg.add(path="./")
include("docs/make.jl")
```
 
Or alternatively, to generate notebooks with more detailed control:

```
using Literate
set1 = ["solid_body_rotation","random_flow_field","detailed_look","particle_cloud"]
[Literate.notebook("examples/basics/"*i*".jl", ".", execute = true, documenter = false) for i in set1]
set2 = ["global_ocean_circulation","three_dimensional_ocean"]
[Literate.notebook("examples/worldwide/"*i*".jl", ".", execute = true, documenter = false) for i in set2]        
```

Or as another alternative one can simply run the julia source file:

```
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(p*"/../examples/worldwide/global_ocean_circulation.jl")
```

## More Examples

-  Examples reproducing trajectories that had been computed earlier in Fortran ([MITgcm/pkg/flt](https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#)) are `detailed_look.jl` and `particle_cloud.jl`. 
- For more examples, see: `example_CyclicArray.jl`, `example123.jl`, `helper_functions.jl`.
