
## Examples

The four examples outlined below form a tutorial of sorts, and thus complement the [User Guide](@ref). They rely on `flow_fields.jl` to define grids and ingest velocity fields. The main goal is to provide a useful jumping off point in order to configure `IndividualDisplacements.jl` for new problems.

Output is in [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format which comes with powerful and convenient analysis methods. Plotting results in space and time can be done as in `recipes_plots.jl`, `recipes_makie.jl`, and `recipes_pyplot.jl` -- see the examples.

To run an example, the recommended method is to copy the corresponding `notebook (code)` link, paste into the [Pluto.jl](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto) prompt, and click `open`.

## Simple Two-Dimensional Flow

[notebook (html)](random_flow_field.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/basics/random_flow_field.jl)

Simulate an ensemble of displacements (and trajectories) in a simple 2D configuration. 

The `convert_to_FlowFields` convenience function defines a simplified gridded domain that matches the velocity array size, adds a time range, and returns a `FlowFields`
data structure `𝐹`. 
All that is left to do at this stage is to define initial conditions, put them together with `𝐹` within the `Individuals` data structure `𝐼`, and call `∫!(𝐼)`.

Exercises include the non-periodic domain case, statistics made easy via `DataFrames.jl`, and replacing the flow field with your own.

## Simple Three-Dimensional Flow

[notebook (html)](solid_body_rotation.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/basics/solid_body_rotation.jl)

Set up a three-dimensional flow field `u,v,w`, initialize a single particle at position `📌`, and wrap everything up within an `Individuals` data structure `𝐼`.

`𝐼` is displaced by integrating the individual velocity, [moving along through space](https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field), over time `𝑇`.  This is the main computation done in this package -- interpolating `u,v,w` to individual positions `𝐼.📌` on the fly, using `𝐼.🚄`, and integrating through time, using `𝐼.∫`.

The flow field consists of [rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), plus a convergent term, plus a sinking term in the vertical direction. This flow field generates a downward, converging spiral -- a idealized version of a relevant case in the Ocean.

## Global Ocean Circulation

[notebook (html)](global_ocean_circulation.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/worldwide/global_ocean_circulation.jl)

A simulation of floating particles over the Global Ocean which illustrates (1) using time variable velocity fields, (2) global connections, (3) particle re-seeding, and (4) output statistics. 

The flow field is based on a data-constrained ocean model solution. The problem is configured in a way to mimic, albeit very crudely, the near-surface tranport of plastics or planktons.

## Three Dimensional Pathways

[notebook (html)](three_dimensional_ocean.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/worldwide/three_dimensional_ocean.jl)

A simulation of particles that follow the three-dimensional ocean circulation. This example illustrates (1) the 3D case in a relatistic configuration, (2) tracking the advent or origin of a water patch, and (3) multifacted visualizations in 3D.

The flow field is based on a data-constrained, realistic, ocean model. The problem configuration mimics, albeit very approximately, ocean tracers / coumpounds transported by water masses .

## Additional Examples

- Interactive UI (Pluto.jl) : [notebook (html)](interactive_UI.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/worldwide/interactive_UI.jl)

- Particle cloud (MITgcm) : [notebook (html)](../particle_cloud/index.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/basics/particle_cloud.jl)

- Detailed look (MITgcm) : [notebook (html)](../detailed_look/index.html) ➭ [notebook (code)](https://github.com/JuliaClimate/IndividualDisplacements.jl/blob/master/examples/basics/detailed_look.jl)

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

[![simulated particle movie (3D)](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)](https://youtu.be/twAAE_WUs_g)
