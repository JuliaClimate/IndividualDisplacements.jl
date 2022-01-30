## User Guide

!!! note
    The [Examples](@ref) section illustrates the features documented here.

### Standard Workflow

As documented in the **examples**, the typical worflow is:

1. set up `FlowFields` data structure
1. set up `Individuals` with initial position `📌`
1. displace `Individuals` (by	`∫🚄dt`) according to `FlowFields`
1. post-process (`🔧`) and record diagnostics in `🔴`
1. go back to `step 2` and continue if needed

The data structures for steps `1` and `2` are documented below. Both steps `3` and `4` normally take place as part of `∫!` which post-processes results, using 🔧, records them in 🔴, and updates individual positions 📌 at the end. As a [DataFrame](https://juliadata.github.io/DataFrames.jl/latest/), 🔴 is easily manipulated, plotted, or saved in step `4` (or after the fact).

## Overview

A central goal of this package is to support scientific analysis of model simulations and observations of materials and particles within the Climate System. The package scope thus includes, for example, drifting plastics in the Ocean or chemical compounds in the Atmosphere. 

As a starting point, the package supports all types of gridded model output (on Arakawa C-grids) from the [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) via the [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) package ([docs found here](https://juliaclimate.github.io/MeshArrays.jl/dev/)). 

By convention, `IndividualDisplacements.jl` expects input flow fields to be provided in a uniform fashion (see [`FlowFields`](@ref)) summarized below: 

1. normalized to grid index units (i.e. in 1/s rather than m/s units)
1. positive towards increasing indices
1. using the Arakawa C-grid, with `u` (resp `v`) staggered by `-0.5` point in direction `1` (resp `2`) from grid cell centers. 

The [Examples](@ref) section documents various, simple methods to prepare and ingest such flow fields (time varying or not; in 2D or 3D) and derive individual displacements / trajectories from them. They cover simple grids often used in process studies, Global Ocean simulations normally done on more complex grids, plotting tools, and data formats. 

For an overview of the examples, please refer to the **example guide**. The rest of this section is focused on the package's **data structures** and **core functions**.

## Data Structures

The `Individuals` struct contains a `FlowFields` struct (incl. e.g. arrays), initial positions for the individuals, and the other elements (e.g. functions) involved in `∫!(𝐼,𝑇)` as documented hereafter.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

## Main Functions

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method ∫, and post-processor 🔧 (all embedded within 𝐼).

```@docs
∫!
```

The velocity interpolation functions (🚄) carry out the central computation of this package -- interpolating gridded flow fields to individual positions 📌. It is normally called via `∫!` to integrate velocity 🚄 over a chosen time period. 

- Velocity interpolation for several array and grid types.
- Preprocessing and postprocessing methods.
- I/O routines to read (write) results from (to) file.

and other functionalities provided in `src/compute.jl` and `src/data_wrangling.jl` are further documented in the _Tool Box_ section. Ingestion of trajectory data which have been collected by the [Ocean Drifting Buoy Program](https://doi.org/10.1002/2016JC011716) ([movie](https://youtu.be/82HPnYBtoVo)) is also supported.

