## Scope / Goals

A central goal of this package is to support scientific analysis of climate model output and observed drifts of e.g. plastics in the Ocean or chemical coumponds in the Atmosphere. To start, the package supports all types of gridded model output [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) by exploiting the [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) package ([docs found here](https://juliaclimate.github.io/MeshArrays.jl/dev/)). The tool box in this package also provides functions to ingest trajectory data which have been collected by the [Ocean Drifting Buoy](https://doi.org/10.1002/2016JC011716) Program over the real Ocean ([movie](https://youtu.be/82HPnYBtoVo)) or computed by MITgcm.


## Typical Workflow

As documented in the **examples**, the typical worflow is:

1. set up the `Individuals` data structure
1. displace them via `∫!`
1. post-process / analyze / plot
1. go back to `2` and continue

The velocity interpolation funtions (🚄 used in step 2; documented below) interpolate gridded output to positions 📌. Steps `3` and `4` are optional. Step `2` also provides the option to post-process results recorded in 🔴 when `∫!` calls 🔧 before updating positions 📌. Since 🔴 is in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, it is easily manipulated or plotted. Ingestion of time varying flow fields, three-dimensional ocean trajectory simulations, process oriented configurations, as well as plotting and data formats are further documented via the **examples**. 

## Core Functions

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method / solver ∫, and post-processing workflow 🔧 (all embedded within 𝐼).

```@docs
∫!
```

## Data Structures

The `Individuals` struct contains velocity fields (arrays), etc, and a record of properties diagnozed along the way.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```
