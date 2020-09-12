## Typical Workflow

As seen in the **examples**, the typical worflow is:

1. set up `Individuals`
1. displace them via `∫!`
1. post-process / analyze / plot
1. go back to `2` and continue

where steps `3` and `4` are optional. 

In fact, step `2` also readily provides the option to post-process results in 🔴 via 🔧 within `∫!` -- 🔧 is called just before updating the individual positions 📌 and exiting `∫!`. Since 🔴 is in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, it is easily manipulated or plotted. 

Please refer to the **examples** section for more on e.g. the use of time varying flow fields, three-dimensional Ocean trajectory simulations, process oriented configurations, as well as plotting and data formats.

A central goal of this package is to support scientific analysis of climate model output and observed drifts of e.g. plastics in the Ocean or chemical coumponds in the Atmosphere. 

To start, the package supports all types of gridded model output [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) by exploiting the [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) package ([docs](https://juliaclimate.github.io/MeshArrays.jl/dev/)). The 🚄 funtions documented as part of the **Tool Box** perform the interpolation from gridded output to positions 📌.

The **Tool Box** section also provides functions to ingest trajectory data which have been collected by the [Ocean Drifting Buoy](https://doi.org/10.1002/2016JC011716) Program over the real Ocean ([movie](https://youtu.be/82HPnYBtoVo)).

## Data Structures

The `Individuals` struct contains velocity fields (arrays), etc, and a record of properties diagnozed along the way.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

## Core Functions

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method / solver ∫, and post-processing workflow 🔧 (all embedded within 𝐼).

```@docs
∫!
```
