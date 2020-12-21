## Scope / Goals

A central goal of this package is to support scientific analysis of climate model output and observed drifts of e.g. plastics in the Ocean or chemical coumponds in the Atmosphere. To start, the package supports all types of gridded model output [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) by exploiting the [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) package ([docs found here](https://juliaclimate.github.io/MeshArrays.jl/dev/)). The tool box in this package also provides functions to ingest trajectory data which have been collected by the [Ocean Drifting Buoy](https://doi.org/10.1002/2016JC011716) Program over the real Ocean ([movie](https://youtu.be/82HPnYBtoVo)) or computed by MITgcm.


## Typical Workflow

As documented in the **examples**, the typical worflow is:

1. set up `FlowFields`
1. set up `Individuals`
1. displace them (`∫!`, 🚄, 📌)
1. record and post-process (🔴, 🔧)
1. go back to `2` and continue 

Steps `1` and `2` are done via constructors documented below in the `Data Structures` section. 
Step `3` and step `4` typically take place within the call to `∫!`. The latter can indeed readily post-process results recorded in 🔴 when `∫!` calls 🔧 before updating positions 📌. Since 🔴 is in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, it is easily manipulated or plotted after the fact (step `4` per se). Step `5` is optional. 

**See the examples** for more documentation regarding the ingestion of time varying flow fields, three-dimensional ocean trajectory simulations, process oriented configurations, as well as plotting and data formats. 

## Core Functions

The velocity interpolation funtions (🚄 used in step `3`; documented in the `Tool Box` section) interpolate flow fields to positions 📌, and `∫!` integrates the result forward in time. 

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method ∫, and post-processor 🔧 (all embedded within 𝐼).

```@docs
∫!
```

## Data Structures

The `Individuals` struct contains a `FlowFields` struct (incl. e.g. arrays), initial positions for the individuals, and the other elements (see below) involved in `∫!(𝐼,𝑇)`.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```
