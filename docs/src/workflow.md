## Scope / Goals

A central goal of this package is to support scientific analysis of climate model output and observed drifts of e.g. plastics in the Ocean or chemical coumponds in the Atmosphere. 

To start, the package supports all types of gridded model output from the [MIT General Circulation Model](https://mitgcm.readthedocs.io/en/latest/?badge=latest) via the [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) package ([docs found here](https://juliaclimate.github.io/MeshArrays.jl/dev/)) and ingestion of trajectory data which have been collected by the [Ocean Drifting Buoy Program](https://doi.org/10.1002/2016JC011716) ([movie](https://youtu.be/82HPnYBtoVo)).


## Typical Workflow

As documented in the **examples**, the typical worflow is:

1. set up `FlowFields` data structure
1. set up `Individuals` with initial position 📌
1. displace `Individuals` according to `FlowFields` (via	∫! and 🚄)
1. post-process and record results (via 🔧 and 🔴)
1. go back to `2` and continue if needed

Steps `1` and `2` are done via data structures documented below. Both steps `3` and step `4` often take place within the call to `∫!`which can readily post-process results via 🔧 before recording them in 🔴 and finally updating positions 📌. Since 🔴 is in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, it is easily manipulated, plotted, or saved after the fact (step `4` per se).

**The examples** document simple methods to ingest time varying flow fields, three-dimensional ocean simulations, process oriented configurations, plotting tools, and data formats. For an overview of the examples, please refer to the **example guide**. The rest of this section is focused on the package's **data structures** and **core functions**.

## Data Structures

The `Individuals` struct contains a `FlowFields` struct (incl. e.g. arrays), initial positions for the individuals, and the other elements (e.g. functions) involved in `∫!(𝐼,𝑇)` as documented hereafter.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

## Core Functions

The velocity interpolation funtions (🚄 used in step `3`; documented in the `Tool Box` section) interpolate flow fields to positions 📌, and `∫!` integrates the result forward in time. 

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method ∫, and post-processor 🔧 (all embedded within 𝐼).

```@docs
∫!
```
