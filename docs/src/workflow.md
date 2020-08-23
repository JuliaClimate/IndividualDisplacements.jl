## Basic Workflow

- set up `Individuals`
- displace them via `∫!`
- post-process and / or repeat

(e.g. see `global_ocean_circulation.jl`)

## Data Structures

The `Individuals` struct contains velocity fields (arrays), etc, and a record of properties diagnozed along the way.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

## Core Functions

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period `𝑇`:

```@docs
∫!
```
