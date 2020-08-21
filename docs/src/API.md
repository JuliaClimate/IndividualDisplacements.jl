
## API Guide

The typical workflow is:

- set up `Individuals`
- displace them via `∫!`
- post-process and / or repeat

(e.g. see `global_ocean_circulation.jl`)

### Core Functionalities

`∫!(𝐼,𝑇)` displaces individuals 𝐼 continuously over time period `𝑇`:

```@docs
∫!
```

### Data Structures

The `Individuals` struct contains velocity fields (arrays), etc, and a record of properties diagnozed along the way.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

### Toolbox

- Velocity functions, interpolating from gridded fields, for different array types.
- Preprocessing and postprocessing methods.
- I/O routines to read / write results from / to file.

#### Velocity Interpolation

The `dxy_dt`, `dxy_dt!`, `dxyz_dt`, etc functions compute the tracked individual velocity. 

```@docs
dxy_dt
dxy_dt!
dxyz_dt
dxy_dt_replay
dxy_dt_CyclicArray
```

#### Setup And Postprocessing 

Convenience functions to initialize a simulation and posprocess the output are provided. 

```@docs
initialize_lonlat
randn_lonlat
postprocess_lonlat
postprocess_xy
```

#### Read Output From File 

Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file using, respectively `MITgcmTools.read_flt` or  `read_drifters`.

```@docs
read_drifters
```
