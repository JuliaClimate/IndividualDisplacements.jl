
## API Guide

See `examples/worldwide/global_ocean_circulation.jl` for an example:

```@docs
start!
displace!
reset!
```

### Data Structures

The main data type used is `Individuals` which contains arrays and a dataframe to store the output diagnostics.

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```

### Velocity Interpolation

The `⬡` and `⬡!` functions compute the tracked point / individual / agent velocities. 

```@docs
⬡!
⬡
```

### Setup And Postprocessing 

Convenience functions to initialize a simulation and posprocess the output are provided. 

```@docs
initialize_lonlat
randn_lonlat
postprocess_lonlat
postprocess_xy
```

### Read Output From File 

Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file.

```@docs
read_flt
read_drifters
```
