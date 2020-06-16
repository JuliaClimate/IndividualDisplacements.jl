
## API Guide

The `⬡` and `⬡!` functions compute the tracked point / individual / agent velocities. 

```@docs
⬡!
⬡
```

### Setup And Postprocessing 

Convenience functions to initialize a simulation and posprocess the output are provided. 

```@docs
setup_periodic_domain
initialize_randn
initialize_gridded
postprocess_lonlat
postprocess_xy
```

### Read Output From File 

Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file.

```@docs
read_flt
read_drifters
```

### Types

```@autodocs
Modules = [IndividualDisplacements]
Order   = [:type]
```
