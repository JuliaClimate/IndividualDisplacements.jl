## Velocity Interpolation

The `dxdt!` etc functions compute the tracked individual velocity. 

```@docs
dxdt!
dxy_dt_replay
dxy_dt_CyclicArray
```

## Setup And Postprocessing 

Convenience functions to initialize a simulation and post-process the output are provided. 

```@docs
postprocess_xy
postprocess_MeshArray
randn_lonlat
interp_to_lonlat
interp_to_xy
nearest_to_xy
gcdist
```

## Read Output From File 

Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file using, respectively `MITgcmTools.read_flt` or  `read_drifters`.

```@docs
read_drifters
read_velocities
read_mds
```
