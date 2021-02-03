## Velocity Interpolation

The `dxy_dt`, `dxy_dt!`, `dxyz_dt`, `dxyz_dt!`, etc functions compute the tracked individual velocity. 

```@docs
dxy_dt
dxy_dt!
dxyz_dt
dxyz_dt!
dxy_dt_replay
dxy_dt_CyclicArray
```

## Setup And Postprocessing 

Convenience functions to initialize a simulation and posprocess the output are provided. 

```@docs
initialize_lonlat
randn_lonlat
postprocess_MeshArray
postprocess_xy
```

## Read Output From File 

Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file using, respectively `MITgcmTools.read_flt` or  `read_drifters`.

```@docs
read_drifters
```
