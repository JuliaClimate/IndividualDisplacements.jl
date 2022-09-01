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
add_lonlat!
```

Basic geography support:

```@docs
gcdist
diff
stproj
stproj_inv
randn_lonlat
interp_to_lonlat
interp_to_xy
nearest_to_xy
```

## Toy Problems

These are used to demonstrate and test the package functionalities:

```@docs
random_flow_field
vortex_flow_field
```

## Read External Files

Trajectory simulated by the [MITgcm](https://mitgcm.readthedocs.io/en/latest/?badge=latest) or observed by the [global drifter program](https://www.aoml.noaa.gov/phod/gdp/index.php) can be read from file using, respectively `MITgcmTools.read_flt` (from [MITgcmTools.jl](https://gaelforget.github.io/MITgcmTools.jl/dev/)) or  `OceanRobots.drifters_hourly_read` (from [OceanRobots.jl](https://gaelforget.github.io/OceanRobots.jl/dev/)).
