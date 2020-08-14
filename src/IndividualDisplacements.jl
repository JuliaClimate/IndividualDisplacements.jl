module IndividualDisplacements

using MeshArrays, OrdinaryDiffEq, StatsBase, DataFrames, Random
using NetCDF, Dates, CFTime, CSV, UnPack, CyclicArrays

include("compute.jl")
include("read.jl")
include("update_locations.jl")
include("data_wrangling.jl")

⬡! = dxy_dt!
⬡ = dxy_dt

export ⬡!, ⬡, dxyz_dt, dxy_dt_CyclicArray, dxy_dt_replay
export initialize_gridded, initialize_lonlat, randn_lonlat
export postprocess_lonlat, add_lonlat!, postprocess_xy
export read_flt, read_mds, read_drifters, read_uvetc
export Individuals, set_up_𝑃, update_𝑃!

day=86400.0
mon=365/12*day
solver_default(prob) = solve(prob,Euler(),dt=2*day)
𝑃_default = ( 𝑇 = [-0.5*mon,0.5*mon] , 🔄 = update_𝑃!,
              u0=[] , u1=[] , v0=[] , v1=[] )
tr_default = DataFrame( ID=[], x=[], y=[], t = [], lon=[], lat=[], fid=[])

"""
    struct Individuals{T}

Contains: xy, id, tr, etc
```
i=Individuals{Float32}(xy=zeros(3,2),id=1:2)
```
"""
Base.@kwdef struct Individuals{T}
   xy  ::Array{T,2} = Array{T,2}(undef, Tuple(Int.(zeros(1,2))))
   id  ::Array{Int,1} = Array{Int,1}(undef, 0)
   tr  ::DataFrame = tr_default
   ⎔  ::Function = dxy_dt
   ⎔! ::Function = dxy_dt!
   □   ::Function = solver_default
   𝑃   ::NamedTuple = 𝑃_default
   𝐷   ::NamedTuple = NamedTuple()
   𝑀  ::NamedTuple = NamedTuple()
end

end # module
