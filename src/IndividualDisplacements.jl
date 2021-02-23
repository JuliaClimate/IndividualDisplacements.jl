module IndividualDisplacements

using MeshArrays, OrdinaryDiffEq, DataFrames
using NetCDF, Dates, CFTime, CSV, UnPack, Random, Pkg.Artifacts
using CyclicArrays, MITgcmTools, OceanStateEstimation

p=dirname(pathof(IndividualDisplacements))
artifact_toml = joinpath(p, "../Artifacts.toml")
flt_example_hash = artifact_hash("flt_example", artifact_toml)
flt_example = joinpath(artifact_path(flt_example_hash)*"/","flt_example-1.0/")

include("API.jl")
include("compute.jl")
include("data_wrangling.jl")
include("read.jl")
include("update_locations.jl")

export Individuals, ∫!
export FlowFields, convert_to_FlowFields
export 𝐹_Array3D, 𝐹_Array2D, 𝐹_MeshArray3D, 𝐹_MeshArray2D
export dxy_dt!, dxy_dt, dxyz_dt!, dxyz_dt, dxy_dt_CyclicArray, dxy_dt_replay
export postprocess_MeshArray, add_lonlat!, postprocess_xy, interp_to_xy
export nearest_to_xy, randn_lonlat, interp_to_lonlat
export gcdist, stproj, stproj_inv
export read_drifters, read_mds, read_velocities

end # module
