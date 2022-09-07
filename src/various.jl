

p=dirname(pathof(IndividualDisplacements))
artifact_toml = joinpath(p, "../Artifacts.toml")
flt_example_hash = artifact_hash("flt_example", artifact_toml)
flt_example_path = joinpath(artifact_path(flt_example_hash)*"/","flt_example-1.0/")
flt_example_download() = artifact"flt_example"

#needed within OrdinaryDiffEq somehow -- begin
#import Base: zero
#zero(tmp::Array) = zero.(tmp)

#import OrdinaryDiffEq.DiffEqBase: abs2_and_sum, UNITLESS_ABS2
#abs2_and_sum(x::Matrix{Vector{Float64}},y::Matrix{Vector{Float64}}) = 
#    reduce(Base.add_sum,x,init=0*x[1]) .+ reduce(Base.add_sum,y,init=0*x[1])
#@inline UNITLESS_ABS2(x::Matrix{Vector{Float64}}) = mapreduce(UNITLESS_ABS2,abs2_and_sum, x, init=0*x[1])
#needed within OrdinaryDiffEq somehow -- end