
## Flow field parameters

"""
    abstract type FlowFields

Data structure that provide access to flow fields (gridded as arrays) which will be 
used to interpolate velocities to individual locations later on (once embedded in
an `Individuals` struct). 

Following the C-grid convention also used in `MITgcm` (https://mitgcm.readthedocs.io) 
flow fields are expected to be staggered as follows: grid cell i,j has its center located at i-1/2,j-1/2 while the
corresponding `u[i,j]` (resp. `v[i,j]) is located at i-1,j-1/2 (resp. i-1/2,j-1). 

Also by convention, velocity fields are expected to have been normalized to grid units (e.g. 1/s rather than m/s)
before sending them to one of the supported `FlowFields` constructors (using either `Array` or `MeshArray`):

```
𝐹_Array2D (u0,u1,v0,v1,𝑇)
𝐹_Array3D (u0,u1,v0,v1,w0,w1,𝑇)
𝐹_MeshArray2D (u0,u1,v0,v1,𝑇,update_location!)
𝐹_MeshArray3D (u0,u1,v0,v1,w0,w1,𝑇,update_location!)
```

For example, constructor calls may look like

```
𝐹=𝐹_Array3D{eltype(u)}(u,u,v,v,0*w,1*w,[0.0,10.0])
or
𝐹=𝐹_MeshArray2D{eltype(u)}(u,u,v,v,[0.0,10.0],func)
```

as shown in the online documentation examples.

"""
abstract type FlowFields end

struct 𝐹_Array2D{T} <: FlowFields
    u0::Array{T,2}
    u1::Array{T,2}
    v0::Array{T,2}
    v1::Array{T,2}
    𝑇::Array{T}
end

struct 𝐹_Array3D{T} <: FlowFields
    u0::Array{T,3}
    u1::Array{T,3}
    v0::Array{T,3}
    v1::Array{T,3}
    w0::Array{T,3}
    w1::Array{T,3}
    𝑇::Array{T}
end

struct 𝐹_MeshArray2D{T} <: FlowFields
    u0::AbstractMeshArray{T,1}
    u1::AbstractMeshArray{T,1}
    v0::AbstractMeshArray{T,1}
    v1::AbstractMeshArray{T,1}
    𝑇::Array{T}
    update_location!::Function
end

struct 𝐹_MeshArray3D{T} <: FlowFields
    u0::AbstractMeshArray{T,2}
    u1::AbstractMeshArray{T,2}
    v0::AbstractMeshArray{T,2}
    v1::AbstractMeshArray{T,2}
    w0::AbstractMeshArray{T,2}
    w1::AbstractMeshArray{T,2}
    𝑇::Array{T}
    update_location!::Function
end


"""
    defaults for Individuals constructor
"""

default_solver(prob) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
a=fill(0.0,1,1)
default_flowfields = 𝐹_Array2D{Float64}(a,a,a,a,[0. 1.])
default_recorder = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
default_postproc = (x->x)

"""
    struct Individuals{T,N}

- Data:           📌 (position),   🔴(record), 🆔 (ID), 𝑃 (`FlowFields`)
- Functions:      🚄 (velocity),   ∫ (integration), 🔧(postprocessing)
- NamedTuples:    𝐷 (diagnostics),      𝑀 (metadata)

The velocity function 🚄 typically computes velocity at individual positions (📌 to start) within the 
specified space-time domain by interpolating gridded variables (provided via 𝑃). Individual trajectories 
are computed by integrating (∫) interpolated velocities through time. Normally, integration is done by 
calling ∫! which updates 📌 at the end and records results in 🔴 via 🔧. Unicode cheatsheet:

- 📌=`\\:pushpin:<tab>`,          🔴=`\\:red_circle:<tab>`, 🆔=`\\:id:<tab>`
- 🚄=`\\:bullettrain_side:<tab>`, ∫=`\\int<tab>`,          🔧=`\\:wrench:<tab>`
- 𝑃=`\\itP<tab>`,                 𝐷=`\\itD<tab>`,           𝑀=`\\itM<tab>`

Simple constructors that use `FlowFields` to choose adequate defaults:

- Individuals(𝐹::𝐹_Array2D,x,y)
- Individuals(𝐹::𝐹_Array3D,x,y,z)
- Individuals(𝐹::𝐹_MeshArray2D,x,y,fid)
- Individuals(𝐹::𝐹_MeshArray3D,x,y,z,fid)

Further customization is achievable via keyword constructors:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
𝐼=Individuals{Float64,2}(📌=zeros(3,10),🆔=1:10,🔴=deepcopy(df))
𝐼=Individuals(📌=zeros(3,2),🆔=collect(1:2),🔴=deepcopy(df))
```

Or via the plain text (or no-unicode) constructors:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
I=(position=zeros(3,2),ID=1:2,record=deepcopy(df))
I=Individuals(I)
```
"""
Base.@kwdef struct Individuals{T,N}
   📌  ::Array{T,N} = Array{T,N}(undef, Tuple(Int.(zeros(1,N)))) #\:pushpin:<tab>
   🔴  ::DataFrame = similar(default_recorder) #\:red_circle:<tab>
   🆔   ::Array{Int,1} = Array{Int,1}(undef, 0) #\:id:<tab>
   🚄  ::Function = dxy_dt #\:bullettrain_side:<tab>
   ∫   ::Function = default_solver #\int<tab>
   🔧  ::Function = default_postproc #\:wrench:<tab>
   𝑃   ::FlowFields = default_flowfields #\itP<tab>
   𝐷   ::NamedTuple = NamedTuple() #\itD<tab>
   𝑀   ::NamedTuple = NamedTuple() #\itM<tab>
end

function Individuals(NT::NamedTuple)

    haskey(NT,:position) ? 📌=NT.position : 📌=Array{Float64,2}(undef, Tuple(Int.(zeros(1,2))))
    haskey(NT,:record) ? 🔴=NT.record : 🔴=similar(default_recorder)
    haskey(NT,:ID) ? 🆔=NT.ID : 🆔=collect(1:size(📌,2))    
    haskey(NT,:velocity) ? 🚄=NT.velocity : 🚄=dxy_dt
    haskey(NT,:integration) ? ∫=NT.integration : ∫=default_solver
    haskey(NT,:postprocessing) ? 🔧=NT.postprocessing : 🔧=default_postproc
    haskey(NT,:parameters) ? 𝑃=NT.parameters : 𝑃=default_flowfields
    haskey(NT,:diagnostics) ? 𝐷=NT.diagnostics : 𝐷=NamedTuple()
    haskey(NT,:metadata) ? 𝑀=NT.metadata : 𝑀=NamedTuple()
    isa(📌,UnitRange) ? 📌=collect(📌) : nothing
    haskey(NT,:type) ? T=NT.type : T=eltype(📌)

    Individuals{T,ndims(📌)}(📌=📌,🔴=🔴,🆔=🆔,🚄=🚄,∫=∫,🔧=🔧,𝑃=𝑃,𝐷=𝐷,𝑀=𝑀)    
end

function Individuals(𝐹::𝐹_Array2D,x,y)
    📌=permutedims([[x[i];y[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
    🔧 = postprocess_xy
    T=eltype(📌)
    🆔=collect(1:size(📌,2))
    
    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxy_dt,∫=default_solver,🔧=🔧)    
end

function Individuals(𝐹::𝐹_Array3D,x,y,z)
    📌=permutedims([[x[i];y[i];z[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], t=Float64[])
    function 🔧(sol,𝑄::FlowFields;id=missing,𝑇=missing)
        df=postprocess_xy(sol,𝐹,id=id,𝑇=𝑇)
        z=sol[3,:]
        df.z=z[:]
        return df
    end
    T=eltype(📌)
    🆔=collect(1:size(📌,2))
    
    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxyz_dt,∫=default_solver,🔧=🔧)    
end

function Individuals(𝐹::𝐹_MeshArray2D,x,y,fid)
    📌=permutedims([[x[i];y[i];fid[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], fid=Int64[], t=Float64[])
    🔧 = postprocess_MeshArray

    T=eltype(📌)
    🆔=collect(1:size(📌,2))

    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxy_dt!,∫=default_solver,🔧=🔧)    
end

function Individuals(𝐹::𝐹_MeshArray3D,x,y,fid)
    📌=permutedims([[x[i];y[i];fid[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], fid=Int64[], t=Float64[])
    function 🔧(sol,𝑄::FlowFields;id=missing,𝑇=missing)
        df=postprocess_MeshArray(sol,𝐹,id=id,𝑇=𝑇)
        z=sol[3,:]
        df.z=z[:]
        return df
    end

    T=eltype(📌)
    🆔=collect(1:size(📌,2))
    ∫=default_solver

    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxyz_dt!,∫=∫,🔧=🔧)    
end

"""
    ∫!(𝐼::Individuals,𝑇::Tuple)

Displace simulated individuals continuously through space over time period 𝑇 starting from position 📌. 

- This is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).
- The current default is `solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)` but all solver options from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package are available.
- After this, `∫!` is also equipped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in `∫!` consists in updating 📌 to be ready for continuing in a subsequent call to `∫!`.
"""
function ∫!(𝐼::Individuals,𝑇::Tuple)
    @unpack 🚄,📌,𝑃, 🔧, 🆔, 🔴, ∫ = 𝐼

    prob = ODEProblem(🚄,📌, 𝑇 ,𝑃)
    sol = ∫(prob)

    tmp = 🔧(sol,𝑃, id=🆔, 𝑇=𝑇)

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])

    nd=length(size(sol))
    nd==3 ? 📌[:,:] = deepcopy(sol[:,:,end]) : 📌[:] = deepcopy(sol[:,end])

end

"""
    ∫!(𝐼::Individuals)

Call ∫!(𝐼::Individuals,𝐼.𝑃.𝑇)
"""
∫!(𝐼::Individuals) = ∫!(𝐼::Individuals,(𝐼.𝑃.𝑇[1],𝐼.𝑃.𝑇[2]))

## Convenience Methods (size,show,similar)

Base.size(A::Individuals) = size(A.📌)

function Base.show(io::IO, 𝐼::Individuals) where {T}
    @unpack 🚄,📌,𝑃, 𝐷, 𝑀, 🔧, 🆔, 🔴, ∫ = 𝐼
    printstyled(io, "  📌 details     = ",color=:normal)
    printstyled(io, "$(size(📌)) $(typeof(𝐼).parameters[1])\n",color=:blue)
    printstyled(io, "  🔴 details     = ",color=:normal)
    printstyled(io, "$(size(🔴)) $(names(🔴))\n",color=:blue)
    printstyled(io, "  🆔 range       = ",color=:normal)
    printstyled(io, "$(extrema(🆔))\n",color=:blue)
    printstyled(io, "  🚄 function    = ",color=:normal)
    printstyled(io, "$(🚄)\n",color=:blue)
    printstyled(io, "  ∫  function    = ",color=:normal)
    printstyled(io, "$(∫)\n",color=:blue)
    printstyled(io, "  🔧 function    = ",color=:normal)
    printstyled(io, "$(🔧)\n",color=:blue)
    printstyled(io, "  𝑃  details     = ",color=:normal)
    printstyled(io, "$(fieldnames(typeof(𝑃)))\n",color=:blue)
  return
end

function Base.similar(𝐼::Individuals)
    @unpack 🚄,📌,𝑃, 𝐷, 𝑀, 🔧, 🆔, 🔴, ∫ = 𝐼
    T = typeof(𝐼).parameters[1]
    return Individuals{T}(📌=similar(📌),🔴=similar(🔴),🆔=similar(🆔),
                          🚄=🚄, ∫=∫, 🔧=🔧, 𝑃=𝑃, 𝐷=𝐷, 𝑀=𝑀)
end

function Base.diff(𝐼::Individuals)
    f(x)=last(x).-first(x)
    🔴_by_ID = groupby(𝐼.🔴, :ID)
    return combine(🔴_by_ID,nrow,:lat => f => :dlat,:lon => f => :dlon)
end
