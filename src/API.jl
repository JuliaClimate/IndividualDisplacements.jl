
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
𝐹_Array2D(u0,u1,v0,v1,𝑇)
𝐹_Array3D(u0,u1,v0,v1,w0,w1,𝑇)
𝐹_MeshArray2D(u0,u1,v0,v1,𝑇,update_location!)
𝐹_MeshArray3D(u0,u1,v0,v1,w0,w1,𝑇,update_location!)
```

Using the `FlowFields` constructor which gets selected by the type of `u0` etc. For example :

```
𝐹=FlowFields(u,u,v,v,0*w,1*w,[0.0,10.0])
𝐹=FlowFields(u,u,v,v,[0.0,10.0],func)
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

function FlowFields(u0::Array{T,2},u1::Array{T,2},
    v0::Array{T,2},v1::Array{T,2},𝑇::Union{Array,Tuple}) where T
    #test for type of 𝑇 and fix if needed
    isa(𝑇,Tuple) ? 𝑇=convert(Array{T},[𝑇...]) : 𝑇=convert(Array{T},𝑇)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp)) for tmp in (u1,v0,v1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    𝐹_Array2D(u0,u1,v0,v1,𝑇)
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

function FlowFields(u0::Array{T,3},u1::Array{T,3},v0::Array{T,3},v1::Array{T,3},
    w0::Array{T,3},w1::Array{T,3},𝑇::Union{Array,Tuple}) where T
    #test for type of 𝑇 and fix if needed
    isa(𝑇,Tuple) ? 𝑇=convert(Array{T},[𝑇...]) : 𝑇=convert(Array{T},𝑇)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp)) for tmp in (u1,v0,v1)])
    tst=tst*prod([(size(u0)==size(tmp).-(0,0,1)) for tmp in (w0,w1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    𝐹_Array3D(u0,u1,v0,v1,w0,w1,𝑇)
end

struct 𝐹_MeshArray2D{T} <: FlowFields
    u0::AbstractMeshArray{T,1}
    u1::AbstractMeshArray{T,1}
    v0::AbstractMeshArray{T,1}
    v1::AbstractMeshArray{T,1}
    𝑇::Array{T}
    update_location!::Function
end

function FlowFields(u0::AbstractMeshArray{T,1},u1::AbstractMeshArray{T,1},
    v0::AbstractMeshArray{T,1},v1::AbstractMeshArray{T,1},
    𝑇::Union{Array,Tuple},update_location!::Function) where T
    #test for type of 𝑇 and fix if needed
    isa(𝑇,Tuple) ? 𝑇=convert(Array{T},[𝑇...]) : 𝑇=convert(Array{T},𝑇)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp))*(u0.fSize==tmp.fSize) for tmp in (u1,v0,v1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    𝐹_MeshArray2D(u0,u1,v0,v1,𝑇,update_location!)
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

function FlowFields(u0::AbstractMeshArray{T,2},u1::AbstractMeshArray{T,2},
    v0::AbstractMeshArray{T,2},v1::AbstractMeshArray{T,2},
    w0::AbstractMeshArray{T,2},w1::AbstractMeshArray{T,2},
    𝑇::Union{Array,Tuple},update_location!::Function) where T
    #test for type of 𝑇 and fix if needed
    isa(𝑇,Tuple) ? 𝑇=convert(Array{T},[𝑇...]) : 𝑇=convert(Array{T},𝑇)
    #check array size consistency
    tst=prod([(size(u0)==size(tmp))*(u0.fSize==tmp.fSize) for tmp in (u1,v0,v1)])
    tst=tst*prod([(size(u0)==size(tmp).-(0,1))*(u0.fSize==tmp.fSize) for tmp in (w0,w1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    𝐹_MeshArray3D(u0,u1,v0,v1,w0,w1,𝑇,update_location!)
end

"""
    defaults for Individuals constructor
"""

default_solver(prob) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

function ensemble_solver(prob;solver=Tsit5(),reltol=1e-8,abstol=1e-8)
	u0 = prob.u0
	prob_func(prob,i,repeat) = remake(prob,u0=u0[i])
	indiv_prob = ODEProblem(prob.f,u0[1],prob.tspan,prob.p)
	ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func)
	solve(ensemble_prob, solver, reltol=reltol, abstol=abstol, trajectories=length(u0))
end

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
calling ∫! which updates 📌 at the end and records results in 🔴 via 🔧. Ancillary data, for use in 
🔧 for example, can be provided in 𝐷 and metadata stored in 𝑀.

Unicode cheatsheet:

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
   🚄  ::Function = dxdt! #\:bullettrain_side:<tab>
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
    haskey(NT,:velocity) ? 🚄=NT.velocity : 🚄=dxdt!
    haskey(NT,:integration) ? ∫=NT.integration : ∫=default_solver
    haskey(NT,:postprocessing) ? 🔧=NT.postprocessing : 🔧=default_postproc
    haskey(NT,:parameters) ? 𝑃=NT.parameters : 𝑃=default_flowfields
    haskey(NT,:diagnostics) ? 𝐷=NT.diagnostics : 𝐷=NamedTuple()
    haskey(NT,:metadata) ? 𝑀=NT.metadata : 𝑀=NamedTuple()
    isa(📌,UnitRange) ? 📌=collect(📌) : nothing
    haskey(NT,:type) ? T=NT.type : T=eltype(📌)

    Individuals{T,ndims(📌)}(📌=📌,🔴=🔴,🆔=🆔,🚄=🚄,∫=∫,🔧=🔧,𝑃=𝑃,𝐷=𝐷,𝑀=𝑀)    
end

function Individuals(𝐹::𝐹_Array2D,x,y, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    🔧 = postprocess_xy
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    ∫=ensemble_solver
    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    𝐷=NamedTuple()
    haskey(NT,:𝐷) ? 𝐷=NT.𝐷 : nothing
    
    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,𝐷=𝐷)
end

function Individuals(𝐹::𝐹_Array3D,x,y,z, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];z[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    function 🔧(sol,𝐹::𝐹_Array3D,𝐷::NamedTuple;id=missing,𝑇=missing)
        df=postprocess_xy(sol,𝐹,𝐷,id=id,𝑇=𝑇)
        z=sol[3,:]
        df.z=z[:]
        return df
    end
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    ∫=ensemble_solver
    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    𝐷=NamedTuple()
    haskey(NT,:𝐷) ? 𝐷=NT.𝐷 : nothing
    
    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,𝐷=𝐷)
end

function Individuals(𝐹::𝐹_MeshArray2D,x,y,fid, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];fid[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], fid=Int64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    🔧 = postprocess_MeshArray
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    ∫=ensemble_solver
    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    𝐷=NamedTuple()
    haskey(NT,:𝐷) ? 𝐷=NT.𝐷 : nothing

    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,𝐷=𝐷)
end

function Individuals(𝐹::𝐹_MeshArray3D,x,y,z,fid, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];z[i];fid[i]] for i in eachindex(x)])
    length(📌)==1 ? 📌=📌[1] : nothing
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], fid=Int64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    function 🔧(sol,𝐹::𝐹_MeshArray3D,𝐷::NamedTuple;id=missing,𝑇=missing)
        df=postprocess_MeshArray(sol,𝐹,𝐷,id=id,𝑇=𝑇)
        z=[sol[1,i,j][1] for i in 1:size(sol,2), j in 1:size(sol,3)]
        df.z=z[:]
        return df
    end
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    ∫=ensemble_solver
    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    𝐷=NamedTuple()
    haskey(NT,:𝐷) ? 𝐷=NT.𝐷 : nothing

    Individuals{T,ndims(📌)}(𝑃=𝐹,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,𝐷=𝐷)
end

"""
    ∫!(𝐼::Individuals,𝑇::Tuple)

Displace simulated individuals continuously through space over time period 𝑇 starting from position 📌. 

- This is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).
- The current default is `solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)` but all solver options from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package are available.
- After this, `∫!` is also equipped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in `∫!` consists in updating 📌 to be ready for continuing in a subsequent call to `∫!`.
"""
function ∫!(𝐼::Individuals,𝑇::Tuple)
    @unpack 🚄,📌,𝑃, 𝐷, 🔧, 🆔, 🔴, ∫ = 𝐼

    prob = ODEProblem(🚄,📌, 𝑇 ,𝑃)
    sol = ∫(prob)

    tmp = 🔧(sol,𝑃,𝐷, id=🆔, 𝑇=𝑇)

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])

    if isa(sol,EnsembleSolution)
        np=size(sol,3)
        📌[:] = deepcopy([sol[i].u[end] for i in 1:np])
    else
        nd=length(size(sol))
        nd==3 ? 📌[:,:] = deepcopy(sol[:,:,end]) : 📌[:] = deepcopy(sol[:,end])
    end

end

∫!(𝐼::Individuals,𝑇::Array) = ∫!(𝐼::Individuals,(𝑇[1],𝑇[2]))

"""
    ∫!(𝐼::Individuals)

Call ∫!(𝐼::Individuals,𝐼.𝑃.𝑇)
"""
∫!(𝐼::Individuals) = ∫!(𝐼::Individuals,𝐼.𝑃.𝑇)

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
    N = ndims(𝐼.📌)
    return Individuals{T,N}(📌=similar(📌),🔴=similar(🔴),🆔=similar(🆔),
                          🚄=🚄, ∫=∫, 🔧=🔧, 𝑃=𝑃, 𝐷=𝐷, 𝑀=𝑀)
end

"""
    Base.diff(𝐼::Individuals)

Difference in grid unit coordinates (dx,dy) between final and initial positions.
"""
function Base.diff(𝐼::Individuals)
    f(x)=last(x).-first(x)
    🔴_by_ID = groupby(𝐼.🔴, :ID)
    return combine(🔴_by_ID,nrow,:x => f => :dx,:y => f => :dy)
end

"""
    gcdist(𝐼::Individuals)

Great circle distance (gcd in radians) between final and initial positions.
"""
function gcdist(𝐼::Individuals)
    🔴_by_ID = groupby(𝐼.🔴, :ID)
    tmp = combine(🔴_by_ID, 
    :lon => first => :lo1,:lon => last => :lo2,
    :lat => first => :la1,:lat => last => :la2)

    gcdist(lo1,lo2,la1,la2) = acos(sind(la1)*sind(la2)+cosd(la1)*cosd(la2)*cosd(lo1-lo2))
    tmp.gcd=[gcdist(tmp.lo1[i],tmp.lo2[i],tmp.la1[i],tmp.la2[i]) for i in 1:size(tmp,1)]
    return tmp
end

