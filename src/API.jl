
day=86400.0
mon=365/12*day
solver_default(prob) = solve(prob,Euler(),dt=2*day)
𝑃_default = ( 𝑇 = [-0.5*mon,0.5*mon] , 🔄 = update_𝑃!,
              u0=[] , u1=[] , v0=[] , v1=[] )
rec_default = DataFrame(fill(Float64, 7),[:ID, :x, :y, :t, :lon, :lat, :fid])
postprocess_default = postprocess_lonlat

"""
    struct Individuals{T}

- Data:           📌 (position),   🔴(record),           🆔 (ID)
- Functions:      🚄 (velocity),   ∫ (integration), 🔧(postprocessing)
- NamedTuples:    𝑃  (parameters), 𝐷 (diagnostics),      𝑀 (metadata)

Keyword constructor -- for example:

```
𝐼=Individuals{Float64}(📌=zeros(3,2),🆔=1:2,🔴=DataFrame( ID=[], x=[], y=[], z=[], t = []))
```

Or alternatively, without unicode:

```
𝐼=Individuals(position=zeros(3,2),ID=1:2,record=DataFrame( ID=[], x=[], y=[], z=[], t = []))
```

Keyword cheatsheet:

- 📌=`\\:pushpin:<tab>`,        🔴=`\\:red_circle:<tab>`, 🆔=`\\:id:<tab>`
- 🚄=`\\bullettrain_side<tab>`, ∫=`\\int<tab>`,           🔧=`\\wrench<tab>`
- 𝑃=`\\itP<tab>`,               𝐷=`\\itD<tab>`,            𝑀 =`\\itM<tab>`
"""
Base.@kwdef struct Individuals{T}
   📌  ::Array{T,2} = Array{T,2}(undef, Tuple(Int.(zeros(1,2)))) #\:pushpin:<tab>
   🔴  ::DataFrame = rec_default #\:red_circle:<tab>
   🆔   ::Array{Int,1} = Array{Int,1}(undef, 0) #\:id:<tab>
   🚄  ::Function = dxy_dt #\bullettrain_side<tab>
   ∫   ::Function = solver_default #\int<tab>
   🔧  ::Function = postprocess_default #\wrench<tab>
   𝑃   ::NamedTuple = 𝑃_default #\itP<tab>
   𝐷   ::NamedTuple = NamedTuple() #\itD<tab>
   𝑀   ::NamedTuple = NamedTuple() #\itM<tab>
end

# constructor that uses plain text keywords:
function Individuals(;
    position::Array{T,2} = Array{T,2}(undef, Tuple(Int.(zeros(1,2)))),
    record::DataFrame = rec_default,
    ID::Union{Array{Int,1},UnitRange{Int}} = Array{Int,1}(undef, 0),
    velocity::Function = dxy_dt,
    integration::Function = solver_default,
    postprocessing::Function = postprocess_default,
    parameters::NamedTuple = 𝑃_default,
    diagnostics::NamedTuple = NamedTuple(),
    metadata::NamedTuple = NamedTuple()
    ) where T
    Individuals{T}(📌=position,🔴=record,🆔=ID,
    🚄=velocity,   ∫=integration, 🔧=postprocessing,
    𝑃=parameters, 𝐷=diagnostics, 𝑀=metadata)    
end

"""
    ∫!(𝐼::Individuals,𝑇::Tuple)

Displace simulated individuals continuously through space over time period 𝑇 starting from position 📌. 

- This is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).
- The current default is `solve(prob,Euler(),dt=2*day)` but all solver options from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package are available.
- After this, `∫!` is also equiped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in `∫!` consiste in updating 📌 to be ready for continuing in a subsequent call to `∫!`.
"""
function ∫!(𝐼::Individuals,𝑇::Tuple)
    @unpack 🚄,📌,𝑃, 🔧, 🆔, 🔴, ∫ = 𝐼

    prob = ODEProblem(🚄,📌, 𝑇 ,𝑃)
    sol = ∫(prob)

    tmp = 🔧(sol,𝑃, id=🆔, 𝑇=𝑇)

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])

    📌[:,:] = deepcopy(sol[:,:,end])
end

"""
    reset_lonlat!(𝐼::Individuals)

Randomly select a fraction (𝐼.𝑃.frac) of the particles and reset their positions.
"""
function reset_lonlat!(𝐼::Individuals)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐼.𝑃.frac*np))
    (lon, lat) = randn_lonlat(2*n_reset)
    (v0, _) = initialize_lonlat(𝐼.𝑃.Γ, lon, lat; msk = 𝐼.𝑃.msk)
    k_reset = rand(1:np, n_reset)
    𝐼.📌[:,k_reset].=v0[:,1:n_reset]
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

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
    printstyled(io, "$(keys(𝑃))\n",color=:blue)
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

