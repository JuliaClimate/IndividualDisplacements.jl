
day=86400.0
mon=365/12*day
solver_default(prob) = solve(prob,Euler(),dt=2*day)
𝑃_default = ( 𝑇 = [-0.5*mon,0.5*mon] , 🔄 = update_𝑃!,
              u0=[] , u1=[] , v0=[] , v1=[] )
tr_default = DataFrame( ID=[], x=[], y=[], t = [], lon=[], lat=[], fid=[])
postprocess_default = postprocess_lonlat

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
   ▽   ::Function = postprocess_default
   𝑃   ::NamedTuple = 𝑃_default
   𝐷   ::NamedTuple = NamedTuple()
   𝑀  ::NamedTuple = NamedTuple()
end

"""
    start!(𝐼::Individuals)

Set up ODE problem over `(0.0,𝐼.𝑃.𝑇[2])`, solve, postprocess, & update `𝐼.xy[:,:]`
"""
function start!(𝐼::Individuals)
    prob = ODEProblem(𝐼.⎔!,𝐼.xy,(0.0,𝐼.𝑃.𝑇[2]),𝐼.𝑃)
    sol = 𝐼.□(prob)
    tmp = 𝐼.▽(sol,𝐼.𝑃,𝐼.id)
    append!(𝐼.tr,tmp)
    𝐼.xy[:,:] = deepcopy(sol[:,:,end])
end

"""
    displace!(𝐼::Individuals)

Update 𝐼.𝑃, set up ODE problem over 𝐼.𝑃.𝑇, solve, postprocess, & update `𝐼.xy[:,:]`
"""
function displace!(𝐼::Individuals)
    𝐼.𝑃.🔄(𝐼.𝑃.k,𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2]),𝐼.𝑃)
    prob = ODEProblem(𝐼.⎔!,𝐼.xy,𝐼.𝑃.𝑇,𝐼.𝑃)
    sol = 𝐼.□(prob)
    tmp = 𝐼.▽(sol,𝐼.𝑃,𝐼.id)
    np=length(𝐼.id)
    append!(𝐼.tr,tmp[np+1:end,:])
    𝐼.xy[:,:] = deepcopy(sol[:,:,end])
end

"""
    reset!(𝐼::Individuals)

Randomly select a fraction (𝐼.𝑃.frac) of the particles and reset their positions.
"""
function reset!(𝐼::Individuals)
    np=length(𝐼.id)
    n_reset = Int(round(𝐼.𝑃.frac*np))
    (lon, lat) = randn_lonlat(2*n_reset)
    (v0, _) = initialize_lonlat(𝐼.𝑃.Γ, lon, lat; msk = 𝐼.𝑃.msk)
    k_reset = rand(1:np, n_reset)
    𝐼.xy[:,k_reset].=v0[:,1:n_reset]
    isempty(𝐼.tr.ID) ? m=maximum(𝐼.id) : m=max(maximum(𝐼.tr.ID),maximum(𝐼.id))
    𝐼.id[k_reset]=collect(1:n_reset) .+ m
end
