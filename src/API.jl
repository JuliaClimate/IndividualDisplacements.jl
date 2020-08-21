
day=86400.0
mon=365/12*day
solver_default(prob) = solve(prob,Euler(),dt=2*day)
𝑃_default = ( 𝑇 = [-0.5*mon,0.5*mon] , 🔄 = update_𝑃!,
              u0=[] , u1=[] , v0=[] , v1=[] )
rec_default = DataFrame( ID=[], x=[], y=[], t = [], lon=[], lat=[], fid=[])
postprocess_default = postprocess_lonlat

"""
    struct Individuals{T}

Contains: 📌 (position), 🔴 (recording), 🆔 (ID), etc
```
i=Individuals{Float32}(📌=zeros(3,2),🆔=1:2)
```
"""
Base.@kwdef struct Individuals{T}
    📌  ::Array{T,2} = Array{T,2}(undef, Tuple(Int.(zeros(1,2)))) #\:pushpin:<tab>
   🔴  ::DataFrame = rec_default #\:red_circle:<tab>
   🆔  ::Array{Int,1} = Array{Int,1}(undef, 0) #\:id:<tab>
   ⎔   ::Function = dxy_dt #\hexagon<tab>
   ∫   ::Function = solver_default #\int<tab>
   ⟁   ::Function = postprocess_default #\whiteinwhitetriangle<tab>
   𝑃   ::NamedTuple = 𝑃_default #\itP<tab>
   𝐷   ::NamedTuple = NamedTuple() #\itD<tab>
   𝑀  ::NamedTuple = NamedTuple() #\itM<tab>
end

#alternative symbol choices?
#⏩  ::Function = dxy_dt #\:fast_forward:<tab>
#🔧  ::Function = postprocess_default #\:wrench:<tab>

"""
    ∫!(𝐼::Individuals,𝑇::Tuple)

Continuously displace individuals (∫! of ⎔), starting from 📌, over time period 𝑇. Then postprocess with ⟁, record data into 🔴, & update 📌
"""
function ∫!(𝐼::Individuals,𝑇::Tuple)
    @unpack ⎔,📌,𝑃, ⟁, 🆔, 🔴, ∫ = 𝐼

    prob = ODEProblem(⎔,📌, 𝑇 ,𝑃)
    sol = ∫(prob)

    tmp = ⟁(sol,𝑃, id=🆔, 𝑇=𝑇)

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])

    📌[:,:] = deepcopy(sol[:,:,end])
end

"""
    start!(𝐼::Individuals)

Set up ODE problem over `(0.0,𝐼.𝑃.𝑇[2])`, solve, postprocess, & update `𝐼.📌[:,:]`
"""
function start!(𝐼::Individuals)
    𝑇=(0.0,𝐼.𝑃.𝑇[2])
    prob = ODEProblem(𝐼.⎔,𝐼.📌, 𝑇 ,𝐼.𝑃)
    sol = 𝐼.∫(prob)
    tmp = 𝐼.⟁(sol,𝐼.𝑃, id=𝐼.🆔, 𝑇=𝑇)
    #tmp.t=0.0 .+ 𝐼.𝑃.𝑇[2] / diff(𝐼.𝑃.𝑇)[1] * tmp.t
    append!(𝐼.🔴,tmp)
    𝐼.📌[:,:] = deepcopy(sol[:,:,end])
end

"""
    displace!(𝐼::Individuals)

Set up ODE problem over 𝐼.𝑃.𝑇, solve, postprocess, & update `𝐼.📌[:,:]`
"""
function displace!(𝐼::Individuals)
    prob = ODEProblem(𝐼.⎔,𝐼.📌,𝐼.𝑃.𝑇,𝐼.𝑃)
    sol = 𝐼.∫(prob)
    tmp = 𝐼.⟁(sol,𝐼.𝑃,id=𝐼.🆔)
    np=length(𝐼.🆔)
    append!(𝐼.🔴,tmp[np+1:end,:])
    𝐼.📌[:,:] = deepcopy(sol[:,:,end])
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
