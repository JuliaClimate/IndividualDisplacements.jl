
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/flow_fields.jl"));

"""
    init_global_range(lons::Tuple = (-160.0, -150.0),lats::Tuple = (35.0, 45.0))

Randomly distribute `np` points over a lon,la region, and 
return position in grid index space (`i,j,subdomain`).
"""
function init_global_range(lons::Tuple = (-160.0, -150.0),lats::Tuple = (35.0, 45.0))
    lo0, lo1 = lons #(-160.0, -150.0)
    la0, la1 = lats #(35.0, 45.0)
    np = 100
    lon = lo0 .+ (lo1 - lo0) .* rand(np)
    lat = la0 .+ (la1 - la0) .* rand(np)
    (u0, _) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
    id=collect(1:np)
    return u0
end

"""
    init_global_randn(np ::Int , 𝑃::NamedTuple)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(2*np)
    (u0, _) = initialize_lonlat(𝑃.Γ, lon, lat; msk = 𝑃.msk)
    u0[:,1:np]
end

"""
    reset_lonlat!(𝐼::Individuals)

Randomly select a fraction (𝐼.𝑃.frac) of the particles and reset their positions.
"""
function reset_lonlat!(𝐼::Individuals,𝐷::NamedTuple)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐷.frac*np))
    (lon, lat) = randn_lonlat(2*n_reset)
    (v0, _) = initialize_lonlat(𝐷.Γ, lon, lat; msk = 𝐷.msk)
    n_reset=min(n_reset,size(v0,2))
    k_reset = rand(1:np, n_reset)
    v0 = permutedims([v0[:,i] for i in 1:size(v0,2)])
    𝐼.📌[k_reset].=v0[1:n_reset]
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

##

"""
    isosurface(θ,T,z)

```
isosurface(𝐼.𝑃.θ0,15,Γ["RC"])
```    
"""
function isosurface(θ,T,z)
    d=NaN*similar(θ[:,1])
    nr=size(θ,2)
    for j=1:size(d,1)
        for k=1:nr-1
            i=findall(isnan.(d[j]).&(θ[j,k].>T).&(θ[j,k+1].<=T))
            a=(θ[j,k][i] .- T)./(θ[j,k][i] .- θ[j,k+1][i])
            d[j][i]=(1 .- a).*Γ["RC"][k] + a.*Γ["RC"][k+1]
            i=findall(isnan.(d[j]).&(θ[j,k].<=T).&(θ[j,k+1].>T))
            a=(θ[j,k+1][i] .- T)./(θ[j,k+1][i] .- θ[j,k][i])
            d[j][i]=(1 .- a).*Γ["RC"][k+1] + a.*Γ["RC"][k]
        end
    end
    return d
end
