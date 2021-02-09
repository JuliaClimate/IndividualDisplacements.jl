
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/flow_fields.jl"));

"""
    init_global_randn(np ::Int , 𝑃::NamedTuple)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(maximum([2*np 10]))
    (_,_,_,_,f,x,y)=InterpolationFactors(𝑃.Γ,lon,lat)
    m=findall(f.!==0)
    m=findall(nearest_to_xy(𝑃.msk,x[m],y[m],f[m]).==1.0)[1:np]
    return permutedims([x[m] y[m] f[m]])
end

"""
    reset_lonlat!(𝐼::Individuals)

Randomly select a fraction (𝐼.𝑃.frac) of the particles and reset their positions.
"""
function reset_lonlat!(𝐼::Individuals,𝐷::NamedTuple)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐷.frac*np))
    v0=init_global_randn(n_reset , 𝐷)
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

"""
    set_up_𝑃(k::Int,t::Float64,Γ::Dict,pth::String)

Define the `𝑃` _parameter_ tuple given grid variables `Γ`, vertical level
choice `k`, time `t` in `seconds`, and velocity fields obtained from
files in `pth`.

The two climatological months (`m0`,`m1`) that bracket time `t` are
read to memory (e.g. months 12 & 1 then 1 & 2 and so on).

_Note: the initial implementation approximates every month duration
to 365 days / 12 months for simplicity._
"""
function set_up_𝑃(k::Int,t::Float64,Γ::Dict,pth::String)
    XC=exchange(Γ["XC"]) #add 1 lon point at each edge
    YC=exchange(Γ["YC"]) #add 1 lat point at each edge
    iDXC=1. ./Γ["DXC"]
    iDYC=1. ./Γ["DYC"]
    γ=Γ["XC"].grid
    mon=86400.0*365.0/12.0
    func=Γ["update_location!"]
    
    𝐷 = (🔄 = update_𝑃!, pth=pth,
         XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC,
         k=k, msk=Γ["hFacC"][:, k])

    tmp = IndividualDisplacements.dict_to_nt(IndividualDisplacements.NeighborTileIndices_cs(Γ))
    𝐷 = merge(𝐷 , tmp)

    𝑃=𝐹_MeshArray2D{Float64}(MeshArray(γ,Float64),MeshArray(γ,Float64),
    MeshArray(γ,Float64),MeshArray(γ,Float64),[-mon/2,mon/2],func)

    𝐷.🔄(𝑃,𝐷,0.0)

    return 𝑃,𝐷
end

"""
    update_𝑃!(𝑃::Union{NamedTuple,FlowFields},t::Float64)

Update input data (velocity arrays) and time period (array) inside 𝑃 (𝑃.u0[:], etc, and 𝑃.𝑇[:])
based on the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) input 𝑃.𝑇 is used to infer `dt` between consecutive velocity fields,
(2) periodicity of 12 monthly records, (3) vertical 𝑃.k is selected -- but this could easily be generalized._ 
"""
function update_𝑃!(𝑃::Union{NamedTuple,FlowFields},𝐷::NamedTuple,t::Float64)
    dt=𝑃.𝑇[2]-𝑃.𝑇[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
    u0=U[:,𝐷.k]; v0=V[:,𝐷.k]
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    u0=u0.*𝐷.iDXC; v0=v0.*𝐷.iDYC; #normalize to grid units
    (u0,v0)=exchange(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=U[:,𝐷.k]; v1=V[:,𝐷.k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1.*𝐷.iDXC; v1=v1.*𝐷.iDYC; #normalize to grid units
    (u1,v1)=exchange(u1,v1,1) #add 1 point at each edge for u and v

    𝑃.u0[:]=u0[:]
    𝑃.u1[:]=u1[:]
    𝑃.v0[:]=v0[:]
    𝑃.v1[:]=v1[:]
    𝑃.𝑇[:]=[t0,t1]

end

