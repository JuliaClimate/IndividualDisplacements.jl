
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/jupyter/flow_fields.jl"));

"""
    init_global_randn(np ::Int , 𝑃::NamedTuple)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(maximum([2*np 10]))
    (_,_,_,_,f,x,y)=InterpolationFactors(𝐷.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(nearest_to_xy(𝐷.msk,x[m],y[m],f[m]).==1.0)[1:np]
    return permutedims([x[m[n]] y[m[n]] f[m[n]]])
end

"""
    reset_📌!(𝐼::Individuals,frac::Number,📌::Array)

Randomly select a fraction (frac) of the particles and reset 
their positions (𝐼.📌) to a random subset of the specificed 📌.
"""
function reset_📌!(𝐼::Individuals,frac::Number,📌::Array)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐷.frac*np))
    k_reset = rand(1:np, n_reset)
    l_reset = rand(1:np, n_reset)
    𝐼.📌[k_reset]=deepcopy(📌[l_reset])
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

##

"""
    isosurface(θ,T,z)

```
isosurface(𝐼.𝑃.θ0,15,Γ.RC)
```    
"""
function isosurface(θ,T,z)
    d=NaN*similar(θ[:,1])
    nr=size(θ,2)
    for j=1:size(d,1)
        for k=1:nr-1
            i=findall(isnan.(d[j]).&(θ[j,k].>T).&(θ[j,k+1].<=T))
            a=(θ[j,k][i] .- T)./(θ[j,k][i] .- θ[j,k+1][i])
            d[j][i]=(1 .- a).*Γ.RC[k] + a.*Γ.RC[k+1]
            i=findall(isnan.(d[j]).&(θ[j,k].<=T).&(θ[j,k+1].>T))
            a=(θ[j,k+1][i] .- T)./(θ[j,k+1][i] .- θ[j,k][i])
            d[j][i]=(1 .- a).*Γ.RC[k+1] + a.*Γ.RC[k]
        end
    end
    return d
end

"""
    set_up_𝑃(k::Int,t::Float64,Γ::NamedTuple,pth::String)

Define `FlowFields` data structure (𝑃) along with ancillary variables (𝐷)
for the specified grid (`Γ` dictionnary), vertical level (`k`), and 
file location (`pth`).
    
_Note: the initial implementation approximates month durations to 
365 days / 12 months for simplicity and sets 𝑃.𝑇 to [-mon/2,mon/2]_
"""
function set_up_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String)
    XC=exchange(Γ.XC) #add 1 lon point at each edge
    YC=exchange(Γ.YC) #add 1 lat point at each edge
    iDXC=1. ./Γ.DXC
    iDYC=1. ./Γ.DYC
    γ=Γ.XC.grid
    mon=86400.0*365.0/12.0
    
    if k==0
        msk=Γ.hFacC
        (_,nr)=size(msk)
        𝑃=FlowFields(MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr+1),MeshArray(γ,Float32,nr+1),
        [-mon/2,mon/2],func)
    else
        msk=Γ.hFacC[:, k]
        𝑃=FlowFields(MeshArray(γ,Float32),MeshArray(γ,Float32),
        MeshArray(γ,Float32),MeshArray(γ,Float32),[-mon/2,mon/2],func)    
    end
    
    𝐷 = (🔄 = update_FlowFields!, pth=pth,
         XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC,
         k=k, msk=msk, θ0=similar(msk), θ1=similar(msk))

    𝐷 = merge(𝐷 , MeshArrays.NeighborTileIndices_cs(Γ))

    return 𝑃,𝐷
end

"""
    update_FlowFields!(𝑃::𝐹_MeshArray2D,𝐷::NamedTuple,t::Float64)

Update flow field arrays (in 𝑃), 𝑃.𝑇, and ancillary variables (in 𝐷) 
according to the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(𝑃.𝑇), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical 𝑃.k is selected_
"""
function update_FlowFields!(𝑃::𝐹_MeshArray2D,𝐷::NamedTuple,t::AbstractFloat)
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

"""
    update_FlowFields!(𝑃::𝐹_MeshArray3D,𝐷::NamedTuple,t::Float64)

Update flow field arrays (in 𝑃), 𝑃.𝑇, and ancillary variables (in 𝐷) 
according to the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(𝑃.𝑇), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical 𝑃.k is selected_
"""
function update_FlowFields!(𝑃::𝐹_MeshArray3D,𝐷::NamedTuple,t::Float64)
    dt=𝑃.𝑇[2]-𝑃.𝑇[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    (_,nr)=size(𝐷.Γ.hFacC)

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
    u0=U; v0=V
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u0[:,k]=u0[:,k].*𝐷.iDXC; v0[:,k]=v0[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u0[:,k],v0[:,k],1) #add 1 point at each edge for u and v
        u0[:,k]=tmpu
        v0[:,k]=tmpv
    end
    w0=IndividualDisplacements.read_nctiles(𝐷.pth*"WVELMASS/WVELMASS","WVELMASS",𝑃.u0.grid,I=(:,:,:,m0))
    w0[findall(isnan.(w0))]=0.0 #mask with 0s rather than NaNs

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=U; v1=V
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u1[:,k]=u1[:,k].*𝐷.iDXC; v1[:,k]=v1[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u1[:,k],v1[:,k],1) #add 1 point at each edge for u and v
        u1[:,k]=tmpu
        v1[:,k]=tmpv
    end
    w1=IndividualDisplacements.read_nctiles(𝐷.pth*"WVELMASS/WVELMASS","WVELMASS",𝑃.u0.grid,I=(:,:,:,m1))
    w1[findall(isnan.(w1))]=0.0 #mask with 0s rather than NaNs

    𝑃.u0[:,:]=u0[:,:]
    𝑃.u1[:,:]=u1[:,:]
    𝑃.v0[:,:]=v0[:,:]
    𝑃.v1[:,:]=v1[:,:]
    for k=1:nr
        tmpw=exchange(-w0[:,k],1)
        𝑃.w0[:,k]=tmpw./𝐷.Γ.DRC[k]
        tmpw=exchange(-w1[:,k],1)
        𝑃.w1[:,k]=tmpw./𝐷.Γ.DRC[k]
    end
    𝑃.w0[:,1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,1]=0*exchange(-w1[:,1],1)
    𝑃.w0[:,nr+1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,nr+1]=0*exchange(-w1[:,1],1)

    θ0=IndividualDisplacements.read_nctiles(𝐷.pth*"THETA/THETA","THETA",𝑃.u0.grid,I=(:,:,:,m0))
    θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ0[:,:]=float32.(θ0[:,:])

    θ1=IndividualDisplacements.read_nctiles(𝐷.pth*"THETA/THETA","THETA",𝑃.u0.grid,I=(:,:,:,m1))
    θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ1[:,:]=float32.(θ1[:,:])

    𝑃.𝑇[:]=[t0,t1]
end
