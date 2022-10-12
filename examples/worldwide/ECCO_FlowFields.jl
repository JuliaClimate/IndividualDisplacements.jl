module ECCO_FlowFields

using IndividualDisplacements, OceanStateEstimation, MITgcmTools

import IndividualDisplacements.DataFrames: DataFrame
import IndividualDisplacements.MeshArrays as MeshArrays
import IndividualDisplacements.MeshArrays: gcmgrid, MeshArray
import IndividualDisplacements.CSV as CSV

import OceanStateEstimation.ECCO_helpers.JLD2 as JLD2

"""
    init_from_file(np ::Int)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_from_file(np ::Int)
    p=dirname(pathof(IndividualDisplacements))
    fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    return DataFrame(CSV.File(fil))[1:np,:]
end

"""
    init_global_randn(np ::Int , 𝐷::NamedTuple)

Randomly distribute `np` points over the Earth, within `𝐷.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , 𝐷::NamedTuple)
    (lon, lat) = randn_lonlat(maximum([2*np 10]))
    (_,_,_,_,f,x,y)=InterpolationFactors(𝐷.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(nearest_to_xy(𝐷.msk,x[m],y[m],f[m]).==1.0)[1:np]
    xyf=permutedims([x[m[n]] y[m[n]] f[m[n]]])
    return DataFrame(x=xyf[1,:],y=xyf[2,:],f=xyf[3,:])
end

"""
    reset_📌!(𝐼::Individuals,frac::Number,📌::Array)

Randomly select a fraction (frac) of the particles and reset 
their positions (𝐼.📌) to a random subset of the specified 📌.
"""
function reset_📌!(𝐼::Individuals,frac::Number,📌::Array)
    np=length(𝐼.🆔)
    n_reset = Int(round(frac*np))
    k_reset = rand(1:np, n_reset)
    l_reset = rand(1:np, n_reset)
    𝐼.📌[k_reset]=deepcopy(📌[l_reset])
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

"""
    setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String)

Define `FlowFields` data structure along with specified grid (`Γ` NamedTuple), 
function `func` (e.g., `(u -> MeshArrays.update_location_llc!(u,Γ)))`, 
and file location (`pth`).
    
_Note: the initial implementation approximates month durations to 
365 days / 12 months for simplicity and sets 𝑃.𝑇 to [-mon/2,mon/2]_
"""
function setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String)
    XC=MeshArrays.exchange(Γ.XC) #add 1 lon point at each edge
    YC=MeshArrays.exchange(Γ.YC) #add 1 lat point at each edge
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
    (u0,v0)=MeshArrays.exchange(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=U[:,𝐷.k]; v1=V[:,𝐷.k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1.*𝐷.iDXC; v1=v1.*𝐷.iDYC; #normalize to grid units
    (u1,v1)=MeshArrays.exchange(u1,v1,1) #add 1 point at each edge for u and v

    𝑃.u0[:]=Float32.(u0[:])
    𝑃.u1[:]=Float32.(u1[:])
    𝑃.v0[:]=Float32.(v0[:])
    𝑃.v1[:]=Float32.(v1[:])
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

"""
    read_velocities(γ::gcmgrid,t::Int,pth::String)

Read velocity components `u,v` from files in `pth`for time `t`
"""
function read_velocities(γ::gcmgrid,t::Int,pth::String)
    u=read_nctiles(joinpath(pth,"UVELMASS/UVELMASS"),"UVELMASS",γ,I=(:,:,:,t))
    v=read_nctiles(joinpath(pth,"VVELMASS/VVELMASS"),"VVELMASS",γ,I=(:,:,:,t))
    return u,v
end

"""
    global_ocean_circulation(;k=1)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.
"""
function global_ocean_circulation(;k=1)

  OceanStateEstimation.get_ecco_velocity_if_needed()

  #read grid and set up connections between subdomains
  γ=MeshArrays.GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
  Γ=MeshArrays.GridLoad(γ)
  f(x,y)=Float32.(MeshArrays.GridLoadVar(x,y))
  tmp=(DXC=f("DXC",γ),DYC=f("DYC",γ),hFacC=f("hFacC",γ),Depth=f("Depth",γ))
  Γ=merge(Γ,tmp)
  Γ=merge(Γ,MeshArrays.NeighborTileIndices_cs(Γ))
  func=(u -> MeshArrays.update_location_llc!(u,Γ))

  #initialize u0,u1 etc
  𝑃,𝐷=setup_FlowFields(k,Γ,func,ScratchSpaces.ECCO)
  𝐷.🔄(𝑃,𝐷,0.0)

  #add background map for plotting
  λ=ECCO_FlowFields.get_interp_coefficients(Γ)
  ODL=ECCO_FlowFields.OceanDepthLog(λ,Γ)
  
  #(optional) fraction of the particles reset per month (e.g., 0.05 for k<=10)
  r_reset = 0.01 

  #add parameters for use in reset!
  tmp=(frac=r_reset, Γ=Γ, ODL=ODL)
  𝐷=merge(𝐷,tmp)

  return 𝑃,𝐷

end

function get_interp_coefficients(Γ)
    fil=joinpath(ScratchSpaces.ECCO,"interp_coeffs_halfdeg.jld2")
    if !isfile(fil)
        url="https://zenodo.org/record/5784905/files/interp_coeffs_halfdeg.jld2"
        OceanStateEstimation.ScratchSpaces.Downloads.download(url,fil;timeout=60000.0)
        #OceanStateEstimation.ECCOdiags_add("interp_coeffs") : nothing
    end
    λ=JLD2.load(fil)
    λ=MeshArrays.Dict_to_NamedTuple(λ)
end

function OceanDepthLog(λ,Γ)
    DL=MeshArrays.Interpolate(λ.μ*Γ.Depth,λ.f,λ.i,λ.j,λ.w)
    DL=reshape(DL,size(λ.lon))
    DL[findall(DL.<0)].=0
    DL=transpose(log10.(DL))
    DL[findall((!isfinite).(DL))].=NaN
    (lon=λ.lon[:,1],lat=λ.lat[1,:],fld=DL,rng=(1.5,5))
end

end #module ECCO_FlowFields
