module ECCO_FlowFields

using IndividualDisplacements, OceanStateEstimation, MITgcmTools

import IndividualDisplacements.DataFrames: DataFrame
import IndividualDisplacements.MeshArrays as MeshArrays
import IndividualDisplacements.MeshArrays: gcmgrid, MeshArray, exchange
import IndividualDisplacements.CSV as CSV

import OceanStateEstimation.ECCO_helpers.JLD2 as JLD2

np=10000 #number of particles
nn=100 #chunk size
backward_time=false

"""
    init_from_file(np ::Int)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_from_file(np ::Int)
    #p=dirname(pathof(IndividualDisplacements))
    #fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    fil="global_ocean_circulation_runs/initial_8_6.csv"
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
    XC=exchange(Γ.XC) #add 1 lon point at each edge
    YC=exchange(Γ.YC) #add 1 lat point at each edge
    iDXC=1. ./Γ.DXC
    iDYC=1. ./Γ.DYC
    γ=Γ.XC.grid
    mon=86400.0*365.0/12.0
    
    if k==0
        msk=Γ.hFacC
        msk=1.0*(msk .> 0.0)
        (_,nr)=size(msk)
        exmsk=similar(msk)
        for k=1:nr
            exmsk[:,k]=exchange(msk[:,k])
        end
        𝑃=FlowFields(MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr+1),MeshArray(γ,Float32,nr+1),
        [-mon/2,mon/2],func)
    else
        msk=Γ.hFacC[:, k]
        msk=1.0*(msk .> 0.0)
        exmsk=exchange(msk)
        𝑃=FlowFields(MeshArray(γ,Float32),MeshArray(γ,Float32),
        MeshArray(γ,Float32),MeshArray(γ,Float32),[-mon/2,mon/2],func)    
    end
    
    𝐷 = (🔄 = update_FlowFields!, pth=pth,
         XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC,
         k=k, msk=msk, exmsk=exmsk, 
         θ0=similar(msk), θ1=similar(msk),
         S0=similar(msk), S1=similar(msk))

    #add parameters related to gridded domain decomposition
    𝐷 = merge(𝐷 , MeshArrays.NeighborTileIndices_cs(Γ))

    frac=0.01 #fraction of the particles reset per month (0.05 for k<=10)
    tmp=(frac=frac, Γ=Γ)
    𝐷=merge(𝐷,tmp)

    k==0 ? nr=50 : nr=1
    𝐷=merge(𝐷, (prof_T=NaN*zeros(np,nr,50),batch_T=zeros(2*nn,nr),local_T=zeros(2*nn)) )
    𝐷=merge(𝐷, (prof_S=NaN*zeros(np,nr,50),batch_S=zeros(2*nn,nr),local_S=zeros(2*nn)) )

    #initialize flow field etc arrays
    #𝐷.🔄(𝑃,𝐷,0.0)

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

    velocity_factor=1.0
    if backward_time
        velocity_factor=-1.0
        m0=13-m0
        m1=13-m1
    end

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
    u0=velocity_factor*U[:,𝐷.k]; v0=velocity_factor*V[:,𝐷.k]
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    u0=u0.*𝐷.iDXC; v0=v0.*𝐷.iDYC; #normalize to grid units
    (u0,v0)=exchange(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=velocity_factor*U[:,𝐷.k]; v1=velocity_factor*V[:,𝐷.k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1.*𝐷.iDXC; v1=v1.*𝐷.iDYC; #normalize to grid units
    (u1,v1)=exchange(u1,v1,1) #add 1 point at each edge for u and v

    𝑃.u0[:]=Float32.(u0[:])
    𝑃.u1[:]=Float32.(u1[:])
    𝑃.v0[:]=Float32.(v0[:])
    𝑃.v1[:]=Float32.(v1[:])

    θ0=read_nctiles(joinpath(𝐷.pth,"THETA/THETA"),"THETA",𝑃.u0.grid,I=(:,:,𝐷.k,m0))
    θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ0[:]=Float32.(θ0[:,1])

    θ1=read_nctiles(joinpath(𝐷.pth,"THETA/THETA"),"THETA",𝑃.u0.grid,I=(:,:,𝐷.k,m1))
    θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ1[:]=Float32.(θ1[:,1])

    S0=read_nctiles(joinpath(𝐷.pth,"SALT/SALT"),"SALT",𝑃.u0.grid,I=(:,:,𝐷.k,m0))
    S0[findall(isnan.(S0))]=0.0 #mask with 0s rather than NaNs
    𝐷.S0[:]=Float32.(S0[:,1])

    S1=read_nctiles(joinpath(𝐷.pth,"SALT/SALT"),"SALT",𝑃.u0.grid,I=(:,:,𝐷.k,m1))
    S1[findall(isnan.(S1))]=0.0 #mask with 0s rather than NaNs
    𝐷.S1[:]=Float32.(S1[:,1])

    𝐷.θ0[:]=exchange(𝐷.θ0)
    𝐷.θ1[:]=exchange(𝐷.θ1)
    𝐷.S0[:]=exchange(𝐷.S0)
    𝐷.S1[:]=exchange(𝐷.S1)

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
function update_FlowFields!(𝑃::𝐹_MeshArray3D,𝐷::NamedTuple,t::AbstractFloat)
    dt=𝑃.𝑇[2]-𝑃.𝑇[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    velocity_factor=1.0
    if backward_time
        velocity_factor=-1.0
        m0=13-m0
        m1=13-m1
    end

    (_,nr)=size(𝐷.Γ.hFacC)

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
    u0=velocity_factor*U; v0=velocity_factor*V
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u0[:,k]=u0[:,k].*𝐷.iDXC; v0[:,k]=v0[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u0[:,k],v0[:,k],1) #add 1 point at each edge for u and v
        u0[:,k]=tmpu
        v0[:,k]=tmpv
    end
    w0=velocity_factor*read_nctiles(joinpath(𝐷.pth,"WVELMASS/WVELMASS"),"WVELMASS",𝑃.u0.grid,I=(:,:,:,m0))
    w0[findall(isnan.(w0))]=0.0 #mask with 0s rather than NaNs

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=velocity_factor*U; v1=velocity_factor*V
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u1[:,k]=u1[:,k].*𝐷.iDXC; v1[:,k]=v1[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u1[:,k],v1[:,k],1) #add 1 point at each edge for u and v
        u1[:,k]=tmpu
        v1[:,k]=tmpv
    end
    w1=velocity_factor*read_nctiles(joinpath(𝐷.pth,"WVELMASS/WVELMASS"),"WVELMASS",𝑃.u0.grid,I=(:,:,:,m1))
    w1[findall(isnan.(w1))]=0.0 #mask with 0s rather than NaNs

    𝑃.u0[:,:]=Float32.(u0[:,:])
    𝑃.u1[:,:]=Float32.(u1[:,:])
    𝑃.v0[:,:]=Float32.(v0[:,:])
    𝑃.v1[:,:]=Float32.(v1[:,:])
    for k=1:nr
        tmpw=exchange(-w0[:,k],1)
        𝑃.w0[:,k]=Float32.(tmpw./𝐷.Γ.DRC[k])
        tmpw=exchange(-w1[:,k],1)
        𝑃.w1[:,k]=Float32.(tmpw./𝐷.Γ.DRC[k])
    end
    𝑃.w0[:,1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,1]=0*exchange(-w1[:,1],1)
    𝑃.w0[:,nr+1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,nr+1]=0*exchange(-w1[:,1],1)

    θ0=read_nctiles(joinpath(𝐷.pth,"THETA/THETA"),"THETA",𝑃.u0.grid,I=(:,:,:,m0))
    θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ0[:,:]=Float32.(θ0[:,:])

    θ1=read_nctiles(joinpath(𝐷.pth,"THETA/THETA"),"THETA",𝑃.u0.grid,I=(:,:,:,m1))
    θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    𝐷.θ1[:,:]=Float32.(θ1[:,:])

    S0=read_nctiles(joinpath(𝐷.pth,"SALT/SALT"),"SALT",𝑃.u0.grid,I=(:,:,:,m0))
    S0[findall(isnan.(S0))]=0.0 #mask with 0s rather than NaNs
    𝐷.S0[:,:]=Float32.(S0[:,:])

    S1=read_nctiles(joinpath(𝐷.pth,"SALT/SALT"),"SALT",𝑃.u0.grid,I=(:,:,:,m1))
    S1[findall(isnan.(S1))]=0.0 #mask with 0s rather than NaNs
    𝐷.S1[:,:]=Float32.(S1[:,:])

    for k=1:nr
        𝐷.θ0[:,k]=exchange(𝐷.θ0[:,k])
        𝐷.θ1[:,k]=exchange(𝐷.θ1[:,k])
        𝐷.S0[:,k]=exchange(𝐷.S0[:,k])
        𝐷.S1[:,k]=exchange(𝐷.S1[:,k])
    end

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
  tmp=( DXC=f("DXC",γ),DYC=f("DYC",γ),hFacC=f("hFacC",γ),
        Depth=f("Depth",γ),RC=f("RC",γ),DRC=f("DRC",γ))
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
	!isfile(fil) ? OceanStateEstimation.ECCOdiags_add("interp_coeffs") : nothing
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
