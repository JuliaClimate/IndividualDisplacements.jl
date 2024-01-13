module ECCO_FlowFields

using IndividualDisplacements, OceanStateEstimation, MITgcmTools, CSV, JLD2

import IndividualDisplacements.OrdinaryDiffEq: solve, Tsit5, ODEProblem
import IndividualDisplacements: update_location!
import IndividualDisplacements.DataFrames: DataFrame
import IndividualDisplacements.MeshArrays as MeshArrays
import IndividualDisplacements.MeshArrays: gcmgrid, MeshArray, exchange

export init_FlowFields, init_positions, init_storage
export custom∫, custom🔧, custom🔴, custom∫!
#export reset_📌!, init_z_if_needed

"""
    init_positions(np ::Int)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_positions(np ::Int; filename="global_ocean_circulation.csv")
    if filename=="global_ocean_circulation.csv"
        p=dirname(pathof(IndividualDisplacements))
        fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    else
        fil=filename
    end
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
    init_gulf_stream(np ::Int , 𝐷::NamedTuple)

Randomly distribute `np` points in the Florida Strait region, within 
`𝐷.msk` region, and return position in grid index space (`i,j,subdomain`).
"""
function init_gulf_stream(np ::Int , 𝐷::NamedTuple)
	lons=[-81,-79]
	lats=[26,28]
	lon=rand(2*np)*diff(lons)[1].+lons[1]
	lat=rand(2*np)*diff(lats)[1].+lats[1]
	
	(_,_,_,_,f,x,y)=IndividualDisplacements.InterpolationFactors(𝐷.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(IndividualDisplacements.nearest_to_xy(𝐷.msk,x[m],y[m],f[m]).==1.0)[1:np]
    xyf=permutedims([x[m[n]] y[m[n]] f[m[n]]])

	#z=rand(np)*27
	z=15 .+rand(np)*12
    return DataFrame(x=xyf[1,:],y=xyf[2,:],z=z,f=xyf[3,:])
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
    #isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    #𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

"""
    init_storage(np)

- np=number of individuals to store
- nn=number of individuals per chunk
- nl=number of vertical levels to store
- nt=number of time steps to store
"""
function init_storage(np,nn,nl,nt=2)
    (   prof_T=NaN*zeros(np,nl,nt),batch_T=zeros(2*nn,nl),local_T=zeros(2*nn),
        prof_S=NaN*zeros(np,nl,nt),batch_S=zeros(2*nn,nl),local_S=zeros(2*nn)
        ) 
end

"""
    setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String)

Define `FlowFields` data structure along with specified grid (`Γ` NamedTuple), 
function `func` (e.g., `(u -> MeshArrays.update_location_llc!(u,Γ)))`, 
and file location (`pth`).
    
_Note: the initial implementation approximates month durations to 
365 days / 12 months for simplicity and sets 𝑃.𝑇 to [-mon/2,mon/2]_
"""
function setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String,backward_time=false)
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

    tmp=(Γ=Γ, backward_time=backward_time)
    𝐷=merge(𝐷,tmp)

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
    if 𝐷.backward_time
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
    if 𝐷.backward_time
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
    update_FlowFields!(𝐼::Individuals)

Update flow field arrays in 𝐼.𝑃, update 𝐼.𝑃.𝑇, and ancillary variables in 𝐼.𝐷
such that 𝐼.𝑃.𝑇 includes time `t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2])`.

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(𝑃.𝑇), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical 𝑃.k is selected_
"""
function update_FlowFields!(𝐼::Individuals)
    t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2])
    𝐼.𝐷.🔄(𝐼.𝑃,𝐼.𝐷,t_ϵ)
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
    init_FlowFields(;k=1)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.
"""
function init_FlowFields(; k=1, backward_time=false)

  OceanStateEstimation.get_ecco_velocity_if_needed()

  #read grid and set up connections between subdomains
  γ=MeshArrays.GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
  Γ=MeshArrays.GridLoad(γ,option="full")
  f(x,y)=Float32.(MeshArrays.GridLoadVar(x,y))
  tmp=( DXC=f("DXC",γ),DYC=f("DYC",γ),hFacC=f("hFacC",γ),
        Depth=f("Depth",γ),RC=f("RC",γ),DRC=f("DRC",γ))
  Γ=merge(Γ,tmp)
  Γ=merge(Γ,MeshArrays.NeighborTileIndices_cs(Γ))
  func=(u -> MeshArrays.update_location_llc!(u,Γ))

  #initialize u0,u1 etc
  𝑃,𝐷=setup_FlowFields(k,Γ,func,ScratchSpaces.ECCO,backward_time)
  𝐷.🔄(𝑃,𝐷,0.0)

  #add background map for plotting
  λ=get_interp_coefficients(Γ)
  ODL=OceanDepthLog(λ,Γ)
  
  #(optional) fraction of the particles reset per month (e.g., 0.05 for k<=10)
  r_reset = 0.05

  #add parameters for use in reset!
  tmp=(frac=r_reset, ODL=ODL)
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

custom∫(prob) = IndividualDisplacements.ensemble_solver(prob,solver=Tsit5(),reltol=1e-5,abstol=1e-5)

custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[], 
lon=Float64[], lat=Float64[], z=Float64[], θ=Float64[], SSθ=Float64[],
S=Float64[], SSS=Float64[], year=Float64[], t=Float64[])

function custom🔧(sol,𝐹::𝐹_MeshArray3D,𝐷::NamedTuple;id=missing,𝑇=missing)

    df=postprocess_MeshArray(sol,𝐹,𝐷,id=id,𝑇=𝑇)
    np=length(sol.u)
    z=[[sol.u[i][1][3] for i in 1:np];[sol.u[i][end][3] for i in 1:np]]
    df.z=z[:]
    df.year=df.t ./86400/365
    add_lonlat!(df,𝐷.XC,𝐷.YC)

    #for k in 1:nr
    # 𝐷.batch_T[:,k]=interp_to_xy(df,𝐷.θ1[:,k])./interp_to_xy(df,𝐷.exmsk[:,k])
    #end

    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;
    
    nr=size(𝐷.exmsk,2)

    #need time interpolation (df.t)
    for k in 1:nr, jj in 1:length(i_c)
        tmp0=(1.0-dx[jj])*(1.0-dy[jj])*𝐷.exmsk[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*𝐷.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*𝐷.exmsk[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*𝐷.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        #
        tmp1=(1.0-dx[jj])*(1.0-dy[jj])*𝐷.θ1[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*𝐷.θ1[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*𝐷.θ1[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*𝐷.θ1[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        𝐷.batch_T[jj,k]=tmp1/tmp0
        #
        tmp1=(1.0-dx[jj])*(1.0-dy[jj])*𝐷.S1[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*𝐷.S1[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*𝐷.S1[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*𝐷.S1[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        𝐷.batch_S[jj,k]=tmp1/tmp0
    end

    #need time interpolation (df.t)
    for p=1:size(df,1)
        #k=max(Int(floor(z[p])),1)
        #local_T[p]=batch_T[p,k]
        k1=floor(z[p]+0.5)
        a2=(z[p]+0.5)-k1
        k2=Int(min(max(k1+1,1),nr))
        k1=Int(min(max(k1,1),nr))
        𝐷.local_T[p]=(1-a2)*𝐷.batch_T[p,k1]+a2*𝐷.batch_T[p,k2]
        𝐷.local_S[p]=(1-a2)*𝐷.batch_S[p,k1]+a2*𝐷.batch_S[p,k2]
    end

    df.SSθ=𝐷.batch_T[:,1]
    df.θ=𝐷.local_T[:]
    df.SSS=𝐷.batch_S[:,1]
    df.S=𝐷.local_S[:]

    if 𝐷.frac==0.0
        t=Int(round(0.5+df.t[end]/(𝑇[2]-𝑇[1])))
        nn=Int(size(𝐷.batch_T,1)/2)
        #df.ID[1]==1 ? println(t) : nothing        
        𝐷.prof_T[df.ID[nn+1:end],:,t].=𝐷.batch_T[nn+1:end,:]
        𝐷.prof_S[df.ID[nn+1:end],:,t].=𝐷.batch_S[nn+1:end,:]
    end

    return df
end

function custom∫!(𝐼::Individuals,𝑇)
    (; 🚄,📌,𝑃,𝐷,🔧,🆔,🔴,∫) = 𝐼

    vel=0*vec(📌)
    [🚄(vel[i],📌[i],𝑃,𝑇[1]) for i in 1:length(vel)]
    nd=ndims(vel[1])-1
    vel=[sqrt(sum(vel[ii][1:nd].^2)) for ii in eachindex(vel)]
    vel=[(ii,vel[ii]) for ii=1:length(vel)]
    sort!(vel, by = x -> x[2])
    ii=[vel[ii][1] for ii=1:length(vel)]
    nn=Int(size(𝐷.batch_T,1)/2)
    ni=Int(ceil(length(ii)/nn))
    
    tmp=deepcopy(custom🔴)
    for i=1:ni
        mm=min(nn,length(ii)-nn*(i-1))
#        println("i="*string(i))
        jj=ii[nn*(i-1) .+ collect(1:mm)]
        prob = ODEProblem(🚄,permutedims(📌[jj]), 𝑇 ,𝑃)
        sol = ∫(prob)
        append!(tmp, 🔧(sol,𝑃,𝐷,id=🆔[jj], 𝑇=𝑇))
        📌[jj] = deepcopy([sol[i].u[end] for i in 1:mm])
        if isa(𝑃,𝐹_MeshArray3D)||isa(𝑃,𝐹_MeshArray2D)
            [update_location!(i,𝑃) for i in 📌[jj]]
        end
    end
    sort!(tmp, IndividualDisplacements.order(:t))

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])
#    isempty(🔴) ? append!(🔴,tmp) : 🔴[:,:]=tmp[:,:]
    return true
end

custom∫!(x) = custom∫!(x,(x.𝑃.𝑇[1],x.𝑃.𝑇[2]))

end #module ECCO_FlowFields
