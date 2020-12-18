
"""
    postprocess_lonlat(sol,𝑃::NamedTuple; id=missing, 𝑇=missing)

Copy `sol` to a `DataFrame` & map position to lon,lat coordinates
using "exchanged" 𝑃.XC, 𝑃.YC via `add_lonlat!`
"""
function postprocess_lonlat(sol,𝑃::Union{NamedTuple,FlowParameters}; id=missing, 𝑇=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(𝑇) ? 𝑇=𝑃.𝑇 : nothing

    nd=length(size(sol))
    nt=size(sol,nd)
    nf=size(sol,nd-1)
    id=id*ones(1,size(sol,nd))
    if (size(sol,1)>1)&&(nd>2)
        x=sol[1,:,:]
        y=sol[2,:,:]
        fIndex=sol[end,:,:]
    elseif (nd>2)
        x=[sol[1,i,j][1] for i in 1:nf, j in 1:nt]
        y=[sol[1,i,j][2] for i in 1:nf, j in 1:nt]
        fIndex=[sol[1,i,j][end] for i in 1:nf, j in 1:nt]
    else
        x=sol[1,:]
        y=sol[2,:]
        fIndex=sol[end,:]
        nf=1
    end

    𝑃.XC.grid.nFaces==1 ? fIndex=ones(size(x)) : nothing

    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    t=𝑇[1] .+ (𝑇[2]-𝑇[1])/t[end].*t
    
    df = DataFrame(ID=Int.(id[:]), x=x[:], y=y[:], fid=Int.(fIndex[:]), t=t[:])
    add_lonlat!(df,𝑃.XC,𝑃.YC)
    return df
end

"""
    add_lonlat!(df::DataFrame,XC,YC)

Add lon & lat to dataframe using "exchanged" XC, YC
"""
function add_lonlat!(df::DataFrame,XC,YC)
    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;

    lon=zeros(length(x),4)
    [lon[k,:]=XC[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]
    lat=zeros(length(x),4)
    [lat[k,:]=YC[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]

    k=findall(vec(maximum(lon,dims=2)-minimum(lon,dims=2)) .> 180.0)
    tmp=view(lon,k,:)
    tmp[findall(tmp.<0.0)]=tmp[findall(tmp.<0.0)] .+ 360.0

    df.lon=(1.0 .-dx).*(1.0 .-dy).*lon[:,1]+dx.*(1.0 .-dy).*lon[:,2] +
         (1.0 .-dx).*dy.*lon[:,3]+dx.*dy.*lon[:,4]
    df.lat=(1.0 .-dx).*(1.0 .-dy).*lat[:,1]+dx.*(1.0 .-dy).*lat[:,2] +
         (1.0 .-dx).*dy.*lat[:,3]+dx.*dy.*lat[:,4]

    return df
end

"""
    postprocess_xy()

Copy `sol` to a `DataFrame` & map position to x,y coordinates,
and define time axis for a simple doubly periodic domain
"""
function postprocess_xy(sol,𝑃::Union{NamedTuple,FlowParameters}; id=missing, 𝑇=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(𝑇) ? 𝑇=𝑃.𝑇 : nothing

    nf=size(sol,2)
    nt=size(sol,3)

    isa(𝑃.u0,MeshArray) ? (nx,ny)=𝑃.u0.grid.ioSize[1:2] : (nx,ny)=size(𝑃.u0)[1:2]
    nd=length(size(sol))

    id=id*ones(1,size(sol,nd))
    if (size(sol,1)>1)&&(nd>2)
        x=mod.(sol[1,:,:],Ref(nx))
        y=mod.(sol[2,:,:],Ref(ny))
    elseif (nd>2)
        x=[mod(sol[1,i,j][1],nx) for i in 1:nf, j in 1:nt]
        y=[mod(sol[1,i,j][2],ny) for i in 1:nf, j in 1:nt]
    else
        x=mod.(sol[1,:],Ref(nx))
        y=mod.(sol[2,:],Ref(ny))
    end
    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    #size(𝑃.XC,1)>1 ? fIndex=sol[3,:,:] : fIndex=fill(1.0,size(x))
    t=𝑇[1] .+ (𝑇[2]-𝑇[1])/t[end].*t

    return DataFrame(ID=Int.(id[:]), t=t[:], x=x[:], y=y[:])
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

    tmp = dict_to_nt(IndividualDisplacements.NeighborTileIndices_cs(Γ))
    𝐷 = merge(𝐷 , tmp)

    𝑃=𝑃_MeshArray2D{Float64}(MeshArray(γ,Float64),MeshArray(γ,Float64),
    MeshArray(γ,Float64),MeshArray(γ,Float64),[-mon/2,mon/2],func)

    𝐷.🔄(𝑃,𝐷,0.0)

    return 𝑃,𝐷
end

"""
    update_𝑃!(𝑃::Union{NamedTuple,FlowParameters},t::Float64)

Update input data (velocity arrays) and time period (array) inside 𝑃 (𝑃.u0[:], etc, and 𝑃.𝑇[:])
based on the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) input 𝑃.𝑇 is used to infer `dt` between consecutive velocity fields,
(2) periodicity of 12 monthly records, (3) vertical 𝑃.k is selected -- but this could easily be generalized._ 
"""
function update_𝑃!(𝑃::Union{NamedTuple,FlowParameters},𝐷::NamedTuple,t::Float64)
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
    initialize_gridded(𝑃::NamedTuple,n_subset::Int=1)

Define initial condition (u0,du) as a subset of grid points
"""
function initialize_gridded(𝑃::NamedTuple,n_subset::Int=1)
    msk=𝑃.msk
    uInitS = Array{Float64,2}(undef, 3, prod(msk.grid.ioSize))

    kk = 0
    for fIndex = 1:length(msk.fSize)
        nx, ny = msk.fSize[fIndex]
        ii1 = 0.5:1.0:nx
        ii2 = 0.5:1.0:ny
        n1 = length(ii1)
        n2 = length(ii2)
        for i1 in eachindex(ii1)
            for i2 in eachindex(ii2)
                if msk[fIndex][Int(round(i1+0.5)),Int(round(i2+0.5))]
                    kk += 1
                    let kk = kk
                        uInitS[1, kk] = ii1[i1]
                        uInitS[2, kk] = ii2[i2]
                        uInitS[3, kk] = fIndex
                    end
                end
            end
        end
    end

    uInitS=uInitS[:,1:kk]
    du=fill(0.0,size(uInitS));

    uInitS=uInitS[:,1:n_subset:end]
    du=du[:,1:n_subset:end]
    return uInitS,du
end

"""
    randn_lonlat(nn=1,seed=missing)

Randomly distributed longitude, latitude positions on the sphere.
"""
function randn_lonlat(nn=1;seed=missing)
    !ismissing(seed) ? rng = MersenneTwister(1234) : rng = MersenneTwister()
    tmp = randn(rng, Float64, (nn, 3))
    tmpn = tmp ./ sqrt.(sum(tmp.^2, dims=2))
    lon = rad2deg.(atan.(tmpn[:,2], tmpn[:,1]))
    lat = 90.0 .- rad2deg.(acos.(tmpn[:,3]))
    return lon, lat
end

"""
    initialize_lonlat(Γ,lon,lat ; msk=missing)

Define initial condition (u0,du) in grid coordinates (Γ) from longitude
& latitude vectors (lon,lat) optionally with a land mask (msk).
"""
function initialize_lonlat(Γ::Dict,lon::Array{Float64,1},lat::Array{Float64,1};msk=missing)
    (f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,lon,lat)
    ii=findall( ((!isnan).(j_x)).&(j_f.!==0) )
    if !ismissing(msk)
        jj=[msk[Int(j_f[i]),1][ Int(round(j_x[i] .+ 0.5)), Int(round(j_y[i] .+ 0.5)) ] for i in ii]
        ii=ii[findall(jj.>0.0)]
    end
    u0=Array(transpose([j_x[ii] j_y[ii] j_f[ii]]))
    du=similar(u0)
    return u0,du
end

initialize_lonlat(Γ::Dict,lon::Float64,lat::Float64;msk=missing) = initialize_lonlat(Γ,[lon],[lat];msk=msk)

"""
    interp_to_lonlat

Use MeshArrays.Interpolate() to interpolate to e.g. a regular grid (e.g. maps for plotting purposes).

```
using MeshArrays, IndividualDisplacements

lon=[i for i=19.5:1.0:379.5, j=-78.5:1.0:78.5]
lat=[j for i=19.5:1.0:379.5, j=-78.5:1.0:78.5]
(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
IntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)

D=Γ["Depth"]
tmp1=interp_to_lonlat(D,Γ,lon,lat)
tmp2=interp_to_lonlat(D,IntFac)
```
"""
function interp_to_lonlat(X::MeshArray,Γ::Dict,lon,lat)
    (f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

function interp_to_lonlat(X::MeshArray,IntFac::NamedTuple)
    @unpack f,i,j,w,lon,lat = IntFac
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

"""
    interp_to_xy(df::DataFrame,Zin)

Interpolate "exchanged" / "hallo-included" Zin to df[!,:x], df[!,:y] on df[!,:fid]
"""
function interp_to_xy(df::DataFrame,Zin)
    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;

    Z=zeros(length(x),4)
    [Z[k,:]=Zin[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]

    return (1.0 .-dx).*(1.0 .-dy).*Z[:,1]+dx.*(1.0 .-dy).*Z[:,2] +
           (1.0 .-dx).*dy.*Z[:,3]+dx.*dy.*Z[:,4]
end
