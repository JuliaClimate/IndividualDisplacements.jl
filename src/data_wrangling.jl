
"""
    postprocess_lonlat(sol,𝑃::NamedTuple)

Copy `sol` to a `DataFrame` & map position to lon,lat coordinates
using "exchanged" 𝑃.XC, 𝑃.YC via `add_lonlat!`
"""
function postprocess_lonlat(sol,𝑃::NamedTuple,id=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    𝐼=id*ones(1,size(sol,3))

    x=sol[1,:,:]
    y=sol[2,:,:]
    𝑃.XC.grid.nFaces>1 ? fIndex=sol[end,:,:] : fIndex=ones(size(x))

    nf=size(sol,2)
    nt=size(sol,3)
    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    t=𝑃.𝑇[1] .+ (𝑃.𝑇[2]-𝑃.𝑇[1])/t[end].*t

    df = DataFrame(ID=Int.(𝐼[:]), x=x[:], y=y[:], fid=Int.(fIndex[:]), t=t[:])
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
function postprocess_xy(sol,𝑃::NamedTuple)
    nf=size(sol,2)
    nt=size(sol,3)
    nx,ny=𝑃.ioSize[1:2]

    x=sol[1,:,:]
    y=sol[2,:,:]
    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    #size(𝑃.XC,1)>1 ? fIndex=sol[3,:,:] : fIndex=fill(1.0,size(x))
    ID=collect(1:size(sol,2))*ones(1,size(sol,3))

    df = DataFrame(ID=Int.(ID[:]), t=𝑃.𝑇[1] .+ (𝑃.𝑇[2]-𝑃.𝑇[1])/t[end].*t,
                   x=mod.(x[:],Ref(nx)), y=mod.(y[:],Ref(ny)))

    return df
end

postprocess_xy(sol,𝑃,id) = postprocess_xy(sol,𝑃)

"""
    read_uvetc(k::Int,Γ::Dict,pth::String)

Define `uvetc` given the grid variables `Γ` and a vertical level choice `k`
including velocities obtained from files in `pth`
"""
function read_uvetc(k::Int,Γ::Dict,pth::String)
    𝑃 = dict_to_nt(IndividualDisplacements.NeighborTileIndices_cs(Γ))
    Γ = dict_to_nt( Γ )
    nt=12; msk=(Γ.hFacC[:,k] .> 0.) #select depth

    u=0. ./Γ.DXC; v=0. ./Γ.DYC;
    for t=1:nt
        (U,V)=read_velocities(Γ.XC.grid,t,pth)
        for i=1:size(u,1)
            u[i]=u[i] + U[i,k]
            v[i]=v[i] + V[i,k]
        end
    end
    u=u ./ nt
    v=v ./ nt #time average

    u[findall(isnan.(u))]=0.0; v[findall(isnan.(v))]=0.0 #mask with 0s rather than NaNs
    u=u./Γ.DXC; v=v./Γ.DYC; #normalize to grid units

    (u,v)=exchange(u,v,1) #add 1 point at each edge for u and v
    XC=exchange(Γ.XC) #add 1 lon point at each edge
    YC=exchange(Γ.YC) #add 1 lat point at each edge

    t0=0.0; t1=86400*366*10.0; dt=10*86400.0;
    tmp = (u0=u, u1=u, v0=v, v1=v, t0=t0, t1=t1, dt=dt, msk=msk, XC=XC, YC=YC)

    return merge(𝑃,tmp)
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

    𝑃 = (u0=MeshArray(γ,Float64), u1=MeshArray(γ,Float64),
         v0=MeshArray(γ,Float64), v1=MeshArray(γ,Float64),
         𝑇=[-mon/2,mon/2], 🔄 = update_𝑃!, pth=pth,
         XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC,
         k=k, msk=Γ["hFacC"][:, k])

    tmp = dict_to_nt(IndividualDisplacements.NeighborTileIndices_cs(Γ))
    𝑃 = merge(𝑃 , tmp)

    𝑃.🔄(k,0.0,𝑃)
    return 𝑃
end

"""
    update_𝑃!(k::Int,t::Float64,𝑃::NamedTuple)

Update velocity and time arrays inside 𝑃 (e.g. 𝑃.u0[:], etc, and 𝑃.𝑇[:])
based on the chosen vertical level `k` and time `t` (in `seconds`).
"""
function update_𝑃!(k::Int,t::Float64,𝑃::NamedTuple)
    dt=𝑃.𝑇[2]-𝑃.𝑇[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝑃.pth)
    u0=U[:,k]; v0=V[:,k]
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    u0=u0.*𝑃.iDXC; v0=v0.*𝑃.iDYC; #normalize to grid units
    (u0,v0)=exchange(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝑃.pth)
    u1=U[:,k]; v1=V[:,k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1.*𝑃.iDXC; v1=v1.*𝑃.iDYC; #normalize to grid units
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
    read_velocities(γ::gcmgrid,t::Int,pth::String)

Read velocity components `u,v` from files in `pth`for time `t`
"""
function read_velocities(γ::gcmgrid,t::Int,pth::String)
    u=read_nctiles("$pth"*"UVELMASS/UVELMASS","UVELMASS",γ,I=(:,:,:,t))
    v=read_nctiles("$pth"*"VVELMASS/VVELMASS","VVELMASS",γ,I=(:,:,:,t))
    return u,v
end

"""
    read_mds(filRoot::String,x::MeshArray)

Read a gridded variable from 2x2 tile files. This is used
in `example2_setup()` with `flt_example/`
"""
function read_mds(filRoot::String,x::MeshArray)
   prec=eltype(x)
   prec==Float64 ? reclen=8 : reclen=4;

   (n1,n2)=Int64.(x.grid.ioSize ./ 2);
   fil=filRoot*".001.001.data"
   tmp1=stat(fil);
   n3=Int64(tmp1.size/n1/n2/reclen);

   v00=x.grid.write(x)
   for ii=1:2; for jj=1:2;
      fid = open(filRoot*".00$ii.00$jj.data")
      fld = Array{prec,1}(undef,(n1*n2*n3))
      read!(fid,fld)
      fld = hton.(fld)

      n3>1 ? s=(n1,n2,n3) : s=(n1,n2)
      v00[1+(ii-1)*n1:ii*n1,1+(jj-1)*n2:jj*n2,:]=reshape(fld,s)
   end; end;

   return x.grid.read(v00,x)
end
