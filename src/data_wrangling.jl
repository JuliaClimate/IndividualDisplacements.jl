
"""
    postprocess_lonlat()

Copy `sol` to a `DataFrame` & map position to lon,lat coordinates
"""
function postprocess_lonlat(sol,XC,YC)
    ID=collect(1:size(sol,2))*ones(1,size(sol,3))
    x=sol[1,:,:]
    y=sol[2,:,:]
    size(sol,1)==3 ? fIndex=sol[3,:,:] : fIndex=ones(size(x))
    df = DataFrame(ID=Int.(ID[:]), x=x[:], y=y[:], fIndex=fIndex[:])

    lon=Array{Float64,1}(undef,size(df,1)); lat=similar(lon)

    for ii=1:length(lon)
        #get location in grid index space
        x=df[ii,:x]; y=df[ii,:y]; fIndex=Int(df[ii,:fIndex])
        dx,dy=[x - floor(x) .+ 0.5,y - floor(y) .+ 0.5]
        i_c,j_c = Int32.(floor.([x y])) .+ 1
        #interpolate lon and lat to position
        tmp=view(YC[fIndex],i_c:i_c+1,j_c:j_c+1)
        lat[ii]=(1.0-dx)*(1.0-dy)*tmp[1,1]+dx*(1.0-dy)*tmp[2,1]+(1.0-dx)*dy*tmp[1,2]+dx*dy*tmp[2,2]

        tmp=view(XC[fIndex],i_c:i_c+1,j_c:j_c+1)
        if (maximum(tmp)>minimum(tmp)+180)
            tmp1=deepcopy(tmp)
            tmp1[findall(tmp.<maximum(tmp)-180)] .+= 360.
            tmp=tmp1
        end
        #kk=findall(tmp.<maximum(tmp)-180); tmp[kk].=tmp[kk].+360.0
        lon[ii]=(1.0-dx)*(1.0-dy)*tmp[1,1]+dx*(1.0-dy)*tmp[2,1]+(1.0-dx)*dy*tmp[1,2]+dx*dy*tmp[2,2]
    end

    df.lon=lon; df.lat=lat; #show(df[end-3:end,:])
    return df
end

postprocess_lonlat(sol,uvetc::Dict) = postprocess_lonlat(sol,uvetc["XC"],uvetc["YC"])

"""
    postprocess_xy()

Copy `sol` to a `DataFrame` & map position to x,y coordinates,
and define time axis for a simple doubly periodic domain
"""
function postprocess_xy(sol,𝑃)
    nf=size(sol,2)
    nt=size(sol,3)
    nx,ny=size(𝑃["XC"][1])

    x=sol[1,:,:]
    y=sol[2,:,:]
    fIndex=sol[3,:,:]
    ID=collect(1:size(sol,2))*ones(1,size(sol,3))
    df = DataFrame(ID=Int.(ID[:]), fIndex=fIndex[:],
    x=mod.(x[:],Ref(nx)), y=mod.(y[:],Ref(ny)))

    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    df[!,:t]=𝑃["t0"] .+ (𝑃["t1"]-𝑃["t0"])/t[end].*t
    return df
end

"""
    read_uvetc(k::Int,γ::Dict,pth::String)

Define `uvetc` given the grid variables `γ` and a vertical level choice `k`
including velocities obtained from files in `pth`
"""
function read_uvetc(k::Int,γ::Dict,pth::String)
    nt=12; msk=(γ["hFacC"][:,k] .> 0.) #select depth

    u=0. *γ["XC"]; v=0. *γ["XC"];
    for t=1:nt
        (U,V)=read_velocities(γ["XC"].grid,t,pth)
        for i=1:size(u,1)
            u[i]=u[i] + U[i,k]
            v[i]=v[i] + V[i,k]
        end
    end
    u=u ./ nt
    v=v ./ nt #time average

    u[findall(isnan.(u))]=0.0; v[findall(isnan.(v))]=0.0 #mask with 0s rather than NaNs
    u=u./γ["DXC"]; v=v./γ["DYC"]; #normalize to grid units

    (u,v)=exchange(u,v,1) #add 1 point at each edge for u and v
    XC=exchange(γ["XC"]) #add 1 lon point at each edge
    YC=exchange(γ["YC"]) #add 1 lat point at each edge

    t0=0.0; t1=86400*366*10.0; dt=10*86400.0;
    uvetc = Dict("u0" => u, "u1" => u, "v0" => v, "v1" => v,
    "t0" => t0, "t1" => t1, "dt" => dt, "msk" => msk, "XC" => XC, "YC" => YC)
    uvetc=merge(uvetc,IndividualDisplacements.NeighborTileIndices_cs(γ));

    return uvetc
end

"""
    read_uvetc(k::Int,t::Float64,γ::Dict,pth::String)

Define `uvetc` given the grid variables `γ`, a vertical level choice `k`, the
time `t` in `seconds` (Float64), and velocities obtained from files in `pth`.

The two climatological months (`m0`,`m1`) that bracket time `t` will be
extracted (e.g. months 12 & 1 then 1 & 2 and so on).

_Note: the nitial implementation does this only approximately by setting
every months duration to 1 year / 12 for simplicity; should be improved..._
"""
function read_uvetc(k::Int,t::Float64,γ::Dict,pth::String)
    dt=86400.0*365.0/12.0
    t<0.0 ? error("time needs to be positive") : nothing

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    #println([t0/dt,t1/dt,m0,m1])

    (U,V)=read_velocities(γ["XC"].grid,m0,pth)
    u0=U[:,k]; v0=V[:,k]
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    u0=u0./γ["DXC"]; v0=v0./γ["DYC"]; #normalize to grid units
    (u0,v0)=exchange(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(γ["XC"].grid,m1,pth)
    u1=U[:,k]; v1=V[:,k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1./γ["DXC"]; v1=v1./γ["DYC"]; #normalize to grid units
    (u1,v1)=exchange(u1,v1,1) #add 1 point at each edge for u and v

    msk=(γ["hFacC"][:,k] .> 0.) #select depth
    XC=exchange(γ["XC"]) #add 1 lon point at each edge
    YC=exchange(γ["YC"]) #add 1 lat point at each edge

    uvetc = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1,
    "t0" => t0, "t1" => t1, "dt" => dt, "msk" => msk, "XC" => XC, "YC" => YC)
    uvetc=merge(uvetc,IndividualDisplacements.NeighborTileIndices_cs(γ));

    return uvetc
end

"""
    initialize_gridded(uvetc::Dict,n_subset::Int=1)

Define initial condition (u0,du) as a subset of grid points
"""
function initialize_gridded(uvetc::Dict,n_subset::Int=1)
    msk=uvetc["msk"]
    uInitS = Array{Float64,2}(undef, 3, prod(msk.grid.ioSize))

    kk = 0
    for fIndex = 1:5
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
    u=Main.read_nctiles("$pth"*"UVELMASS/UVELMASS","UVELMASS",γ,I=(:,:,:,t))
    v=Main.read_nctiles("$pth"*"VVELMASS/VVELMASS","VVELMASS",γ,I=(:,:,:,t))
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
