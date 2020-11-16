"""
    read_drifters(pth,lst;chnk=Inf,rng=(missing,missing))

Read near-surface [drifter data](https://doi.org/10.1002/2016JC011716) from the
[Global Drifter Program](https://doi.org/10.25921/7ntx-z961) into a DataFrame.

Note: need to use NetCDF.jl as NCDatasets.jl errors when TIME = Inf
```
pth,list=drifter_files()
df=read_drifters( pth*lst[end],chnk=1000,rng=(2014.1,2014.2) )
#sort!(df, [:t, :lat])
#CSV.write(pth*"Drifter_hourly_2005_2019.csv", df)
```
"""
function read_drifters(fil::String;chnk=Inf,rng=(missing,missing))
   t=ncread(fil,"TIME")
   t_u=ncgetatt(fil,"TIME","units")
   lo=ncread(fil,"LON")
   la=ncread(fil,"LAT")
   ID=ncread(fil,"ID")

   ii=findall(isfinite.(lo.*la.*t))
   t_ii=t[ii]
   t_ii=timedecode.(t_ii, t_u)
   tmp=dayofyear.(t_ii)+(hour.(t_ii) + minute.(t_ii)/60 ) /24
   t_ii=year.(t_ii)+tmp./daysinyear.(t_ii)
   isfinite(rng[1]) ? jj=findall( (t_ii.>rng[1]).&(t_ii.<=rng[2])) : jj=1:length(ii)

   ii=ii[jj]
   t=t_ii[jj]
   lo=lo[ii]
   la=la[ii]
   ID=ID[ii]

   df = DataFrame(ID=Int[], lon=Float64[], lat=Float64[], t=Float64[])
   
   !isinf(chnk) ? nn=Int(ceil(length(ii)/chnk)) : nn=1
   for jj=1:nn
      #println([jj nn])
      !isinf(chnk) ? i=(jj-1)*chnk.+(1:chnk) : i=(1:length(ii))
      i=i[findall(i.<length(ii))]
      append!(df,DataFrame(lon=lo[i], lat=la[i], t=t[i], ID=Int.(ID[i])))
   end

   return df
end

"""
    read_drifters( pth, lst )

Read near-surface [hourly drifter data](https://doi.org/10.1002/2016JC011716) from the
[Global Drifter Program](https://doi.org/10.25921/7ntx-z961) into a DataFrame.

Note: need to use NetCDF.jl as NCDatasets.jl errors when TIME = Inf

```
pth,list=drifter_files()
df=read_drifters( pth, lst)
```
"""
function read_drifters( pth, lst )
   df = DataFrame([fill(Int, 1) ; fill(Float64, 3)], [:ID, :lon, :lat, :t])
   for fil in lst
      println(fil)
      append!(df,read_drifters( pth*fil,chnk=10000,rng=(2005.0,2020.0) ))
      #println(size(df))
   end
   return df
end
   
function drifter_files()
   pth="Drifter_hourly_v013/"
   lst=["driftertrajGPS_1.03.nc","driftertrajWMLE_1.02_block1.nc","driftertrajWMLE_1.02_block2.nc",
      "driftertrajWMLE_1.02_block3.nc","driftertrajWMLE_1.02_block4.nc","driftertrajWMLE_1.02_block5.nc",
      "driftertrajWMLE_1.02_block6.nc","driftertrajWMLE_1.03_block7.nc"]
   return pth,list
end   

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

"""
    read_uvetc(k::Int,Γ::Dict,pth::String)

Define `uvetc` given the grid variables `Γ` and a vertical level choice `k`
including velocities obtained from files in `pth`

**deprecated: use set_up_𝑃 and update_𝑃! instead -- only still used in GlobalOceanNotebooks ?**
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
    get_ecco_velocity_if_needed()

Download `MITgcm` transport output to `examples/nctiles_climatology` if needed
"""
function get_ecco_velocity_if_needed()
    p=dirname(pathof(OceanStateEstimation))
    lst=joinpath(p,"../examples/nctiles_climatology.csv")
    pth=ECCOclim_path
    !isdir(pth*"UVELMASS") ? get_from_dataverse(lst,"UVELMASS",pth) : nothing
    !isdir(pth*"VVELMASS") ? get_from_dataverse(lst,"VVELMASS",pth) : nothing
    !isdir(pth*"WVELMASS") ? get_from_dataverse(lst,"WVELMASS",pth) : nothing
end

"""
    get_occa_velocity_if_needed()

Download `MITgcm` transport output to `examples/OCCA_climatology` if needed
"""
function get_occa_velocity_if_needed()
    p=dirname(pathof(OceanStateEstimation))
    lst=joinpath(p,"../examples/OCCA_climatology.csv")
    pth=OCCAclim_path
    nams = ("DDuvel.0406clim.nc","DDvvel.0406clim.nc","DDwvel.0406clim.nc","DDtheta.0406clim.nc","DDsalt.0406clim.nc")
    if !isfile("$pth"*"DDuvel.0406clim.nc") 
        tmp=joinpath(pth,"tmp/")
        !isdir(tmp) ? mkdir(tmp) : nothing
        [get_from_dataverse(lst,nam,tmp) for nam in nams]
        [mv(joinpath(tmp,nam,nam),joinpath(pth,nam)) for nam in nams]
    end
end
