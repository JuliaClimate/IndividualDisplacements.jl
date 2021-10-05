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

