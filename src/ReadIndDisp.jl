
using StatsBase
using MeshArrays
using DataFrames

"""
    ReadGriddedFields()

Read gridded variables from file using MeshArrays and
return result in uvetc Dictionary.
"""
function ReadGriddedFields()

###### 1) Get gridded variables via MeshArrays.jl


mygrid=gcmgrid("flt_example/","ll",1,[(80,42)], [80 42], Float32, read, write)
nr=8

## Put grid variables in a dictionary:

GridVariables=Dict("XC" => myread(mygrid.path*"XC",MeshArray(mygrid,Float32)),
"YC" => myread(mygrid.path*"YC",MeshArray(mygrid,Float32)),
"XG" => myread(mygrid.path*"XG",MeshArray(mygrid,Float32)),
"YG" => myread(mygrid.path*"YG",MeshArray(mygrid,Float32)),
"dx" => 5000.0)

## Put velocity fields in a dictionary:

t0=0.0 #approximation / simplification
t1=18001.0*3600.0

u0=myread(mygrid.path*"U.0000000001",MeshArray(mygrid,Float32,nr))
u1=myread(mygrid.path*"U.0000018001",MeshArray(mygrid,Float32,nr))
v0=myread(mygrid.path*"V.0000000001",MeshArray(mygrid,Float32,nr))
v1=myread(mygrid.path*"V.0000018001",MeshArray(mygrid,Float32,nr))

kk=3 #3 to match -1406.25 in pkg/flt output
u0=u0[:,kk]; u1=u1[:,kk];
v0=v0[:,kk]; v1=v1[:,kk];

u0=u0./GridVariables["dx"]
u1=u1./GridVariables["dx"]
v0=v0./GridVariables["dx"]
v1=v1./GridVariables["dx"]

## Merge the two dictionaries:

uvetc=Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1)

uvetc=merge(uvetc,GridVariables)

## Visualize velocity fields

mskW=myread(mygrid.path*"hFacW",MeshArray(mygrid,Float32,nr))
mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
mskS=myread(mygrid.path*"hFacS",MeshArray(mygrid,Float32,nr))
mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)

msk=Dict("mskW" => mskW, "mskS" => mskS)

uvetc=merge(uvetc,msk)

end

"""
    ReadDisplacements(dirIn::String,prec::DataType)

Read displacements from MITgcm output file using MeshArrays
and return as a DataFrame.
"""
function ReadDisplacements(dirIn::String,prec::DataType)

   #load the data into one array
   prec==Float64 ? reclen=8 : reclen=4
   n1=13

   filIn="float_trajectories"
   tmp1=readdir(dirIn)
   tmp1=filter(x -> occursin(filIn,x),tmp1)
   filList=filter(x -> occursin(".data",x),tmp1)
   #hack:
   #filList=filter(x -> occursin("002.002.data",x),tmp1)
   nf=length(filList)

   n2=Array{Int,1}(undef,nf)
   for ff=1:nf
      fil=dirIn*filList[ff]
      #println(fil)
      tmp=stat(fil)
      n2[ff]=Int64(tmp.size/n1/reclen)-1
   end

   arr = Array{prec,2}(undef,(n1+1,sum(n2)));
   ii=0;
   #@softscope for ff=1:nf
   for ff=1:nf
      fil=dirIn*filList[ff]
      fid = open(fil)
      tmp = Array{prec,2}(undef,(n1,n2[ff]+1))
      read!(fid,tmp)
      arr[1:n1,ii+1:ii+n2[ff]] = hton.(tmp[:,2:n2[ff]+1])
      arr[n1+1,ii+1:ii+n2[ff]] .= ff
      ii=ii+n2[ff]
   end

   #sort the whole dataset by time
   jj = sort!([1:ii;], by=i->arr[2,i]); arr=arr[:,jj];
   #arr = sort!(arr, dims=2, by=i->arr[2,i]);

   #nfloats=Int(maximum(arr[1,:]))
   #npoints=counts(Int.(arr[1,:]))

   #reformat data as a DataFrame
   df=DataFrame()
   df.ID=Int.(arr[1,:])
   df.time=Int.(arr[2,:])
   df.lon=arr[3,:]
   df.lat=arr[4,:]
   df.dep=arr[5,:]
   if true
      df.i=arr[6,:]
      df.j=arr[7,:]
      df.k=arr[8,:]
      df.etaN=arr[9,:]
      df.uVel=arr[10,:]
      df.vVel=arr[11,:]
      df.theta=arr[12,:]
      df.salt=arr[13,:]
      df.tile=Int.(arr[14,:])
   end

   nfloats=maximum(df.ID);
   nsteps=maximum(counts(df.ID));

   println("# floats=$nfloats")
   println("# steps=$nsteps")

   return df
end

"""
    myread()

Read a gridded variable from 2x2 tile files. This is used
in `ReadGriddedFields()` with `flt_example/`
"""
function myread(filRoot::String,x::MeshArray)
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
