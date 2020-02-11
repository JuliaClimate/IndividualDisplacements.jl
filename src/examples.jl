
"""
    example1()

Global ocean case -- just reading from file for now.

```
df=IndividualDisplacements.example1()

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"plot_pyplot.jl"))
PyPlot.figure(); PlotMapProj(df,300); gcf()
```
"""
function example1()
   dirIn="run_offflt/"
   prec=Float32
   df=ReadDisplacements(dirIn,prec)
end

"""
    example2()

Reproducing `MITgcm/verification/flt_example/` case. This is based on an
extended and modified configuration of the standard MITgcm test case.

```
(df,ref,sol)=IndividualDisplacements.example2();

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"plot_pyplot.jl"))
PyPlot.figure(); PlotBasic(df,300,100000.0); gcf()

using Plots
Plots.plot(sol[1,:],sol[2,:],linewidth=5,lc=:black, title="One Trajectory Example",
xaxis="x",yaxis="y",label="Julia Solution") # legend=false
pl=Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,lc=:red,label="MITgcm Solution")
```
"""
function example2()
   dirIn="flt_example/"
   prec=Float32
   df=ReadDisplacements(dirIn,prec)
   uvetc=IndividualDisplacements.example2_setup()
   #
   tmp=df[df.ID .== 200, :]
   nSteps=Int32(tmp[end,:time]/3600)-2
   ref=transpose([tmp[1:nSteps,:lon] tmp[1:nSteps,:lat]])
   maxLon=80*5.e3
   maxLat=42*5.e3
   for i=1:nSteps-1
       ref[1,i+1]-ref[1,i]>maxLon/2 ? ref[1,i+1:end]-=fill(maxLon,(nSteps-i)) : nothing
       ref[1,i+1]-ref[1,i]<-maxLon/2 ? ref[1,i+1:end]+=fill(maxLon,(nSteps-i)) : nothing
       ref[2,i+1]-ref[2,i]>maxLat/2 ? ref[2,i+1:end]-=fill(maxLat,(nSteps-i)) : nothing
       ref[2,i+1]-ref[2,i]<-maxLat/2 ? ref[2,i+1:end]+=fill(maxLat,(nSteps-i)) : nothing
   end
   ref=ref./uvetc["dx"]
   #
   comp_vel=IndividualDisplacements.VelComp
   get_vel=IndividualDisplacements.VelCopy
   uInit=[tmp[1,:lon];tmp[1,:lat]]./uvetc["dx"]
   du=fill(0.0,2)
   #
   tspan = (0.0,nSteps*3600.0)
   #prob = ODEProblem(get_vel,uInit,tspan,tmp)
   prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
   sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
   #
   return df,ref,sol
end

"""
example2_setup()

Read gridded variables from file using MeshArrays and
return result in uvetc Dictionary.
"""
function example2_setup()

   ###### 1) Get gridded variables via MeshArrays.jl


   mygrid=gcmgrid("flt_example/","PeriodicChanel",1,[(80,42)], [80 42], Float32, read, write)
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
myread()

Read a gridded variable from 2x2 tile files. This is used
in `example2_setup()` with `flt_example/`
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
