
using Random, Makie, MeshArrays, DataFrames, ColorSchemes, Statistics

"""
    PlotMakie(df::DataFrame,nn::Integer)

Plot random subset of size nn trajectories.
"""
function PlotMakie(df::DataFrame,nn::Integer,dMax::Float64=0.)
   IDs = randperm(maximum(df.ID))
   COs=[:gray76 :yellow2 :limegreen :black]

   #scene=Scene(limits=FRect(0, 0, 40, 40),show_axis = false)
   #scene=Scene(limits=FRect(-185, -95, 370, 190),show_axis = false)
   scene=Scene(show_axis = false)
   for ii=1:nn
      tmp=df[df.ID .== IDs[ii], :]
      if dMax > 0.
         d=abs.(diff(tmp[!,:lon]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
         d=abs.(diff(tmp[!,:lat]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
      end
      CO=COs[mod(ii,4)+1]
      Makie.lines!(scene,tmp[!,:lon],tmp[!,:lat],color=CO,linewidth=0.5)
   end

   return scene
end

"""
MakieScatterMovie(scene,df,tt,fil::String)

Animate positions, according to time vector tt, and save movie to mp4 file.

```
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_Makie.jl"))
include(joinpath(p,"../examples/helper_functions.jl"))
module ex3
    fil="../examples/worldwide/three_dimensional_ocean.jl"
    include(joinpath(Main.p,fil))
    export set_up_individuals, Γ
end

using .ex3
𝐼=set_up_individuals(ex3.𝐼,nf=10000);
𝑇=(0.0,2*𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

θ=0.5*(𝐼.𝑃.θ0+𝐼.𝑃.θ1)
scene = MakieBase(θ,2:2:28,Tiso=15)
MakieScatterMovie(scene,𝐼.🔴,0:0.02:4,"tmp.mp4")
```
"""
function MakieScatterMovie(scene::Scene,df,tt,fil::String)

   🔴_by_t = groupby(𝐼.🔴, :t)
   _, threeD, twoD = MakieScatter(scene,🔴_by_t[end])

   np,nmax=length(🔴_by_t[end][:lon]),100000
   (xs,ys,zs)=fill(NaN,nmax),fill(NaN,nmax),fill(NaN,nmax)
   cs=fill(:blue,nmax)
   zmul=1/5
   
   ye=[🔴_by_t[i][1,:year] for i in 1:length(🔴_by_t)]
   tt,dt=collect(tt),0.25

   scene
   record(scene, fil, 1:length(tt); framerate=12) do i
       jj = findall( (ye.>tt[i]-dt).&(ye.<=tt[i]) )
       [xs[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:lon] for j in jj]
       [ys[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:lat] for j in jj]
       [zs[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:z] for j in jj]
       [cs[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:col] for j in jj]
       
       #xs[xs.>180.0] .-= 360.0
       xs[xs.<20.0] .+= 360.0
       #
       threeD[1] = xs
       threeD[2] = ys
       threeD[3] = zmul*zs
       threeD[:color] = zmul*zs
       #
       twoD[1] = xs
       twoD[2] = ys
       twoD[:color] = cs
       #
       rotate_cam!(scene,-0.05, 0., 0.)
   end

   return scene
end

"""
    MakieScatter(scene::Scene,df)

Add a scatter plot of e.g. x,y,z

```
scene = MakieBase(θ,2:2:28;Tiso=15)
_, threeD, twoD = MakieScatter(scene,🔴_by_t[end])
```
"""
function MakieScatter(scene::Scene,df)

    zmul=1/5
    nmax=100000
 
    xs=fill(NaN,nmax)
    ys=fill(NaN,nmax)
    zs=fill(NaN,nmax)
    cs=fill(:black,nmax)
    nt=length(df[!, :lon])
    xs[1:nt] = deepcopy(df[!, :lon])
    #xs[xs.>180.0] .-= 360.0
    xs[xs.<20.0] .+= 360.0  
    ys[1:nt] = deepcopy(df[!, :lat])
    zs[1:nt] = deepcopy(df[!, :z])
    z0=0*zs .- 200.0
    cs[1:nt]=deepcopy(df[!, :col])

    Makie.scatter!(scene, xs, ys, zmul*zs, markersize = 500.0, 
    show_axis = false, color=zs, colormap = :turbo, strokewidth=0.0)[end]
    threeD = scene[end]
 
    Makie.scatter!(scene, xs, ys, zmul*z0, markersize = 500.0, 
    show_axis = false, color=cs, strokewidth=0.0)[end]
    twoD = scene[end]

    return scene, threeD, twoD
end

"""
    MakieBase(θ,T; LONin=140.:0.5:250.,LATin=10.:0.5:50.,DEPin=0.:10:200.)

Contour plot of a gridded 2D array projected onto e.g. 200m depth plane.

```
scene = MakieBase(θ,2:2:28)
```
"""
function MakieBase(θ,T; Tiso=12, LONin=140.:0.5:250.,LATin=10.:0.5:50.,DEPin=0.:10:200.)

    sds=(x1=false,xn=true,y1=true,yn=false)
    zmul=1/5

    isa(T,AbstractRange) ? T=collect(T) : nothing

    lo=[i for i=LONin, j=LATin]
    la=[j for i=LONin, j=LATin]
    lon=deepcopy(lo); lat=deepcopy(la);  

    xlims=extrema(Γ["XC"][1])
    lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    lo[findall(lo.>xlims[2])]=lo[findall(lo.>xlims[2])].-360
    (f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
    IntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)
    
    dMin,dMax=extrema(DEPin)
    d=isosurface(θ,Tiso,Γ["RC"])
    d[1][findall(isnan.(d[1]))].=0.
    dd=interp_to_lonlat(d,IntFac)
    dd[findall(dd.>=-dMin)].=NaN #-dMin
    dd[findall(dd.<=-dMax)].=NaN #-dMax

    kMax=maximum(findall(Γ["RC"].>-dMax))
    θbox=fill(NaN,(size(lo)...,kMax));
    [θbox[:,:,k].=interp_to_lonlat(θ[:,k],IntFac) for k=1:kMax];

    lim=FRect3D([minimum(lon) minimum(lat) -dMax*zmul],
            [maximum(lon)-minimum(lon) maximum(lat)-minimum(lat) (dMax-dMin)*zmul])
    #eyepos=Vec3f0(minimum(lon),mean(lat),-1.1*dMax*zmul)
    #lookat=Vec3f0(maximum(lon)+60,mean(lat)-40,-0.1*dMax*zmul)

    scene=Scene(camera = cam3d!,colormap = :blues)
    #scene.center = false # prevent scene from recentering on display
    #update_cam!(scene, lookat, eyepos)

    #bottom
    scene = Makie.contour!(scene,vec(lon[:,1]), vec(lat[1,:]), θbox[:,:,kMax], levels = T,
    transformation = (:xy, -dMax*zmul), limits=lim, color=:black, linewidth = 2)#,show_axis = false)

    #sides
    sds.x1 ? Makie.contour!(scene,vec(lat[1,:]),vec(Γ["RC"][1:kMax])*zmul,θbox[1,:,:],
    levels=T, transformation = (:yz, lon[1,1]), color=:black, linewidth = 2) : nothing

    sds.xn ? Makie.contour!(scene,vec(lat[end,:]),vec(Γ["RC"][1:kMax])*zmul,θbox[end,:,:],
    levels=T, transformation = (:yz, lon[end,1]), color=:black, linewidth = 2) : nothing

    sds.y1 ? Makie.contour!(scene,vec(lon[:,1]),vec(Γ["RC"][1:kMax])*zmul,θbox[:,1,:],
    levels=T, transformation = (:xz, lat[1,1]), color=:black, linewidth = 2) : nothing

    sds.yn ? Makie.contour!(scene,vec(lon[:,end]),vec(Γ["RC"][1:kMax])*zmul,θbox[:,end,:],
    levels=T, transformation = (:xz, lat[1,end]), color=:black, linewidth = 2) : nothing

    #isotherm
    bw(i::Float64) = ColorSchemes.RGB{Float64}(i, i, i)
    wireframe!(scene,vec(lon[:,1]),vec(lat[1,:]), zmul*dd, 
    linewidth=0.25, color=bw(0.7))

    #xlabel!("lon"); ylabel!("lat"); zlabel!("$zmul x depth")
    xlabel!(""); ylabel!(""); zlabel!("")

    return scene
end
