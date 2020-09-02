
using Random, Makie, DataFrames, ColorSchemes

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
𝐼,Γ=example3("OCCA");

lon,lat=Float64.(Γ["XC"][1]),Float64.(Γ["YC"][1])
DepthLog=Float64.(log10.(Γ["Depth"][1]))
DepthLog[(!isfinite).(DepthLog)].=minimum(DepthLog[isfinite.(DepthLog)])
(lon,lat,DepthLog)=circshift.((lon,lat,DepthLog),Ref((-200,0)));
lon[findall(lon.<20)] .+= 360.0;

𝐼.🔴.year=𝐼.🔴.t ./86400/365; 
🔴_by_t = groupby(𝐼.🔴, :t);

scene = MakieMap(DepthLog,colorrange=(3.,4.))
MakieScatterMovie(scene,𝐼.🔴,0:0.05:2,"tmp.mp4")
```
"""
function MakieScatterMovie(scene::Scene,df,tt,fil::String)

   🔴_by_t = groupby(𝐼.🔴, :t)
   _, threeD, twoD = MakieScatter(scene,🔴_by_t[end])

   np,nmax=length(🔴_by_t[end][:lon]),100000
   xs,ys,zs=fill(NaN,nmax),fill(NaN,nmax),fill(NaN,nmax)
   zmul=1/5
   
   ye=[🔴_by_t[i][1,:year] for i in 1:length(🔴_by_t)]
   tt,dt=collect(tt),0.25

   scene
   record(scene, fil, 1:length(tt); framerate=12) do i
       jj = findall( (ye.>tt[i]-dt).&(ye.<=tt[i]) )
       [xs[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:lon] for j in jj]
       [ys[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:lat] for j in jj]
       [zs[collect((1:np).+(j-jj[1])*np)]=🔴_by_t[j][:,:z] for j in jj]
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
   end

   return scene
end

"""
    MakieScatter(scene::Scene,df)

Add a scatter plot of e.g. x,y,z

```
scene = MakieMap(DepthLog,colorrange=(3.,4.))
_, threeD, twoD = MakieScatter(scene,🔴_by_t[end])
```
"""
function MakieScatter(scene::Scene,df)

    zmul=1/5
    nmax=100000
 
    xs=fill(NaN,nmax)
    ys=fill(NaN,nmax)
    zs=fill(NaN,nmax)
    nt=length(df[!, :lon])
    xs[1:nt] = deepcopy(df[!, :lon])
    #xs[xs.>180.0] .-= 360.0
    xs[xs.<20.0] .+= 360.0
    ys[1:nt] = deepcopy(df[!, :lat])
    zs[1:nt] = deepcopy(df[!, :z])
    z0=0*zs .- 200.0
  
    Makie.scatter!(scene, xs, ys, zmul*zs, markersize = 2.0, 
    show_axis = false, color=zs)[end]
    threeD = scene[end]
 
    Makie.scatter!(scene, xs, ys, zmul*z0, markersize = 1.0, 
    show_axis = false, color=:black)[end]
    twoD = scene[end]

    return scene, threeD, twoD
end

"""
    MakieMap(col;colorrange=AbstractPlotting.Automatic())

Contour plot of a gridded 2D array, `col`, projected onto e.g. 200m depth plane.

```
scene = MakieMap(OceanDepth,colorrange=(0.,6000.))
```
"""
function MakieMap(col;colorrange=AbstractPlotting.Automatic())
    zmul=1/5
    xs = [x for y in lat[1,91:130], x in lon[121:230,1]]
    ys = [y for y in lat[1,91:130], x in lon[121:230,1]]
    cs = zmul*col[121:230,91:130]

#    xs = [x for y in lat[1,:], x in lon[:,1]]
#    ys = [y for y in lat[1,:], x in lon[:,1]]
#    cs = col[:,:]

    lim=FRect3D([minimum(xs) minimum(ys) -200.0*zmul],
                [maximum(xs)-minimum(xs) maximum(ys)-minimum(ys) 200.0*zmul])

    scene = Makie.contour(vec(xs[1,:]), vec(ys[:,1]), cs, levels = 15, linewidth = 2, 
    transformation = (:xy, -200.0*zmul), color=:black, limits=lim)#,show_axis = false)

    scene.center = false # prevent scene from recentering on display
    update_cam!(scene, Vec3f0(maximum(xs)+60,mean(ys),-20.0*zmul), 
                       Vec3f0(minimum(xs),mean(ys),-220.0*zmul))
    xlabel!("lon"); ylabel!("lat"); zlabel!("$zmul x depth")

    return scene
end
