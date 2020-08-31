
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

Plot positions at time tt[t] in a movie

```
𝐼,Γ=example3("OCCA")

lon,lat=Float64.(Γ["XC"][1]),Float64.(Γ["YC"][1])
DepthLog=Float64.(log10.(Γ["Depth"][1]))
DepthLog[(!isfinite).(DepthLog)].=minimum(DepthLog[isfinite.(DepthLog)])
𝐼.🔴.year=2000 .+ 𝐼.🔴.t ./86400/365; 

scene = MakieMap(DepthLog,colorrange=(3.,4.))
MakieScatterMovie(scene,𝐼.🔴,2000:0.1:2002,"tmp.mp4")
```
"""
function MakieScatterMovie(scene::Scene,df,tt,fil::String)

   xmul,ymul=2,4
   nmax=100000

   tt=collect(tt)
   dt=0.25
   df_t(i) = df[ (df.year.>tt[i]-dt).&(df.year.<=tt[i]) , :]

   xs=fill(NaN,nmax)
   ys=fill(NaN,nmax)
   zs=fill(NaN,nmax)
   nt=length(df_t(1)[!, :lon])
   xs[1:nt] = deepcopy(df_t(1)[!, :lon])
   xs[xs.>180.0]=xs[xs.>180.0] .-360.0
   ys[1:nt] = deepcopy(df_t(1)[!, :lat])
   zs[1:nt] = deepcopy(df_t(1)[!, :z])
   z0=0*zs .- 200.0
 
   Makie.scatter!(scene, xmul*xs, ymul*ys, zs, markersize = 2.0, 
   show_axis = false, color=zs)[end]
   threeD = scene[end]

   Makie.scatter!(scene, xmul*xs, ymul*ys, z0, markersize = 1.0, 
   show_axis = false, color=:black)[end]
   twoD = scene[end]

   scene
   record(scene, fil, 1:length(tt); framerate=12) do i
       xs=fill(NaN,nmax)
       ys=fill(NaN,nmax)
       zs=fill(NaN,nmax)
       nt=length(df_t(i)[!, :lon])
       xs[1:nt] = deepcopy(df_t(i)[!, :lon])
       xs[xs.>180.0]=xs[xs.>180.0] .-360.0
       ys[1:nt] = deepcopy(df_t(i)[!, :lat])
       zs[1:nt] = deepcopy(df_t(i)[!, :z])
       #
       threeD[1] = xmul*xs
       threeD[2] = ymul*ys
       threeD[3] = zs
       threeD[:color] = zs
       #
       twoD[1] = xmul*xs
       twoD[2] = ymul*ys
   end

   return scene
end

"""
    MakieMap(c;colorrange=AbstractPlotting.Automatic())

Map a gridded 2D array, `col`, using Makie

```
scene = MakieMap(OceanDepth,colorrange=(0.,6000.))
```
"""
function MakieMap(col;colorrange=AbstractPlotting.Automatic())
    xmul,ymul=2,4
    xs = xmul*[x for y in lat[1,80:140], x in lon[1:70,1]]
    ys = ymul*[y for y in lat[1,80:140], x in lon[1:70,1]]
    cs = col[1:70,80:140]

    #lim=FRect(-180.0, -90.0, 360.0, 180.0)
    lim=FRect3D([minimum(xs) minimum(ys) -200.0],
                [maximum(xs)-minimum(xs) maximum(ys)-minimum(ys) 200.0])

    scene = Makie.contour(vec(xs[1,:]), vec(ys[:,1]), cs, levels = 15, linewidth = 2, 
    transformation = (:xy, -200.0), color=:black, limits=lim)#,show_axis = false)

    scene.center = false # prevent scene from recentering on display
    update_cam!(scene, Vec3f0(maximum(xs)+200,mean(ys),80.0), 
                       Vec3f0(minimum(xs)-100,mean(ys),-200.0))
                       
    return scene
end
