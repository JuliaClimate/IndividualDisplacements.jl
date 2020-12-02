using Random, Plots, DataFrames, ColorSchemes

"""
    PlotBasic(df::DataFrame,nn::Integer,dMax::Float64=0.)

Plot random subset of size nn trajectories.
"""
function PlotBasic(df::DataFrame,nn::Integer,dMax::Float64=0.)
   IDs = randperm(maximum(df.ID))
   COs=["w" "y" "g" "k"]

   plt=plot(leg=false)
   df_by_ID = groupby(df, :ID)
   ID_list=[df_by_ID[i][1,:ID] for i in 1:length(df_by_ID)]
   for ii=1:nn
      jj=findall(ID_list.==IDs[ii])[1]
      tmp=df_by_ID[jj]
      if dMax > 0.
         d=abs.(diff(tmp[!,:lon]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
         d=abs.(diff(tmp[!,:lat]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
      end
      CO=COs[mod(ii,4)+1]
      
      hasproperty(df,:z) ? plot!(tmp[!,:lon],tmp[!,:lat],tmp[!,:z],linewidth=0.3) : nothing
      !hasproperty(df,:z) ? plot!(tmp[!,:lon],tmp[!,:lat],linewidth=0.3) : nothing
      #plot!(tmp[!,:lon],tmp[!,:lat],tmp[!,:z],linewidth=0.3)
   end
   return plt
end

"""
    scatter_subset(df,t)

```
nf=size(u0,2)
t=[ceil(i/nf) for i in 1:367*nf]
df[!,:t]=2000 .+ 10/365 * t

@gif for t in 2000:0.1:2016
   scatter_subset(df,t)
end
```
"""
function scatter_subset(df,t)
    dt=0.25
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    scatter(df_t.lon,df_t.lat,markersize=2,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0))
end

"""
    scatter_zcolor(df,t,zc; plt=plot(),cam=(0, 90))

```
t=maximum(df[!,:t])
scatter_zcolor(df,t,df.z)
```
"""
function scatter_zcolor(df,t,zc; plt=plot(),cam=(0, 90), dt=1.0)
    lo=extrema(df.lon); lo=round.(lo).+(-5.0,5.0)
    la=extrema(df.lat); la=round.(la).+(-5.0,5.0)
    de=extrema(zc); de=round.(de).+(-5.0,5.0)

    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    zc_t = Float64.(zc[ (df.t.>t-dt).&(df.t.<=t)])

    #fig=deepcopy(plt)
    scatter(df_t.lon,df_t.lat,zc_t,zcolor = zc_t,
    markersize=2,markerstrokewidth=0.1,camera = cam,
    xlims=lo,ylims=la,zlims=de,clims=de)
end

"""
    scatter_movie(𝐼; cam=(0, 90))

Animation using `scatter_zcolor()
```
𝐼,Γ=example3("OCCA", lon_rng=(-165.0,-155.0),lat_rng=(25.0,35.0), z_init=5.5,)
scatter_movie(𝐼,cam=(70, 70))
```
"""
function scatter_movie(𝐼; cam=(0, 90))
   df=𝐼.🔴
   nf=maximum(df.ID)
   nt=min(size(df,1)/nf,100)
   dt=maximum(df.t)/(nt-1)
   #println("nt="*"$nt"*"dt="*"$dt")
   return @gif for t in 0:nt-1
        scatter_zcolor(df,t*dt,df.z;cam=cam,dt=dt)
   end
end

"""
    phi_scatter(ϕ,df)

```
phi_scatter(ϕ,df)
```
"""
function phi_scatter(ϕ,df)
    nx,ny=size(ϕ)
    contourf(-0.5 .+ (1:nx),-0.5 .+ (1:ny),
             transpose(ϕ),c = :blues,linewidth = 0.1)
    scatter!(df.x,df.y,markersize=4.0,c=:red,marker = (:circle, stroke(0)),
             xlims=(0,nx),ylims=(0,ny),leg=:none)
end

"""
    DL()

Compute Ocean depth logarithm.
"""
function DL(Γ)
    lon=[i for i=19.5:1.0:379.5, j=-78.5:1.0:78.5]
    lat=[j for i=19.5:1.0:379.5, j=-78.5:1.0:78.5]
    DL=interp_to_lonlat(Γ["Depth"],Γ,lon,lat)
    DL[findall(DL.<0)].=0
    DL=transpose(log10.(DL))
    DL[findall((!isfinite).(DL))].=NaN
    return lon[:,1],lat[1,:],DL
end

"""
    plot_end_points(𝐼::Individuals,Γ)

Plot initial and final positions, superimposed on a map of ocean depth log.
"""
function plot_end_points(𝐼::Individuals,Γ)
    lo,la,dl=DL(Γ)
    xlims=extrema(lo)
    ylims=extrema(la)
    plt=contourf(lo,la,dl,clims=(1.5,5),c = :ice, colorbar=false, xlims=xlims,ylims=ylims)

    t=𝑃.𝑇[2]
    df = 𝐼.🔴[ (𝐼.🔴.t.>t-1.0).&(𝐼.🔴.t.<=t) , :]
    lo=deepcopy(df.lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,df.lat,markersize=1.5,c=:red,leg=:none,marker = (:circle, stroke(0)))

    t=0.0
    df = 𝐼.🔴[ (𝐼.🔴.t.>t-1.0).&(𝐼.🔴.t.<=t) , :]
    lo=deepcopy(df.lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,df.lat,markersize=1.5,c=:yellow,leg=:none,marker = (:dot, stroke(0)))

    return plt
end
