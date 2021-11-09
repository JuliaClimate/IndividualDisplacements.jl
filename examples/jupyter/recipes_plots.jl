using Random, Plots, DataFrames, ColorSchemes

import Plots: plot

"""
    plot(𝐼::Individuals)

Plot the initial and final positions as scatter plot in x,y plane.
"""
function plot(𝐼::Individuals)
    🔴_by_t = groupby(𝐼.🔴, :t)
    if (sum(names(🔴_by_t).=="lon")==0)
        fig=scatter(🔴_by_t[1].x,🔴_by_t[1].y,c=:red,label="t0",marker = (:circle, stroke(0)))
        scatter!(🔴_by_t[end].x,🔴_by_t[end].y,c=:blue,label="t1",marker = (:circle, stroke(0)))
    else
        fig=scatter(🔴_by_t[1].lon,🔴_by_t[1].lat,c=:red,label="t0",marker = (:circle, stroke(0)))
        scatter!(🔴_by_t[end].lon,🔴_by_t[end].lat,c=:blue,label="t1",marker = (:circle, stroke(0)))
    end
    return fig
end

"""
    plot_paths(df::DataFrame,nn::Integer,dMax::Float64=0.)

Plot random subset of size nn trajectories / paths.
"""
function plot_paths(df::DataFrame,nn::Integer,dMax::Float64=0.)
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
    scatter_zcolor(df,zc; cam=(0, 90))

```
df=groupby(𝐼.🔴, :t)
scatter_zcolor(df[end],df[end].z)
```
"""
function scatter_zcolor(df,zc; cam=(0, 90))
    lo=extrema(df.lon); lo=round.(lo).+(-5.0,5.0)
    la=extrema(df.lat); la=round.(la).+(-5.0,5.0)
    de=extrema(zc); de=round.(de).+(-5.0,5.0)

    scatter(df.lon,df.lat,zc,zcolor = zc,
    markersize=2,markerstrokewidth=0.1,camera = cam,
    xlims=lo,ylims=la,zlims=de,clims=de,leg=false)
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
    df=groupby(𝐼.🔴, :t)
    return @gif for t in 1:length(df)
        scatter_zcolor(df[t],df[t].z;cam=cam)
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
    OceanDepthLog()

Compute Ocean depth logarithm on regular grid.
"""
function OceanDepthLog(Γ)
    lon=[i for i=20.:2.0:380., j=-79.:2.0:89.]
    lat=[j for i=20.:2.0:380., j=-79.:2.0:89.]
    DL=interp_to_lonlat(Γ.Depth,Γ,lon,lat)
    DL[findall(DL.<0)].=0
    DL=transpose(log10.(DL))
    DL[findall((!isfinite).(DL))].=NaN
    return (lon=lon[:,1],lat=lat[1,:],fld=DL,rng=(1.5,5))
end

"""
    map(𝐼::Individuals,background::NamedTuple)

Plot initial and final positions, superimposed on a map of ocean depth log.
"""
function map(𝐼::Individuals,𝐵::NamedTuple)
    xlims=extrema(𝐵.lon)
    ylims=extrema(𝐵.lat)
    plt=contourf(𝐵.lon,𝐵.lat,𝐵.fld,clims=𝐵.rng,c = :ice, 
    colorbar=false, xlims=xlims,ylims=ylims)

    🔴_by_t = groupby(𝐼.🔴, :t)
    lo=deepcopy(🔴_by_t[1].lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,🔴_by_t[1].lat,markersize=1.5,c=:red,leg=:none,marker = (:circle, stroke(0)))
    lo=deepcopy(🔴_by_t[end].lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,🔴_by_t[end].lat,markersize=1.5,c=:yellow,leg=:none,marker = (:dot, stroke(0)))

    return plt
end
