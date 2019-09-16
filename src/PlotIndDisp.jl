
#PlotIndDisp.jl

using DataFrames, Random
#using PyPlot, PyCall

"""
    PlotBasic(df::DataFrame,nn::Integer)

Plot random subset of size nn trajectories.
"""
function PlotBasic(df::DataFrame,nn::Integer)

   #PyPlot.figure()

   IDs = randperm(maximum(df.ID))
   COs=["k" "r" "b" "m" "c"]
   #for global ocean case: dMax=90.
   #for flt_example case:
   dMax=100000.
   println("dMax=$dMax")

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
      CO=COs[mod(ii,3)+1]
      PyPlot.plot(tmp[!,:lon],tmp[!,:lat],color=CO,linewidth=0.5)
   end

   #to display the figure in JUNO (not needed in REPL or Jupyter):
   #gcf()

end

"""
    PlotMapProj(df::DataFrame,nn::Integer)

Plot random subset of size nn trajectories using PyPlot & basemap.
"""
function PlotMapProj(df::DataFrame,nn::Integer)

   #PyPlot.figure()

   # Set up Equidistant cylindrical map projection. Use low resolution coastlines.
   basemap=pyimport("mpl_toolkits.basemap")
   map = basemap.Basemap(projection="cyl", llcrnrlat=-80, urcrnrlat=90,
   llcrnrlon=20,urcrnrlon=380, resolution="l")

   # Draw coastlines, country boundaries, fill continents.
   map.drawcoastlines(linewidth=0.25)
   map.drawcountries(linewidth=0.25)
   map.fillcontinents(color="grey")

   # Draw the edge of the map projection region (the projection limb)
   #map[:drawmapboundary](fill_color="aqua")

   # Draw lat/lon grid lines every 30 degrees.
   map.drawmeridians(collect(0:30:360))
   map.drawparallels(collect(-90:30:90))

   # Draw trajectories
   IDs = randperm(maximum(df.ID))
   COs=["k" "r" "b" "m" "c"]

   #for global ocean case:
   dMax=90.
   println("dMax=$dMax")

   for ii=1:nn
      tmp=df[df.ID .== IDs[ii], :]
      jj=findall(tmp[!,:lon] .< 20)
      tmp[jj,:lon]=tmp[jj,:lon] .+ 360.0;
      if dMax > 0.
         d=abs.(diff(tmp[!,:lon]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
         d=abs.(diff(tmp[!,:lat]))
         jj=findall(d .> dMax)
         tmp[jj,:lon].=NaN; tmp[jj,:lat].=NaN
      end
      #
      CO=COs[mod(ii,3)+1]
      x, y = map(tmp[!,:lon], tmp[!,:lat])
      map.plot(x, y,color=CO,linewidth=0.5)
   end

   #to display the figure in JUNO (not needed in REPL or Jupyter):
   #gcf()

end
