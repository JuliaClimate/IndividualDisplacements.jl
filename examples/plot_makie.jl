
using Random, Makie

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
