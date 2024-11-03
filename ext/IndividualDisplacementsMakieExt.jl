module IndividualDisplacementsMakieExt

	using Makie, IndividualDisplacements
	import IndividualDisplacements: InDiPlot
	import Makie: plot

	function plot(x::InDiPlot)
		if !isempty(x.options)
			o=x.options
			if string(o.plot_type)=="simple_plot1"
				simple_plot1(x.data[:I],x.data[:ϕ])
            else
				println("unknown option (b)")	
			end
		else
			println("unknown option (a)")
		end
	end

##

"""
    simple_plot1(I,ϕ)

```
using IndividualDisplacements, CairoMakie
include("basics/random_flow_field.jl")
x=InDiPlot( data=(I=I,ϕ=ϕ), options=(plot_type=:simple_plot1,) )
plot(x)
```
"""
function simple_plot1(I,ϕ)
	I_t = groupby(I, :t)
	nt=length(I_t)

	time = Observable(nt)
	xp=@lift( I_t[$time].x )
	yp=@lift( I_t[$time].y )

	siz=size(ϕ)
	xx=-0.5 .+ collect(1:siz[1])
	yy=-0.5 .+ collect(1:siz[2])
	ll=collect(-1.0:0.2:1.0)*maximum(ϕ)
	
	fig=Figure(size = (900, 600)); set_theme!(theme_light())
	a = Axis(fig[1, 1],xlabel="x",ylabel="y", title="Positions and Streamfunction")		
	contourf!(a, xx,yy, ϕ, levels=ll, colormap=:grayC)
	scatter!(a,I_t[1].x,I_t[1].y,color=:green2)
	scatter!(a,xp,yp,color=:red)

	fig
end
##

end
