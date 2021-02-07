
"""
    RelocationFunctions_cs_check(xmpl,RF,trgt)

Visualize that RelocationFunctions_cs behaves as expected
"""
function RelocationFunctions_cs_check(
    xmpl::MeshArray,
    RF::Array{Function,2},
    trgt::Int,
)

    s = size.(xmpl.f)
    nFaces = length(s)
    nFaces == 5 ? s = vcat(s, s[3]) : nothing

    (aW, aE, aS, aN, iW, iE, iS, iN) = MeshArrays.exch_cs_sources(trgt, s, 1)
    nx, ny = s[trgt]
    p = plot([0.0, nx], [0.0, ny], color = :black)
    plot!([0.0, nx], [ny, 0.0], color = :black)
    for i = 1:nFaces
        (nx, ny) = s[i]
        x = [i - 0.5 for i = 1:nx, j = 1:ny]
        y = [j - 0.5 for i = 1:nx, j = 1:ny]
        c = missing
        if aW == i
            println("source West=$i")
            c = :red
        end
        if aE == i
            println("source East=$i")
            c = :orange
        end
        if aS == i
            println("source South=$i")
            c = :blue
        end
        if aN == i
            println("source North=$i")
            c = :cyan
        end
        if !ismissing(c)
            (x, y) = RF[trgt, i](x, y)
            p = scatter!(
                x,
                y,
                color = c,
                legend = false,
                marker = :rect,
                markerstrokewidth = 0.0,
                markersize = 1.0,
            )
        end
    end
    return p
end
