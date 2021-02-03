
"""
    scatter_subset(𝑃,df,t=missing,dt=1.0)

```
@gif for t in 𝑃.t0:1.0:𝑃.t1
   scatter_subset(𝑃,df,t)
end
```
"""
function scatter_subset(𝑃,df,t=missing,dt=1.0)
    ismissing(t) ? t=maximum(df[!,:t]) : nothing
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    nx,ny=size(𝑃.XC[1])
    scatter(df_t.x,df_t.y,markersize=2.0,c=:red,
    xlims=(0,nx),ylims=(0,ny),leg=:none,marker = (:circle, stroke(0)))
end
