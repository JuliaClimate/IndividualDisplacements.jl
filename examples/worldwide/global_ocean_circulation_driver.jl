
using Distributed, IndividualDisplacements

if !isdefined(Main,:ECCO_FlowFields) && myid()==1
    @everywhere include("ECCO_FlowFields.jl")
    @everywhere using Main.ECCO_FlowFields, IndividualDisplacements,  CSV

    @everywhere path="global_ocean_circulation_inputs"
    @everywhere list=readdir(path)
    @everywhere ii=findall(occursin.(Ref("initial_"),list) .& occursin.(Ref(".csv"),list))
    @everywhere n_per_worker=Int(ceil(length(ii)/nworkers()))

    @everywhere output_path=joinpath(tempdir(),"global_ocean_tmp")
end

"""
main_comp(kk)

To run in parallel, see global_ocean_circulation_support.jl and example below

```
include("global_ocean_circulation_driver.jl")
@everywhere include("global_ocean_circulation_driver.jl")

output_path=joinpath(tempdir(),"global_ocean_tmp")
!ispath(output_path) ? mkdir(output_path) : nothing

@sync @distributed for i in 1:nworkers()
for j in 1:n_per_worker
    k=(i-1)*n_per_worker+j
    main_comp(k)
end
end
```

and for plotting

```
fil=joinpath("global_ocean_circulation_outputs","initial_5_4_▶▶.csv")
df=CSV.read(fil,DataFrame)

include("global_ocean_plotting.jl")
fig,tt=PlottingFunctions.plot([],df)

file_output_mp4=tempname()*".mp4"
PlottingFunctions.record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end
```    
"""
function main_comp(kk)

println(kk)

"""## Initial Position Files"""

list=readdir(path)
ii=findall(occursin.(Ref("initial_"),list) .& occursin.(Ref(".csv"),list))
list=list[ii]

file_input=joinpath(path,list[kk])
file_base = list[kk][1:end-4]
backward_time = false
backward_time ? file_output=file_base*"_◀◀" : file_base=file_base*"_▶▶"

k=0
np=10000
ny=3 #number of years
nm=12 #number of months

P,D=init_FlowFields(k=k,backward_time=backward_time)

df = init_positions(np,filename=file_input)
#"z" in names(df) ? nothing : df.z=10.0 .+ 0.0*df.x

S = init_storage(np,100,length(D.Γ.RC),50)
I = Individuals(P,df.x,df.y,df.z,df.f,
    (D=merge(D,S),∫=custom∫,🔧=custom🔧,🔴=deepcopy(custom🔴)))

T=(0.0,I.P.T[2])
custom∫!(I,T)

[step!(I) for y=1:ny, m=1:nm]

file_output=joinpath(output_path,file_base*".csv")
CSV.write(file_output, Float32.(I.🔴))

I

end

function step!(I::Individuals)
    I.D.🔄(I)
    custom∫!(I)
end
