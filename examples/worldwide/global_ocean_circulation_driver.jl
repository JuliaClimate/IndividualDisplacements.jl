
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

𝑃,𝐷=init_FlowFields(k=k,backward_time=backward_time)

df = init_positions(np,filename=file_input)
#"z" in names(df) ? nothing : df.z=10.0 .+ 0.0*df.x

𝑆 = init_storage(np,100,length(𝐷.Γ.RC),50)
𝐼 = Individuals(𝑃,df.x,df.y,df.z,df.f,
    (𝐷=merge(𝐷,𝑆),∫=custom∫,🔧=custom🔧,🔴=deepcopy(custom🔴)))

𝑇=(0.0,𝐼.𝑃.𝑇[2])
custom∫!(𝐼,𝑇)

[step!(𝐼) for y=1:ny, m=1:nm]

file_output=joinpath(output_path,file_base*".csv")
CSV.write(file_output, Float32.(𝐼.🔴))

𝐼

end

function step!(𝐼::Individuals)
    𝐼.𝐷.🔄(𝐼)
    custom∫!(𝐼)
end
