
"""
    ReadDisplacements(dirIn::String,prec::DataType)

Read displacements from MITgcm output file using MeshArrays
and return as a DataFrame.
"""
function ReadDisplacements(dirIn::String,prec::DataType)

   #load the data into one array
   prec==Float64 ? reclen=8 : reclen=4
   n1=13

   filIn="float_trajectories"
   tmp1=readdir(dirIn)
   tmp1=filter(x -> occursin(filIn,x),tmp1)
   filList=filter(x -> occursin(".data",x),tmp1)
   #hack:
   #filList=filter(x -> occursin("002.002.data",x),tmp1)
   nf=length(filList)

   n2=Array{Int,1}(undef,nf)
   for ff=1:nf
      fil=dirIn*filList[ff]
      #println(fil)
      tmp=stat(fil)
      n2[ff]=Int64(tmp.size/n1/reclen)-1
   end

   arr = Array{prec,2}(undef,(n1+1,sum(n2)));
   ii=0;
   #@softscope for ff=1:nf
   for ff=1:nf
      fil=dirIn*filList[ff]
      fid = open(fil)
      tmp = Array{prec,2}(undef,(n1,n2[ff]+1))
      read!(fid,tmp)
      arr[1:n1,ii+1:ii+n2[ff]] = hton.(tmp[:,2:n2[ff]+1])
      arr[n1+1,ii+1:ii+n2[ff]] .= ff
      ii=ii+n2[ff]
   end

   #sort the whole dataset by time
   jj = sort!([1:ii;], by=i->arr[2,i]); arr=arr[:,jj];
   #arr = sort!(arr, dims=2, by=i->arr[2,i]);

   #nfloats=Int(maximum(arr[1,:]))
   #npoints=counts(Int.(arr[1,:]))

   #reformat data as a DataFrame
   df=DataFrame()
   df.ID=Int.(arr[1,:])
   df.time=Int.(arr[2,:])
   df.lon=arr[3,:]
   df.lat=arr[4,:]
   df.dep=arr[5,:]
   if true
      df.i=arr[6,:]
      df.j=arr[7,:]
      df.k=arr[8,:]
      df.etaN=arr[9,:]
      df.uVel=arr[10,:]
      df.vVel=arr[11,:]
      df.theta=arr[12,:]
      df.salt=arr[13,:]
      df.tile=Int.(arr[14,:])
   end

   nfloats=maximum(df.ID);
   nsteps=maximum(counts(df.ID));

   println("# floats=$nfloats")
   println("# steps=$nsteps")

   return df
end
