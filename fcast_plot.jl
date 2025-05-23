n = size(Yfor3D,2);
Yfor_low1 = percentile_mat(Yfor3D,0.05,dims=3);
Yfor_low = percentile_mat(Yfor3D,0.16,dims=3);
Yfor_med = median(Yfor3D,dims=3)
Yfor_hih = percentile_mat(Yfor3D,0.84,dims=3);
Yfor_hih1 = percentile_mat(Yfor3D,0.95,dims=3);

if n < 4
    lout = (n,1)
elseif n == 4
    lout = (2,2)
elseif n>4 
    lout = (convert(Int,ceil(n/4)),4)
# elseif n>11
#     lout = (convert(Int,ceil(n/4)),4)
end

p = plot(layout=lout)
for ik in 1:n
    plot!(p,Yfor_med[:,ik],w=0;ribbon=(Yfor_med[:,ik]-Yfor_low1[:,ik],Yfor_hih1[:,ik]-Yfor_med[:,ik]),fillalpha = 0.1,color=1,legend=false,subplot=ik,title=varList[ik])
    plot!(p,Yfor_med[:,ik],w=2;ribbon = (Yfor_med[:,ik]-Yfor_low[:,ik],Yfor_hih[:,ik]-Yfor_med[:,ik]),fillalpha=0.05,color=1,subplot=ik)
    
end
display(p)