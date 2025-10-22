function fanChart(toPlot_mat)
    prct_dim = length(size(toPlot_mat)); # the dimension with the draws (always the last)
    if prct_dim == 3
        n=size(toPlot_mat,2);
    elseif prct_dim ==2
        n = 1;
    end
    Yfor_low1 = percentile_mat(toPlot_mat,0.05,dims=prct_dim);
    Yfor_low = percentile_mat(toPlot_mat,0.16,dims=prct_dim);
    Yfor_med = percentile_mat(toPlot_mat,0.5,dims=prct_dim);
    Yfor_hih = percentile_mat(toPlot_mat,0.84,dims=prct_dim);
    Yfor_hih1 = percentile_mat(toPlot_mat,0.95,dims=prct_dim);

    if n < 4
        lout = (n,1)
    elseif n == 4
        lout = (2,2)
    elseif n>4 
        lout = (convert(Int,ceil(n/4)),4)
    # elseif n>11
    #     lout = (convert(Int,ceil(n/4)),4)
    end

    p_hand = plot(layout=lout)
    for ik in 1:n
        plot!(p_hand,Yfor_med[:,ik],w=0;ribbon=(Yfor_med[:,ik]-Yfor_low1[:,ik],Yfor_hih1[:,ik]-Yfor_med[:,ik]),fillalpha = 0.1,color=1,legend=false,subplot=ik) #,title=varList[ik]
        plot!(p_hand,Yfor_med[:,ik],w=2;ribbon = (Yfor_med[:,ik]-Yfor_low[:,ik],Yfor_hih[:,ik]-Yfor_med[:,ik]),fillalpha=0.05,color=1,subplot=ik)
        
    end
    display(p_hand)
    return Yfor_low1, Yfor_low, Yfor_med, Yfor_hih, Yfor_hih1
end



function fanChart(toPlot_mat,dates::Vector{DateTime})
    prct_dim = length(size(toPlot_mat)); # the dimension with the draws (always theh last)
    if prct_dim == 3
        n=size(toPlot_mat,2);
    elseif prct_dim ==2
        n = 1;
    end
    Yfor_low1 = percentile_mat(toPlot_mat,0.05,dims=prct_dim);
    Yfor_low = percentile_mat(toPlot_mat,0.16,dims=prct_dim);
    Yfor_med = percentile_mat(toPlot_mat,0.5,dims=prct_dim);
    Yfor_hih = percentile_mat(toPlot_mat,0.84,dims=prct_dim);
    Yfor_hih1 = percentile_mat(toPlot_mat,0.95,dims=prct_dim);

    if n < 4
        lout = (n,1)
    elseif n == 4
        lout = (2,2)
    elseif n>4 
        lout = (convert(Int,ceil(n/4)),4)
    # elseif n>11
    #     lout = (convert(Int,ceil(n/4)),4)
    end

    p_hand = plot(layout=lout)
    for ik in 1:n
        plot!(p_hand,dates,Yfor_med[:,ik],w=0;ribbon=(Yfor_med[:,ik]-Yfor_low1[:,ik],Yfor_hih1[:,ik]-Yfor_med[:,ik]),fillalpha = 0.1,color=1,legend=false,subplot=ik) #,title=varList[ik]
        plot!(p_hand,dates,Yfor_med[:,ik],w=2;ribbon = (Yfor_med[:,ik]-Yfor_low[:,ik],Yfor_hih[:,ik]-Yfor_med[:,ik]),fillalpha=0.05,color=1,subplot=ik)
        
    end
    display(p_hand)
    return Yfor_low1, Yfor_low, Yfor_med, Yfor_hih, Yfor_hih1
end