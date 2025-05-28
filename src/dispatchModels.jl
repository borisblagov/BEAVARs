

function dispatchModel(::Chan2020minn_type,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello Minn")
    intercept = 1;
    if isa(YY_tup[1],Array{})
        YY = YY_tup[1];
    elseif isa(YY_tup[1],TimeArray{})
        YY_TA = YY_tup[1];
        YY = values(YY_TA)
        varList = colnames(YY_TA)
    end
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_β, store_Σ = Chan2020minn(YY,set_strct,hyper_str);
    out_strct = VAROutput_Chan2020minn(store_β,store_Σ,YY)
    return out_strct, set_strct
end

function dispatchModel(::Chan2020iniw_type,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello Independent Normal Inverse Wishart")
    intercept = 1;
    if isa(YY_tup[1],Array{})
        YY = YY_tup[1];
    elseif isa(YY_tup[1],TimeArray{})
        YY_TA = YY_tup[1];
        YY = values(YY_TA)
        varList = colnames(YY_TA)
    end
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_β, store_Σ = Chan2020iniw(YY,set_strct,hyper_str);
    out_strct = VAROutput_Chan2020iniw(store_β,store_Σ,YY)
    return out_strct, set_strct
end

function dispatchModel(::Chan2020iniw_type2,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello Independent Normal Inverse Wishart")
    intercept = 1;
    if isa(YY_tup[1],Array{})
        YY = YY_tup[1];
    elseif isa(YY_tup[1],TimeArray{})
        YY_TA = YY_tup[1];
        YY = values(YY_TA)
        varList = colnames(YY_TA)
    end
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_β, store_Σ = Chan2020iniw(YY,set_strct,hyper_str);
    out_strct = VAROutput_Chan2020iniw(store_β,store_Σ,YY)
    return out_strct, set_strct
end



function dispatchModel(::Chan2020csv_type,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello csv")
    intercept = 1;
    if isa(YY_tup[1],Array{})
        YY = YY_tup[1];
    elseif isa(YY_tup[1],TimeArray{})
        YY_TA = YY_tup[1];
        YY = values(YY_TA)
        varList = colnames(YY_TA)
    end
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_β, store_h, store_Σ, s2_h_store, store_ρ, store_σ_h2, eh_store = Chan2020csv(YY,set_strct,hyper_str);
    out_strct = VAROutput_Chan2020csv(store_β,store_Σ,store_h,s2_h_store, store_ρ, store_σ_h2, eh_store,YY)
    return out_strct, set_strct
end

function dispatchModel(::BGR2010_type,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello BGR2010")
    intercept = 0;
    if isa(YY_tup[1],Array{})
        YY = YY_tup[1];
    elseif isa(YY_tup[1],TimeArray{})
        YY_TA = YY_tup[1];
        YY = values(YY_TA)
        varList = colnames(YY_TA)
    end
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_β, store_Σ = BGR2010(YY,set_strct,hyper_str);
    out_strct = VAROutput(store_β,store_Σ)
    return out_strct, set_strct
end

function dispatchModel(::CPZ2024_type,YY_tup, hyp_strct, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello CPZ2024")
    intercept = 1;
    dataHF_tab  = YY_tup[1]
    dataLF_tab  = YY_tup[2]
    varList     = YY_tup[3]
    trans       = YY_tup[4] # transformation of the LF variables (0: growth rates or 1: log-levels)
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_YY,store_β, store_Σt_inv, M_zsp, z_vec, Sm_bit,store_Σt = CPZ2024(dataHF_tab,dataLF_tab,varList,set_strct,hyp_strct,trans)    
    out_strct = VAROutput_CPZ2024(store_β,store_Σt_inv,store_YY,M_zsp, z_vec, Sm_bit,store_Σt)
    return out_strct, set_strct
end
