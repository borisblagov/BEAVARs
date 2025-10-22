# ------------------------------------------------------------------
# Minnesota prior with independent Normal inverse Wishart prior

function makeHypSetup(::Chan2020iniw_type)
    return hypChan2020()
end

@doc raw"""
    makeDataSetup(::Chan2020iniw_type,data_tab::TimeArray; var_list =  colnames(data_tab))

Generate a dataset strcture for use with the single-frequency models.
"""
function makeDataSetup(::Chan2020iniw_type,data_tab::TimeArray; var_list =  colnames(data_tab))
    return dataBVAR_TA(data_tab, var_list)
end


function Chan2020iniw(YY,VARSetup::BVARmodelSetup,hypSetup::BVARmodelHypSetup)
    @unpack p,nburn,nsave = VARSetup
    ndraws  = nsave+nburn;

    Y, X, T, n, sigmaP, S_0, Σt_inv, Vβminn_inv, Vβminn_inv_elview, Σ_invsp, Σt_LI, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept = BEAVARs.initMinn(YY,p);

    (idx_kappa1,idx_kappa2, Vβminn, βminn) = prior_Minn(n,p,sigmaP,hypSetup)

    Vβminn_inv_elview[:] = 1.0./Vβminn;             # update the diagonal of Vβminn_inv
    Xsur_den[Xsur_CI] = X[X_CI];                    # update Xsur  

    # allocate output for saving
    store_β = zeros(n^2*p+n,nsave);
    store_Σt = zeros(n,n,nsave);
    for ii = 1:ndraws 
        beta = BEAVARs.Chan2020_drawβ(Σ_invsp,Xsur_den,XtΣ_inv_den,XtΣ_inv_X,Vβminn_inv,βminn,K_β,Y,n,k);
        Σt, Σt_inv = Chan2020_drawΣt(Y,Xsur_den,beta,n,T,S_0,hypSetup.nu0);

        Σ_invsp.nzval[:] = Σt_inv[Σt_LI];               # update ( I(T) ⊗ Σ^{-1} )

        if ii>nburn
            store_β[:,ii-nburn] = beta;
            store_Σt[:,:,ii-nburn] = Σt;
        end
    end

    return store_β, store_Σt
end



@with_kw struct VAROutput_Chan2020iniw <: BVARmodelOutput
    store_β::Array{}      # 
    store_Σ::Array{}      # 
    YY::Array{}           #
end


function forecast(VAROutput::VAROutput_Chan2020iniw,VARSetup)
    @unpack store_β, store_Σ, YY = VAROutput
    @unpack n_fcst,p,nsave = VARSetup
    n = size(YY,2);

    @unpack store_β, store_Σ, YY = VAROutput
    @unpack n_fcst,p,nsave = VARSetup
    n = size(YY,2);

    Yfor3D    = fill(NaN,(p+n_fcst,n,nsave))
    Yfor3D[1:p,:,:] .= @views YY[end-p+1:end,:];
    
    for i_draw = 1:nsave
        Yfor = @views Yfor3D[:,:,i_draw];
        A_draw = @views reshape(store_β[:,i_draw],n*p+1,n);
        Σ_draw = @views store_Σ[:,:,i_draw];
        # Σ_draw =  dropdims(median(store_Σ,dims=3),dims=3);
                
        for i_for = 1:n_fcst
            tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
            tclass = [1;tclass];
            Yfor[p+i_for,:]=tclass'*A_draw  .+ (cholesky(Hermitian(Σ_draw)).U*randn(n,1))';    
        end
    end
    return Yfor3D


end # end function forecast Chan2020iniw()

