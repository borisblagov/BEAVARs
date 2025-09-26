# ------------------------------------------------------------------
# Classical Minnesota prior with fixed variance-covariance matrix Σ 


function makeHypSetup(::Chan2020minn_type)
    return hypChan2020()
end

@doc raw"""
    Prepares the structure containg the data for Bayesian VAR using the Chann2020 type. Uses Time Arrays from the TimeSeries package
"""
function makeDataSetup(::Chan2020minn_type,data_tab::TimeArray; var_list =  colnames(data_tab))
    return dataBVAR_TA(data_tab, var_list)
end



@doc raw"""
    Chan2020minn(YY,VARSetup,hypSetup)

Implements the classic homoscedastic Minnesota prior with a SUR form following Chan (2020)

"""
function Chan2020minn(YY,VARSetup::BVARmodelSetup,hypSetup::BVARmodelHypSetup)
    @unpack p,nburn,nsave = VARSetup
    
    Y, X, T, n, sigmaP, S_0, Σt_inv, Vβminn_inv, Vβminn_inv_elview, Σ_invsp, Σt_LI, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept = BEAVARs.initMinn(YY,p);

    (idx_kappa1,idx_kappa2, Vβminn, βminn) = prior_Minn(n,p,sigmaP,hypSetup)

    Vβminn_inv_elview[:] = 1.0./Vβminn;             # update the diagonal of Vβminn_inv
    Xsur_den[Xsur_CI] = X[X_CI];                    # update Xsur  
    mul!(XtΣ_inv_den,Xsur_den',Σ_invsp);            #  X'*( I(T) ⊗ Σ^{-1} )
    mul!(XtΣ_inv_X,XtΣ_inv_den,Xsur_den);           #  X'*( I(T) ⊗ Σ^{-1} )*X
    K_β[:,:] .= Vβminn_inv .+ XtΣ_inv_X;            #  K_β = V^{-1} + X'*( I(T) ⊗ Σ^{-1} )*X
    prior_mean = Vβminn_inv*βminn;                  #  V^-1 * βminn 
    mul!(prior_mean,XtΣ_inv_den, vec(Y'),1.0,1.0);  # (V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y
    cholK_β = cholesky(Hermitian(K_β));             # C is lower triangular, C' is upper triangular
    beta_hat = ldiv!(cholK_β.U,ldiv!(cholK_β.L,prior_mean));    # C'\(C*(V^-1_Minn * beta_Minn + X' ( I(T) ⊗ Σ-1 ) y)
    

    ndraws = nsave+nburn;
    store_β=zeros(n^2*p+n,nsave)
    for ii = 1:ndraws 
        beta = beta_hat + ldiv!(cholK_β.U,randn(k*n,)); # draw for β
        if ii>nburn
            store_β[:,ii-nburn] = beta;
        end
    end

    store_Σ = repeat(vec(S_0),1,nsave);
    return store_β, store_Σ
end



# types for output export
@with_kw struct VAROutput_Chan2020minn <: BVARmodelOutput
    store_β::Array{}      # 
    store_Σ::Array{}      # 
    YY::Array{}             #
end



#------------------------------
# Forecasting block
#------------------------------

function forecast(VAROutput::VAROutput_Chan2020minn,VARSetup)
    @unpack store_β, store_Σ, YY = VAROutput
    @unpack n_fcst,p,nsave = VARSetup
    n = size(YY,2);

    Yfor3D    = fill(NaN,(p+n_fcst,n,nsave))
    Yfor3D[1:p,:,:] .= @views YY[end-p+1:end,:];
    
    for i_draw = 1:nsave
        Yfor = @views Yfor3D[:,:,i_draw];
        A_draw = @views reshape(store_β[:,i_draw],n*p+1,n);
        Σ_draw = @views reshape(store_Σ[:,i_draw],n,n);
                
        for i_for = 1:n_fcst
            tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
            tclass = [1;tclass];
            Yfor[p+i_for,:]=tclass'*A_draw  .+ (cholesky(Σ_draw).U*randn(n,1))';    
        end
    end
    return Yfor3D

end # end function fcastChan2020minn()


# depreciated after moving to the syntax beavar(strcts)
#------------------------------
# dispatchModel block
#------------------------------
# function dispatchModel(::Chan2020minn_type,YY_tup, hyper_str, p,n_burn,n_save,n_irf,n_fcst)
#     println("Hello Minn")
#     intercept = 1;
#     if isa(YY_tup[1],Array{})
#         YY = YY_tup[1];
#     elseif isa(YY_tup[1],TimeArray{})
#         YY_TA = YY_tup[1];
#         YY = values(YY_TA)
#         varList = colnames(YY_TA)
#     end
#     set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
#     store_β, store_Σ = Chan2020minn(YY,set_strct,hyper_str);
#     out_strct = VAROutput_Chan2020minn(store_β,store_Σ,YY)
#     return out_strct, set_strct
# end





"""

"""
function initMinn(YY,p)
    Y, X, T, n, intercept       = mlagL(YY,p);
    k                           = n*p+intercept
    sigmaP                      = ar4!(YY,zeros(n,));                       # do OLS to initialize priors
    S_0                         = Diagonal(sigmaP);              
    Σt_inv                      = S_0\I;                                    # initialize Σ^-1              
    Vβminn_inv                  = 1.0*Matrix(I,n*k,n*k);                    # prior matrix
    Vβminn_inv_elview           = @view(Vβminn_inv[diagind(Vβminn_inv)]);   # will be used to update the diagonal    
    Σ_invsp, Σt_LI              = BEAVARs.makeBlkDiag(T*n,n,0,Σt_inv);      # I(T) ⊗ Σ-1 and its indices for update
    XtΣ_inv_den                 = zeros(k*n,T*n);                           # will be X' ( I(T) ⊗ Σ-1 )   from page 6 in Chan 2020 LBA
    XtΣ_inv_X                   = zeros(n*k,n*k);                           # will be X' ( I(T) ⊗ Σ-1 ) X from page 6 in Chan 2020 LBA    
    Xsur_den, Xsur_CI, X_CI     = BEAVARs.SUR_form_dense(X,n);              # prepares the SUR form and the indices of the parameters for updating
    K_β                         = zeros(n*k,n*k);                           # Variance covariance matrix of the parameters
    beta                        = zeros(n*k,);                              # the parameters in a vector
    
    return Y, X, T, n, sigmaP, S_0, Σt_inv, Vβminn_inv, Vβminn_inv_elview, Σ_invsp, Σt_LI, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept
end