
@doc raw"""
# Chan2020iniw(YY,VARSetup,hypSetup)

Implements the classic homoscedastic Minnesota prior with a SUR form following Chan (2020)

"""
function Chan2020iniw(YY,VARSetup::modelSetup,hypSetup::modelHypSetup)
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
        mul!(XtΣ_inv_den,Xsur_den',Σ_invsp);            #  X'*( I(T) ⊗ Σ^{-1} )
        mul!(XtΣ_inv_X,XtΣ_inv_den,Xsur_den);           #  X'*( I(T) ⊗ Σ^{-1} )*X
        K_β[:,:] .= Vβminn_inv .+ XtΣ_inv_X;            #  K_β = V^{-1} + X'*( I(T) ⊗ Σ^{-1} )*X
        prior_mean = Vβminn_inv*βminn;                  #  V^-1 * βminn 
        mul!(prior_mean,XtΣ_inv_den, vec(Y'),1.0,1.0);  # (V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y
        cholK_β = cholesky(Hermitian(K_β));             # C is lower triangular, C' is upper triangular
        beta_hat = ldiv!(cholK_β.U,ldiv!(cholK_β.L,prior_mean));    # C'\(C*(V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y)
        beta = beta_hat + ldiv!(cholK_β.U,randn(k*n,)); # draw for β

        U = reshape(vec(Y') - Xsur_den*beta,n,T);       
        Σt = rand(InverseWishart(hypSetup.nu0+n+T,S_0+U*U'));    # draw for Σ
        Σt_inv = Σt\I;
        Σ_invsp.nzval[:] = Σt_inv[Σt_LI];               # update ( I(T) ⊗ Σ^{-1} )

        if ii>nburn
            store_β[:,ii-nburn] = beta;
            store_Σt[:,:,ii-nburn] = Σt;
        end
    end

    return store_β, store_Σt
end

function Chan2020iniw2(YY,VARSetup::modelSetup,hypSetup::modelHypSetup)
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
        beta = BEAVARs.Chan2020_drawβ(Σ_invsp,Xsur_den,XtΣ_inv_den,XtΣ_inv_X,Vβminn_inv,K_β,Y);
        Σt, Σt_inv = Chan2020_drawΣt(Y,Xsur_den,beta);

        Σ_invsp.nzval[:] = Σt_inv[Σt_LI];               # update ( I(T) ⊗ Σ^{-1} )

        if ii>nburn
            store_β[:,ii-nburn] = beta;
            store_Σt[:,:,ii-nburn] = Σt;
        end
    end

    return store_β, store_Σt
end

"""
    
"""
function Chan2020_drawβ(Σ_invsp,Xsur_den,XtΣ_inv_den,XtΣ_inv_X,Vβminn_inv,K_β,Y)
        mul!(XtΣ_inv_den,Xsur_den',Σ_invsp);            #  X'*( I(T) ⊗ Σ^{-1} )
        mul!(XtΣ_inv_X,XtΣ_inv_den,Xsur_den);           #  X'*( I(T) ⊗ Σ^{-1} )*X
        K_β[:,:] .= Vβminn_inv .+ XtΣ_inv_X;            #  K_β = V^{-1} + X'*( I(T) ⊗ Σ^{-1} )*X
        prior_mean = Vβminn_inv*βminn;                  #  V^-1 * βminn 
        mul!(prior_mean,XtΣ_inv_den, vec(Y'),1.0,1.0);  # (V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y
        cholK_β = cholesky(Hermitian(K_β));             # C is lower triangular, C' is upper triangular
        beta_hat = ldiv!(cholK_β.U,ldiv!(cholK_β.L,prior_mean));    # C'\(C*(V^-1_Minn * beta_Minn + X' ( I(T) ⊗ Σ-1 ) y)
        beta = beta_hat + ldiv!(cholK_β.U,randn(k*n,)); # draw for β
        return beta
end

"""
"""
function Chan2020_drawΣt(Y,Xsur_den,beta)
    U = reshape(vec(Y') - Xsur_den*beta,n,T);       
    Σt = rand(InverseWishart(hypSetup.nu0+n+T,S_0+U*U'));    # draw for Σ
    Σt_inv = Σt\I;
    return Σt, Σt_inv
end



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