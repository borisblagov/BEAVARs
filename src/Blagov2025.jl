function makeHypSetup(::Blagov2025_type)
    return hypChan2020()
end


@doc raw"""
    Updates parameters using an independennt Normal-Wishart prior
"""
function CPZcsv!(YY,p,hypSetup,n,k,b0,B_draw,Σt_inv,structB_draw,Y,X,T,beta,h,dg_ind_Ωinv)
 
end




@doc raw"""
     
"""
function Blagov2025(dataHF_tab,dataLF_tab,varList,varSetup,hypSetup,trans)
    @unpack p, nburn,nsave, const_loc = varSetup
    ndraws = nsave+nburn;
    nmdraws = 10;               # given a draw from the parameters to draw multiple time from the distribution of the missing data for better confidence intervals

    fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TimeArrays(dataLF_tab,dataHF_tab,varList,trans)

    YYwNA = values(fdataHF_tab);
    YY = deepcopy(YYwNA);
    Tf,n = size(YY);
    
    B_draw, structB_draw, Σt_inv, b0 = BEAVARs.initParamMatrices(n,p,const_loc) 

    YYt, Y0, longyo, nm, H_B, H_B_CI, strctBdraw_LI, Σ_invsp, Σt_LI, Σp_invsp, Σpt_ind, Xb, cB, cB_b0_LI, Smsp, Sosp, Sm_bit, Gm, Go, GΣ, Kym = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);
    
    M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);

    # YY has missing values so we need to draw them once to be able to initialize matrices and prior values
    YYt = BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI,nmdraws);
    
    
    # Initialize matrices for updating the parameter draws from CPZ_iniv  
    # ------------------------------------
    Y, X, T, deltaP, sigmaP, mu_prior, V_Minn_inv, V_Minn_inv_elview, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept = CPZ_initMinn(YY,p)
    

    # prepare matrices for storage
    store_YY    = zeros(Tf,n,nsave);
    # store_β     = zeros(n^2*p+n,nsave);
    # store_Σt_inv= zeros(n,n,nsave);
    # store_Σt    = zeros(n,n,nsave);

    @unpack ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ = hypSetup
    @unpack p, nsave, nburn = VARSetup

    Y, X, T, n = mlagL(YY,p)
    (deltaP, sigmaP, mu_prior) = trainPriors(YY,p)
    np1 = n*p+1; # number of parameters per equation
    
    (idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_NonConj(n,p,sigmaP,hypSetup);
    
    A_0 = reshape(beta_Minn,np1,n);
    V_Ainv = sparse(1:np1,1:np1,1.0./V_Minn)
    
    S_0 = Diagonal(sigmaP);
    Σ = S_0;
    v_0 = n+3;
    h = zeros(T,)
    H_ρ = sparse(Matrix(1.0I, T, T)) - sparse(ρ*diagm(-1=>repeat(1:1,T-1)));
    # dv = ones(T,); ev = -ρ.*ones(T-1,);
    # H_ρ = Bidiagonal(dv,ev,:L)
    
    # This part follows page 19 of Chan, J. (2020)
    ndraws = nsave+nburn;
    store_β = zeros(np1*n,nsave);
    store_h = zeros(T,nsave);
    store_Σ = zeros(n,n,nsave);
    store_s2_h = zeros(T,nsave);
    ρ_store = zeros(nsave,);
    σ_h2_store = zeros(nsave,); 
    eh_store = zeros(T,nsave);
    
    eh = similar(h);

    Ωinv = sparse(1:T,1:T,exp.(-h));
    dg_ind_Ωinv = diagind(Ωinv);

    @showprogress for ii in 1:ndraws
        # draw of the missing values
        BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI,nmdraws);
        
        # draw of the parameters
        Y, X = mlagL!(YY,Y,X,p,n)
        Ωinv[dg_ind_Ωinv]=exp.(-h);
        ZtΩinv = X'*Ωinv;
        
        K_A = V_Ainv + ZtΩinv*X;
        A_hat = K_A\(V_Ainv\A_0 + ZtΩinv*Y);
        S_hat = S_0 + A_0'*V_Ainv*A_0 + Y'*Ωinv*Y - A_hat'*K_A*A_hat;
        S_hat = (S_hat+S_hat')/2;
    
        Σ = rand(InverseWishart(hypSetup.nu0+T,S_hat));
        cholΣ  = cholesky(Σ).U; # if we get the upper we don't need constant transpose
        A = A_hat + (cholesky(Hermitian(K_A)).U\randn(np1,n))*cholΣ;
    
        # Errors
        U = Y - X*A
        s2_h = sum((U/cholΣ).^2,dims=2)
    
        BEAVARs.draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ)
        eh[1,] = h[1]*sqrt(1-ρ^2);
        eh[2:end,] = h[2:end,].-ρ.*h[1:end-1,]
        σ_h2 = 1.0./rand(Gamma(v_h0+T/2.0, 1/( S_h0 + sum(eh.^2)./2.0 ) ) )
    
        K_rho = 1.0/V_ρ + sum(h[1:T-1,].^2)./σ_h2;
        ρ_hat = K_rho\(ρ_0/V_ρ + h[1:T-1,]'*h[2:T,]./σ_h2);
        ρ_c = ρ_hat + sqrt(K_rho)\randn();
    
        g_ρ(x) = -0.5*log(σ_h2./(1-x.^2))-0.5*(1.0-x.^2)/σ_h2*h[1]^2;
        if abs(ρ_c)<.999
            alpMH = exp(g_ρ(ρ_c)-g_ρ(ρ));
            if alpMH>rand()
                ρ = ρ_c;            
            end
        end    
        H_ρ[diagind(H_ρ,-1)]=fill(-ρ,T-1);


        if ii>nburn
            store_β[:,ii-nburn]  = beta;
            store_YY[:,:,ii-nburn]  = YY;
            # store_Σt_inv[:,:,ii-nburn]    = Σt_inv;
            # store_Σt[:,:,ii-nburn] = Σt;
            store_β[:,ii-nburn] = vec(A);
            store_h[:,ii-nburn] = h;
            store_Σ[:,:,ii-nburn] = Σ;
            store_s2_h[:,ii-nburn] = s2_h;
            ρ_store[ii-nburn,] = ρ;
            σ_h2_store[ii-nburn,] = σ_h2;
            eh_store[:,ii-nburn] = eh;
        end
    end

    return store_YY,store_β, store_Σt_inv, M_zsp, z_vec, Sm_bit, store_Σt
end





#------------------------------
# Output structure
@with_kw struct VAROutput_Blagov2025 <: modelOutput
    store_β::Array{}        # 
    store_Σt_inv::Array{}        # 
    store_YY::Array{}
    M_zsp::Array{} 
    z_vec::Array{} 
    Sm_bit::Array{}
    store_Σt::Array{}        # 
end
# end of output strcutres
#------------------------------

#--------------------------------------
# Forecast CPZ2024
function forecast(VAROutput::VAROutput_Blagov2025,VARSetup)
    @unpack store_β, store_Σt, store_YY = VAROutput
    @unpack n_fcst,p,nsave = VARSetup

    YY = median(store_YY,dims=3)

    n = size(YY,2);

    Yfor3D    = fill(NaN,(p+n_fcst,n,nsave))
    Yfor3D[1:p,:,:] .= @views YY[end-p+1:end,:];
    
    for i_draw = 1:nsave
        # Yfor3D[1:p,:,i_draw] .= @views store_YY[end-p+1:end,:,i_draw];
        Yfor3D[1:p,:,i_draw] .= @views YY[end-p+1:end,:];
        Yfor = @views Yfor3D[:,:,i_draw];
        A_draw = @views reshape(store_β[:,i_draw],n*p+1,n);
        Σ_draw = @views store_Σt[:,:,i_draw];
                
        for i_for = 1:n_fcst
            tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
            tclass = [1;tclass];
            Yfor[p+i_for,:]=tclass'*A_draw  .+ (cholesky(Hermitian(Σ_draw)).U*randn(n,1))';    
        end
    end
    return Yfor3D

end # end function fcastCPZ2024()




function dispatchModel(::Blagov2025_type,YY_tup, hyp_strct, p,n_burn,n_save,n_irf,n_fcst)
    println("Hello CPZ2024")
    intercept = 1;
    dataHF_tab  = YY_tup[1]
    dataLF_tab  = YY_tup[2]
    varList     = YY_tup[3]
    trans       = YY_tup[4] # transformation of the LF variables (0: growth rates or 1: log-levels)
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,intercept);
    store_YY,store_β, store_Σt_inv, M_zsp, z_vec, Sm_bit,store_Σt = Blagov2025(dataHF_tab,dataLF_tab,varList,set_strct,hyp_strct,trans)    
    out_strct = VAROutput_Blagov2025(store_β,store_Σt_inv,store_YY,M_zsp, z_vec, Sm_bit,store_Σt)
    return out_strct, set_strct
end
