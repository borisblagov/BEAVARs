"""
    irf_chol(beta_vec,sigma_vec,n,p,intercept,n_irf)

    Calculates impulse responses for a set of parameters in a vector beta_vec and variance-covariance matrix in a vector sigma_vec
    Uses the Cholesky decomposition and the companion form.
        - beta_vec needs to have the parameters for each equation stacked on top of each other
        - intercept: 
            - the function needs to know whether the intercept is at the left or the right to drop it correctly
            - use 1 if it is on the right, 0 if there is none, and -1 if it is on the left
"""
function irf_chol_opt(beta_vec,sigma_vec,n::Integer,p::Integer,intercept,n_irf::Integer,IRF_mat;shSize::String = "stdev")

k_nc = n*p;                     # number of perameters with no constant
# reshape the B matrix and drop the constant
B_draw = reshape(beta_vec,k_nc+1,n)'
if intercept == 1
    Bnc_draw = B_draw[:,1:k_nc]
end

# reshape the Sigma matrix
Σ_draw = reshape(sigma_vec,n,n)

F = [Bnc_draw; 1.0I(n*(p-1)) zeros(n*(p-1),n)]; # Companion form

A_chol = cholesky(Σ_draw);                  # structural identification via Cholesky

# Shock size (if nothing is selected it will default to one s.d.)
if shSize == "unity"        
    shocks_mat    = (1.0./diag(A_chol.L)).*1.0I(n);    
else 
    shocks_mat    = 1.0I(n);           
end
impulse_full = zeros(k_nc,);
resp_all = similar(impulse_full);
responses = zeros(n,n_irf);


# now cycle over the shocks
for i_var = 1:n
    responses[:,1] = A_chol.L*shocks_mat[:,i_var];
    impulse_full[1:n,1] = responses[:,1];

    for i_irf = 2:n_irf
        mul!(resp_all,F^(i_irf-1), impulse_full);
        responses[:,i_irf] = resp_all[1:n,]
    end

    IRF_mat[:,:,i_var] = transpose(responses)

end
    
return IRF_mat
end


function irf_chol1(beta_vec,sigma_vec,n,p,intercept,n_irf)

k_nc = n*p;
IRF_mat = zeros(n_irf,n,n);
for i_var = 1:n
    B_draw = reshape(beta_vec,k_nc+1,n)'
    if intercept == 1
        Bnc_draw = B_draw[:,1:k_nc]
    end

    Σ_draw = reshape(sigma_vec,n,n)

    F = [Bnc_draw; 1.0I(n*(p-1)) zeros(n*(p-1),n)]

    A_chol = cholesky(Σ_draw)

    responses = zeros(n,n_irf);
    shocks    = zeros(n,);     # shocks vector (either with unity shocks or std.dev shocks)

    shocks[i_var,] = 1.0

    responses[:,1] = A_chol.L*shocks;
    impulses = [responses[:,1]; zeros(n*(p-1))]

    for i_irf = 2:n_irf
        F_n = F^(i_irf-1)
        resp_all = F_n*impulses
        responses[:,i_irf] = resp_all[1:n]
    end

    IRF_mat[:,:,i_var] = transpose(responses)

end
    
return IRF_mat
end


function irf_chol_Bayesian(store_beta,store_sigma,n,p,intercept,n_irf,nsave)
    
    IRF_4d = zeros(n_irf,n,n,nsave);
    for i_draw = 1:nsave
        beta_vec=store_beta[:,i_draw];
        sigma_vec = store_sigma[:,i_draw];
    
        IRF_4d[:,:,:,i_draw] = irf_chol(beta_vec,sigma_vec,n,p,intercept,n_irf);
    end
    return IRF_4d
end


function irf_chol_Bayesian_opt(store_beta,store_sigma,n,p,intercept,n_irf,nsave;shSize = "stdev")
    IRF_4d = zeros(n_irf,n,n,nsave);
    IRF_mat = zeros(n_irf,n,n);
    for i_draw = 1:nsave
        beta_vec = store_beta[:,i_draw];
        sigma_vec = store_sigma[:,i_draw];
    
        IRF_4d[:,:,:,i_draw] = irf_chol_opt(beta_vec,sigma_vec,n,p,intercept,n_irf,IRF_mat,shSize=shSize);
    end
    return IRF_4d
end


function irf_chol_BayesianMulti(store_beta,store_sigma,n,p,intercept,n_irf,nsave,IRF_4d)
    a = 1:nsave;
    chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads());
    for i_draw in chunks
        Threads.@spawn begin
        beta_vec=store_beta[:,i_draw];
        sigma_vec = store_sigma[:,i_draw];
    
        IRF_4d[:,:,:,i_draw] = irf_chol(beta_vec,sigma_vec,n,p,intercept,n_irf);
        end
    end
    return IRF_4d
end