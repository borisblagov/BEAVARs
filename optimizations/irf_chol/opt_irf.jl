include("devPkgs.jl")
using DelimitedFiles
using LinearAlgebra

p = 2;
n = 3;

intercept = 1;
n_irf = 16;
beta_vec = readdlm("bvec.txt");
sigma_vec = readdlm("svec.txt");

includet("optimizations/irfs.jl")
IRF_mat = zeros(n_irf,n,n);
one_irf1 = irf_chol1(beta_vec,sigma_vec,n,p,intercept,n_irf);
one_irf2 = irf_chol_opt(beta_vec,sigma_vec,n,p,intercept,n_irf,IRF_mat,shSize="unity")
one_irf1==one_irf2


# @time one_irf1 = irf_chol1(beta_vec,sigma_vec,n,p,intercept,n_irf);
# @time one_irf2 = irf_chol_opt(beta_vec,sigma_vec,n,p,intercept,n_irf,IRF_mat,shSize="unity");

# @btime one_irf1 = irf_chol1(beta_vec,sigma_vec,n,p,intercept,n_irf);
# @btime one_irf2 = irf_chol_opt(beta_vec,sigma_vec,n,p,intercept,n_irf,IRF_mat);
