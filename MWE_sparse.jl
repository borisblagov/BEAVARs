using LinearAlgebra
using SparseArrays
using Revise
using BenchmarkTools


# -------------------------------------------
# this code simply generates the necessary matrices
Tf = 242; p = 5; T=Tf-p; n = 2;
X = rand(T,n*p+1);
repX = kron(X,ones(n,1));
r,c = size(X);
idi = repeat(1:r*n,inner=c);
idj=repeat(1:c*n,r);
Σhatt_inv_sp = randn(n,n);
# Σhatt_inv_sp = 1.0I(n)
Σ_invsp_full = kron(sparse(Matrix(1.0I, Tf, Tf)),Σhatt_inv_sp);
# -------------------------------------------



# main matrix No. 1
Xsur_sp = sparse(idi,idj,vec(repX'));

# main matrix No. 2
Σ_invsp = Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];

# alternative as a view
Σ_invsp_view = @view Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];

# the line in question using sparse matrices
XtΣ_inv = Xsur_sp'*Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];  

# timing 
print("n = 2, sparse array impementation:")
@btime $XtΣ_inv = $Xsur_sp'*$Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];


# -------------------------------------------
# alternative with dense matrices
# -------------------------------------------
XsurD = Matrix(Xsur_sp);
Σ_inv_full = Matrix(Σ_invsp_full);
XtΣ_invD = Matrix(XtΣ_inv);
Σ_invD   = Matrix(Σ_invsp);
Σ_inv_viewD = @view Σ_inv_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];


# the line in question with dense matrices
mul!(XtΣ_invD,XsurD',Σ_invD);

# timing 
print("n = 2, dense matrix impementation:")
@btime mul!($XtΣ_invD,$XsurD',$Σ_invD);
# timing with the view
print("n = 2, dense matrix with view:")
@btime mul!($XtΣ_invD,$XsurD',$Σ_inv_viewD);


# -------------------------------------------
# Side note: mul! with sparse matrices just hangs forever for n=17

print("n = 2, timing with sparse view:")
@btime $XtΣ_inv = $Xsur_sp'*$Σ_invsp_view;

print("n = 2, timing with sparse mul!:")
@btime mul!($XtΣ_inv,$Xsur_sp',$Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end]);





# n = 17



# -------------------------------------------
# this code simply generates the necessary matrices
Tf = 242; p = 5; T=Tf-p; n = 17;
X = rand(T,n*p+1);
repX = kron(X,ones(n,1));
r,c = size(X);
idi = repeat(1:r*n,inner=c);
idj=repeat(1:c*n,r);
Σhatt_inv_sp = randn(n,n);
Σ_invsp_full = kron(sparse(Matrix(1.0I, Tf, Tf)),Σhatt_inv_sp);
# -------------------------------------------



# main matrix No. 1
Xsur_sp = sparse(idi,idj,vec(repX'));

# main matrix No. 2
Σ_invsp = Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];

# alternative as a view
Σ_invsp_view = @view Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];

# the line in question using sparse matrices
XtΣ_inv = Xsur_sp'*Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];  

# timing 
print("n = 17, sparse array impementation:")
@btime $XtΣ_inv = $Xsur_sp'*$Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];


# -------------------------------------------
# alternative with dense matrices
# -------------------------------------------
XsurD = Matrix(Xsur_sp);
Σ_inv_full = Matrix(Σ_invsp_full);
XtΣ_invD = Matrix(XtΣ_inv);
Σ_invD   = Matrix(Σ_invsp);
Σ_inv_viewD = @view Σ_inv_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end];


# the line in question with dense matrices
mul!(XtΣ_invD,XsurD',Σ_invD);

# timing 
println("n = 17, dense matrix impementation:")
@btime mul!($XtΣ_invD,$XsurD',$Σ_invD);

# timing with the view
print("n = 17, dense matrix with view:")
@btime mul!($XtΣ_invD,$XsurD',$Σ_inv_viewD);


# Side note: mul! with sparse matrices just hangs forever for n=17

print("n = 17, timing with sparse view never completes on my machine\n")
# @btime $XtΣ_inv = $Xsur_sp'*$Σ_invsp_view;

print("n = 17, timing with sparse mul!() never completes on my machine")
# @btime mul!($XtΣ_inv,$Xsur_sp',$Σ_invsp_full[(Tf-T)*n+1:end, (Tf-T)*n+1:end]);
