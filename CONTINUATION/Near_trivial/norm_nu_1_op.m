function nrm = norm_nu_1_op(A,nu)
N = length(A);
omega = [1,2*nu.^(1:1:N-1)];
OMEGA = repmat(omega.',[1,N]);
nrm = max(1./omega.*sum(abs(A).*OMEGA,1));
end