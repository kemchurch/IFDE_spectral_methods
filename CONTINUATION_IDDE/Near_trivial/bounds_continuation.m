function [Y0sup,Z0sup,Z1sup,Z2_times_r_sup,R,pint,psup] = bounds_continuation(a0,a1,eps0,eps1,beta0,beta1,N,alpha,u2,delta,u1,nu)
Q = u1-1;
nu = intval(nu);
alpha = intval(alpha);
a0 = intval(a0);
a1 = intval(a1);
beta0 = intval(beta0);
beta1 = intval(beta1);
eps0 = intval(eps0);
eps1 = intval(eps1);
A_dagger = intval(D_scriptF_int(a0,eps0,N,alpha,beta0,u2,delta,u1));
Abar = inv(A_dagger);
abars = a0 + infsup(0,1)*(a1-a0);
epsbars = eps0 + infsup(0,1)*(eps1-eps0);
betabars = beta0 + infsup(0,1)*(beta1-beta0);
Fzero = scriptF_int(a0,eps0,N,alpha,beta0,u2,delta,u1);

abars_q = reshape(abars,[N+1,Q+1]);
a0_q = reshape(a0,[N+1,Q+1]);
na0_q = intval(zeros(1,u1));
nabars_q = intval(zeros(1,u1));
a0_q_ext = intval([a0_q;zeros(N+1,Q+1)]);
abars_q_ext = [abars_q;intval(zeros(N+1,Q+1))];
abars_qa_qdelta = intval(zeros(2*N+2,Q+1));
a0_qa_qdelta = intval(zeros(2*N+2,Q+1));
na0_qa_qdelta = intval(zeros(1,Q+1));
nabars_qa_qdelta = intval(zeros(1,Q+1));
for n = 0:Q
    a0_qa_qdelta(:,n+1) = quadratic_cheb_intval(a0_q_ext(:,n+1),a0_q_ext(:,mod(n+delta,Q+1)+1));
    abars_qa_qdelta(:,n+1) = quadratic_cheb_intval(abars_q_ext(:,n+1),abars_q_ext(:,mod(n+delta,Q+1)+1));
    na0_qa_qdelta(n+1) = norm_1_nu(intval(a0_qa_qdelta(:,n+1)),nu);
    nabars_qa_qdelta(:,n+1) = norm_1_nu(abars_qa_qdelta(:,n+1),nu);
    na0_q(n+1) = norm_1_nu(intval(a0_q(:,n+1)),nu);
    nabars_q(n+1) = norm_1_nu(abars_q(:,n+1),nu);
end
n = ones(1,2*N+1);
nn = 1:1:2*N+1;
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;

% Precomputation: Y0
gammabar = alpha/2/u2*T*(abars_q_ext-eps0*abars_qa_qdelta);
gammabar_tail = gammabar;
gammabar_tail(1:N,:)=0;
Linv = intval(zeros(2*N+2));
Linv(2:end,2:end) = diag(nu.^(nn)/2./nn);
gammabar_tail = Linv*gammabar_tail;
% Bound Y0
AFzero = Abar*Fzero;
DF_bars = function_D_xbeta_F_int(abars,epsbars,N,alpha,betabars,u2,delta,u1);
ADF_barsDelta_xbeta = DF_bars*[a1-a0;eps1-eps0;beta1-beta0];
Y0_10 = 0;
Y0_1_MVT = 0;
for i=0:u1-1
    Y0_10 = max(Y0_10,norm_1_nu(AFzero(1+i*(N+1):(i+1)*(N+1)),nu));
    Y0_1_MVT = max(Y0_1_MVT,norm_1_nu(ADF_barsDelta_xbeta(1+i*(N+1):(i+1)*(N+1)),nu));
end
Y0_10 = max(Y0_10,abs(AFzero(end)));
Y0_1_MVT = max(Y0_1_MVT,abs(ADF_barsDelta_xbeta(end)));
Y02 = max(sum(abs(gammabar_tail),1));
Y0 = Y0_10 + Y0_1_MVT + Y02;
Y0sup = sup(Y0);

% Bound Z0
Z0 = norm_nu_1_XR_mat_int(intval(eye((N+1)*u1+1)) - Abar*A_dagger,N,u1,nu);
Z0sup = sup(Z0);

% Precomputation: Z1
nT = 2*nu + 1/nu;
I = intval(eye((N+1)));
w_bar = infsup(-1,1)*[2*max(1,beta0)*repmat(I(:,1),[u1,1]);1];
h_bar = infsup(-1,1)*alpha*nT/4/u2*[repmat(I(:,N+1),[u1,1]);0];
Ahbar = Abar*h_bar;
Awbar = Abar*w_bar;
BAXX = compute_bold_B(Abar,N,u1,nu);
BAXR = compute_bold_B_star(Abar,N,u1,nu);
aconv_tail = intval(zeros(u1,1));
for n=0:u1-1
    aconv_tail(n+1) = norm_1_nu(convcheb_tailbound(a0_q(:,n+1),nu),nu);
end

% Bound Z1
Z1_inf = nT/(4*(N+1))*alpha/u2*(1+eps0*2*max(na0_q) + max(na0_qa_qdelta));
Z1_bars = 1/nu^(N+1)*(max([norm_1_nu(reshape(Ahbar(1:u1*(N+1)),[N+1,u1]),nu),abs(Ahbar(end))])...
    + max([norm_1_nu(reshape(Awbar(1:u1*(N+1)),[N+1,u1]),nu),abs(Awbar(end))]));
Z1_bold_XX = BAXX*aconv_tail;
Z1_bold_XR = BAXR'*aconv_tail;
Z1 = Z1_inf + Z1_bars + alpha*eps0/2/u2*max([Z1_bold_XX;Z1_bold_XR]);
Z1sup = sup(Z1);

% Precompuitation: Z2
mu_bar = intval([zeros((u1-1)*(N+1),1);1;zeros(N+1,1)]);
Z_2q_2 = intval(zeros(u1,3));
Z_2q_3 = intval(zeros(u1,3));
Z_2_bold = intval(zeros(u1,3));
da = a1-a0;
da_q = reshape(da,[N+1,Q+1]);
for q=0:u1-1
   Z_2q_2(q+1,1) = alpha/2/u2;
   Z_2q_2(q+1,2) = alpha/2/u2*(abs(epsbars) + nabars_q(q+1));
   Z_2q_2(q+1,3) = alpha/2/u2*(abs(eps0)*da_q(q+1) + ...
       abs(eps1-eps0)*nabars_q(q+1));
   Z_2q_3(q+1,1) = alpha/2/u2;
   Z_2q_3(q+1,2) = alpha/2/u2*2*nabars_q(q+1);
   Z_2q_3(q+1,3) = alpha/2/u2*norm_1_nu(...
       quadratic_cheb_intval(da_q(:,q+1),a0_q(:,mod(q+delta,Q+1)+1)) + ...
       quadratic_cheb_intval(da_q(:,mod(q+delta,Q+1)+1),a0_q(:,q+1)) + ...
       quadratic_cheb_intval(da_q(:,mod(q+delta,Q+1)+1),da_q(:,q+1)), nu);
end
for q=0:u1-1
    Z_2_bold(q+1,:) = nT*(Z_2q_2(q+1,:) + Z_2q_2(mod(q+delta,Q+1)+1,:)...
        + Z_2q_3(q+1,:) );
end
% nAbar = norm_nu_1_XR_mat_int(Abar,N,u1,nu);
% nA = max(nAbar,1/(2*(N+1)));
% dbeta = abs(beta1-beta0);
% deps = abs(eps1-eps0);
% da = a1-a0;
% da_q = reshape(da,[N+1,Q+1]);
% norm_a = max(na0_q);
% n_da_q = intval(zeros(1,u1));
% for n = 0:Q
%     n_da_q(n+1) = norm_1_nu(da_q(:,n+1),nu);
% end
% norm_da = max(n_da_q);
% Bound Z2
Amubar = Abar*mu_bar;
n_Amubar = max([norm_1_nu(reshape(Amubar(1:u1*(N+1)),[N+1,u1]),nu),abs(Amubar(end))]);
Z2_X = 1/(2*intval(N+1))*Z_2_bold + BAXX*Z_2_bold;
Z2_R = BAXR'*Z_2_bold;
Z2 = max([Z2_X;Z2_R]) + [0,0,abs(beta1-beta0)*n_Amubar];
Z2_sup = sup(Z2);
Z2_times_r = [Z2,0];
Z2_times_r_sup = sup(Z2_times_r);

% R
pint = Z2_times_r + [0,0,-(1-Z1-Z0),Y0];
psup = Z2_times_r_sup + [0,0,-(1-Z1sup-Z0sup),Y0sup];
R = sort(roots(psup));
end