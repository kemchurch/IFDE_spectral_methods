

function [Y0,Y0_hat,Z0,Z0_hat,Z1,Z2,R] = Bounds_para_cont(a_0_not_int,a_1_not_int,beta_0_not_int,beta_1_not_int,a_0,a_s,N,alpha,beta_0,beta_s,u2,delta,Q,nu)

A_N = intval(inv(Derivative_DF(a_0_not_int,N,alpha,beta_0_not_int,u2,delta,Q)));
alpha = intval(alpha);
u2 = intval(u2);
nu = intval(nu);

a_0_ext = [reshape(a_0,N+1,Q+1);zeros(N+1,Q+1)];
norm_a_0 = norm_1_nu(reshape(a_0,N+1,Q+1),nu);
a_s_ext = [reshape(a_s,N+1,Q+1);zeros(N+1,Q+1)];
% Convolution
a_0_convo = intval(zeros(size(a_0_ext)));
a_s_convo = intval(zeros(size(a_s_ext)));
for i = 0:Q
    a_0_convo(:,i+1) = quadratic_cheb_intval(a_0_ext(:,i+1),a_0_ext(:,mod(i+delta,Q)+1));
    a_s_convo(:,i+1) = quadratic_cheb_intval(a_s_ext(:,i+1),a_s_ext(:,mod(i+delta,Q)+1));
end

F_0 = function_F_intval(a_0,N,alpha,beta_0,u2,delta,Q);
DF_0 = Derivative_DF_intval(a_0,N,alpha,beta_0,u2,delta,Q);
F_s = function_F_intval(a_s,N,alpha,beta_s,u2,delta,Q);
DF_s = Derivative_DF_intval(a_s,N,alpha,beta_s,u2,delta,Q);
DF_para = [DF_s, Derivative_D_beta(a_s,N,Q)];

% Bound Y0 and Y0_hat
delta_a = intval(a_1_not_int-a_0_not_int);
delta_beta = intval(beta_1_not_int-beta_0_not_int);
z_para = DF_para*[delta_a;delta_beta];
Y0_fin = intval(zeros(Q+1));
Y0_fin_hat = intval(zeros(Q+1));
Psi_0 = alpha/2/u2*(a_0_ext-a_0_convo);
Psi_s = alpha/2/u2*(a_s_ext-a_s_convo);
n = ones(1,2*N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
TPsi_0 = intval(zeros(size(Psi_0)));
TPsi_s = intval(zeros(size(Psi_s)));
for i = 0:Q
    for j=0:Q
        Y0_fin(i+1,j+1) = norm_1_nu(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1))*F_0(j*(N+1)+1:(j+1)*(N+1)),nu);
        Y0_fin_hat(i+1,j+1) = norm_1_nu(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1))*z_para(j*(N+1)+1:(j+1)*(N+1)),nu);
    end
    TPsi_0(:,i+1) = T*Psi_0(:,i+1);  
    TPsi_s(:,i+1) = T*Psi_s(:,i+1);  
end
TPsi_0(1:N+1,:) = 0;
TPsi_s(1:N+1,:) = 0;
Y0_inf = eye(Q+1).*norm_1_nu(TPsi_s,nu)/(2*(intval(N)+1));
Y0_inf_hat = eye(Q+1).*norm_1_nu(TPsi_s-TPsi_0,nu)/(2*(intval(N)+1));
Y0 = max(sup(sum(Y0_inf+Y0_fin,2)));
Y0_hat = max(sup(sum(Y0_inf_hat+Y0_fin_hat,2)));


% Bound Z0
Z0 = intval(zeros(Q+1));
I = intval(eye(size(DF_0)));
B = I-A_N*DF_0;
for i = 0:Q
    for j=0:Q
        Z0(i+1,j+1) = norm_nu_N(B(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N);
    end
end
Z0 = max(sup(sum(Z0,2)));

% Xi vector
a_0_reshape = reshape(a_0,N+1,Q+1);
a_s_reshape = reshape(a_s,N+1,Q+1);
Delta_a_reshape = a_s_reshape - a_0_reshape;
norm_delta_a = intval(zeros(1,Q+1));
Xi_a = intval(zeros(N+2,Q+1));
Xi_Delta_a = intval(zeros(N+2,Q+1));
for i = 0:Q
    norm_delta_a(i+1) = norm_1_nu(Delta_a_reshape(:,i+1),nu);
    Xi_a(1:N+1,i+1) = intval(xi_intval(a_0_reshape(:,i+1),nu));
    Xi_Delta_a(1:N+1,i+1) = intval(xi_intval(Delta_a_reshape (:,i+1),nu));
end
n = ones(1,N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
Txi = abs(T)*Xi_a;
Txi = Txi(1:N+1,1:Q+1);
Tdeltaxi = abs(T)*Xi_Delta_a;
Tdeltaxi = Tdeltaxi(1:N+1,1:Q+1);

% Bound Z1
Z1_head = intval(zeros(Q+1));
Z1_body = intval(zeros(Q+1));
Z1_tails = intval(zeros(Q+1));

for i = 0:Q
    Z1_tails(i+1,i+1) = alpha/2/u2*(nu+1/nu)/(2*(intval(N)+1))*(norm_a_0(i+1)+norm_a_0(mod(i+delta,Q+1)+1)+1);
    for j=0:Q
        Z1_head(i+1,j+1) = 2*norm_1_nu(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1),nu)/nu^(N+1);
        Z1_body(i+1,j+1) = alpha/2/u2*norm_1_nu(abs(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)))*(Txi(:,i+1)+Txi(:,mod(i+delta,Q+1)+1)),nu);
    end
end
Z1 = max(sup(sum(Z1_head + Z1_body + Z1_tails,2)));

% Bound Z2
norm_A = intval(zeros(Q+1));
norm_DPsi = intval(zeros(Q+1));
z_1_2 = intval(zeros(1,Q+1));
for i = 0:Q
    z_1_2(i+1) = norm_delta_a(i+1)+norm_delta_a(mod(i+delta,Q+1)+1);
    for j=0:Q
        if i==j
            norm_A(i+1,j+1) = max([norm_nu_N(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N),1/(2*(intval(N)))]);
            norm_DPsi(i+1,j+1) = abs(alpha/u2);
        else
            norm_A(i+1,j+1) = norm_nu_N(A_N(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N);
        end
        if j == mod(i + delta,Q+1) 
            norm_DPsi(i+1,j+1) = abs(alpha/u2);
        end
    end
end
z_1_2 = sup(sum(z_1_2));
norm_T = 1/nu+2*nu;
Z0_hat = max(sup(sum(norm_T*norm_A.*norm_DPsi,2)))*z_1_2;
Z2 = max(sup(sum(norm_T*norm_A.*norm_DPsi,2)));

p = [Z2,-(1-Z1-Z0-Z0_hat),Y0+Y0_hat];
R = sort(roots(p));

end

