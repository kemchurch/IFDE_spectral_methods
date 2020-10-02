function [Y0,Z0,Z1,Z2,R] = Bounds(a,N,alpha,beta,u2,delta,Q,nu)
A = intval(inv(Derivative_DF(a,N,alpha,beta,u2,delta,Q)));

a = intval(a);
alpha = intval(alpha);
beta  = intval(beta);
u2 = intval(u2);
nu = intval(nu);

a_ext = [reshape(a,N+1,Q+1);zeros(N+1,Q+1)];
norm_a = norm_1_nu(reshape(a,N+1,Q+1),nu);
% Convolution
a_convo = intval(zeros(size(a_ext)));
for i = 0:Q
    a_convo(:,i+1) = quadratic_cheb_intval(a_ext(:,i+1),a_ext(:,mod(i+delta,Q)+1));
end

F = function_F_intval(a,N,alpha,beta,u2,delta,Q);
DF = Derivative_DF_intval(a,N,alpha,beta,u2,delta,Q);

% Bound Y0
Y0_fin = intval(zeros(Q+1));
Psi = alpha/2/u2*(a_ext-a_convo);
n = ones(1,2*N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
TPsi = intval(zeros(size(Psi)));
for i = 0:Q
    for j=0:Q
        Y0_fin(i+1,j+1) = norm_1_nu(A(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1))*F(j*(N+1)+1:(j+1)*(N+1)),nu);
    end
    TPsi(:,i+1) = T*Psi(:,i+1);  
end
TPsi(1:N+1,:) = 0;
Y0_inf = eye(Q+1).*norm_1_nu(TPsi,nu)/(2*(intval(N)+1));
Y0 = max(sup(sum(Y0_inf+Y0_fin,2)));

% Bound Z0
Z0 = intval(zeros(Q+1));
I = intval(eye(size(DF)));
B = I-A*DF;
for i = 0:Q
    for j=0:Q
        Z0(i+1,j+1) = norm_nu_N(B(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N);
    end
end
Z0 = max(sup(sum(Z0,2)));

% Bound Z1
Z1_head = intval(zeros(Q+1));
Z1_body = intval(zeros(Q+1));
Z1_tails = intval(zeros(Q+1));
beta_j = intval(ones(Q+1,1));
beta_j(Q+1) = beta;

% Xi vector
a_reshape = reshape(a,N+1,Q+1);
Xi_a = intval(zeros(N+2,Q+1));
for i = 0:Q
    Xi_a(1:N+1,i+1) = intval(xi_intval(a_reshape(:,i+1),nu));
end
n = ones(1,N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
Txi = abs(T)*Xi_a;
Txi = Txi(1:N+1,1:Q+1);

for i = 0:Q
    Z1_tails(i+1,i+1) = alpha/2/u2*(nu+1/nu)/(2*(intval(N)+1))*(norm_a(i+1)+norm_a(mod(i+delta,Q+1)+1)+1);
    for j=0:Q
        Z1_head(i+1,j+1) = 2*norm_1_nu(A(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1),nu)/nu^(N+1);
        Z1_body(i+1,j+1) = alpha/2/u2*norm_1_nu(abs(A(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)))*(Txi(:,i+1)+Txi(:,mod(i+delta,Q+1)+1)),nu);
    end
end
Z1 = max(sup(sum(Z1_head + Z1_body + Z1_tails,2)));

% Bound Z2
norm_A = intval(zeros(Q+1));
norm_DPsi = intval(zeros(Q+1));
for i = 0:Q
    for j=0:Q
        if i==j
            norm_A(i+1,j+1) = max([norm_nu_N(A(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N),1/(2*(intval(N)))]);
            norm_DPsi(i+1,j+1) = abs(alpha/u2);
        else
            norm_A(i+1,j+1) = norm_nu_N(A(i*(N+1)+1:(i+1)*(N+1),j*(N+1)+1:(j+1)*(N+1)),nu,N);
        end
        if j == mod(i + delta,Q+1) 
            norm_DPsi(i+1,j+1) = abs(alpha/u2);
        end
    end
end
norm_T = 1/nu+2*nu;
Z2 = max(sup(sum(norm_T*norm_A.*norm_DPsi,2)));

p = [Z2,-(1-Z1-Z0),Y0];
R = sort(roots(p));
end

