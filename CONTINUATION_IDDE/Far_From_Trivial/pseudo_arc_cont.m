function [] = pseudo_arc_cont(N,alpha,beta,u2,delta,Q,a,r,T,tau,delta_s,nu)
j = 1;
while beta > exp(-r*T)
    a_q = reshape(a,[N+1,Q+1]);
    D_beta_F_q = zeros(size(a_q));
    D_beta_F_q(1,Q+1) = sum([a_q(1,Q+1);2*a_q(2:end,Q+1)]);
    D_beta_F = reshape(D_beta_F_q,[],1);
    
    DF_D_beta_F = [D_beta_F,Derivative_DF(a,N,alpha,beta,u2,delta,Q)];
    
    X0_dot = null(DF_D_beta_F);
    if X0_dot(1) > 0
        X0_dot = -X0_dot;
    end
    X0 = [beta;a];
    X1_hat = X0 + delta_s*X0_dot;
    
    E = @(a,beta) ([beta;a]-X1_hat)'*X0_dot;
    
    F = @(a,beta) [E(a,beta);function_F(a,N,alpha,beta,u2,delta,Q)];
    
    DF = @(a,beta) [ transpose(X0_dot) ;D_beta_F,Derivative_DF(a,N,alpha,beta,u2,delta,Q)];
    
    Fn = F(a,beta);
    Xn = X1_hat;
    beta_n = Xn(1);
    a_n = Xn(2:end);
    k=0;
    while norm(Fn) > 10^(-14) && k < 25
        k = k+1;
        Xn = Xn - DF(a_n,beta_n)\Fn;
        beta_n = Xn(1);
        a_n = Xn(2:end);
        Fn = F(a_n,beta_n);
    end
    if k == 25 || norm(a) < 10^(-12)
        break
    end
    a = a_n;
    beta = beta_n;
    A(:,j) = a;
    B(:,j) = beta;
    [Y0,Z0,Z1,Z2,R] = Bounds(a,N,alpha,beta,u2,delta,Q,nu);
    Y0_MAT(:,j) = Y0;
    Z0_MAT(:,j) = Z0;
    Z1_MAT(:,j) = Z1;
    Z2_MAT(:,j) = Z2;
    R_MAT(:,j) = R;
    j = j+1;
end

A = A(:,1:end);
B = B(:,1:end);
Y0_MAT = Y0_MAT(:,1:end);
Z0_MAT = Z0_MAT(:,1:end);
Z1_MAT = Z1_MAT(:,1:end);
Z2_MAT = Z2_MAT(:,1:end);
R_MAT = R_MAT(:,1:end);

for i = 1:length(B)
    if R_MAT(1,i) > 0 && Z1_MAT(1,i) < 1
        plot(B(i),abs(A(1,i)),'r*');hold on;
    else
        plot(B(i),abs(A(1,i)),'b*');hold on;
    end
end
xlabel('\beta');
ylabel('|a_0|');
h = beta-1;
% Proven Solution
save('data_h0_proven','A','B','Y0_MAT','Z0_MAT','Z1_MAT','Z2_MAT','R_MAT','N','Q','u2','delta','alpha','r','tau','T','h','nu')
% Begining of next step
save('data_h1','a','N','Q','u2','delta','alpha','beta','r','tau','T','h')
end

