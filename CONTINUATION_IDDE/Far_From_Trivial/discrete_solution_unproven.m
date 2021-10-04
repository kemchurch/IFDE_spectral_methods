function [] = discrete_solution_unproven(N,alpha,beta,u2,delta,Q,a,r,T,tau,nu,delta_s)
j = 1;
while beta-delta_s > exp(-r*T)
    a_q = reshape(a,[N+1,Q+1]);
    D_beta_F_q = zeros(size(a_q));
    D_beta_F_q(1,Q+1) = sum([a_q(1,Q+1);2*a_q(2:end,Q+1)]);
    D_beta_F = reshape(D_beta_F_q,[],1);
    
    a0_dot = -inv(Derivative_DF(a,N,alpha,beta,u2,delta,Q))*D_beta_F;
    if j == 1
        beta_n = beta;
        a_n = a;
    else
        beta_n = beta-delta_s;
        a_n = a-delta_s*a0_dot;
    end
    
    Fn = function_F(a_n,N,alpha,beta_n,u2,delta,Q);
    
    k=0;
    while norm(Fn) > 10^(-14) && k < 25
        k = k+1;
        a_n = a_n - Derivative_DF(a_n,N,alpha,beta_n,u2,delta,Q)\Fn;
        Fn = function_F(a_n,N,alpha,beta_n,u2,delta,Q);
    end
    if k == 25 || norm(a_n) < 10^(-10)
        break
    end
    a = a_n;
    beta = beta_n;
    A(:,j) = a;
    B(:,j) = beta;
    j = j+1;
end
h = beta-1;
% Proven Solution
save('Data\h0_2_h1_discrete_solution_unproven','A','B','h')
% Begining of next step
save('Data\data_h1','a','N','beta','r','h')
end

