function [] = para_cont(N,alpha,beta,u2,delta,Q,a,r,T,tau,nu)
delta_s = 0.0001;
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
    [Y0,Z0,Z1,Z2,R] = Bounds(a,N,alpha,beta,u2,delta,Q,nu);
    if R(1) < 0 || Z1+Z0> 1
        break
    end
    A(:,j) = a;
    B(:,j) = beta;
    Y0_MAT(:,j) = Y0;
    Z0_MAT(:,j) = Z0;
    Z1_MAT(:,j) = Z1;
    Z2_MAT(:,j) = Z2;
    R_MAT(:,j) = R;
    if R_MAT(1,j) >= 0 && Z1_MAT(1,j) < 1
        subplot(1,2,1);
        plot(B(j),abs(A(1,j)),'g.');hold on;
        title('Away from \beta^*')
        set(gca,'Color','k')
        xlabel('\beta');
        ylabel('|a_0|');
        xlim([exp(-r*T) 1])
        ylim([0 1])
        subplot(1,2,2);
        set(gca,'Color','k')
        title('Near \beta^*')
        xlabel('\beta');
        ylabel('|a_0|');
        drawnow
    else
        subplot(1,2,1);
        plot(B(j),abs(A(1,j)),'b*');hold on;
        title('Away from \beta^*')
        set(gca,'Color','k')
        xlabel('\beta');
        ylabel('|a_0|');
        xlim([exp(-r*T) 1])
        ylim([0 1.1])
                subplot(1,2,2);
        set(gca,'Color','k')
        title('Near \beta^*')
        xlabel('\beta');
        ylabel('|a_0|');
        drawnow
    end
    j = j+1;
end
Y0_MAT_existence = exist('Y0_MAT');
if Y0_MAT_existence == 0
    error('Unable to do a step foward in parameter beta.')
end
subplot(1,2,1);
plot(exp(-r*T),0,'ro')
set(gca,'Color','k')
xlabel('\beta');
ylabel('|a_0|');
xlim([exp(-r*T) 1.1])
ylim([0 1.1])
h = beta-1;
drawnow
% Proven Solution
save('Data\h0_2_h1_proven','A','B','Y0_MAT','Z0_MAT','Z1_MAT','Z2_MAT','R_MAT','N','Q','u2','delta','alpha','r','tau','T','h','nu')
% Begining of next step
save('Data\data_h1','a','N','Q','u2','delta','alpha','beta','r','tau','T','h')
end



