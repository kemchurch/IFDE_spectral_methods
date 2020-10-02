function [] = proof_para_cont(A,N,alpha,B,u2,delta,Q,nu,r,T,tau)
k=0;
for i = 1:length(B)-1
    x_0 = intval(A(:,i));
    x_1 = intval(A(:,i+1));
    s = midrad(0.5,0.5);
    x_s = (x_1-x_0)*s+x_0;
    beta_0 = intval(B(i));
    beta_1 = intval(B(i+1));
    beta_s = (beta_1-beta_0)*s+beta_0;
    [Y0,Y0_hat,Z0,Z0_hat,Z1,Z2,R] = Bounds_para_cont(A(:,i),A(:,i+1),B(i),B(i+1),x_0,x_s,N,alpha,beta_0,beta_s,u2,delta,Q,nu);
    
    if R(1) < 0 || Z1+Z0+Z0 > 1 || imag(R(1)) ~= 0
        break
    else
        k = k+1;
    end
    Y0_MAT(:,i) = Y0;
    Z0_MAT(:,i) = Z0;
    Y0_hat_MAT(:,i) = Y0_hat;
    Z0_hat_MAT(:,i) = Z0_hat;
    Z1_MAT(:,i) = Z1;
    Z2_MAT(:,i) = Z2;
    R_MAT(:,i) = R;
    subplot(1,2,1);
    plot(B(i:i+1),abs(A(1,i:i+1)),'g.-');hold on;
    title('Away from \beta^*')
    xlabel('\beta');
    ylabel('|a_0|');
    set(gca,'Color','k')
    xlim([exp(-r*T) 1])
    ylim([0 1])
    subplot(1,2,2);
    set(gca,'Color','k')
    title('Near \beta^*')
    xlabel('\beta');
    ylabel('|a_0|');
    drawnow
    
end
if k == 0
    error("The branch fail to compute")
else
    A = A(:,1:k);
    B = B(:,1:k);
end
% Proven Solution
save('Data\data_h0_2_h1_proven','A','B','Y0_MAT','Z0_MAT','Y0_hat_MAT','Z0_hat_MAT','Z1_MAT','Z2_MAT','R_MAT','N','Q','u2','delta','alpha','r','tau','T','nu')
end
