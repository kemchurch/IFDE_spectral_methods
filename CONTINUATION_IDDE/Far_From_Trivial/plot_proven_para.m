function [] = plot_proven_para()
close all
figure
counter = 0;
load('Proof_cont\data_h0_2_h1_proven');
R1 = R_MAT_cont;
for i = 1:length(B)-1
    A0 = reshape(A(:,i),[N+1,Q+1]);
    a0_0 = A0(1,:);
    A1 = reshape(A(:,i+1),[N+1,Q+1]);
    a0_1 = A1(1,:);
    norm_a0 = [norm(a0_0,inf),norm(a0_1,inf)];
    if R_MAT(1,i) >= 0 && Z1_MAT_cont(1,i)+Z0_MAT_cont(1,i)+Z0_hat_MAT_cont(1,i) < 1
        plot(B(i:i+1),norm_a0,'bs-','LineWidth',1.5,'Markerfacecolor','b');hold on;
    else
        plot(B(i),norm(a0_0),'white*');hold on;
    end
end
load('Proof_cont\data_h1_2_h2_proven');
for i = 1:length(B)-1
    A0 = reshape(A(:,i),[N+1,Q+1]);
    a0_0 = A0(1,:);
    A1 = reshape(A(:,i+1),[N+1,Q+1]);
    a0_1 = A1(1,:);
    norm_a0 = [norm(a0_0,inf),norm(a0_1,inf)];
    if R_MAT(1,i) >= 0 && Z1_MAT_cont(1,i)+Z0_MAT_cont(1,i)+Z0_hat_MAT_cont(1,i) < 1
        plot(B(i:i+1),norm_a0,'cs-','LineWidth',1.5,'Markerfacecolor','c');hold on;
    else
        plot(B(i),norm(a0_0),'g*');hold on;
    end
end
load('Proof_cont\data_h1_2_h2_proven');
R2 = R_MAT_cont;
r_max = max([R1(1,:),R2(1,:)]);
plot(exp(-r*T),0,'ro','Markerfacecolor','r')
xlabel('\beta');
ylabel('||a_0||_\infty');
set(gca,'Color','k')
axis([exp(-r*T) 1 0 1])
title(['r_0 = ' num2str(r_max)])
% figure
% load('movie_of_sol.mat')
% movie(F,1,24)





end
