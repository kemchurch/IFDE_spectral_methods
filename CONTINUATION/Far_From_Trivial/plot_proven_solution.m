function [] = plot_proven_solution()
close all
fig = figure('Position',  [100, 150, 1700, 600]);
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
        subplot(1,2,1);
        plot(B(i:i+1),norm_a0,'bs-','LineWidth',1.5,'Markerfacecolor','b');hold on;
    else
        subplot(1,2,1);
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
        subplot(1,2,1);
        plot(B(i:i+1),norm_a0,'cs-','LineWidth',1.5,'Markerfacecolor','c');hold on;
    else
        subplot(1,2,1);
        plot(B(i),norm(a0_0),'g*');hold on;
    end
end
load('Proof_cont\data_h1_2_h2_proven');
R2 = R_MAT_cont;
r_max = max([R1(1,:),R2(1,:)]);
subplot(1,2,1);
plot(exp(-r*T),0,'ro','Markerfacecolor','r')
xlabel('\beta');
ylabel('||a_0||_\infty');
set(gca,'Color','k')
axis([exp(-r*T) 1 0 1])
title(['r_0 = ' num2str(r_max)])
load('Proof_Discrete\h1_2_h2_proven');
B2 = B;
load('Proof_Discrete\h0_2_h1_proven');
B1 = B;
m = 1;
for j = 1:length(B)
    subplot(1,2,1);
    A0 = reshape(A(:,j),[N+1,Q+1]);
    a0_0 = A0(1,:);
    plot(B(j),norm(a0_0,inf),'gs-','LineWidth',1.5,'Markerfacecolor','g');hold on;
    a_q = reshape(A(:,j),[N+1,Q+1]);
    a_q(2:end,:) = 2*a_q(2:end,:);
    phi_q = zeros(1,Q+1);
    delta_t = 2/(Q+1);
    tt = -1:0.1:1;
    for q = 0:Q
        i = 1;
        for t = tt
            phi_q(i,q+1) = sum(a_q(:,q+1).*transpose(chebyshevT(0:N,t)));
            i = i+1;
        end
        subplot(1,2,2);
        plot(tt+2*q-Q,phi_q(:,q+1),'LineWidth',1.5);hold on;
        axis([-10 10 -0.1 1.1])
        set(gca,'Color','k')
        xlabel('t');
        ylabel('\phi(t)')
        drawnow
        F(m) = getframe(fig);
        m=m+1;
    end
    subplot(1,2,1);
    plot(B(j),norm(a0_0,inf),'bs-','LineWidth',1.5,'Markerfacecolor','b');hold on;
end
load('Proof_Discrete\h1_2_h2_proven');
subplot(1,2,1);
A0 = reshape(A(:,1),[N+1,Q+1]);
a0_0 = A0(1,:);
axis([exp(-r*T) B(1) 0 norm(a0_0,inf)])
for j = 2:length(B)
    subplot(1,2,1);
    A0 = reshape(A(:,j),[N+1,Q+1]);
    a0_0 = A0(1,:);
    plot(B(j),norm(a0_0,inf),'gs-','LineWidth',1.5,'Markerfacecolor','g');hold on;
    a_q = reshape(A(:,j),[N+1,Q+1]);
    a_q(2:end,:) = 2*a_q(2:end,:);
    phi_q = zeros(1,Q+1);
    delta_t = 2/(Q+1);
    tt = -1:0.1:1;
    for q = 0:Q
        i = 1;
        for t = tt
            phi_q(i,q+1) = sum(a_q(:,q+1).*transpose(chebyshevT(0:N,t)));
            i = i+1;
        end
        subplot(1,2,2);
        plot(tt+2*q-Q,phi_q(:,q+1),'LineWidth',1.5);hold on;
        axis([-10 10 -0.1 1.1])
        drawnow
        F(m) = getframe(fig);
        m=m+1;
    end
    subplot(1,2,1);
    plot(B(j),norm(a0_0,inf),'cs-','LineWidth',1.5,'Markerfacecolor','c');hold on;
end
close all
[h0, w, p] = size(F(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h0]);
axis off
movie(hf,F);
save('movie_of_sol.mat','F')
end
