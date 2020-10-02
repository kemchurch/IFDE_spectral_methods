function [a,epsilon,beta,nF] = scriptF_continue_beta(a0,epsilon0,steps,N,alpha,beta0,beta1,u2,delta,u1)
if beta0 < beta1
    beta = linspace(beta0,beta1,steps);
else
    beta = flip(linspace(beta1,beta0,steps));
end
a = zeros(u1*(N+1),steps);
epsilon = zeros(1,steps);
nF = zeros(1,steps);
[a(:,1),epsilon(1),~,nF(1)] = Newton_scriptF(a0,epsilon0,N,alpha,beta(1),u2,delta,u1);
for i=2:steps
    [a(:,i),epsilon(i),~,nF(i)] = Newton_scriptF(a(:,i-1),epsilon(i-1),N,alpha,beta(i),u2,delta,u1);
end
end