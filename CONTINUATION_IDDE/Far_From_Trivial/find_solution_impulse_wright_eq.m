function [] = find_solution_impulse_wright_eq(N,r,tau,T,h)
% Description:
%       Description of the function goes here
% Input: 
%       N = length(a_q)-1, Truncation number
%       r = 
%       h = 
%       tau = 
%       T = 
% Output:
%
% singularitey at h = 0.259181779318282
% How close can we get to that : 0.258



% Parameters
delta_s = 0.001;
h0=0;
alpha = r*tau;
beta = 1-h0;
u  = T/tau;
[u1,u2] = rat(u);
if u1 >= u2
    delta = u1-u2;
else
    k=2;
    while k*u1-u2<0
        k = k+1;
    end
    delta = k*u1 - u2;
end
Q = u1-1;
a = zeros(N+1,1);
a(1) =1;
a = repmat(a,Q+1,1);
while h0 < h
    beta_temp = beta-delta_s;
    [anplus1,k] = Newton_iwe(a,N,alpha,beta_temp,u2,delta,Q);
    if norm(anplus1,1) < 10^(-10)
        delta_s = delta_s/2;
        if delta_s < 10^(-6)
            error('delta_s too small and converge to trivial solution')
        end
    else
        a = anplus1;
        beta = beta_temp;
        h0 = 1-beta;
        delta_s = delta_s*2^((4-k)/3);
        %         beta
        %         error('You reach the trivial solution')
    end
end
if h0 ~= h
    a = Newton_iwe(a,N,alpha,1-h,u2,delta,Q);
    beta = 1-h;
end
save('data_h0','a','N','Q','u2','delta','alpha','beta','r','tau','T','h')
end


















