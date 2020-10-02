function [F] = function_D_beta_F_int(a,epsilon,N,alpha,beta,u2,delta,u1)
% To do this partial derivative we basically recycle the function_F and
% delete everything that should not be there.
Q = u1-1;
% Preallocation
a_q = reshape(a,[N+1,Q+1]);
F_q = intval(zeros(N+1,Q+1));

% Initial condition
F_q(1,Q+1) = sum([beta*a_q(1,Q+1);2*(beta*a_q(2:end,Q+1))]);

% Function F
F = reshape(F_q,[],1);
end