function [scriptF] = scriptF(a,epsilon,N,alpha,beta,u2,delta,u1)
a_q = reshape(a,[N+1,u1]);
F1 = function_F(a,epsilon,N,alpha,beta,u2,delta,u1);
F2 = 1 - a_q(1,1) - 2*sum(a_q(2:end,1));
scriptF = [F1;F2];
end