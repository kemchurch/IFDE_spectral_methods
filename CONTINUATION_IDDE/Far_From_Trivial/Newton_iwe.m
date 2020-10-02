function [an,k] = Newton_iwe(a,N,alpha,beta,u2,delta,Q)
an = a;
k = 0;
while norm(function_F(an,N,alpha,beta,u2,delta,Q)) > 10^(-14) && k <200
    k = k+1;
    anplus1 = an - Derivative_DF(an,N,alpha,beta,u2,delta,Q)\function_F(an,N,alpha,beta,u2,delta,Q);
    an = anplus1;
end
end