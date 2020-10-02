function [an,epsn,k,nF] = Newton_scriptF(a,epsilon,N,alpha,beta,u2,delta,u1)
an = a;
epsn = epsilon;
k = 0;
while norm(scriptF(an,epsn,N,alpha,beta,u2,delta,u1)) > 10^(-15) && k <200
    k = k+1;
    iterplus1 = [an;epsn] - D_scriptF(an,epsn,N,alpha,beta,u2,delta,u1)\scriptF(an,epsn,N,alpha,beta,u2,delta,u1);
    an = iterplus1(1:end-1);
    epsn = iterplus1(end);
end
nF = norm(scriptF(an,epsn,N,alpha,beta,u2,delta,u1));
end