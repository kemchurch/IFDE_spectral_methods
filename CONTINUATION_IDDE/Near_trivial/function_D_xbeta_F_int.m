function DF = function_D_xbeta_F_int(a,epsilon,N,alpha,beta,u2,delta,u1)
D_a_eps_F = D_scriptF_int(a,epsilon,N,alpha,beta,u2,delta,u1);
D_beta_F = function_D_beta_F_int(a,epsilon,N,alpha,beta,u2,delta,u1);
DF = intval(zeros(u1*(N+1)+1,u1*(N+1)+2));
DF(1:end,1:end-1) = D_a_eps_F;
DF(1:end-1,end) = D_beta_F;
end