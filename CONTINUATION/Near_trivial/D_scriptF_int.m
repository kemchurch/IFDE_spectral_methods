function [D_scriptF] = D_scriptF_int(a,epsilon,N,alpha,beta,u2,delta,u1)
% Preallocate
D_scriptF = intval(zeros(u1*(N+1)+1,u1*(N+1)+1));
% Build
[DF_ak,DF_eps] = Derivative_DF_int(a,epsilon,N,alpha,beta,u2,delta,u1);
D11 = DF_ak;
D12 = DF_eps;
D21 = intval([-[1,2*ones(1,N)],zeros(1,(u1-1)*(N+1))]);
% Place entries
D_scriptF(1:u1*(N+1),1:u1*(N+1)) = D11;
D_scriptF(1:end-1,end) = D12;
D_scriptF(end,1:u1*(N+1)) = D21;
end