function [DF_beta] = Derivative_D_beta(a,N,Q)

a_reshape = reshape(a,N+1,Q+1);
D_beta_F_q = intval(zeros(N+1,1));
D_beta_F_q(1) = sum([a_reshape(1,end);2*a_reshape(2:end,end) ]);
DF_beta = [zeros(Q*(N+1),1);D_beta_F_q];

end