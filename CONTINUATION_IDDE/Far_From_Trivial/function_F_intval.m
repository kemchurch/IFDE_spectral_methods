function [F] = function_F_intval(a,N,alpha,beta,u2,delta,Q)

% Preallocation
a_q = reshape(a,[N+1,Q+1]);
a_q_ext = [a_q;zeros(1,Q+1)];
a_qa_qdelta = intval(zeros(N+2,Q+1));
F_q = 2*(0:N)'.*a_q;
n = ones(1,N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;

% Convolution
for n = 0:Q
    a_qa_qdelta(:,n+1) = quadratic_cheb_intval(a_q_ext(:,n+1),a_q_ext(:,mod(n+delta,Q+1)+1));
end

% Psi
Psi_q = alpha/(2*u2)*(a_q_ext - a_qa_qdelta);
Tpsi = T*Psi_q;

% Initial condition
F_q(1,1:Q) = sum([a_q(1,1:Q)-a_q(1,2:Q+1);2*(a_q(2:end,1:Q)-((-1).^(1:N))'.*a_q(2:end,2:Q+1))]);
F_q(1,Q+1) = sum([beta*a_q(1,Q+1)-a_q(1,1);2*(beta*a_q(2:end,Q+1)-((-1).^(1:N))'.*a_q(2:end,1))]);

% Function F_q
F_q(2:end,:) = F_q(2:end,:) + Tpsi(2:N+1,:);

% Function F
F = reshape(F_q,[],1);
end