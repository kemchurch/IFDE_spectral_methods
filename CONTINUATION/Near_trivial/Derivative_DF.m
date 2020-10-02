function [DF_ak,DF_eps] = Derivative_DF(a,epsilon,N,alpha,beta,u2,delta,u1)
Q = u1-1;
% Preallocation, precomputation
a_q = reshape(a,[N+1,Q+1]);
DF_ak = zeros((Q+1)*(N+1),(Q+1)*(N+1));
n = ones(1,N+1);
T = -diag(n,-1) + diag(n,1);
T(1,:) = 0;
a_qa_qdelta = zeros(N+2,Q+1);
a_q_ext = [a_q;zeros(1,Q+1)];
for n = 0:Q
    a_qa_qdelta(:,n+1) = quadratic_cheb(a_q_ext(:,n+1),a_q_ext(:,mod(n+delta,Q+1)+1));
end

% DF_ak
for q = 0:Q
    for k = 0:Q
        DF_ak(q*(N+1)+1:(q+1)*(N+1),k*(N+1)+1:(k+1)*(N+1)) = Derivative_D_a_kF_q(a_q(:,q+1),a_q(:,mod(q+delta,Q+1)+1),epsilon,q,k,N,alpha,beta,u2,delta,Q);
    end
end

% DF_eps
TPsi = T*a_qa_qdelta;
DF_eps = -alpha/(2*u2)*reshape(TPsi(1:N+1,:),[],1);
% for q=0:Q
%     DF_eps(1+q*(N+1):(q+1)*(N+1)) = -alpha/(2*u2)*reshape(T*a_qa_qdelta,[],1);
% end

end