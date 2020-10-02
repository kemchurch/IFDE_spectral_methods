function [df] = finite_diff(a,N,alpha,beta,u2,delta,Q)

% Preallocation
h = 1e-6;
E = eye(N+1);
df = zeros((Q+1)*(N+1));
a_k = reshape(a,[N+1,Q+1]);
dftemp = zeros(N+1);

% Finite Difference
F = function_F(a,N,alpha,beta,u2,delta,Q);

for q=0:Q
    for k = 0:Q
        for j=1:N+1
            a_kh=a_k(:,k+1)+h*E(:,j);
            ah = a_k;
            ah(:,k+1) = a_kh;
            ah = reshape(ah,[],1);
            Fah = function_F(ah,N,alpha,beta,u2,delta,Q);
            dftemp(:,j) = (Fah(q*(N+1)+1:(q+1)*(N+1)) - F(q*(N+1)+1:(q+1)*(N+1)) )/h;
        end
        df((q*(N+1)+1:(q+1)*(N+1)),(k*(N+1)+1:(k+1)*(N+1))) = dftemp;
    end
end
end

% repmat(F(q*(N+1)+1:(q+1)*(N+1)),Q+1,1)