function [D_a_kF_j] = Derivative_D_a_kF_q_intval(a_k,a_delta,q,k,N,alpha,beta,u2,delta,Q)

% Case k = q
if k == q
    % Parameters
    n = ones(1,N+1);
    T = -diag(n,-1) + diag(n,1);
    T(1,:) = 0;
    n = (1:N+2)';
    J = (1:N+1).*ones(N+2,N+1);
    nminusJ = abs(n-J)+1;
    nplusJ = abs(n+J)-1;
    
    % Preallocation
    D_a_kF_j = intval(2*diag(0:N));
    D_a_kPsi_q = intval(zeros(N+2));
    a_delta_ext = [a_delta;zeros(N+2,1)];
    
    % First row
    if q == Q
        D_a_kF_j(1,:) = beta*[1,2*ones(1,N)];
    else
        D_a_kF_j(1,:) = [1,2*ones(1,N)];
    end
    
    % D_a_kPsi_q
    Id = intval(eye(N+2));
    D_a_kPsi_q(:,1) = alpha/(2*u2)*(Id(:,1)-a_delta_ext(nminusJ(:,1)));
    D_a_kPsi_q(:,2:end) = alpha/(2*u2)*(Id(:,2:end)-[(a_delta_ext(nminusJ(:,2:end) )+a_delta_ext(nplusJ(:,2:end))),zeros(N+2,1)]);
    
    % D_a_kF_j
    TD_a_kPsi_q = T*D_a_kPsi_q;
    D_a_kF_j = D_a_kF_j + TD_a_kPsi_q(1:N+1,1:N+1);
else
    % Preallocation
    D_a_kF_j = intval(zeros(N+1));
end

% Case k = q+1
if k == q+1 || (k == 0 && q == Q)
    % First row
    D_a_kF_j(1,:) = intval([1,2*ones(1,N)].*(-1).^(1:N+1));
end

% Case k = mod(q + delta,Q)
if k == mod(q + delta,Q+1) 
    % Parameters
    n = ones(1,N+1);
    T = -diag(n,-1) + diag(n,1);
    T(1,:) = 0;
    n = (1:N+2)';
    J = (1:N+1).*ones(N+2,N+1);
    nminusJ = abs(n-J)+1;
    nplusJ = abs(n+J)-1;
    
    % Preallocation
    D_a_kPsi_q = intval(zeros(N+2));
    a_k_ext = [a_k;zeros(N+2,1)];
    
    % D_a_kPsi_q
    D_a_kPsi_q(:,1) = -alpha/(2*u2)*a_k_ext(nminusJ(:,1));
    D_a_kPsi_q(:,2:end) = -alpha/(2*u2)*[(a_k_ext(nminusJ(:,2:end) )+a_k_ext(nplusJ(:,2:end))),zeros(N+2,1)];
    
    % D_a_kF_j
    TD_a_kPsi_q = T*D_a_kPsi_q;
    D_a_kF_j = D_a_kF_j + TD_a_kPsi_q(1:N+1,1:N+1);
end
end

