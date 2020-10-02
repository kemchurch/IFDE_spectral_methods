function BXX = compute_bold_B(A,N,u1,nu)
omega = [1,2*nu.^(1:N)];
BXX = intval(zeros(u1,u1));
for m=0:u1-1
    for j=0:u1-1
        BXX(m+1,j+1) = max(1./omega.*norm_1_nu(A(1+m*(N+1):(m+1)*(N+1), 1+j*(N+1):(j+1)*(N+1)),nu));
    end
end
end

