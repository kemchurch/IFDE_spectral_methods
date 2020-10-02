function BXR = compute_bold_B_star(A,N,u1,nu)
omega = [1,2*nu.^(1:N)];
BXR = intval(zeros(u1,1));
for j=0:u1-1
    BXR(j+1) = max(1./omega.*abs(A(end, 1+j*(N+1):(j+1)*(N+1))));
end
end