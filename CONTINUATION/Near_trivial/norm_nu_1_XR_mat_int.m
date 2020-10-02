function nrm = norm_nu_1_XR_mat_int(A,N,u1,nu)
nu = intval(nu);
A11 = A(1:end-1,1:end-1);
A12 = A(1:end-1,end);
A21 = A(end,1:end-1);
A22 = A(end,end);
omega = [1,2*nu.^(1:1:N-1)];
normsA11 = intval(zeros(u1,1));
normsA12 = intval(zeros(u1,1));
normsA21 = intval(zeros(u1,1));

for i=0:u1-1
    normsA12(i+1) = norm_1_nu(A12(1+i*N:(i+1)*N),nu);
    normsA21(i+1) = sum(abs(A21(1+i*N:(i+1)*N))./omega);
    for q=0:u1-1
        normsA11(i+1) = normsA11(i+1) + norm_nu_1_op(A11(1+i*N:(i+1)*N, 1+q*N:(q+1)*N),nu);
    end
end
nA11 = max(normsA11);
nA22 = abs(A22);
nA12 = max(normsA12);
nA21 = max(normsA21);
nrm = nA11 + nA12 + nA21 + nA22;
end