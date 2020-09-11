function W = W_truncation_pq1(A,C1,C2,N2)
W = E_truncation_pq1(A,C1,C2,N2);
d = size(C1)*[1;0];
W(:,1:d) = zeros(d*(N2+1),d);
end