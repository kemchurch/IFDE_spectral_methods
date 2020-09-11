function [M,E,W] = M_truncation_pq1(A,B,C1,C2,N2)
%inputs:
% B: Chebyshev coefficient tensor for the matrix function B
% C1,C2: Matrices C1 and C2, both constant
% N2: Number of (non-constant) Chebyshev modes to take in discretization
%outputs:
% M: discretization of monodromy operator given input datum
% E,W: individual matrices used to define M; needed for rigorous error
% bounds and eigenvalue validation
W = W_truncation_pq1(A,C1,C2,N2+1);
E = E_truncation_pq1(B,C1,C2,N2+1);
M = (eye(length(C1)*(N2+2),length(C1)*(N2+2))-W)\E;
end