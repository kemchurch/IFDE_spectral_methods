function [M,E,W] = M_truncation_pLq(A,B,C1,C2,p,q,N2)
%inputs:
% A,B: Cell arrays of Chebyshev coefficient tensors for the matrix
% functions A and B
% C1,C2: Matrices C1 and C2, both constant
% N2: Number of (non-constant) Chebyshev modes to take in discretization
%outputs:
% M: discretization of monodromy operator given input datum
% E,W: individual matrices used to define M; needed for rigorous error
% bounds and eigenvalue validation
d=size(C1)*[1;0];
W = W_truncation_pLq(A,C1,N2,p,q);
E = E_truncation_pLq(B,C1,C2,N2,p,q);
M = (eye(d*(q*(N2+1) + 1))-W)\E;
end