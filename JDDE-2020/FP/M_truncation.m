function M=M_truncation(A,B,C1,C2,p,q,N2)

% Compute discretized monodromy operator for an impulsive delay delay 
% differential equation (IDDE) with specified order Chebyshev spectral 
% method.
% M = M_truncation(A,B,C1,C2,p,q,N2) computes discretized monodromy 
% operator M for the IDDE
% x'=A(t)*x(t) + B(t)*x(t-q), t/=/pk
% Delta(x) = C1*x(t-) + C2*x(t-q), t=pk
% for integers k. The delay (q) and period (p) must be positive integers.
% C1 and C2 are real 2-dimensional arrays of dimension n. If A(t) is a 
% constant matrix, A can be input as such. Otherwise, it must 
% be input as a cell array A=cell{N1,p} such that A{m,n} is the m'th matrix 
% Chebyshev coefficient of the function t->A_n(t)=A((t-1)/2 + n). That is, 
% A*(t)=A{1,n} + 2*sum_{m>1} A{m,n}T_m(t) for T_m the m'th Chebyshev
% polynomial of the first kind. The input B is similar. The order of the
% method is N2, a positive integer.

if size(C1)~=size(C2)
    error('Mismatch: size of C1 and C2')
end
d = length(C1);

% Convert matrix to 1x1 cell if input as such and p=1
if p==1
    if isa(A,'double')
        A=mat2cell(A,rows(A),columns(A));
    end
    if isa(B,'double')
        B=mat2cell(B,rows(B),columns(B));
    end
end

% Convert matrix to 1xp cell if the input A and/or B is constant.
if p>1 && length(A)==d
    A = mat2cell(repmat(A,[1,p]),d,d*ones(1,p));
end
if p>1 && length(B)==d
    B = mat2cell(repmat(B,[1,p]),d,d*ones(1,p));
end


% Simple mismatch errors
if p>1
    if ~isa(A,'cell')
        error('p>1 requires input A to be a cell array')
    elseif ~isa(B,'cell')
        error('p>1 requires input B to be a cell array')
    elseif length(A)~=p
        error('Mismatch: length of A and period p')
    elseif length(B)~=p
        error('Mismatch: length of B and period p')
    end
end


% The matrix computations
if p<q
    [M,~,~]=M_truncation_pLq(A,B,C1,C2,p,q,N2);
end

if p==q
    if p==1
        [M,~,~]=M_truncation_pLq(A,B,C1,C2,p,q,N2);
    else
        M=eye(d*(q*(N2+1) + 1));
        for k=1:p-1
            Ak=mat2cell(A{k},rows(A{k}),columns(A{k}));
            Bk=mat2cell(B{k},rows(B{k}),columns(B{k}));
            [Mk,~,~]=M_truncation_pLq(Ak,Bk,zeros(size(C1)),zeros(size(C2)),1,q,N2);
            M=Mk*M;
        end
        Ap=mat2cell(A{p},rows(A{p}),columns(A{p}));
        Bp=mat2cell(B{p},rows(B{p}),columns(B{p}));
        [Mp,~,~]=M_truncation_pLq(Ap,Bp,C1,C2,1,q,N2);
        M=Mp*M;
    end    
end

if p>q
    M=eye(d*(q*(N2+1) + 1));
    for m=1:p-1
        Am=mat2cell(A{1},rows(A{1}),columns(A{1}));
        Bm=mat2cell(B{1},rows(B{1}),columns(B{1}));
        [Mm,~,~]=M_truncation_pLq(Am,Bm,zeros(size(C1)),zeros(size(C2)),1,q,N2);
        M=Mm*M;
        A(1)=[];
        B(1)=[];
    end
    [Mp,~,~]=M_truncation_pLq(A,B,C1,C2,1,q,N2);
    M=Mp*M;
        
end

end