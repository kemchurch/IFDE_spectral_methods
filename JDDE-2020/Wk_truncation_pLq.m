function [W,W_extra] = Wk_truncation_pLq(A,C1,N2,k)
W = Ek_truncation_pLq(A,C1,N2);
d = size(C1)*[1;0];
I = eye(d,d);
W_extra = zeros(d,N2*d);
if k~=1
    W_extra(:,1:d)=I;
    for n=1:N2-1
        W_extra(:,1+n*d:(n+1)*d) = 2*I;
    end
end
end