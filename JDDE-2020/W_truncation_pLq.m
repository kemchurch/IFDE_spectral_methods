function W = W_truncation_pLq(A,C1,N2,p,q)
d = size(C1)*[1;0];
W = zeros(d*(q*(N2+1) + 1),d*(q*(N2+1) + 1));
I=eye(d);

% Wp0 stuff
[Wp0,Wp0special1] = Wp0_truncation_pLq(A{p},C1,N2+1,p);
W(1:d+d*(N2+1),1:d+d*(N2+1)) = Wp0;
if p~=1
    W(1:d,1+d+d*(N2+1):d+2*d*(N2+1)) = (I+C1)*Wp0special1;
    W(d+1:2*d,1+d+d*(N2+1):d+2*d*(N2+1)) = Wp0special1;
end

%Remainder (extra terms to handle shifted indexes k-p-1)
if p>1
    for k=1:p-1
        nk = p-1 -(k-1);
        [Wnk,Wnk_extra] = Wk_truncation_pLq(A{nk},C1,N2+1,nk);
        W(1+d+d*(N2+1)+(k-1)*d*(N2+1):d+d*(N2+1)+(k)*d*(N2+1), ...
            1+d+d*(N2+1)+(k-1)*d*(N2+1):d+d*(N2+1)+(k)*d*(N2+1)) = Wnk;
        if nk~=1
            W(1+d+d*(N2+1)+(k-1)*d*(N2+1) : d+d*(N2+1)+(k-1)*d*(N2+1) + d, ...
                1+d+d*(N2+1)+(k)*d*(N2+1):d+d*(N2+1)+(k)*d*(N2+1) + (N2+1)*d)=Wnk_extra;
        end
    end
end

end