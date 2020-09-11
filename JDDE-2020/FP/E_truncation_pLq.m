function E = E_truncation_pLq(B,C1,C2,N2,p,q)
d = size(C1)*[1;0];
E = zeros(d*(q*(N2+1) + 1),d*(q*(N2+1) + 1));

% Identity part:
% Leftovers from E1
E(1+d+(p-1)*d*(N2+1):2*d+(p-1)*d*(N2+1),1:d)=eye(d);
% Everything else
E(1+d+p*d*(N2+1):end,1+d:d+(q-p)*d*(N2+1))=eye(d*(q-p)*(N2+1));

%Handling Ep and Ep0
[Ep,Ep01,Ep02] = Ep0_truncation_pLq(B{p},C1,C2,N2+1,p,q);
E(1:d+d*(N2+1), 1+d+(q-p-1)*d*(N2+1) + d*(N2+1) : d+(q-p-1)*d*(N2+1) + 2*d*(N2+1)) = Ep;
E(1:d+d*(N2+1),1:d)=Ep02;
if p~=q
    E(1:d, 1+d+(q-p-1)*d*(N2+1):d+(q-p-1)*d*(N2+1) + d*(N2+1)) = Ep01;
end

%Remainder
if p>1
    Eblock = zeros((p-1)*d*(N2+1),(p-1)*d*(N2+1));
    for k=1:p-1
        nk = p-1 -(k-1);
        Enk = Ek_truncation_pLq(B{nk},C1,N2+1);
        Eblock(1+d*(N2+1)*(k-1):d*(N2+1)*(k),1+d*(N2+1)*(k-1):d*(N2+1)*(k)) = Enk;
    end
    if p~=q
        E(1+d+d*(N2+1):d+d*(N2+1)+(p-1)*d*(N2+1),...
            d+(q-p-1)*d*(N2+1) + 2*d*(N2+1) + 1:end) = Eblock;
    else
        E(1+d+d*(N2+1):end ,1+d+(q-p-1)*d*(N2+1) + 2*d*(N2+1):end) = Eblock;
    end
end
end