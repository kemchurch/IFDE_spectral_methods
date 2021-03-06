function E = Ek_truncation_pLq(B,C1,N2)
d = size(C1)*[1;0];
E = zeros(d*(N2+1),d*(N2+1));
I = eye(d,d);
B(length(B)+1:d*(2*(2*N2)+4),:) = zeros(d*(2*(2*N2)+4)-length(B),d); %padding 

% 0-Chebyshev row (exclude 1st column)
d_0sum = (1/2)*B(1:d,:) - (1/4)*B(d+1:2*d,:);
for n=2:N2
    d_0sum = d_0sum - ((-1)^n/(n^2-1))*B(n*d+1:(n+1)*d,:);
end
E(d+1:2*d,d+1:2*d) = d_0sum;
for k=1:N2-1
    d_ksum = B(d*k+1:d*(k+1),:) - (1/4)*(B(1+(1+k)*d:(2+k)*d,:) + ...
        B(1+(k-1)*d:k*d,:));
    for j = 2:N2+k
        d_ksum = d_ksum - ((-1)^j/(j^2-1))*( B(1+d*abs(j-k):d*(abs(j-k)+1),:)...
            + B(1+d*abs(j+k):d*(abs(j+k)+1),:) );
    end
    E(d+1:2*d,1+(1+k)*d:(2+k)*d) = d_ksum;
end

% f1-column (column 2, rows 3-end)
for n=2:N2
   E(1+n*d:(n+1)*d,1+d:2*d) = (1/(4*(n-1)))*(B(1+d*(n-2):d*(n-1),:) - ...
       B(1+d*n:d*(n+1),:));
end

% tridiagonal operator
T = zeros(d*(N2+1),d*(N2+1));
if N2==2
    T = 0;
else
    T(1+d:2*d,1:d) = eye(d,d);
    for i=2:N2
        T(1+d*(i-2):d*(i-1),1+d*(i-1):d*i) = -eye(d,d);
        T(1+d*i:d*(i+1),1+d*(i-1):d*i) = eye(d,d);
    end
    if N2>=2
    T(1+d*(N2-1):d*(N2),1+d*(N2):d*(N2+1)) = -eye(d,d);
    end
end

% B-convolution
Bconv = zeros(d*(N2+1),d*(N2+1));
for j=1:N2+1
    Bconv(1:d,1+d*(j-1):d*j) = 2*B(1+d*(j-1):d*j,:);
    for i=2:N2+1
        Bconv(1+d*(i-1):d*i,1+d*(j-1):d*j) = ...
            B(1+d*(i-1+j-1):d*(1+(i-1+j-1)),:) + ...
            B(1+d*abs(i-1-(j-1)):d*(1+abs(i-1-(j-1))),:);
    end
end
% Now repair the first column
Bconv(1:d*(N2+1),1:d) = B(1:d*(N2+1),1:d);

% Bottom Right tridiagonal weighted product with convolution
scale = zeros(d*N2,d*N2);
for j=1:N2
    scale(1+(j-1)*d:j*d,1+(j-1)*d:j*d) = eye(d,d)/(4*j);
end
Ebar = blkdiag(zeros(d,d),scale)*T*Bconv;
E_block = Ebar(1+d:(N2)*d,1+d:(N2)*d);
E(1+2*d:d*(N2+1),1+2*d:d*(N2+1)) = E_block;

% Final cleanup
E = E(d+1:end,d+1:end);

end