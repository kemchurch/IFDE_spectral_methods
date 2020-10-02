function [P] = xi(a,nu)
% a = transpose(abs(a));
% N = length(a)-1;
% P = zeros(m,1);
% for k = 0:m-1
%     if k==0 && N == m-1
%         P(1) = 0;
%     else
%         j1 = m:(N+k);
%         max1 =  max(a(abs(k-j1)+1)./(nu.^abs(j1)));
%         if k-N <= -m
%             j2 = k-N:-m;
%             max2 =  max(a(abs(k-j2)+1)./(nu.^abs(j2)));
%         else
%             max2 = 0;
%         end
%         P(k+1) = max(max1,max2);
%         %         infmax = max(inf(max1),inf(max2));
%         %         supmax = max(sup(max1),sup(max2));
%         %         P(k+1) = infsup(infmax,supmax);
%     end
% end

N = length(a)-1;
PL = zeros(2*N+1,1);
PR = zeros(2*N+1,1);
a_ext = [zeros(N,1);flip(a(2:end));a;zeros(N,1)];
for k = -N:N
    % Left
    if k-N <= -N-1
        i = 1;
        tempL = zeros(k+2,1);
        for j = k-N:-N-1
            tempL(i) = abs(a_ext(k-j+2*N+1))/(nu^(abs(j)));
            i = i+1;
        end
        PL(k+N+1) = max(tempL);
    end
    
    % Right
    if N+1 <= k+N
        i = 1;
        tempR = zeros(k+2,1);
        for j = k-N:-N-1
            tempR(i) = abs(a_ext(k-j+2*N+1))/(nu^(abs(j)));
            i = i+1;
        end
        PR(k+N+1) = max(tempL);
    end
end
P = max([PL, PR],[],2);
P = [P(N+1);max([P(N+2:end),flip(P(1:N))],[],2)];
end

