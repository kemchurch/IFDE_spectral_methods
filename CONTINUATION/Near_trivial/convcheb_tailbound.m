function ab = convcheb_tailbound(a,nu)
N = length(a)-1;
a = [a;intval(zeros(N,1))];
ab = intval(zeros(N+1,1));
for k=1:N
    ab(k) = max(a(N+2-k:N+1)./(2*nu.^(N+1:1:N+k)'));
end
end