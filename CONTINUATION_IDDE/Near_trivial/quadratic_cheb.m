function [s]=quadratic_cheb(a,b)
a = [flip(a);a(2:end)];
b = [flip(b);b(2:end)];
m=(length(a)+1)/2;
B=zeros(2*m-1);
tb=[zeros(m-1,1);flipdim(b,1);zeros(m-1,1)].';
for k=-m+1:m-1
    B(k+m,:)=tb((-m+1:m-1)-k+(2*m-1));
end
s=B*a;
N = (length(s)+1)/2; 
s = s(N:end);
end