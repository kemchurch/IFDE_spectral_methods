function [s] = conv_cheb(a1,a2)
N=length(a1);

a1=[flip(a1(2:N));a1];
a2=[flip(a2(2:N));a2];

ta1=[zeros(N,1);a1;zeros(N,1)]; tu1=ifft(ifftshift(ta1));
ta2=[zeros(N,1);a2;zeros(N,1)]; tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=real((4*N-1)*F(2*N:3*N-1));
end

