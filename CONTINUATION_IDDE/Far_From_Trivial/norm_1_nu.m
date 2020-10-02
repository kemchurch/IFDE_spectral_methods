function [Norm_nu] = norm_1_nu(a,nu)
N = length(a(:,1))-1;
vect_nu = nu.^((0:N)');
vect_nu(2:end) = 2*vect_nu(2:end);
Norm_nu = sum(abs(a).*vect_nu);
end

