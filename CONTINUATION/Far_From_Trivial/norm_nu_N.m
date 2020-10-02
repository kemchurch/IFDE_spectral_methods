function [norm_B] = norm_nu_N(B,nu,N)
vec_nu = nu.^(0:N)';
vec_nu(2:end) = 2*vec_nu(2:end);
norm_B = max(sum(vec_nu.*abs(B))./(vec_nu)');
end

