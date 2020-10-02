% Script to compute and validate the branch on the parameter interval 
% [beta*-eps,beta0] for beta0 the point where continuation with the base 
% map fails.
load('data_h0_2_h1_proven.mat')
a = A(:,end);
beta = B(:,end);
fail = 0;
u1 = Q+1;
nu = 1.4;
pts = 30;
disp('Data loaded. Transforming data and computing the branch.')
[a_tr,eps_tr]=desingularize_transform(a,N,u1); % Convert data
[acont,epscont,betacont,~] = scriptF_continue_beta(a_tr,eps_tr,pts,N,...
    alpha,beta,inf(exp(-intval(alpha*u1/u2)))-10^-3,u2,delta,u1); %Compute approximate branch in double
disp('Branch computed.')
if epscont(1)*epscont(pts)<0
    disp('Quasi-amplitude crosses through zero along the numerical branch. To be verified later via interval methods.')
else
    disp('Quasi-amplitude does not seem to cross through zero along the numerical branch.')
    return
end
disp(['Attempting validation at nu = ',num2str(nu)]);
R = zeros(3,pts-1);
for n=1:pts-1
   [~,~,~,~,R(:,n),~,~] = bounds_continuation(acont(:,n),acont(:,n+1),...
       epscont(n),epscont(n+1),betacont(n),betacont(n+1),...
       N,alpha,u2,delta,u1,nu);
   if R(3,n)>0
       disp(['Success at branch segment n = ',num2str(n)])
   else
       disp(['Failure at branch segment n = ',num2str(n)])
       fail = 1;
       break
   end
end
if epscont(pts)+R(2,pts-1)>=0 || epscont(1)-R(2,1) <= 0
    disp('Can not prove that quasi-amplitude crossed through zero.')
    fail = 1;
elseif fail == 1
    disp('Can not prove that quasi-amplitude crossed through zero.')
else
    disp('Successfully proved the quasi-ampliutude crossed through zero.')
end
if fail == 0
    disp('Proof successful!')
    disp(['The branch has been validated for beta in the interval [',num2str(inf(exp(-intval(alpha*u1/u2)))-10^-4),',',num2str(beta+10^-4),'].'])
end
if fail == 0
    disp('Computing existence interval for the branch connection.')
    [~,~,~,~,R_connection,~,~] = bounds_continuation(intval(acont(:,1)),intval(acont(:,1)),...
    epscont(1),epscont(1),betacont(1),betacont(1),N,alpha,u2,delta,u1,nu);
    disp(['Smallest existence interval: ', num2str(R_connection(2))])
    center = acont(:,1)*epscont(1);
    na = 0;
    for k = 0:u1-1
        na = max(na,norm_1_nu(intval(acont(1+k*(N+1):(k+1)*(N+1))),nu));
    end
    radius = R_connection(2)^2 + R_connection(2)*(na+abs(intval(epscont(1))));
    disp('Center and radius of ball in original (X) coordinates stored as variables (center) and (radius)')
end
