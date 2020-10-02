function [] = script_compute_validate_desingularized_function(A,N,alpha,B,u2,delta,Q,nu,Y0_MAT,Z0_MAT,Z1_MAT,Z2_MAT,R_MAT,r,T,tau,delta_s_near,bif_tol)
% Script to compute and validate the branch on the parameter interval 
% [beta*-tol,beta0] for beta0 the point where continuation with the base 
% map fails. Stock tolerance: 1E-3
a = A(:,end);
beta = B(:,end);
fail = 0;
u1 = Q+1;
disp('Transforming data and computing the branch near beta^star.')
pts = ceil(abs(beta-(exp(-(alpha*u1/u2))-bif_tol))/delta_s_near);
[a_tr,eps_tr]=desingularize_transform(a,N,u1); % Convert data
[a_cont,eps_cont,beta_cont,~] = scriptF_continue_beta(a_tr,eps_tr,pts,N,...
    alpha,beta,inf(exp(-intval(alpha*u1/u2)))-bif_tol,u2,delta,u1); %Compute approximate branch in double
disp('Branch computed.')
if eps_cont(1)*eps_cont(pts)<0
    disp('Quasi-amplitude crosses through zero along the numerical branch. To be verified later via interval methods.')
else
    disp('Quasi-amplitude does not seem to cross through zero along the numerical branch.')
    return
end
disp(['Attempting to validate the branch.']);
R_MAT_cont = zeros(3,pts-1);
Y0_MAT_cont = zeros(1,pts-1);
Z0_MAT_cont = zeros(1,pts-1);
Z1_MAT_cont = zeros(1,pts-1);
Z2_MAT_cont = zeros(pts-1,4);
for n=1:pts-1
    [Y0_MAT_cont(n),Z0_MAT_cont(n),Z1_MAT_cont(n),Z2_MAT_cont(n,:),R_MAT_cont(:,n),~,~] = bounds_continuation(a_cont(:,n),a_cont(:,n+1),...
        eps_cont(n),eps_cont(n+1),beta_cont(n),beta_cont(n+1),...
        N,alpha,u2,delta,u1,nu);
    if R_MAT_cont(3,n)>0
        diary myDiaryFile
        %disp(['Success at branch segment n = ',num2str(n)])
        diary off
        subplot(1,2,2);
        plot(beta_cont(n:n+1),abs(a_cont(1,n:n+1)),'b.-');hold on;
        xlabel('\beta');
        ylabel('|a_0|');
        set(gca,'Color','k')
        xlim([beta_cont(1,end),beta_cont(1,1)])
        ylim([abs(a_cont(1,end)),abs(a_cont(1,1))])
        drawnow
    else
        disp(['Failure at branch segment n = ',num2str(n)])
        fail = 1;
        break
    end
end
if eps_cont(pts)+R_MAT_cont(2,pts-1)>=0 || eps_cont(1)-R_MAT_cont(2,1) <= 0
    disp('The enclosure of the branch was not tight enough to prove a zero-crossing. Trying to prove the crossing by validating first and last points on the branch as isolated solutions.')
    [~,~,~,~,R_start,~] = bounds_continuation(a_cont(:,1),a_cont(:,1),...
        eps_cont(1),eps_cont(1),beta_cont(1),beta_cont(1),...
        N,alpha,u2,delta,u1,nu);
    [~,~,~,~,R_end,~] = bounds_continuation(a_cont(:,pts),a_cont(:,pts),...
        eps_cont(pts),eps_cont(pts),beta_cont(pts),beta_cont(pts),...
        N,alpha,u2,delta,u1,nu);
    if eps_cont(pts)+R_end(2)>=0 || eps_cont(1)-R_start(2)<=0
        disp('Failed to prove that the branch crossed through zero. Try a increasing bif_tol.');
        fail = 1;
    else
    end
elseif fail == 1
    disp('The branch validation failed.')
else
    disp('Successfully proved the quasi-ampliutude crossed through zero.')
end
if fail == 0
    disp('Proof successful!')
    disp(['The branch has been validated for beta in the interval [',num2str(inf(exp(-intval(alpha*u1/u2)))-10^-4),',',num2str(beta+10^-4),'].'])
end
if fail == 0
    disp('Computing existence interval for the branch connection.')
    [~,~,~,~,R_connection,~,~] = bounds_continuation(intval(a_cont(:,1)),intval(a_cont(:,1)),...
    eps_cont(1),eps_cont(1),beta_cont(1),beta_cont(1),N,alpha,u2,delta,u1,nu);
    disp(['Smallest existence interval: ', num2str(R_connection(2))])
    center = a_cont(:,1)*eps_cont(1);
    na = 0;
    for k = 0:u1-1
        na = max(na,norm_1_nu(intval(a_cont(1+k*(N+1):(k+1)*(N+1))),nu));
    end
    radius = R_connection(2)^2 + R_connection(2)*(na+abs(intval(eps_cont(1))));
    disp('Center and radius of ball in original (X) coordinates stored as variables (center) and (radius)')
else
    center = a_cont(:,1)*eps_cont(1);
    radius = NaN;
end
Z2_MAT_cont = Z2_MAT_cont';
save('Data\data_h1_2_beta_star_proven','a_cont','beta_cont','eps_cont','N','Q','u2','delta','alpha','r','tau','T','nu',...
     'Y0_MAT_cont','Z0_MAT_cont','Z1_MAT_cont','Z2_MAT_cont','R_MAT_cont','center','radius');
end