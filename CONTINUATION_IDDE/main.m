function runtime = main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
close all
% input: 
%       - tau,r,T : model parameters
%       - N :  Number of nodes minus 1
%       - nu : sequence space weight >=1 for the proof.
%       - delta_s_far :  step size far away from bifurcation
%       - delta_s_near :  step size near the bifurcation
%       - bif_tol: amount to cross bifurcation
tic;
beta = 1;
alpha = r*tau;
u  = T/tau;
[u1,u2] = rat(u);
if u1 >= u2
    delta = u1-u2;
else
    k=2;
    while k*u1-u2<0
        k = k+1;
    end
    delta = k*u1 - u2;
end
Q = u1-1;
a = zeros(N+1,1);
a(1) =1;
a = repmat(a,Q+1,1);

%Computation Away from beta^star discrete solution
addpath('.\Far_From_Trivial')
disp('Computing points of the branch away from beta^star.')
discrete_solution_unproven(N,alpha,beta,u2,delta,Q,a,r,T,tau,nu,delta_s_far);
load('Data\h0_2_h1_discrete_solution_unproven.mat');
disp('Proving branches away from beta^star.')
%Computation Away from beta^star branch
proof_para_cont(A,N,alpha,B,u2,delta,Q,nu,r,T,tau)
load('Data\data_h0_2_h1_proven.mat');
disp('For branches aways from \beta^star ')
disp(['Smallest existence interval: ', num2str(max(R_MAT(1,:)))])
rmpath('.\Far_From_Trivial')
%Computation Away from beta^star
addpath('.\Near_Trivial')
script_compute_validate_desingularized_function(A,N,alpha,B,u2,delta,Q,nu,Y0_MAT,Z0_MAT,Z1_MAT,Z2_MAT,R_MAT,r,T,tau,delta_s_near,bif_tol);
rmpath('.\Near_Trivial')
%Branch connection
load('Data\data_h0_2_h1_proven.mat');
R_far_last_largest = R_MAT(2,end);
load('Data\data_h1_2_beta_star_proven.mat');
if isnan(radius)
    disp('Proof failed.');
else
    disp('Checking branch connection.');
    if intval(radius) < intval(R_far_last_largest)
        disp('Success! The branches are connected. Proof successful.');
    else
    disp('I can not prove that these branches are connected.');
    end
end
runtime = toc;
disp(['Proof runtime: ',num2str(runtime),' seconds.']);
end