function [] = run_proofs
% This function completes the proof of Theorem 10. It will prompt you for
% an integer 1,2,...8 for the relevant proof.
PROOFNUMBER = input('Which Proof # do you want to run? Enter an integer 1,2,...,8. : ');
if PROOFNUMBER==1
    load('test_params_1.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==2
    load('test_params_2.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==3
    load('test_params_3.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==4
    load('test_params_4.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==5
    load('test_params_5.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==6
    load('test_params_7.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
elseif PROOFNUMBER==8
    load('test_params_8.mat');
    main(tau,r,T,N,nu,delta_s_far,delta_s_near,bif_tol)
else
    disp('Invalid proof number.');
end
end