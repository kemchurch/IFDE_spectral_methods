%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the monodromy operator for the system
% x' = A*x(t) + B*x(t-2),           t\notin Z
% \Delta x = C1*x(t^-) + C2*x(t-2), t\in Z
% for the matrices A=[0,1;-1,0], B=[0.8,-0.2;0,-0.1];
% C1 = [-0.25,0;-0.35,0]; C2 = [0.1,0;0.1,0.3];
% It then plots the eigenvalues (Floquet multipliers) and a 
% unit circle for scale (all eigenvalues inside the circle ->
% stability). You can adjust the number of modes (N) in the 
% script at the appropriate place.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[0,1;-1,0]; B=[0.8,-0.2;0,-0.1];
C1 = [-0.25,0;-0.35,0]; C2 = [0.1,0;0.1,0.3];
N = 20; %Adjust the number of modes here if you want.
p=1; %period
q=2; %delay
M=M_truncation(A,B,C1,C2,p,q,N);
plot(eig(M),'o');
hold on
u = 0:1e-2:2*pi;
plot(cos(u),sin(u));