%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the monodromy operator for the system
% x' = A(t)*x(t) + B(t)*x(t-2),     t\notin 3Z
% \Delta x = C1*x(t^-) + C2*x(t-2), t\in 3Z
% for A(t)=
% [1,0;0,1] for t\in[0,1]
% [-1,0;-1,0] for t\in[1,2]
% [0,1;1,0] for t\in[2,3]
% and extended periodically, while B(t)=[-1,0;0,-1],
% C1=[-0.5,0;0,0], C2=[-0.2,0;0.2,-0.3];
% Note the period (p) is 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=cell(3,1);
B=cell(3,1);
A{1} = [1,0;0,1];
A{2} = [-1,0;-1,0];
A{3} = [0,1;1,0];
B{1} = [-1,0;0,-1];
B{2} = B{1};
B{3} = B{1};
C1 = [-0.5,0;0,0]; C2 = [-0.2,0;0.2,-0.3];
N = 20; %Adjust the number of modes here if you want.
p=3; %period
q=2; %delay
M=M_truncation(A,B,C1,C2,p,q,N);
plot(eig(M),'o');
hold on
u = 0:1e-2:2*pi;
plot(cos(u),sin(u));
