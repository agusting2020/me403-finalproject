%
% compute_G.m
%
% Computes a scalar transfer function from matrix equations
%
% Written by: Hua Wang
%
% Last Change: November 16, 2004
% 

%
% Compute random matrices for A(4x4) and B(4x2)having the same
% dimension as aircraft matrices
% 

% A = rand(4,4);
% B = rand(4,2);

%
% Define C and D so that y=x
%
% Note that A,B,C,D define a system in state-space form
%

C = eye(4,4);
D = zeros(4,2);

%
% Define input element (row of u)
%

i_v = 1;


%
% Convert to transfer function form % 
%

[NUM, DEN] = ss2tf(A, B, C, D, i_v);

%
% Define output element % 
%

o_v = 1;

%
% Create transfer function from from input iu to output ou % 
%

G = tf(NUM(o_v,:), DEN);



