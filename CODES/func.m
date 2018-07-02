%func.m file for variable input

%Input Variables
%sigma : volitility
%N     : number of time steps
%T     : Trading Horizon Length
%gamma : Index for utility function (exponential in this case)


delta = [1 1 0]';
a = [-1 0 1]';
a0 = 0;
sigma = [0.200 0 0;0.0375 0.1452 0;0.0250 0.0039 0.0967];
T = 1;
N = 2^8;
gamma = 0.6;

%run the main function
simufunct(delta,a,a0,sigma,N,T,gamma)