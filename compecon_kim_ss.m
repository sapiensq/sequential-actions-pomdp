%compecon_kim_ss.m
function [ALFA,BETTA,GAMA,DELTA,RHO,RHOV,c,cp,n,np,k,kp,a,ap,v,vp,C,N,K,A,V]=compecon_kim_ss

%This program computes a second-order  approximation to the policy functions 
%of the simple dynamic stochastic disequilibrium model 
%developed by ''Sequential Action and Beliefs under Partially Observable
%DSGE Environments '' Seong-Hoon Kim (2011)

%Companion files: compecon_kim_model.m;  compecon_kim_ss.m; 
%and computation codes from Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775),
%downloadable from https://www.columbia.edu/~mu2166/2nd_order.htm

%compecon_kim_run.m is the main execution file for the present model.

%For details, refer to Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775).

ALFA=0.36; %capital share
BETTA=0.99; %discount rate
DELTA=0.025; %depreciation rate
RHO=0.99; %persistence of permanent technology shock
RHOV=0; %persistence of noisy, transitory shock

A = 1; %steady-state value of aggregate technology (Z in Section 2 of the paper)
V = 1; %steady-state value of transitory technology  
N = 1/3; %steady-state value of worked hours
K = N*((1/BETTA+DELTA-1)/(ALFA*A*V))^(1/(ALFA-1)); %steady-state value of capital
C = A*V*K^(ALFA)* N^(1-ALFA) - DELTA*K; %steady-state value of consumption 
GAMA = (1-N)*(1-ALFA)*C^(-1)*A*V*K^ALFA*N^(-ALFA); % relative weight on leisure to consumption, with N=1/3 at SS. 


c = log(C);
n = log(N);
k = log(K);
a = log(A); 
v = log(V);


cp = c;
np = n;
kp=k;
ap=a;
vp = v;

