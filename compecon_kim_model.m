%compecon_kim_model.m
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = compecon_kim_model

%This program computes a second-order  approximation to the policy functions 
%of the simple dynamic stochastic disequilibrium model 
%developed by ''Sequential Action and Beliefs under Partially Observable
%DSGE Environments '' Seong-Hoon Kim (2011)

%Companion files: compecon_kim_model.m;  compecon_kim_ss.m; 
%and computation codes from Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775),
%downloadable from https://www.columbia.edu/~mu2166/2nd_order.htm

%compecon_kim_run.m is the main execution file for the present model.

% The function f is of the following form as a collection of dynamic
% equations; E_t f(yp,y,xp,x) =0. 

%For details, refer to Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775).

%Define parameters
syms ALFA BETTA GAMA DELTA RHO RHOV

%Define variables 
syms c cp n np k kp a ap v vp  
% All notations are as in the paper, except a and ap 
% we use a (ap) for z (zp ) in the paper, because zp is a matlab built-in function).

%Write equations fi, i=1:4
f1 = - c^(-1) + BETTA * cp^(-1) * (ap * vp * ALFA * kp^(ALFA-1) * np^(1-ALFA) + 1 - DELTA);
f2 = -GAMA*(1-n)^(-1) + (1-ALFA) * c^(-1) * a * v * k^ALFA * n^(-ALFA);
f3 = c + kp - (1-DELTA) * k - a * v * k^ALFA * n^(1-ALFA);
f4 = - log(ap) + RHO * log(a);
f5 = - log(vp) + RHOV * log(v);

%Create function f
f = [f1;f2;f3;f4;f5];

% Define the vector of controls, y, and states, x
x = [k a v];
y = [c n];
xp = [kp ap vp];
yp = [cp np];

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);

