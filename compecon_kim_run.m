%compecon_kim_run.m

%This program computes a second-order  approximation to the policy
%functions 
%of the simple dynamic stochastic disequilibrium model 
%developed by ''Sequential Action and Beliefs under Partially Observable
%DSGE Environments '' Seong-Hoon Kim (2011)

%Companion files: compecon_kim_model.m;  compecon_kim_ss.m; 
%and computation codes from Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775),
%downloadable from https://www.columbia.edu/~mu2166/2nd_order.htm

%compecon_kim_run.m is the main execution file for the present model.

%For details, refer to Schmitt-Grohe and Uribe (2004, JEDC, pp. 755-775).

clear;

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = compecon_kim_model;

%% Quantitative Example

%% The Fundamentals

%Steady State and Parameter Values
[ALFA,BETTA,GAMA,DELTA,RHO,RHOV,c,cp,n,np,k,kp,a,ap,v,vp,C,N,K,A,V]=compecon_kim_ss;

%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

%Second-order approximation
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 

%System coefficients

Azw = [1 0; 0 1]; % <------ assumed at the first place
nz = length(Azw); % number of exogenous shocks
[nc,nx] = size(gx); % nc, number of controls. nx, number of states including predetermined variables, i.e., x --> theta
nk = nx-nz;

AA = [nfyp, nfxp];
BB = [nfy, nfx];

Axx = AA(nc+1:size(AA,1), nc+1:size(AA,2));
Akk = Axx(1:nk,1:nk);
Azz = Axx(nk+1:nx,nk+1:nx);
Bxc = BB(nc+1:size(BB,1), 1:nc);
Bxlc = -inv(Axx)*Bxc;
Bklc = Bxlc(1:nk,:);
Bxx = BB(nc+1:size(BB,1), nc+1:size(BB,2));
Bkz = Bxx(1:nk,nk+1:nx);
Bzz = Bxx(nk+1:nx,nk+1:nx);
Bxlx = -inv(Axx)*Bxx;
Bklz = -inv(Akk)*Bkz;
Bzlz = -inv(Azz)*Bzz;
Axw = [zeros(nk,nz); Azw];
eta = -inv(Axx)*Axw; % matrix defining driving force

[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);

%Optimal control
Hcx = gx;
Hxx = hx;
Hxw = eta;

%Observation channel
Byy = -eye(1);
ny = length(Byy);
Byx = [ALFA 1 1];
Bylx = -inv(Byy)*Byx;
Byc = [0 1-ALFA];
Bylc = -inv(Byy)*Byc;

%Shock base
std1 =0.05;
std2 = 0.05;
covw = [std1^2 0; 0 std2^2];


%% The POMDP version of the Neoclassical Growth Model

Byx = [ALFA 1 1];
Bylx = -inv(Byy)*Byx;

%coverging setup
ddx = 1;
itx = 1;
maxit = 1e3;

%initial condition
covxhat = 0.1*eye(nx); % covxhat is SIGMA_hat(t+1|t) in the work, covxifr is SIGMA(t|t)
mx = covxhat * Bylx'*inv((Bylx*covxhat*Bylx'));

%iteration
while ddx>1e-7
     itx = itx + 1; 
     last_mx(:,itx) = mx;
     last_covxhat(:,:,itx) = covxhat;
     last_Mx(:,itx) = last_mx(:,itx)*inv(eye(ny) + Bylc*Hcx*last_mx(:,itx));
     mx = covxhat * Bylx'*inv((Bylx*covxhat*Bylx'));
     covxifr = (eye(nx) - mx*Bylx)*covxhat;
     covxhat = Bxlx*covxifr*Bxlx' + Hxw*covw*Hxw';
     
     dd_mx = max(max(abs(mx-last_mx(:,itx))));
     dd_covxhat = max(max(abs(covxhat - last_covxhat(:,:,itx))));
     ddx = max(dd_covxhat,dd_mx);
     
     if itx > maxit
         error('mx failed to converge');
     end
end

text = ['mx iteration converged in ' num2str(itx) ' iterations' ];
disp(text);

% The optimal inference-forecast solution to the POMDP
Mx = mx*inv(eye(ny) + Bylc*Hcx*mx);
mx;
Mx;
covxifr;
covxhat;
Hxx*Mx;


%% The HMM version of the Neoclassical Growth Model

Hcz = Hcx(:,nk+1:nx);
Hck = Hcx(:,1:nk);
Hzz = Hxx(nk+1:nx,nk+1:nx); % = Bzlz
Hkk = Hxx(1:nk,1:nk);
Hzw = Hxw(nk+1:nx,:);

Byk = Byx(:,1:nk);
Bylk = -inv(Byy)*Byk;
Byz = Byx(:,nk+1:nx);
Bylz = -inv(Byy)*Byz;

%coverging setup
ddz = 1;
itz = 1;

%initial condition
covzhat = 0.1*eye(nz); % covxhat is SIGMA_hat(t+1|t) in the work, covxifr is SIGMA(t|t)
mz = covzhat * Bylz'*inv((Bylz*covzhat*Bylz'));

%iteration
while ddz>1e-7
     itz = itz + 1; 
     last_mz(:,itz) = mz;
     last_covzhat(:,:,itz) = covzhat;
     last_Mz(:,itz) = last_mz(:,itz)*inv(eye(ny) + Bylc*Hcz*last_mz(:,itz));
     mz = covzhat * Bylz'*inv((Bylz*covzhat*Bylz'));
     covzifr = (eye(nz) - mz*Bylz)*covzhat;
     covzhat = Bzlz*covzifr*Bzlz' + Hzw*covw*Hzw';
     
     dd_mz = max(max(abs(mz-last_mz(:,itz))));
     dd_covzhat = max(max(abs(covzhat - last_covzhat(:,:,itz))));
     ddz = max(dd_covzhat,dd_mz);
     
     if itz > maxit
         error('mz failed to converge');
     end
     
end

text = ['mz iteration converged in ' num2str(itz) ' iterations' ];
disp(text);

% The optimal inference-forecast solution to the HMM
Mz = mz*inv(eye(ny) + Bylc*Hcz*mz);
mz;
Mz;
covzifr;
covzhat;


%% Impulse Response Analysis

irdate = 101;

%%% initial conditions

% the POMDP version
xtr(:,1) = [0;0;0];
xhat(:,1) = [0;0;0];
yhat(:,1) = 0;
yobs(:,1) = inv(eye(ny)-Bylc*Hcx*Mx)*(Bylx*xtr(:,1) + Bylc*Hcx*xhat(:,1) - Bylc*Hcx*Mx*yhat(:,1));
xifr(:,1)= [1;0;0];%xhat(:,1) + Mx*(yobs(:,1)-yhat(:,1));
cnt(:,1) = Hcx*xifr(:,1);
chat(:,1) = Hcx*xhat(:,1);


%%% impulse response

for t=2:irdate
         
    % POMDP
    xtr(:,t) = Bxlc*Hcx*xifr(:,t-1) + Bxlx*xtr(:,t-1);
    xhat(:,t) = Hxx*xifr(:,t-1);
    yhat(:,t) = (Bylc*Hcx+Bylx)*xhat(:,t);
    yobs(:,t) = inv(eye(ny)-Bylc*Hcx*Mx)*(Bylx*xtr(:,t) + Bylc*Hcx*xhat(:,t) - Bylc*Hcx*Mx*yhat(:,t));
    xifr(:,t)= xhat(:,t) + Mx*(yobs(:,t)-yhat(:,t));
    cnt(:,t) = Hcx*xifr(:,t);
    chat(:,t) = Hcx*xhat(:,t);
        
end

close all;
t = 0:1:irdate-1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

subplot(2,3,1)
plot(t,cnt(1,:),'-b','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('consumption')
ylabel('% deviation from ss')
grid on

subplot(2,3,2)
plot(t,cnt(2,:),'-b','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('worked hours')
ylabel('% deviation from ss')
grid on

subplot(2,3,3)
plot(t,yobs(1,:),'-b','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('production')
ylabel('% deviation from ss')
grid on

subplot(2,3,4)
plot(t,xtr(1,:),'-b', t,xifr(1,:),'--r','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('capital')
ylabel('% deviation from ss')
grid on

subplot(2,3,5)
plot(t,xtr(2,:),'-b', t,xifr(2,:),'--r','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('technology')
ylabel('% deviation from ss')
grid on

subplot(2,3,6)
plot(t,xtr(3,:),'-b', t,xifr(3,:),'--r','LineWidth',1.5)
axis([0 irdate-1 -1.1 2])
xlabel('transitory component')
ylabel('% deviation from ss')
grid on

figure;clf;hold on; % yy, c, n 
plot(t,xtr(1,:),'-r',t,xtr(2,:),'-b','LineWidth',1.5)
plot(t,xifr(1,:),'-.r',t,xifr(2,:),'-.b','LineWidth',2)
set(gcf,'Position',[100 100 400 400])
legend('true capital','true technology','capital-valuation','technology-estimation')

figure;clf;hold on; % yy, c, n 
plot(t,cnt(1,:),'-ob','MarkerSize',2.5)
plot(t,cnt(2,:),'-r',t,yobs(1,:),'-.k','LineWidth',1.5)
set(gcf,'Position',[100 100 400 400])
legend('consumption','worked hours','output')