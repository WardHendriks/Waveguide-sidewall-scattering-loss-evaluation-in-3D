% Please read "README.md" before runing this code

% The example computes scattering loss from TE00 mode within AlN waveguide
% WG geometry: H = 0.5 um, W = 1.2 um
clear; %clc;
FunctionsPath = genpath('Functions');
addpath(FunctionsPath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nd2=4; % Number of width to iterate, setting to 1 to disable the loop
nh=1; % Number of height to iterate, setting to 1 to disable the loop

lambda = 405e-9; % wavelength defined in m

% waveguide
height=70e-9;%linspace(100,100,nh).*1e-9;                                          % Waveguide height 
d2record=linspace(400,700,nd2).*1e-9;                                       % Waveguide width 

% Refractive indices:
n1 = 1.474196;                                                             % lower cladding
n2 = 1.741319;                                                             % Core
n3 = 1.468748;                                                             % upper cladding

ea=n2.^2; % epsilon of core
es=n1.^2; % epsilon of cladding

% Roughness
sigma=4.2*1e-9; % sidewall roughness, usually 3-5 nm
Lc=50*1e-9; % correlation length, usually 50-100 nm

np=100; % number of points used to describe E field 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dB_sidewall=zeros(nh,nd2);                                                  % prealocate result matrix
C=physconst('lightspeed');

for ii=1:nh
for jj=1:nd2
as_ratio=d2record(jj)/(height(ii));
d1=0;
d2=d2record(jj); % d2-d1 equal to the waveguide width, in um.
r=10000*lambda;
omega=2*pi*C/lambda;
u=4*pi*1e-7;
sl=d2/4;
delta_epsilon_dl=1; %It defines the pertubation of epsilon
delta_epsilon_sd=ea-es;
dS_dl=1*10^(-18);
dS_sd=dS_dl;
% To define current density, the mode is calculated in this section
%%%%%%%%%%%%%%%%%%%% mode calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmodes=12; % TE and TM
TE=1;
if as_ratio>=1
    if TE==1
        order=1;
    else
        order=2;
    end
else
    if TE==1
        order=2;
    else
        order=1;
    end
end
%order=1 % NOTE: we can get certain mode by either setting "order" or symmetry from "000S" to "000A"
pole=false; % This term is related to the simplification of modal. In most of the situation, it should be false.
% 1st order, TE
% 2nd order, TM
% Layer heights:
h2 =height(ii)*1e6;           % Core thickness
h3 = 2;           % Upper cladding

% Horizontal dimensions:
rh = h2;           % Ridge height
%rw=1.5;
rw = d2*1e6/2;           % Ridge half-width
%THIS IS RELATED TO d2!
side = 2;         % Space on side
% Grid size:
dx = 0.005;        % grid size (horizontal)
dy = 0.005;        % grid size (vertical)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,xc,yc,nx,ny,eps,edges] = ...
waveguidemeshfull([n1,n2,n3],[side,h2,h3],rh,rw, ...
                  side,dx,dy);
[Hx_mode,Hy_mode,ntemp] = wgmodes(lambda*10^(6),n2,nmodes,dx,dy,eps,'0000');
neff=ntemp(order);
beta=neff*2*pi/lambda;
[Hz_mode,Ex_mode,Ey_mode,Ez_mode] = postprocess (lambda*10^(6), neff, Hx_mode(:,:,order), Hy_mode(:,:,order), dx, dy, eps, '0000');
fprintf(1,'neff = %.6f\n',neff);

%% PLOT mode
% figure(1);
% subplot(131);
% contourmode(x,y,Ex_mode(:,:));
% title('Ex');
% for v = edges, line(v{:}); end
% subplot(132);
% contourmode(x,y,Ey_mode(:,:));
% title('Ey');
% for v = edges, line(v{:}); end
% subplot(133);
% contourmode(x,y,Ez_mode(:,:));
% title('Ez');
% for v = edges, line(v{:}); end

%%%comput the power of mode
[nx_mode,ny_mode] = size(Ex_mode);
S_mode=zeros(nx_mode,ny_mode);
Z=377/2;
for iiii=1:nx_mode
    for jjjj=1:ny_mode
        if TE==1
            E=Ex_mode(iiii,jjjj);
        else
            E=Ey_mode(iiii,jjjj);
        end
        S_temp=0.5*real(E*conj(E/Z));
        S_mode(iiii,jjjj)=S_temp;
    end
end
P=sum(sum(S_mode))*dx*dy*10^(-12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[E_sd,S_sd,Sr_sd,power_p_sd,power_sr_sd,ratio_sd,Er_sd,P_dB ]=farfield_v4(2,Ex_mode,Ey_mode,Ez_mode,sl,dx,dy,h2,side,delta_epsilon_sd,false,dS_sd,np,r,d1,d2,es,ea,lambda,omega,u,pole,sigma,Lc,beta,0);
disp(['power_p_sd=',num2str(power_p_sd)]);
disp(['power_sr_sd=',num2str(power_sr_sd)]);

P_dB=P_dB*10^(27);
power_sd_rc=power_p_sd;
power_sd_rc=power_sd_rc';
dB_sidewall(ii,jj)=-2*10*log10((P-P_dB)/P)*10^(7);

end
end
dB_sidewall_real=dB_sidewall*(ea-es);
disp(['loss=',num2str(dB_sidewall_real),'dB/cm']);

figure;
plot(d2record*1e9,dB_sidewall_real)
xlabel('waveguide width [nm]')
ylabel('loss [dB/cm]')