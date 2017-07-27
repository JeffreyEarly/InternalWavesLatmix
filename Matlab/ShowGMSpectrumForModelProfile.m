Lz = 5000;
latitude = 33;
z = linspace(-Lz,0,1000)';

nEVP = 256;
nGrid = 2^14+1;

f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the Latmix stratification profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes for this profile are here:
% /Users/jearly/Dropbox/Documents/Notes/latmix-glider-profiles/ADifferentAnalyticalProfile/ADifferentAnalyticalProfile.tex

rho0 = 1025;
g = 9.81;

delta_p = 0.9;
L_s = 8;
L_d = 100;
z_p = -17;
D = Lz; % intended to be 5000m

N0 = 2.8e-3; % Surface
Nq = 1.4e-2; % Pycnocline portion to exponentials
Np = 4.7e-2;
Nd = 1.75e-3;

A = Np*Np - Nq*Nq;
B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
C = N0*N0-B*exp(-2*z_p^2/L_s);
E = (Nq*Nq - Nd*Nd)/( 1 - exp(-2*(-D-z_p)^2/L_d^2) );
F = Nd*Nd - E*exp( -2*(-D-z_p)^2/L_d^2 );

rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * sqrt(pi/2) *( erf(sqrt(2)*z_p/L_s) - erf(-sqrt(2)*(z-z_p)/L_s) ) - C*z/g);
rho_deep = @(z) rho_surface(z_p) - (rho0*L_d*E/(2*g)) * sqrt(pi/2) * (erf(sqrt(2)*(z-z_p)/L_d)) - rho0*F*(z-z_p)/g;
rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
rho = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p).*rho_deep(z) + rho_p(z);

% N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
% N2_deep = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
% N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
% N2 = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute measures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,[-Lz 0],latitude, 'nEVP', nEVP, 'nGrid', nGrid);
end
Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.HorizontalIsopycnalVariance(z);
Ew = GM.HorizontalVerticalVelocityVariance(z);
N2 = GM.N2(z);
N = sqrt(N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build exponential profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 5.2e-3;
L_gm = 1.3e3;
rho = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));

if ~exist('GMExp','var')
    GMExp = GarrettMunkSpectrum(rho,[-Lz 0],latitude);
end
Euv_exp = GMExp.HorizontalVelocityVariance(z);
Eeta_exp = GMExp.HorizontalIsopycnalVariance(z);
Ew_exp = GMExp.HorizontalVerticalVelocityVariance(z);
N_exp = sqrt(GMExp.N2(z));

figure
subplot(1,3,1)
plot(1e4*Euv,z), hold on
plot(1e4*Euv_exp,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,3,2)
plot(Eeta,z), hold on
plot(Eeta_exp,z)
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
subplot(1,3,3)
plot(1e4*(Euv + Ew + N2.*Eeta)/2,z), hold on
plot(1e4*(Euv_exp + Ew_exp + GMExp.N2(z).*Eeta_exp)/2,z)
xlabel('total energy (cm^2/s^2)')
ylabel('depth (m)')
legend('Exponential profile (GM reference)', 'Constant stratification')

figure
subplot(1,2,1)
plot(1e4*Euv.*(N0./N),z), hold on
plot(1e4*Euv_exp.*(N0./N_exp),z)
vlines(44,'k--')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
title('wkb scaled')
subplot(1,2,2)
plot(Eeta.*(N/N0),z),hold on
plot(Eeta_exp.*(N_exp/N0),z)
vlines(53,'k--')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
title('wkb scaled')

E = 0.5*(Euv + Ew + N2.*Eeta);
Etotal = trapz(z,E);

return

figure

subplot(2,2,[1 2])
omega = linspace(-N0,N0,200);
S = GM.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,S), ylog
ylim([1e-4 1e2])
xlim(1.05*[-N0 N0])
title('horizontal velocity spectra')
xlabel('radians per second')

subplot(2,2,3)
omega = linspace(0,N0,500);
Siso = GM.HorizontalIsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Siso), ylog, xlog
Sref = omega.^(-2); Sref(omega<f0) = 0; refIndex = find(omega>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
hold on, plot(omega,Sref,'k','LineWidth',2)
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('horizontal isopycnal spectra')
xlabel('radians per second')

subplot(2,2,4)
omega = linspace(0,N0,500);
Sw = GM.HorizontalVerticalVelocitySpectrumAtFrequencies(linspace(-500,0,20),omega);
plot(omega,Sw), ylog, xlog
ylim([1e-6 1e-1])
xlim(1.05*[0 N0])
title('horizontal vertical velocity spectra')
xlabel('radians per second')


return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
